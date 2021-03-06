#include "cfdproblem.h"


CFDProblem::CFDProblem(double LX, double RX, double LY, double RY, int XSize, int YSize, double dtime, double Viscousity)
{
    ThreadNum = 4;

    if(RX>LX)
    {
        lx = LX;
        rx = RX;
    }
    else
    {
        lx = RX;
        rx = LX;
    }

    if(RY>LY)
    {
        ly = LY;
        ry = RY;
    }
    else
    {
        ly = RY;
        ry = LY;
    }

    xSize = XSize;
    ySize = YSize;

    hx = (rx - lx)/(xSize-1);
    hy = (ry - ly)/(ySize-1);

    this->dt = dtime;
    this->nu = Viscousity;

    StreamFunc = new PoissonTaskWDiscreteForce(lx,rx,
                                               ly,ry,
                                               xSize, ySize);
    CurlFunc = new HeatTaskWDiscreteForce(lx,rx,
                                          ly,ry,
                                          xSize, ySize,
                                          dt,nu);

    int i,j;
    vx = new double* [xSize];
    vy = new double* [xSize];
    for(i = 0; i < xSize; i++)
    {
        vx[i] = new double[ySize];
        vy[i] = new double[ySize];
    }
    for(i=0;i<xSize;i++)
        for(j=0;j<ySize;j++)
        {
            vx[i][j] = 0.;
            vy[i][j] = 0.;
        }
}

void CFDProblem::SetInitialConditions(double AngleOfAttack)
{
    int i,j;
    //Setting boundary distributions for velocity and stream function
    for(i=0; i<xSize; i++)
    {
        StreamFunc->u[i][0]     =   StreamBoundaryFunc((StreamFunc->lx+StreamFunc->hx*i), StreamFunc->ly);
        StreamFunc->u[i][ySize-1]  =   StreamBoundaryFunc((StreamFunc->lx+StreamFunc->hx*i), StreamFunc->ry);

        vx[i][0]    = (StreamBoundaryFunc(StreamFunc->lx+i*hx,
                                          StreamFunc->ly+1e-8)
                      -StreamBoundaryFunc(StreamFunc->lx+i*hx,
                                          StreamFunc->ly-1e-8))/2.e-8;
        vx[i][ySize-1] = (StreamBoundaryFunc(StreamFunc->lx+i*hx,
                                          StreamFunc->ry+1e-8)
                      -StreamBoundaryFunc(StreamFunc->lx+i*hx,
                                          StreamFunc->ry-1e-8))/2.e-8;

        vy[i][0]    =-(StreamBoundaryFunc(StreamFunc->lx+i*hx+1e-8,
                                          StreamFunc->ly)
                      -StreamBoundaryFunc(StreamFunc->lx+i*hx-1e-8,
                                          StreamFunc->ly))/2.e-8;
        vy[i][ySize-1] =-(StreamBoundaryFunc(StreamFunc->lx+i*hx+1e-8,
                                          StreamFunc->ry)
                      -StreamBoundaryFunc(StreamFunc->lx+i*hx-1e-8,
                                          StreamFunc->ry))/2.e-8;

    }
    for(i=0; i<ySize; i++)
    {
        StreamFunc->u[0][i]    =   StreamBoundaryFunc(StreamFunc->lx,(StreamFunc->ly+StreamFunc->hy*i));
        StreamFunc->u[xSize-1][i] =   StreamBoundaryFunc(StreamFunc->rx,(StreamFunc->ly+StreamFunc->hy*i));

        vx[0][i]    = (StreamBoundaryFunc(StreamFunc->lx,
                                          StreamFunc->ly+i*hy+1e-8)
                      -StreamBoundaryFunc(StreamFunc->lx,
                                          StreamFunc->ly+i*hy-1e-8))/2.e-8;
        vx[xSize-1][i] = (StreamBoundaryFunc(StreamFunc->rx,
                                          StreamFunc->ly+i*hy+1e-8)
                      -StreamBoundaryFunc(StreamFunc->rx,
                                          StreamFunc->ly+i*hy-1e-8))/2.e-8;

        vy[0][i]    =-(StreamBoundaryFunc(StreamFunc->lx+1e-8,
                                          StreamFunc->ly+i*hy)
                      -StreamBoundaryFunc(StreamFunc->lx-1e-8,
                                          StreamFunc->ly+i*hy))/2.e-8;
        vy[xSize-1][i] =-(StreamBoundaryFunc(StreamFunc->rx+1e-8,
                                             StreamFunc->ly+i*hy)
                        -StreamBoundaryFunc(StreamFunc->rx-1e-8,
                                            StreamFunc->ly+i*hy))/2.e-8;
    }


    //Setting in-flow boundaries (for solid bodies)
    for(i=1;i<xSize-1;i++)
        for(j=1;j<ySize-1;j++)
        {
            StreamFunc->u[i][j] = StreamBoundaryFunc(StreamFunc->lx+i*StreamFunc->hx,
                                                 StreamFunc->ly+j*StreamFunc->hy);
            if(profileRfunc((StreamFunc->lx+i*StreamFunc->hx),
                            (StreamFunc->ly+j*StreamFunc->hy),
                            M_PI/180.*AngleOfAttack) > 0.)
            {
                StreamFunc->NodeState[i][j] = 1;
                StreamFunc->u[i][j] = 0.5;
            }
        }


    //setting initial and boundary distributions for Vorticity
    for(i=1; i<xSize-1; i++)
        for(j=1; j<ySize-1; j++)
            CurlFunc->u[i][j] =   -((StreamFunc->u[i+1][j]-2.*StreamFunc->u[i][j]+StreamFunc->u[i-1][j])/hx/hx+
                                    (StreamFunc->u[i][j+1]-2.*StreamFunc->u[i][j]+StreamFunc->u[i][j-1])/hy/hy);
    for(i=0; i<xSize; i++)
    {
        CurlFunc->u[i][0]    = CurlFunc->u[i][1];
        CurlFunc->u[i][ySize-1] = CurlFunc->u[i][ySize-2];
    }
    for(i=0; i<ySize; i++)
    {
        CurlFunc->u[0][i]    = CurlFunc->u[1][i];
        CurlFunc->u[xSize-1][i] = CurlFunc->u[xSize-2][i];
    }
}

double PoiseuilleStreamFunc(double t, double p)
{
    if(t<=0.)
        return 0.;
    if(t>0 && t<p)
        return t*t/p/p*(3.-2.*t/p);
    else
        return 1.;
}

double CFDProblem::StreamBoundaryFunc(double x, double y)
{
    return PoiseuilleStreamFunc((y+1.)/2.+0.*x,1.);
}

void CFDProblem::StepInTime()
{
    double hy_s = StreamFunc->hy*StreamFunc->hy;

    int i,j;

    for(i=1; i<xSize-1; i++)
    {
        StreamFunc->u[i][1]     = StreamFunc->u[i][0];
        StreamFunc->u[i][ySize-2]  = StreamFunc->u[i][ySize-1];


        CurlFunc->u[i][0]       = -(StreamFunc->u[i][0]       -2.*StreamFunc->u[i][1]        +StreamFunc->u[i][2]       )/hy_s;
        CurlFunc->u[i][ySize-1] = -(StreamFunc->u[i][ySize-1] -2.*StreamFunc->u[i][ySize-2]  +StreamFunc->u[i][ySize-3] )/hy_s;
    } //условия на границе твёрдого тела (стенки трубы)

    for(i=1; i<ySize-1; i++)
    {
        CurlFunc->u[0][i]    = //-(StreamFunc->u[0][i]        -2.*StreamFunc->u[1][i]   +StreamFunc->u[2][i]      )/hx_s
                                 -(StreamFunc->u[0][i+1]         -2.*StreamFunc->u[0][i]      +StreamFunc->u[0][i-1]          )/hy_s;
        CurlFunc->u[xSize-1][i] = //-(StreamFunc->u[Nx-1][i]     -2.*StreamFunc->u[Nx-2][i]+StreamFunc->u[Nx-3][i]   )/hx_s
                                    -(StreamFunc->u[xSize-1][i+1]   -2.*StreamFunc->u[xSize-1][i]+StreamFunc->u[xSize-1][i-1] )/hy_s;
    } //условия на границе "внешнего" течения (на стоки и источнике жидкости)

    UpdateConvectiveForce();

    CurlFunc->StepInTime_Crank();

    for(i=1;i<xSize-1;i++)
        for(j=1;j<ySize-1;j++)
        {
            StreamFunc->Force[i][j] = -CurlFunc->u[i][j];
        }

    StreamFunc->Iterate(4,4);

    //EHD actuator's part
//    Imitate_EHD_Actuator();
}

void CFDProblem::UpdateConvectiveForce()
{
    int i;
    std::vector<std::thread> ThrPool;

    for(i=0;i<ThreadNum;i++)
    {
        ThrPool.push_back(std::thread(UpdateConvectiveForce_Crutch,this,ThreadNum,i));
    }

    for(i=0;i<ThreadNum;i++)
    {
        ThrPool[i].join();
    }
}

void CFDProblem::UpdateConvectiveForce_Thread(int ThreadNum, int ThreadID)
{
    int i, j;
    int H = (xSize-2)/ThreadNum;
    int Istart  = 1+ThreadID*H;
    int Iend    =   Istart + H;

    double hx_s = StreamFunc->hx*StreamFunc->hx,
           hy_s = StreamFunc->hy*StreamFunc->hy,
           hxhy = 0.25/CurlFunc->hx/CurlFunc->hy;



    for(i = Istart; i < Iend; i++)
    {
        for( j = 1; j< ySize-1; j++)
        {
            if(StreamFunc->NodeState[i][j]==1)
            {
                StreamFunc->u[i][j]  = 0.5;
                CurlFunc->u[i][j]    =  -(StreamFunc->u[i-1][j] -2.*StreamFunc->u[i][j] +StreamFunc->u[i+1][j])/hx_s
                                        -(StreamFunc->u[i][j-1] -2.*StreamFunc->u[i][j] +StreamFunc->u[i][j+1])/hy_s;
            }//условия на твёрдые тела внутри потока


            //консервативная форма переносной силы
            CurlFunc->Force[i][j] =
                    -((CurlFunc->u[i+1][j]*(StreamFunc->u[i+1][j+1]-StreamFunc->u[i+1][j-1])
                      -CurlFunc->u[i-1][j]*(StreamFunc->u[i-1][j+1]-StreamFunc->u[i-1][j-1])) -
                      (CurlFunc->u[i][j+1]*(StreamFunc->u[i+1][j+1]-StreamFunc->u[i-1][j+1])
                      -CurlFunc->u[i][j-1]*(StreamFunc->u[i+1][j-1]-StreamFunc->u[i-1][j-1])))*hxhy;
        }
    }

}

void CFDProblem::UpdateConvectiveForce_Crutch(CFDProblem *Task, int ThreadNum, int ThreadID)
{
    Task->UpdateConvectiveForce_Thread(ThreadNum, ThreadID);
}

void CFDProblem::Imitate_EHD_Actuator()
{
    int i,j;

    for(i=1; i<xSize-1; i++)
        for(j=1; j<ySize-1; j++)
            if(StreamFunc->NodeState[i][j-1]==1 && StreamFunc->NodeState[i][j]==0)
                StreamFunc->u[i][j] = -hy/4.5+StreamFunc->u[i][j-1];
}

void CFDProblem::ParaViewOutput_TecPlotASCII(const char *filename)
{
    int i, j;
    FILE *velocityValue;
    velocityValue = fopen(filename,"w");
    fprintf(velocityValue,"TITLE = \"Flow Model\"\n");
    fprintf(velocityValue,"VARIABLES = \"x\", \"y\", \"v_x\", \"v_y\"\n");
    fprintf(velocityValue,"ZONE T=\"Frame 0\", I=%d, J=%d\n", ySize, xSize);

    for(i=0;i<xSize;i++)
        for(j=0;j<ySize;j++)
            fprintf(velocityValue,"%f %f %f %f\n",hx*i+lx, hy*j+ly, vx[i][j], vy[i][j]);

    fclose(velocityValue);
}

void CFDProblem::ParaViewOutput_NetCDF(const char *filename)
{
    try
    {
        NcFile dataFile(filename, NcFile::replace);

        NcDim xDim = dataFile.addDim("x", xSize);
        NcDim yDim = dataFile.addDim("y", ySize);

        double xLine[xSize], yLine[ySize];

        for(int i = 0; i < xSize; i++)
        {
            xLine[i] = lx + i*hx;
        }

        for(int i = 0; i < ySize; i++)
        {
            yLine[i] = ly + i*hy;
        }

        NcVar xVar = dataFile.addVar("x", ncDouble, xDim);
        NcVar yVar = dataFile.addVar("y", ncDouble, yDim);

        xVar.putVar(xLine);
        yVar.putVar(yLine);

        xVar.putAtt("units","m");
        yVar.putAtt("units","m");

        vector<NcDim> dims;
        dims.push_back(yDim);
        dims.push_back(xDim);
        NcVar vx_arr = dataFile.addVar("vx", ncDouble, dims);
        NcVar vy_arr = dataFile.addVar("vy", ncDouble, dims);
        NcVar BodyRf = dataFile.addVar("Body",  ncInt, dims);


        vx_arr.putAtt("units","m/s");
        vy_arr.putAtt("units","m/s");

        double VX[ySize][xSize];
        double VY[ySize][xSize];
        int    body[ySize][xSize];

        for(int i=1;i<xSize-1;i++)
            for(int j=1;j<ySize-1;j++)
            {
                vx[i][j] =  (StreamFunc->u[i][j+1]-StreamFunc->u[i][j-1])/StreamFunc->hy*0.5;
                vy[i][j] = -(StreamFunc->u[i+1][j]-StreamFunc->u[i-1][j])/StreamFunc->hx*0.5;
            }


        for(int i = 0; i<xSize; i++)
        {
            for(int j = 0; j<ySize; j++)
            {
                VX[j][i] = vx[i][j];
                VY[j][i] = vy[i][j];
                body[j][i] = StreamFunc->NodeState[i][j];
            }
        }


        vx_arr.putVar(VX);
        vy_arr.putVar(VY);
        BodyRf.putVar(body);

        return;
    }
    catch(NcException& e)
    {
        e.what();
        return;
    }
}

void CFDProblem::SetThreadNum(int n)
{
    ThreadNum = n;
    StreamFunc->ThreadNum = n;
    CurlFunc->ThreadNum = n;
}
