#include "cfdproblem.h"


CFDProblem::CFDProblem(double LX, double RX, double LY, double RY, int XSize, int YSize, double dtime, double Viscousity)
{
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
//    double StreamBoundaryFunc;
//    int xSize, ySize;
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
    return PoiseuilleStreamFunc((y+1.)/2.,1.);
}

//virtual double CFDProblem::bodyRfunc(double x, double y, double theta)
//{
//    return 0.04-((x+0.5)*(x+0.5)+y*y);
//}

void CFDProblem::StepInTime()
{
    int i,j;
    for(i=1;i<xSize-1;i++)
        for(j=1;j<ySize-1;j++)
        {
            vx[i][j] =  (StreamFunc->u[i][j+1]-StreamFunc->u[i][j-1])/StreamFunc->hy*0.5;
            vy[i][j] = -(StreamFunc->u[i+1][j]-StreamFunc->u[i-1][j])/StreamFunc->hx*0.5;
        }

    for(i=1; i<xSize-1; i++)
    {
        StreamFunc->u[i][1]     = StreamFunc->u[i][0];
        StreamFunc->u[i][ySize-2]  = StreamFunc->u[i][ySize-1];


        CurlFunc->u[i][0]    = -(StreamFunc->u[i][0]    -2.*StreamFunc->u[i][1]     +StreamFunc->u[i][2]    )/StreamFunc->hy/StreamFunc->hy;
        CurlFunc->u[i][ySize-1] = -(StreamFunc->u[i][ySize-1] -2.*StreamFunc->u[i][ySize-2]  +StreamFunc->u[i][ySize-3] )/StreamFunc->hy/StreamFunc->hy;
    } //условия на границе твёрдого тела

    for(i=1; i<ySize-1; i++)
    {
        CurlFunc->u[0][i]    = //-(StreamFunc->u[0][i]        -2.*StreamFunc->u[1][i]   +StreamFunc->u[2][i]      )/StreamFunc->hx/StreamFunc->hx
                -(StreamFunc->u[0][i+1]      -2.*StreamFunc->u[0][i]   +StreamFunc->u[0][i-1]    )/StreamFunc->hy/StreamFunc->hy;
        CurlFunc->u[xSize-1][i] = //-(StreamFunc->u[Nx-1][i]     -2.*StreamFunc->u[Nx-2][i]+StreamFunc->u[Nx-3][i]   )/StreamFunc->hx/StreamFunc->hx
                -(StreamFunc->u[xSize-1][i+1]   -2.*StreamFunc->u[xSize-1][i]+StreamFunc->u[xSize-1][i-1] )/StreamFunc->hy/StreamFunc->hy;
    } //условия на границе "внешнего" течения (на стоки и источнике жидкости)


    for(i = 1; i < xSize-1; i++)
    {
        for( j = 1; j< ySize-1; j++)
        {
            if(StreamFunc->NodeState[i][j]==1)
            {
                StreamFunc->u[i][j] = 0.5;
                CurlFunc->u[i][j]    =  -(StreamFunc->u[i-1][j] -2.*StreamFunc->u[i][j] +StreamFunc->u[i+1][j])/StreamFunc->hx/StreamFunc->hx
                                        -(StreamFunc->u[i][j-1] -2.*StreamFunc->u[i][j] +StreamFunc->u[i][j+1])/StreamFunc->hy/StreamFunc->hy;
            }//условия на твёрдые тела внутри потока


            //консервативная форма переносной силы
            CurlFunc->Force[i][j] =
                    -((CurlFunc->u[i+1][j]*vx[i+1][j]-CurlFunc->u[i-1][j]*vx[i-1][j])/CurlFunc->hx +
                      (CurlFunc->u[i][j+1]*vy[i][j+1]-CurlFunc->u[i][j-1]*vy[i][j-1])/CurlFunc->hy)*0.5;
        }
    }

    if(CurlFunc->time == 0.)
        CurlFunc->StepInTime_Euler();
    else
        CurlFunc->StepInTime_Adams();


    for(i=1;i<xSize-1;i++)
        for(j=1;j<ySize-1;j++)
        {
            StreamFunc->Force[i][j] = -CurlFunc->u[i][j];
        }
    StreamFunc->IterateWAutostop(30,1e-10);


    //EHD actuator's part
//    for(i=1; i<xSize-1; i++)
//    {
//        for(j=1; j<ySize-1; j++)
//        {
//            if(StreamFunc->NodeState[i][j-1]==1 && StreamFunc->NodeState[i][j]==0)
//            {
//                StreamFunc->u[i][j] = -hy/4.5+StreamFunc->u[i][j-1];
////                StreamFunc->NodeState[i][j] = 2;

//            }
//        }
//    }


}

void CFDProblem::ParaViewOutput(const char *filename)
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
                                               //sqrt(vx[i][j]*vx[i][j]+vy[i][j]*vy[i][j]));
    fclose(velocityValue);
}
