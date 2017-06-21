#include "poissontask.h"

PoissonTask::PoissonTask(double LX, double RX, double LY, double RY, int XSize, int YSize)
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

    int i,j;

    u = new double* [xSize];
    for(i = 0; i < xSize; i++)
            u[i] = new double[ySize];

    tmp_u = new double* [xSize];
    for(i = 0; i < xSize; i++)
        tmp_u[i] = new double[ySize];

    u_err = new double* [xSize];
    for(i = 0; i < xSize; i++)
        u_err[i] = new double[ySize];

    for(i=0;i<xSize;i++)
    {
        for(j=0;j<ySize;j++)
        {
            u[i][j] = 0.;
            tmp_u[i][j] = 0.;
            u_err[i][j] = 0.;
        }
    }

    for(i=0;i<xSize;i++)
    {
        u[i][0] = boundary_u(lx+i*hx,ly);
        tmp_u[i][0] = u[i][0];

        u[i][ySize-1] = boundary_u(lx+i*hx,ry);
        tmp_u[i][ySize-1] = u[i][ySize-1];
    }

    for(i=0;i<ySize;i++)
    {
        u[0][i] = boundary_u(lx,ly+i*hy);
        tmp_u[0][i] = u[0][i];

        u[i][ySize-1] = boundary_u(lx,ly+i*hy);
        tmp_u[i][ySize-1] = u[i][ySize-1];
    }

}

void PoissonTask::Iterate(int n)
{
    int i, j, K;
    for(K=0; K<n; K++)
    {

        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                tmp_u[i][j] = ((u[i+1][j]+u[i-1][j])*hy*hy
                              +(u[i][j+1]+u[i][j-1])*hx*hx
                              - rightpart_f(lx+i*hx,ly+j*hy)*hx*hx*hy*hy)/2/(hx*hx+hy*hy);
            }
        }

        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                u[i][j] = tmp_u[i][j];
            }
        }
    }

}

double PoissonTask::IterateWAutostop(int maxIters, double stop_criteria)
{
    int i;
    double current_err = 1.;
    for(i=0;i<maxIters && current_err>stop_criteria && i>=0;i++)
    {
            this->Iterate(30);
            current_err = this->EstimateConvolution();
            printf("Iters\n");
    }
    return current_err;
}

void PoissonTask::Output()
{
    int i, j;

    FILE *output;
    output = fopen("result.dat","w");

    for(i=0;i<xSize;i++)
    {
        for(j=0;j<ySize;j++)
        {
            fprintf(output,"%f %f %f\n",lx+hx*i,ly+hy*j,u[i][j]);
        }
        fprintf(output,"\n");
    }

    fclose(output);
}

void PoissonTask::Plot()
{
    Output();
    if(!Graph.isOpened())
    {
        std::vector<std::string> script;
        script.push_back("set terminal wxt");
        script.push_back("splot \"result.dat\" w surface ");

        Graph.open();
        Graph.execute(script);
    }
    else
    {
        Graph.write("replot");
        Graph.flush();
    }
}

void PoissonTask::ClosePlot()
{
    if(Graph.isOpened())
        Graph.close();
}

double PoissonTask::ExactError()
{
    int i, j;
    double maxErr = 0.;

    for(i=0;i<xSize;i++)
    {
        for(j=0;j<ySize;j++)
        {
//            fprintf(output,"%f %f %f\n",lx+hx*i,ly+hy*j,u[i][j]);
            u_err[i][j] = fabs(u[i][j]-exact_u(lx+hx*i,ly+hy*j));
            if(u_err[i][j]>maxErr)
                maxErr = u_err[i][j];
        }
    }
    return maxErr;
}

double PoissonTask::EstimateConvolution()
{
    int i,j;
    double maxErr;
    for(i = 1; i < xSize-1; i++)
    {
        for(j = 1; j < ySize-1; j++)
        {
            tmp_u[i][j] = ((u[i+1][j]+u[i-1][j])*hy*hy
                          +(u[i][j+1]+u[i][j-1])*hx*hx
                          - rightpart_f(lx+i*hx,ly+j*hy)*hx*hx*hy*hy)/2/(hx*hx+hy*hy);
        }
    }

    for(i = 1; i < xSize-1; i++)
    {
        for(j = 1; j < ySize-1; j++)
        {
            maxErr = fmax(maxErr, fabs(u[i][j]-tmp_u[i][j]));
        }
    }

    for(i = 1; i < xSize-1; i++)
    {
        for(j = 1; j < ySize-1; j++)
        {
            u[i][j] = tmp_u[i][j];
        }
    }

    return maxErr;
}

void PoissonTask::DoubleGrid()
{
    int i,j;
    for(i=0;i<xSize;i++)
    {
        for(j=0;j<ySize;j++)
        {
            tmp_u[i][j] = u[i][j];
        }
    }

    int NxSize = 2*xSize-1;
    int NySize = 2*ySize-1;

    delete [] u;
    delete [] u_err;

    u = new double* [NxSize];
    for(i = 0; i < NxSize; i++)
            u[i] = new double[NySize];


    u_err = new double* [NxSize];
    for(i = 0; i < NxSize; i++)
        u_err[i] = new double[NySize];

    //перенос старых значений
    for(i=0;i<xSize;i++)
        for(j=0;j<ySize;j++)
            u[2*i][2*j] = tmp_u[i][j];

    //межузловая интерполяция
    for(i=1;i<NxSize;i+=2)
        for(j=0;j<ySize;j++)
            u[i][2*j] = (tmp_u[i/2][j]+tmp_u[i/2+1][j])/2;
    for(i=0;i<xSize;i++)
        for(j=1;j<NySize;j+=2)
            u[2*i][j] = (tmp_u[i][j/2]+tmp_u[i][j/2+1])/2;

    //инетрполяция центра ячейки
    //(имеет ли смысл заморачиваться и прописывать это значение из дискретизации уравнения Пуассона?)
    for(i=1;i<NxSize;i+=2)
        for(j=1;j<NySize;j+=2)
            u[i][j] = (u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1])/4;

    xSize = NxSize;
    ySize = NySize;

    hx = (rx-lx)/(xSize-1);
    hy = (ry-ly)/(ySize-1);

    delete [] tmp_u;

    tmp_u = new double* [xSize];
    for(i = 0; i < xSize; i++)
        tmp_u[i] = new double[ySize];

    for(i=0;i<xSize;i++)
        for(j=0;j<ySize;j++)
            tmp_u[i][j] = u[i][j];

}
