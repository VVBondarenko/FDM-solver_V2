#ifndef POISSONTASK_H
#define POISSONTASK_H
#include <iostream>
extern "C"
{
#include <math.h>
#include <string.h>
}
#include <gnuplointerface.h>
#include <omp.h>
#include <thread>
#include <vector>

class PoissonTask
{
public:
    PoissonTask(double LX, double RX, double LY, double RY, int XSize, int YSize);

    double lx, rx;
    double ly, ry;
    int xSize, ySize;

    virtual double rightpart_f(double x, double y)
    {
        return 0.*x*y;
    }

    virtual double boundary_u(double x, double y)
    {
        return 0.*x*y;
    }

    virtual double exact_u(double x, double y)
    {
        return 0.*x*y;
    }

    void Iterate(int n);
    void DoubleGrid();

    double ExactError();
    double EstimateConvolution();

    double IterateWAutostop(int maxIters, double stop_criteria);

    void Output();
    void Plot(int Kind);
    void ClosePlot();

    double hx, hy;
    double **u;
    double **tmp_u;
    double **u_err;

    GnuploInterface Graph;
    int PreviousPlotKind;
};


class PoissonTaskWDiscreteForce : public PoissonTask
{
public:
    PoissonTaskWDiscreteForce(double LX, double RX,
                              double LY, double RY,
                              int XSize, int YSize) : PoissonTask(LX,   RX,
                                                                  LY,   RY,
                                                                  XSize,YSize)
    {
        int i, j;

        Force = new double* [xSize];
        for(i = 0; i < xSize; i++)
            Force[i] = new double[ySize];

        NodeState = new int* [xSize];
        for(i=0;i<xSize;i++)
            NodeState[i] = new int[ySize];

        for(i=0;i<xSize;i++)
            for(j=0;j<ySize;j++)
            {
                Force[i][j] = 0.;
                NodeState[i][j]=0.;
            }
        ThreadNum = 4;
    }

    double **Force;
    int **NodeState;
    int ThreadNum;
    static void IteratorCrutch(PoissonTaskWDiscreteForce* Temp, int ID, int Q)
    {
        Temp->IteratorThread(ID,Q);
    }

    void IteratorThread(int ThreadID, int ThreadNum)
    {
        int i, j;
        int H = (xSize-2)/ThreadNum;
        int Istart  = 1+ThreadID*H;
        int Iend    =   Istart + H;

        double hx_s = hx*hx,
               hy_s = hy*hy;
        double u_temp;
//        int Nactive = 0;
        for(i = Istart; i < Iend; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                if(NodeState[i][j]==0)
                {
                    u_temp= ((u[i+1][j]+u[i-1][j]-Force[i][j]*hx_s)*hy_s
                            +(u[i][j+1]+u[i][j-1])*hx_s
                            /*- Force[i][j]*hx_s*hy_s*/)*0.5/(hx_s+hy_s);
                    if(fabs(u_temp-u[i][j])<1.e-7)
                        NodeState[i][j] = -1;
                    u[i][j] = u_temp;
//                    Nactive++;
                }
            }
        }
//        printf("%d\n",Nactive);
    }

    void Iterate(int n, int ResetActivesNum)
    {

        int i,j,K;
//        int ThreadNum = 6;
        for(K=0; K<n; K++)
        {
            std::vector<std::thread> ThrPool;

            for(i=0;i<ThreadNum;i++)
            {
                ThrPool.push_back(std::thread(IteratorCrutch,this,i,ThreadNum));
            }

            for(i=0;i<ThreadNum;i++)
            {
                ThrPool[i].join();
            }

            if(K%ResetActivesNum==ResetActivesNum-1 || K==n-1)
                for(i=1;i<xSize-1;i++)
                    for(j=1;j<ySize-1;j++)
                        if(NodeState[i][j] == -1)
                            NodeState[i][j] = 0;
        }

    }

    void IterateOnMultyRes(int n, int factor)
    {
        int K;
        for(K=0; K<n; K++)
        {
            int i, j;
//            int H = (xSize-2)/4;
//            int Istart  = 1+ThreadID*H;
//            int Iend    =   Istart + H;
//            if(ThreadID==3)
//                Iend = xSize-1;

            double hx_s = hx*hx,
                   hy_s = hy*hy;

            for(i = factor; i < xSize-factor; i+=factor)
            {
                for(j = factor; j < ySize-factor; j+=factor)
                {
                    if(NodeState[i][j]==0)
                        u[i][j] = ((u[i+factor][j]+u[i-factor][j]-Force[i][j]*hx_s)*hy_s
                                  +(u[i][j+factor]+u[i][j-factor])*hx_s
                                  /*- Force[i][j]*hx_s*hy_s*/)*0.5/(hx_s+hy_s);
                }
            }
        }
    }

    static void EstimatorCrutch(PoissonTaskWDiscreteForce* Temp, int ID, double &maxErr)
    {
        Temp->EstimatorThread(ID, maxErr);
    }

    void EstimatorThread(int ThreadID, double &maxErr)
    {
        int i,j;
        maxErr = 0.;

        int H = (xSize-2)/4;
        int Istart  = 1+ThreadID*H;
        int Iend    =   Istart + H;
        if(ThreadID==3)
            Iend = xSize-1;

        double hx_s = hx*hx,
               hy_s = hy*hy;

        for(i = Istart; i < Iend; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {

                tmp_u[i][j] = ((u[i+1][j]+u[i-1][j]-Force[i][j]*hx_s)*hy_s
                              +(u[i][j+1]+u[i][j-1])*hx_s
                              /*- Force[i][j]*hx_s*hy_s*/)*0.5/(hx_s+hy_s);

            }
        }

        for(i = Istart; i < Iend; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                maxErr = fmax(maxErr, fabs(u[i][j]-tmp_u[i][j]));
            }
        }

        for(i = Istart; i < Iend; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                u[i][j] = tmp_u[i][j];
            }
        }

    }

    double EstimateConvolution()
    {
        double err0, err1, err2, err3;

        std::thread Sol0 (EstimatorCrutch,this,0,std::ref(err0));
        std::thread Sol1 (EstimatorCrutch,this,1,std::ref(err1));
        std::thread Sol2 (EstimatorCrutch,this,2,std::ref(err2));
        std::thread Sol3 (EstimatorCrutch,this,3,std::ref(err3));
//            std::thread Sol4 (EstimatorCrutch,this,4,std::ref(err4));
//            std::thread Sol5 (EstimatorCrutch,this,5,std::ref(err5));
//            std::thread Sol6 (EstimatorCrutch,this,6,std::ref(err6));
//            std::thread Sol7 (EstimatorCrutch,this,7,std::ref(err7));

        Sol0.join();
        Sol1.join();
        Sol2.join();
        Sol3.join();
//            Sol4.join();
//            Sol5.join();
//            Sol6.join();
//            Sol7.join();

        double maxErr = fmax(err0,fmax(err1,fmax(err2,err3)));

        return maxErr;
    }

    double IterateWAutostop(int maxIters, double stop_criteria)
    {
        int i;
        double current_err = 1., prev_err = 2.;
        for(i=0;i<maxIters && current_err>stop_criteria/* && i>=0*/;i++)
        {
            this->Iterate(2,1);
            prev_err = current_err;
//            current_err = this->EstimateConvolution();
            printf("%e\t%e\t%d\n",prev_err,current_err = this->EstimateConvolution(),i);
            if(fabs(prev_err-current_err)/prev_err<1.+1e-7)
            {
                printf("brake\n");
                break;
            }
        }
        return current_err;
    }
};



#endif // POISSONTASK_H
