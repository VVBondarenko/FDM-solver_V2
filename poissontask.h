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
//private:
    double **u_err;

    GnuploInterface Graph;
    int PreviousPlotKind;
};

//extern "C" void CMemCpy(void *dest, void *src)
//{
//    memcpy(dest, src, sizeof(src));
//}


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
    }

    double **Force;
    int **NodeState;

    static void IteratorCrutch(PoissonTaskWDiscreteForce* Temp, int ID)
    {
        Temp->IteratorThread(ID);
    }

    void IteratorThread(int ThreadID)
    {
        int i, j;
        int H = (xSize-2)/4;
        int Istart  = 1+ThreadID*H;
        int Iend    =   Istart + H;
        if(ThreadID==3)
            Iend = xSize-1;

        for(i = Istart; i < Iend; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                if(NodeState[i][j]==0)
                    u[i][j] = ((u[i+1][j]+u[i-1][j])*hy*hy
                              +(u[i][j+1]+u[i][j-1])*hx*hx
                              - Force[i][j]*hx*hx*hy*hy)/2./(hx*hx+hy*hy);
            }
        }
    }

    void Iterate(int n, double SORomega)
    {
//        double SORomega = 1.1;

        int K;
        for(K=0; K<n; K++)
        {
//            double **U, **Tmp_U;
//            U       = this->u;
//            Tmp_U   = this->tmp_u;

//#pragma omp parallel for collapse(2) //shared(U,Tmp_U) private(i,j)
/*
            for(i = 1; i < xSize-1; i++)
            {
                for(j = 1; j < ySize-1; j++)
                {
                    if(NodeState[i][j]==0)
//                        Tmp_U[i][j] = ((U[i+1][j]+U[i-1][j])*hy*hy
//                        U[i][j] = (1-SORomega)*U[i][j] +SORomega*((U[i+1][j]+U[i-1][j])*hy*hy
//                                                                 +(U[i][j+1]+U[i][j-1])*hx*hx
//                                                                 -Force[i][j]*hx*hx*hy*hy)/2./(hx*hx+hy*hy);
                    u[i][j] = ((u[i+1][j]+u[i-1][j])*hy*hy
                              +(u[i][j+1]+u[i][j-1])*hx*hx
                              - Force[i][j]*hx*hx*hy*hy)/2./(hx*hx+hy*hy);
                }
            }
*/
            std::thread Sol0 (IteratorCrutch,this,0);
            std::thread Sol1 (IteratorCrutch,this,1);
            std::thread Sol2 (IteratorCrutch,this,2);
            std::thread Sol3 (IteratorCrutch,this,3);
//            std::thread Sol4 (IteratorCrutch,this,4);
//            std::thread Sol5 (IteratorCrutch,this,5);
//            std::thread Sol6 (IteratorCrutch,this,6);
//            std::thread Sol7 (IteratorCrutch,this,7);

            Sol0.join();
            Sol1.join();
            Sol2.join();
            Sol3.join();
//            Sol4.join();
//            Sol5.join();
//            Sol6.join();
//            Sol7.join();


//#pragma omp parallel for collapse(2)
//            for(i = 1; i < xSize-1; i++)
//            {
//#pragma GCC ivdep
//                for(j = 1; j < ySize-1; j++)
//                {
//                    if(NodeState[i][j]==0)
//                        U[i][j] = Tmp_U[i][j];
//                }
//            }
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
//#pragma omp parallel for collapse(2) //shared(U,Tmp_U) private(i,j)

        int H = (xSize-2)/4;
        int Istart  = 1+ThreadID*H;
        int Iend    =   Istart + H;
        if(ThreadID==3)
            Iend = xSize-1;

        for(i = Istart; i < Iend; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                tmp_u[i][j] = ((u[i+1][j]+u[i-1][j])*hy*hy
                              +(u[i][j+1]+u[i][j-1])*hx*hx
                              - Force[i][j]*hx*hx*hy*hy)/2/(hx*hx+hy*hy);
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
/*
        int i,j;
        double maxErr = 0.;
//#pragma omp parallel for collapse(2) //shared(U,Tmp_U) private(i,j)
        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                tmp_u[i][j] = ((u[i+1][j]+u[i-1][j])*hy*hy
                              +(u[i][j+1]+u[i][j-1])*hx*hx
                              - Force[i][j]*hx*hx*hy*hy)/2/(hx*hx+hy*hy);
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
*/
        double err0, err1, err2, err3;
//        std::ref
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
        for(i=0;i<maxIters && current_err>stop_criteria && i>=0;i++)
        {
            this->Iterate(10,1.0);
            prev_err = current_err;
            current_err = this->EstimateConvolution();
            if(fabs(prev_err-current_err)/prev_err<1.+1e-6)
                break;
    //            printf("Iters\n");
        }
        return current_err;
    }
};



#endif // POISSONTASK_H
