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



    void Iterate(int n)
    {
        int i, j, K;
        for(K=0; K<n; K++)
        {
            double **U/*, **Tmp_U*/;
            U       = this->u;
//            Tmp_U   = this->tmp_u;

#pragma omp parallel for collapse(2) //shared(U,Tmp_U) private(i,j)
            for(i = 1; i < xSize-1; i++)
            {
//#pragma GCC ivdep
                for(j = 1; j < ySize-1; j++)
                {
                    if(NodeState[i][j]==0)
//                        Tmp_U[i][j] = ((U[i+1][j]+U[i-1][j])*hy*hy
                        U[i][j] = ((U[i+1][j]+U[i-1][j])*hy*hy
                                  +(U[i][j+1]+U[i][j-1])*hx*hx
                                  - Force[i][j]*hx*hx*hy*hy)/2./(hx*hx+hy*hy);
                }
            }

//#pragma omp parallel for collapse(2)
//#pragma GCC ivdep
//            for(i = 1; i < xSize-1; i++)
//            {
//                for(j = 1; j < ySize-1; j++)
//                {
//                    if(NodeState[i][j]==0)
//                        U[i][j] = Tmp_U[i][j];
//                }
//            }
        }
    }

    double IterateWAutostop(int maxIters, double stop_criteria)
    {
        int i;
        double current_err = 1., prev_err = 2.;
        for(i=0;i<maxIters && current_err>stop_criteria && i>=0;i++)
        {
            this->Iterate(30);
            prev_err = current_err;
            current_err = this->EstimateConvolution();
            if(fabs(prev_err-current_err)/prev_err<1.+1e-9)
                break;
    //            printf("Iters\n");
        }
        return current_err;
    }
};


#endif // POISSONTASK_H
