#ifndef POISSONTASK_H
#define POISSONTASK_H
#include <iostream>
#include <math.h>
#include <gnuplointerface.h>


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
        for(i=0;i<xSize;i++)
            for(j=0;j<ySize;j++)
                Force[i][j] = 0.;
    }

    double **Force;
    void Iterate(int n)
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
                                  - Force[i][j]*hx*hx*hy*hy)/2/(hx*hx+hy*hy);
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
};


#endif // POISSONTASK_H
