#ifndef POISSONTASK_H
#define POISSONTASK_H
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
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
        return 0.;
    }

    virtual double boundary_u(double x, double y)
    {
        return 0.;
    }

    virtual double exact_u(double x, double y)
    {
        return 0.;
    }

    void Iterate(int n);
    void DoubleGrid();

    double ExactError();
    double EstimateConvolution();

    double IterateWAutostop(int maxIters, double stop_criteria);

    void Output();
    void Plot();
    void ClosePlot();

private:
    double hx, hy;
    double **u;
    double **tmp_u;
    double **u_err;

    GnuploInterface Graph;
};

#endif // POISSONTASK_H
