#ifndef HEATTASK_H
#define HEATTASK_H
#include <math.h>
#include <gnuplointerface.h>

class HeatTask
{
public:
    HeatTask(double LX, double RX, double LY, double RY, int XSize, int YSize, double dt);

    double lx, rx;
    double ly, ry;
    int xSize, ySize;

    double time = 0.;
    double TCC = 1.; // ThermoConductivityCoefficient
    virtual double rightpart_f(double x, double y, double t)
    {
        return 0.*x*y*t;
    }

    virtual double boundary_u(double x, double y, double t)
    {
        return 0.*x*y*t;
    }

    virtual double initial_f(double x, double y)
    {
        return 0.*x*y;
    }

    void StepInTime_Euler();
    void StepInTime_Adams();
    void DoubleGrid(); //is it necessary at that moment?

    void Output();
    void Plot(int Kind);
    void ClosePlot();

private:
    double hx, hy;
    double dt;
    double **u;
    double **Prev_du_dt;
    double **Curr_du_dt;

    GnuploInterface Graph;
    int PreviousPlotKind;
};

#endif // HEATTASK_H
