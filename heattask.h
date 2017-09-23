#ifndef HEATTASK_H
#define HEATTASK_H
#include <unistd.h>
#include <math.h>
#include <gnuplointerface.h>
#include <thread>
#include <vector>

class HeatTask
{
public:
    HeatTask(double LX, double RX,
             double LY, double RY,
             int XSize, int YSize,
             double dt, double TCC);

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
    void Animate(double FinalTime);
    void ClosePlot();

    double hx, hy;
    double **u;
//private:
    double dt;
    double **Prev_du_dt;
    double **Curr_du_dt;

    GnuploInterface Graph;
    int PreviousPlotKind;
};

class HeatTaskWDiscreteForce : public HeatTask
{
public:
    double **Force;
    double **next_u;
    int ThreadNum;
    HeatTaskWDiscreteForce(double LX, double RX,
                           double LY, double RY,
                           int XSize, int YSize,
                           double dt, double TCC);

    void StepInTime_Euler();

    void StepInTime_Adams();

    void CrankIterator(int n);

    void CrankIterator_Thread(int ThreadNum, int ThreadID);

    static void CrankIterator_Crutch(HeatTaskWDiscreteForce *Task, int ThreadNum, int ThreadID);

    void StepInTime_Crank();
};



#endif // HEATTASK_H
