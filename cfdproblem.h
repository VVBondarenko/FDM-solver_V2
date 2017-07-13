#ifndef CFDPROBLEM_H
#define CFDPROBLEM_H

#include <poissontask.h>
#include <heattask.h>


class CFDProblem
{
public:
    CFDProblem(double LX, double RX,
               double LY, double RY,
               int XSize, int YSize,
               double dtime, double Viscousity);

    double lx, rx;
    double ly, ry;
    int xSize, ySize;
    double dt, nu;

    double hx, hy;
    double **vx, **vy;

    PoissonTaskWDiscreteForce *StreamFunc;
    HeatTaskWDiscreteForce *CurlFunc;

    void SetInitialConditions(double AngleOfAttack);

    double StreamBoundaryFunc(double x, double y);
    virtual double profileRfunc(double x, double y, double theta)
    {
//        return 0.04-((x+0.5)*(x+0.5)+y*y);
        return 0.;
    }

    void UpdateBoundaryCond();
    void StepInTime();

    void ParaViewOutput(const char *filename);
};

#endif // CFDPROBLEM_H
