#ifndef CFDPROBLEM_H
#define CFDPROBLEM_H

#include <poissontask.h>
#include <heattask.h>
#include <iostream>
#include <netcdf>
#include <vector>
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;


class CFDProblem
{
public:
    CFDProblem(double LX, double RX,
               double LY, double RY,
               int XSize, int YSize,
               double dtime, double Viscousity);

    int ThreadNum;

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
        return 0.*x*y*theta;
    }

    void UpdateBoundaryCond();
    void StepInTime();

    void UpdateConvectiveForce();
    void UpdateConvectiveForce_Thread(int ThreadNum, int ThreadID);
    static void UpdateConvectiveForce_Crutch(CFDProblem *Task, int ThreadNum, int ThreadID);


    void ParaViewOutput     (const char *filename);
    void ParaViewOutput_v2  (const char *filename);


    void SetThreadNum(int n);
};

#endif // CFDPROBLEM_H
