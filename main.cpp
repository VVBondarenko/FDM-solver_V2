#include <cfdproblem.h>

using namespace std;

/*
 * ToDo:
 *
 *  add normal error control
 *  add selectivity to multy-resolution
 *
 */

class CFDairfoil : public CFDProblem
{

public:
    CFDairfoil(double LX, double RX,
               double LY, double RY,
               int XSize, int YSize,
               double dt, double Nu) : CFDProblem(LX,RX,LY,RY,XSize,YSize,dt,Nu)
    {

    }

    virtual double profileRfunc(double x, double y, double theta)
    {
        double X = cos(theta)*x-sin(theta)*y+0.54,
               Y = sin(theta)*x+cos(theta)*y;

        if(X<0)
            return -1.;
        double temp = 0.2969*sqrt(X)+ X*(-0.126 + X*(-0.3516 + (0.2843 - 0.1015*X)*X));
        return 1.5*temp - sqrt(pow(0.75*temp-Y,2) + pow(0.75*temp + Y,2));
    }
};

int main()
{
    CFDProblem *Test = new CFDairfoil(-1.,3.,
                                      -1.,1.,
                                      256, 128,
                                      1e-4,3./520.);
    Test->SetInitialConditions(25.);
    int k;
    char name[50];
    for(k=0; k<60000; k++)
    {
        Test->StepInTime();
        if(k%500==0)
        {
            sprintf(name, "result_%6.6d.dat", k);
            Test->ParaViewOutput(name);
        }
    }
    return 0;
}
