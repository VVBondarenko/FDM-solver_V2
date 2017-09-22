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

void CallOutput(CFDProblem *Task, const char *name)
{
    Task->ParaViewOutput_v2(name);
}

int main()
{
    CFDProblem *Test = new CFDairfoil(-1.,3.,
                                      -1.,1.,
                                      256, 128,
                                      1e-3,3./520.);
    Test->SetInitialConditions(24.);
    Test->StreamFunc->Iterate(20,4);
    int k;
//    scanf("%d",&k);
    char name[50];

/*

//  part of possible future balancer mechanism

    int nodes = 0,p;
    for(k=1;k<Test->xSize-1;k++)
    {
        if(k%(Test->xSize/8)==0)
        {
            printf("%d\n",nodes);
            nodes = 0;
        }
        for(p=1;p<Test->ySize-1;p++)
        {

            if(Test->StreamFunc->NodeState[k][p] == 0)
                nodes++;
        }
    }
    printf("%d\n",nodes);
*/
    for(k=0; k<6000; k++)
    {
        Test->StepInTime();
        if(k%50==0)
        {
            sprintf(name, "result_%6.6d.nc", k);
            Test->ParaViewOutput_v2(name);
        }
    }
    return 0;
}
