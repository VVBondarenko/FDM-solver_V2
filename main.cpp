//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

#include <poissontask.h>
#include <heattask.h>
//#include <gnuplointerface.h>


using namespace std;

/*
 * ToDo:
 *
 *  add normal error control
 *  add selectivity to multy-resolution
 *  add inner boundaries with conditions
 *
 *
 */

class test1 : public PoissonTask
{

public:
    test1(double LX, double RX, double LY, double RY,
          int XSize, int YSize) : PoissonTask(LX,RX,LY,RY,XSize,YSize)
    {    }

    virtual double rightpart_f(double x, double y)
    {
        return -2.*sin(x)*sin(y);
    }
    virtual double boundary_u(double x, double y)
    {
        return 0.*x*y;
    }
    virtual double exact_u(double x, double y)
    {
        return sin(x)*sin(y);
    }
};

int PoissonTester()
{
    PoissonTask* test = new test1(0,2.*M_PI,0,2.*M_PI,16,16);
    double err;
    err = test->IterateWAutostop(30,1e-9);
    printf("%e\t%e\t%e\n", err, test->EstimateConvolution(), test->ExactError());

    test->Plot(1);
    scanf("%f",&err);

    test->Plot(0);
    scanf("%f",&err);

    test->ClosePlot();

    return 0;
}

class test2 : public HeatTask
{

public:
    test2(double LX, double RX, double LY, double RY,
          int XSize, int YSize, double dt) : HeatTask(LX, RX, LY, RY, XSize, YSize, dt)
    {   }

    virtual double rightpart_f(double x, double y, double t)
    {
        return 0.*x*y*t;
    }

    virtual double boundary_u(double x, double y, double t)
    {
        if(x==0. || y==0.)
            return 1.;

        return 0.*t;
    }

    virtual double initial_f(double x, double y)
    {
        if(x==0. || y==0.)
            return 1.;

        return 0.;
    }

};

int main()
{
    int i, itMax = 2, temp;
    HeatTask *test = new test2(0.,1.,0.,1.,32,32,0.0001);
    test->StepInTime_Euler();
    test->Plot(0);

    for(i = 0; i < itMax; i++)
    {
        if(i == itMax-1)
        {
            test->Plot(0);
            scanf("%d", &temp);
            itMax += temp;
        }
        else
        {
            test->StepInTime_Adams();
        }
    }

    test->ClosePlot();
    return 0;
}
