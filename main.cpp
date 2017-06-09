//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

#include <poissontask.h>

using namespace std;

/*
 * ToDo:
 *
 *  add automatical iteration contol
 *  add selectivity to multy-resolution
 *
 *
 */

class test1 : public PoissonTask {
public:
    test1(double LX, double RX, double LY, double RY, int XSize, int YSize) : PoissonTask(LX,RX,LY,RY,XSize,YSize)
    {    }

    virtual double rightpart_f(double x, double y)
    {
        return -2.*sin(x)*sin(y);
    }
    virtual double boundary_u(double x, double y)
    {
        return 0.;
    }
    virtual double exact_u(double x, double y)
    {
        return sin(x)*sin(y);
    }
};


int main()
{
    PoissonTask* test = new test1(0,M_PI,0,M_PI,16,16);
    test->Iterate(100);
    test->Output();
    test->DoubleGrid();
    test->Iterate(100);
    test->Output();
    test->DoubleGrid();
    test->Iterate(1000);
    test->Output();
    return 0;
}
