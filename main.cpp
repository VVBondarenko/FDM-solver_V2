//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

#include <poissontask.h>

using namespace std;

/*
 * ToDo:
 *  form library from this code:
 *      devide to methods
 *      devide to classes
 *
 *  add automatical iteration contol
 *  add multy-resolution grids
 *  add selectivity to multy-resolution
 *
 *
 *
 */


double f(double x, double y)
{
    return -2.*sin(x)*sin(y);
}

double u_boundary(double x, double y)
{
    return 0.;
}

double u_exact(double x, double y)
{
    return sin(x)*sin(y);
}

int main_old(int argc, char *argv[])
{
    double lx = 0., rx = M_PI;
    double ly = 0., ry = M_PI;
    int xSize = 16, ySize = 16;


    double  hx = (rx-lx)/(xSize-1),
            hy = (ry-ly)/(ySize-1);
    double     u[xSize][ySize];
    double tmp_u[xSize][ySize];

    //init boundary conditions
    int i, j, K;

    for(i=0;i<xSize;i++)
        for(j=0;j<ySize;j++)
        {
            u[i][j] = 0.;
            tmp_u[i][j] = 0.;
        }

    for(i = 0; i < xSize; i++)
    {
        u[i][0]             = u_boundary(lx+i*hx, ly);
        u[i][ySize-1]       = u_boundary(lx+i*hx, ry);
        tmp_u[i][0]         = u[i][0];
        tmp_u[i][ySize-1]   = u[i][ySize-1];
    }

    for(i = 0; i < ySize; i++)
    {
        u[0][i]             = u_boundary(lx, ly+i*hy);
        u[xSize-1][i]       = u_boundary(rx, ly+i*hy);
        tmp_u[0][i]         = u[0][i];
        tmp_u[xSize-1][i]   = u[xSize-1][i];
    }

    //jacoby iterations
    for(K=0; K<100; K++)
    {

        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                tmp_u[i][j] = ((u[i+1][j]+u[i-1][j])*hy*hy
                              +(u[i][j+1]+u[i][j-1])*hx*hx
                              - f(lx+i*hx,ly+j*hy)*hx*hx*hy*hy)/2/(hx*hx+hy*hy);
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

    FILE *output;
    output = fopen("result.dat","w");

    for(i=0;i<xSize;i++)
    {
        for(j=0;j<ySize;j++)
        {
            fprintf(output,"%f %f %f\n",lx+hx*i,ly+hy*j,u[i][j]);
        }
        fprintf(output,"\n");
    }

    return 0;
}

class test1 : public PoissonTask {
public:
    test1(double LX, double RX, double LY, double RY, int XSize, int YSize) : PoissonTask(LX,RX,LY,RY,XSize,YSize)
    {

    }

    virtual double rightpart_f(double x, double y)
    {
        return f(x,y);
    }
    virtual double boundary_u(double x, double y)
    {
        return 0.;
    }
    virtual double exact_u(double x, double y)
    {
        return u_exact(x,y);
    }
};


int main()
{
    PoissonTask* test = new test1(0,M_PI,0,M_PI,16,16);
    test->Iterate(100);
    test->Output();
    return 0;
}
