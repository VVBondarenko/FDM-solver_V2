#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;


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

int main(int argc, char *argv[])
{
    double l_x = 0., r_x = M_PI;
    double l_y = 0., r_y = M_PI;
    int xSize = 16, ySize = 16;


    double  hx = (r_x-l_x)/(xSize-1),
            hy = (r_y-l_y)/(ySize-1);
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
        u[i][0]             = u_boundary(l_x+i*hx, l_y);
        u[i][ySize-1]       = u_boundary(l_x+i*hx, r_y);
        tmp_u[i][0]         = u[i][0];
        tmp_u[i][ySize-1]   = u[i][ySize-1];
    }

    for(i = 0; i < ySize; i++)
    {
        u[0][i]             = u_boundary(l_x, l_y+i*hy);
        u[xSize-1][i]       = u_boundary(r_x, l_y+i*hy);
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
                              - f(l_x+i*hx,l_y+j*hy)*hx*hx*hy*hy)/2/(hx*hx+hy*hy);
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
            fprintf(output,"%f %f %f\n",l_x+hx*i,l_y+hy*j,u[i][j]);
        }
        fprintf(output,"\n");
    }

    return 0;
}
