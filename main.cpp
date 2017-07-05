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
 *  fix bug in animation
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
    test2(double LX, double RX,
          double LY, double RY,
          int XSize, int YSize,
          double dt, double TCC) : HeatTask(LX, RX, LY, RY, XSize, YSize, dt, TCC)
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

int heatEqWAnimation() //animated heat equation solution...
{
    int graphics_mode = 0;
    HeatTask *test = new test2(0.,1.,
                               0.,1.,
                               32,32,
                               0.0001, 0.2);
    test->StepInTime_Euler();
    test->Plot(graphics_mode);
    system("sleep 0.2");

    int i, itMax = 300;
    for(i = 0; i < itMax; i++)
    {
        test->StepInTime_Adams();
        test->Plot(graphics_mode);
        system("sleep 0.1");
    }

    test->ClosePlot();
//    test->Animate(0.1);

    return 0;
}



double f(double t, double p)
{
    if(t<=0.)
        return 0.;
    if(t>0 && t<p)
        return t*t/p/p*(3.-2.*t/p);
    else
        return 1.;
}

double streamBoundary(double x, double y)
{
    return f(x-2./3.*y, 1./3.);
}

int main()
{
    //объявление задачи для Пуассона и Теплопроводности
    //  Пуассона - начальное приближение для задачи обтекания невязкой несжимаемой
    //  Теплопроводность - уравнение переноса вихря
    int i,j;

    int Nx = 64, Ny = 64;
    PoissonTaskWDiscreteForce *StreamFunc = new PoissonTaskWDiscreteForce(-1.,1.,
                                                                          -1.,1.,
                                                                          Nx,Ny);
    HeatTaskWDiscreteForce *CurlFunc = new HeatTaskWDiscreteForce(-1.,1.,
                                                                  -1.,1.,
                                                                  Nx,Ny,
                                                                  0.00001,1./5.);
    double **vx, **vy;
    vx = new double* [Nx];
    vy = new double* [Nx];
    for(i = 0; i < Nx; i++)
    {
        vx[i] = new double[Ny];
        vy[i] = new double[Ny];
    }
    for(i=0;i<Nx;i++)
        for(j=0;j<Ny;j++)
        {
            vx[i][j] = 0.;
            vy[i][j] = 0.;
        }


    //заполнение ГУ для функции тока
    for(i=0; i<Nx; i++)
    {
        StreamFunc->u[i][0]     =   streamBoundary((-1.+StreamFunc->hx*i), -1.);
        StreamFunc->u[i][Ny-1]  =   streamBoundary((-1.+StreamFunc->hx*i),  1.);
    }
    for(i=0; i<Ny; i++)
    {
        StreamFunc->u[0][i]    =   streamBoundary(-1.,(-1.+StreamFunc->hx*i));
        StreamFunc->u[Nx-1][i] =   streamBoundary( 1.,(-1.+StreamFunc->hx*i));
    }

    //получаем начальное приближение для уравнения Пуассона
    StreamFunc->IterateWAutostop(100,1e-8);
    StreamFunc->Plot(1);
    system("sleep 10");

    for(i=1;i<Nx-1;i++)
        for(j=1;j<Ny-1;j++)
        {
            vx[i][j] =  (StreamFunc->u[i][j+1]-StreamFunc->u[i][j-1])/StreamFunc->hy*0.5;
            vy[i][j] = -(StreamFunc->u[i+1][j]-StreamFunc->u[i-1][j])/StreamFunc->hx*0.5;
        }

    //сетаем начальное распределение для уравнения переноса вихря (+ГУ) (+переносная сила)
    for(i=1; i<Nx-1; i++)
    {
        for(j=1; j<Ny-1; j++)
        {
            CurlFunc->u[i][j] = - ((StreamFunc->u[i+1][j]-2.*StreamFunc->u[i][j]+StreamFunc->u[i-1][j])/StreamFunc->hx/StreamFunc->hx+
                                   (StreamFunc->u[i][j+1]-2.*StreamFunc->u[i][j]+StreamFunc->u[i][j-1])/StreamFunc->hy/StreamFunc->hy);
        }
    }
    for(i=0; i<Nx; i++)
    {
        CurlFunc->u[i][0]    = CurlFunc->u[i][1];
        CurlFunc->u[i][Nx-1] = CurlFunc->u[i][Nx-2];
    }
    for(i=0; i<Ny; i++)
    {
        CurlFunc->u[0][i]    = CurlFunc->u[1][i];
        CurlFunc->u[Ny-1][i] = CurlFunc->u[Nx-2][i];
    }

    for(i=1;i<Nx-1;i++)
        for(j=1;j<Ny-1;j++)
        {
            CurlFunc->Force[i][j] = -(vx[i][j]*(CurlFunc->u[i-1][j]-CurlFunc->u[i+1][j])/CurlFunc->hx +
                                      vy[i][j]*(CurlFunc->u[i][j+1]-CurlFunc->u[i][j-1])/CurlFunc->hy);
        }

    //итерируем по времени
    CurlFunc->StepInTime_Euler();
    //обновляем решение уравнения Пуассона
    //  а вот тут возможны нюансы... решено
    for(i=1;i<Nx-1;i++)
        for(j=1;j<Ny-1;j++)
        {
            StreamFunc->Force[i][j] = -CurlFunc->u[i][j];
        }
    StreamFunc->IterateWAutostop(100,1e-8);
    StreamFunc->Plot(1);
    system("sleep 10");

    return 0;
}
