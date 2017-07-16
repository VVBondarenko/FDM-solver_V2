//#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

//#include <poissontask.h>
//#include <heattask.h>
//#include <gnuplointerface.h>
#include <cfdproblem.h>

using namespace std;

/*
 * ToDo:
 *
 *  add normal error control
 *  add selectivity to multy-resolution
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
    int graphics_mode = 0, ret_code;
    HeatTask *test = new test2(0.,1.,
                               0.,1.,
                               32,32,
                               0.0001, 0.2);
    test->StepInTime_Euler();
    test->Plot(graphics_mode);
    ret_code = system("sleep 0.2");

    int i, itMax = 300;
    for(i = 0; i < itMax; i++)
    {
        test->StepInTime_Adams();
        test->Plot(graphics_mode);
        ret_code = system("sleep 0.1");
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
    //течение Пуазейля
    return f((y+1.)/2.,1.)+0.*x;
//    return y;
    //течение в каверне
//    return 0.*x*y;
    //течение Куэтта
//    return 8./3.*f((y+1.)/4.,1.);
//    return 8./3.*f((y+1.)/2.,1.);
}
double bodyRfunc_airfoil(double X, double Y, double theta)
{
    //theta - ange of attack
    double x = cos(theta)*X-sin(theta)*Y+0.333333,
           y = sin(theta)*X+cos(theta)*Y;

    if(x<0)
        return -1.;
    double temp = 0.2969*sqrt(x)+ x*(-0.126 + x*(-0.3516 + (0.2843 - 0.1015*x)*x));
    return 1.5*temp - sqrt(pow(0.75*temp-y,2) + pow(0.75*temp + y,2));
}

double bodyRfunc(double X, double Y, double theta)
{
    return 0.04-((X+0.5)*(X+0.5)+Y*Y)+0.*theta;
}

int main_CFD_test()
{
    //объявление задачи для Пуассона и Теплопроводности
    //  Пуассона - начальное приближение для задачи обтекания невязкой несжимаемой
    //  Теплопроводность - уравнение переноса вихря
    int i,j;

    int ret_code;




    int     Nx = 128,
            Ny = 128;
    PoissonTaskWDiscreteForce *StreamFunc = new PoissonTaskWDiscreteForce(-1.,1.,
                                                                          -1.,1.,
                                                                          Nx,Ny);
    HeatTaskWDiscreteForce *CurlFunc = new HeatTaskWDiscreteForce(-1.,1.,
                                                                  -1.,1.,
                                                                  Nx,Ny,
                                                                  0.0001,/*1./20.*/3./520.);
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
        StreamFunc->u[i][0]     =   streamBoundary((StreamFunc->lx+StreamFunc->hx*i), StreamFunc->ly);
        StreamFunc->u[i][Ny-1]  =   streamBoundary((StreamFunc->lx+StreamFunc->hx*i), StreamFunc->ry);

//        StreamFunc->NodeState[i][1] = 1;
//        StreamFunc->NodeState[i][Ny-2] = 1;

//        vx[i][Ny-1] = 1.;

    }
    for(i=0; i<Ny; i++)
    {
        StreamFunc->u[0][i]    =   streamBoundary(StreamFunc->lx,(StreamFunc->ly+StreamFunc->hy*i));
        StreamFunc->u[Nx-1][i] =   streamBoundary(StreamFunc->rx,(StreamFunc->ly+StreamFunc->hy*i));

//        StreamFunc->NodeState[1][i] = 1;
//        StreamFunc->NodeState[Nx-2][i] = 1;


//        if(i>2*Ny/5 && i<3*Ny/5)    //условия для пластинки в потоке
//        {
//            StreamFunc->NodeState[Nx/4][i] = 1;
//            StreamFunc->u[Nx/4][i] = 0.5;
//        }

        vx[0][i]    = (streamBoundary(StreamFunc->lx,
                                      StreamFunc->ly+i*StreamFunc->hy+1e-8)
                      -streamBoundary(StreamFunc->lx,
                                      StreamFunc->ly+i*StreamFunc->hy-1e-8))/2.e-8;
        vx[Nx-1][i] = (streamBoundary(StreamFunc->rx,
                                      StreamFunc->ly+i*StreamFunc->hy+1e-8)
                      -streamBoundary(StreamFunc->rx,
                                      StreamFunc->ly+i*StreamFunc->hy-1e-8))/2.e-8;
    }

    //получаем начальное приближение для уравнения Пуассона (а надо ли?..)
//    StreamFunc->IterateWAutostop(200,1e-8);
    for(i=1;i<Nx-1;i++)
        for(j=1;j<Ny-1;j++)
        {
            StreamFunc->u[i][j] = streamBoundary(StreamFunc->lx+i*StreamFunc->hx,
                                                 StreamFunc->ly+j*StreamFunc->hy);
            if(bodyRfunc((StreamFunc->lx+i*StreamFunc->hx),
                         (StreamFunc->ly+j*StreamFunc->hy     ),
                         M_PI/180.*30.) > 0.)
            {
                StreamFunc->NodeState[i][j] = 1;
                StreamFunc->u[i][j] = 0.5;
            }
        }
//    StreamFunc->Plot(2);
//    system("sleep 2");
//    StreamFunc->ClosePlot();

    for(i=1; i<Nx-1; i++)
        for(j=1; j<Ny-1; j++)
            CurlFunc->u[i][j] =   -((StreamFunc->u[i+1][j]-2.*StreamFunc->u[i][j]+StreamFunc->u[i-1][j])/StreamFunc->hx/StreamFunc->hx+
                                    (StreamFunc->u[i][j+1]-2.*StreamFunc->u[i][j]+StreamFunc->u[i][j-1])/StreamFunc->hy/StreamFunc->hy);
    for(i=0; i<Nx; i++)
    {
        CurlFunc->u[i][0]    = CurlFunc->u[i][1];
        CurlFunc->u[i][Ny-1] = CurlFunc->u[i][Ny-2];
    }

    for(i=0; i<Ny; i++)
    {
        CurlFunc->u[0][i]    = CurlFunc->u[1][i];
        CurlFunc->u[Nx-1][i] = CurlFunc->u[Nx-2][i];
    }

    int k;
    for(k=0;k<12000;k++)
    {
        //получаем поле скоростей по предыдущему виду функции тока


        for(i=1;i<Nx-1;i++)
            for(j=1;j<Ny-1;j++)
            {
                vx[i][j] =  (StreamFunc->u[i][j+1]-StreamFunc->u[i][j-1])/StreamFunc->hy*0.5;
                vy[i][j] = -(StreamFunc->u[i+1][j]-StreamFunc->u[i-1][j])/StreamFunc->hx*0.5;
            }


        //сетаем начальное распределение для уравнения переноса вихря

/*      //ГУ для течения Пуазейля
        for(i=1; i<Nx-1; i++)
        {
            StreamFunc->u[i][1]     = StreamFunc->u[i][0];
            StreamFunc->u[i][Ny-2]  = StreamFunc->u[i][Ny-1];


            CurlFunc->u[i][0]    = -(StreamFunc->u[i][0]    -2.*StreamFunc->u[i][1]     +StreamFunc->u[i][2]    )/StreamFunc->hy/StreamFunc->hy;
            CurlFunc->u[i][Ny-1] = -(StreamFunc->u[i][Ny-1] -2.*StreamFunc->u[i][Ny-2]  +StreamFunc->u[i][Ny-3] )/StreamFunc->hy/StreamFunc->hy;
        } //условия на границе твёрдого тела

        for(i=1; i<Ny-1; i++)
        {
            CurlFunc->u[0][i]    = //-(StreamFunc->u[0][i]        -2.*StreamFunc->u[1][i]   +StreamFunc->u[2][i]      )/StreamFunc->hx/StreamFunc->hx
                                   -(StreamFunc->u[0][i+1]      -2.*StreamFunc->u[0][i]   +StreamFunc->u[0][i-1]    )/StreamFunc->hy/StreamFunc->hy;
            CurlFunc->u[Nx-1][i] = //-(StreamFunc->u[Nx-1][i]     -2.*StreamFunc->u[Nx-2][i]+StreamFunc->u[Nx-3][i]   )/StreamFunc->hx/StreamFunc->hx
                                   -(StreamFunc->u[Nx-1][i+1]   -2.*StreamFunc->u[Nx-1][i]+StreamFunc->u[Nx-1][i-1] )/StreamFunc->hy/StreamFunc->hy;

//            if(i>2*Ny/5 && i<3*Ny/5)
//            {
//                StreamFunc->NodeState[Nx/4][i] = 1;
//                StreamFunc->u[Nx/4][i] = 0.5;
//                CurlFunc->u[Nx/4][i]    = -(StreamFunc->u[Nx/4-1][i] -2.*StreamFunc->u[Nx/4][i] +StreamFunc->u[Nx/4+1][i])/StreamFunc->hx/StreamFunc->hx;
//            }

        } //условия на границе "внешнего" течения (на стоки и источнике жидкости)
*/

/*      //ГУ для течения в каверне
        //
        for(i=1; i<Nx-1; i++)
        {
            StreamFunc->u[i][1]     = StreamFunc->u[i][0];
            StreamFunc->u[i][Ny-2]  = -StreamFunc->hy;

            CurlFunc->u[i][0]    = -(StreamFunc->u[i][0]    -2.*StreamFunc->u[i][1]     +StreamFunc->u[i][2]    )/StreamFunc->hy/StreamFunc->hy;
            CurlFunc->u[i][Ny-1] = -(StreamFunc->u[i][Ny-1] -2.*StreamFunc->u[i][Ny-2]  +StreamFunc->u[i][Ny-3] )/StreamFunc->hy/StreamFunc->hy;
        }

        for(i=1; i<Ny-1; i++)
        {
            StreamFunc->u[1][i]     = StreamFunc->u[0][i];
            StreamFunc->u[Nx-2][i]  = StreamFunc->u[Nx-1][i];

            CurlFunc->u[0][i]    = -(StreamFunc->u[0][i]        -2.*StreamFunc->u[1][i]   +StreamFunc->u[2][i]      )/StreamFunc->hx/StreamFunc->hx;
            CurlFunc->u[Nx-1][i] = -(StreamFunc->u[Nx-1][i]     -2.*StreamFunc->u[Nx-2][i]+StreamFunc->u[Nx-3][i]   )/StreamFunc->hx/StreamFunc->hx;
        }
*/

/*      //ГУ для течения Куэтта
        //
        for(i=1; i<Nx-1; i++)
        {
            StreamFunc->u[i][1]     = StreamFunc->u[i][0];
            StreamFunc->u[i][Ny-2]  = StreamFunc->u[i][Ny-1]-StreamFunc->hy;

            CurlFunc->u[i][0]    = -(StreamFunc->u[i][0]    -2.*StreamFunc->u[i][1]     +StreamFunc->u[i][2]    )/StreamFunc->hy/StreamFunc->hy;
            CurlFunc->u[i][Ny-1] = -(StreamFunc->u[i][Ny-1] -2.*StreamFunc->u[i][Ny-2]  +StreamFunc->u[i][Ny-3] )/StreamFunc->hy/StreamFunc->hy;
        }//два твёрдых тела, верхнее - подвижное, нижнее - покоится

        for(i=1; i<Ny-1; i++)
        {
            CurlFunc->u[0][i]    = //-(StreamFunc->u[0][i]        -2.*StreamFunc->u[1][i]   +StreamFunc->u[2][i]      )/StreamFunc->hx/StreamFunc->hx
                                   -(StreamFunc->u[0][i+1]      -2.*StreamFunc->u[0][i]   +StreamFunc->u[0][i-1]    )/StreamFunc->hy/StreamFunc->hy;
            CurlFunc->u[Nx-1][i] = //-(StreamFunc->u[Nx-1][i]     -2.*StreamFunc->u[Nx-2][i]+StreamFunc->u[Nx-3][i]   )/StreamFunc->hx/StreamFunc->hx
                                   -(StreamFunc->u[Nx-1][i+1]   -2.*StreamFunc->u[Nx-1][i]+StreamFunc->u[Nx-1][i-1] )/StreamFunc->hy/StreamFunc->hy;
        }//сток и исток, профили заданы...
*/

        //ГУ для пластинки в потоке (вертикальная)

        for(i=1; i<Nx-1; i++)
        {
            StreamFunc->u[i][1]     = StreamFunc->u[i][0];
            StreamFunc->u[i][Ny-2]  = StreamFunc->u[i][Ny-1];


            CurlFunc->u[i][0]    = -(StreamFunc->u[i][0]    -2.*StreamFunc->u[i][1]     +StreamFunc->u[i][2]    )/StreamFunc->hy/StreamFunc->hy;
            CurlFunc->u[i][Ny-1] = -(StreamFunc->u[i][Ny-1] -2.*StreamFunc->u[i][Ny-2]  +StreamFunc->u[i][Ny-3] )/StreamFunc->hy/StreamFunc->hy;
        } //условия на границе твёрдого тела

        for(i=1; i<Ny-1; i++)
        {
            CurlFunc->u[0][i]    = //-(StreamFunc->u[0][i]        -2.*StreamFunc->u[1][i]   +StreamFunc->u[2][i]      )/StreamFunc->hx/StreamFunc->hx
                    -(StreamFunc->u[0][i+1]      -2.*StreamFunc->u[0][i]   +StreamFunc->u[0][i-1]    )/StreamFunc->hy/StreamFunc->hy;
            CurlFunc->u[Nx-1][i] = //-(StreamFunc->u[Nx-1][i]     -2.*StreamFunc->u[Nx-2][i]+StreamFunc->u[Nx-3][i]   )/StreamFunc->hx/StreamFunc->hx
                    -(StreamFunc->u[Nx-1][i+1]   -2.*StreamFunc->u[Nx-1][i]+StreamFunc->u[Nx-1][i-1] )/StreamFunc->hy/StreamFunc->hy;
        } //условия на границе "внешнего" течения (на стоки и источнике жидкости)

        for(i = 1; i < Nx-1; i++)
            for( j = 1; j< Ny-1; j++)
                if(StreamFunc->NodeState[i][j]==1)
                {
                    StreamFunc->u[i][j] = 0.5;
                    CurlFunc->u[i][j]    =  -(StreamFunc->u[i-1][j] -2.*StreamFunc->u[i][j] +StreamFunc->u[i+1][j])/StreamFunc->hx/StreamFunc->hx
                                            -(StreamFunc->u[i][j-1] -2.*StreamFunc->u[i][j] +StreamFunc->u[i][j+1])/StreamFunc->hy/StreamFunc->hy;
                }//условия на твёрдые тела внутри потока


        int CurlTimeIterator = 0;
        for(CurlTimeIterator = 0; CurlTimeIterator < 1; CurlTimeIterator++) //интегрирование уравнения переноса вихря по времени
        {
            //переносная сила для вихря (форма, соотвующая субстанционной производной)
//#pragma omp parallel for collapse(2)
//            for(i=1;i<Nx-1;i++)
//                for(j=1;j<Ny-1;j++)
//                    CurlFunc->Force[i][j] = -(vx[i][j]*(CurlFunc->u[i+1][j]-CurlFunc->u[i-1][j])/CurlFunc->hx +
//                            vy[i][j]*(CurlFunc->u[i][j+1]-CurlFunc->u[i][j-1])/CurlFunc->hy)*0.5;

            //переносная сила для вихря в консервативной форме
            for(i=1;i<Nx-1;i++)
                for(j=1;j<Ny-1;j++)
                    CurlFunc->Force[i][j] =
                            -((CurlFunc->u[i+1][j]*vx[i+1][j]-CurlFunc->u[i-1][j]*vx[i-1][j])/CurlFunc->hx +
                              (CurlFunc->u[i][j+1]*vy[i][j+1]-CurlFunc->u[i][j-1]*vy[i][j-1])/CurlFunc->hy)*0.5;

            if(k==0 && CurlTimeIterator == 0)
                CurlFunc->StepInTime_Euler();
            else
                CurlFunc->StepInTime_Adams();
        }

        //обновляем вихревое возмущение для ФТ, а потом и её саму
        for(i=1;i<Nx-1;i++)
            for(j=1;j<Ny-1;j++)
            {
                StreamFunc->Force[i][j] = -CurlFunc->u[i][j];
            }
        StreamFunc->IterateWAutostop(50,1e-8);
//        StreamFunc->Plot(2);
        if(k%100==0)
        {
            CurlFunc->Plot(1);
            printf("updated\n");
            ret_code = system("sleep 0.1");
        }
        //        StreamFunc->ClosePlot();
    }
    StreamFunc->Plot(1);
    ret_code = system("sleep 3");

    FILE *velocityValue;
    velocityValue = fopen("velocity.dat","w");
    fprintf(velocityValue,"TITLE = \"Flow Model\"\n");
    fprintf(velocityValue,"VARIABLES = \"x\", \"y\", \"v_x\", \"v_y\"\n");
    fprintf(velocityValue,"ZONE T=\"Frame 0\", I=%d, J=%d\n", Ny, Nx);

    for(i=0;i<Nx;i++)
        for(j=0;j<Ny;j++)
            fprintf(velocityValue,"%f %f %f %f\n",StreamFunc->hx*i+StreamFunc->lx,
                                                         StreamFunc->hy*j+StreamFunc->ly,
                                                         vx[i][j], vy[i][j]);
                                               //sqrt(vx[i][j]*vx[i][j]+vy[i][j]*vy[i][j]));
    fclose(velocityValue);

    return 0;
}
/*
    fprintf(velocityValue,"# vtk DataFile Version 2.0\nsimple plot");
    fprintf(velocityValue,"ASCII\nDATASET RECTILINEAR_GRID\n");
    fprintf(velocityValue,"DIMENSIONS %d %d %d\n", Nx, Ny, 1);

    fprintf(velocityValue,"X_COORDINATES %d float\n", Nx);
    for(i=0;i<Nx;i++)
        fprintf(velocityValue,"%f\n",StreamFunc->hx*i+StreamFunc->lx);

    fprintf(velocityValue,"Y_COORDINATES %d float\n", Ny);
    for(i=0;i<Ny;i++)
        fprintf(velocityValue,"%f\n",StreamFunc->hy*i+StreamFunc->ly);

    fprintf(velocityValue,"Z_COORDINATES %d float\n", 1);
    fprintf(velocityValue,"%f\n",0.);

    fprintf(velocityValue,"POINT_DATA %d\n", Nx*Ny);
//    fprintf(velocityValue,"VECTORS velocity float\n");
    fprintf(velocityValue,"SCALARS streamFunc float 1\nLOOKUP_TABLE default");

    for(i=0;i<Nx;i++)
        for(j=0;j<Ny;j++)
//            fprintf(velocityValue,"%f %f 0.0\n", vx[i][j], vy[i][j]);
    fprintf(velocityValue,"%f\n", StreamFunc->u[i][j]);
                                               //sqrt(vx[i][j]*vx[i][j]+vy[i][j]*vy[i][j]));
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
