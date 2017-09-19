#ifndef HEATTASK_H
#define HEATTASK_H
#include <unistd.h>
#include <math.h>
#include <gnuplointerface.h>

class HeatTask
{
public:
    HeatTask(double LX, double RX,
             double LY, double RY,
             int XSize, int YSize,
             double dt, double TCC);

    double lx, rx;
    double ly, ry;
    int xSize, ySize;

    double time = 0.;
    double TCC = 1.; // ThermoConductivityCoefficient
    virtual double rightpart_f(double x, double y, double t)
    {
        return 0.*x*y*t;
    }

    virtual double boundary_u(double x, double y, double t)
    {
        return 0.*x*y*t;
    }

    virtual double initial_f(double x, double y)
    {
        return 0.*x*y;
    }

    void StepInTime_Euler();
    void StepInTime_Adams();
    void DoubleGrid(); //is it necessary at that moment?

    void Output();
    void Plot(int Kind);
    void Animate(double FinalTime);
    void ClosePlot();

    double hx, hy;
    double **u;
//private:
    double dt;
    double **Prev_du_dt;
    double **Curr_du_dt;

    GnuploInterface Graph;
    int PreviousPlotKind;
};

class HeatTaskWDiscreteForce : public HeatTask
{
public:
    double **Force;
    double **next_u;
    HeatTaskWDiscreteForce(double LX, double RX,
                           double LY, double RY,
                           int XSize, int YSize,
                           double dt, double TCC) : HeatTask(LX, RX,
                                                             LY, RY,
                                                             XSize, YSize,
                                                             dt, TCC)
    {
        int i, j;
        Force = new double* [xSize];
        next_u= new double* [xSize];
        for(i = 0; i < xSize; i++)
        {
            Force[i] = new double[ySize];
            next_u[i]= new double[ySize];
        }
        for(i=0;i<xSize;i++)
            for(j=0;j<ySize;j++)
            {
                Force[i][j] = 0.;
                next_u[i][j]= 0.;
            }
    }

    void StepInTime_Euler()
    {
        int i,j;
        // du/dt = a \Delta u + f
        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                Prev_du_dt[i][j] = TCC*( (u[i-1][j]-2.*u[i][j]+u[i+1][j])/hx/hx  +
                                         (u[i][j-1]-2.*u[i][j]+u[i][j+1])/hy/hy) +
                                   Force[i][j];
            }
        }

        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                u[i][j] += dt*Prev_du_dt[i][j];
            }
        }

        time += dt;

/*
        for(i=0;i<xSize;i++)
        {
            u[i][0]         = boundary_u(lx+i*hx,ly,time);
            u[i][ySize-1]   = boundary_u(lx+i*hx,ry,time);
        }

        for(i=0;i<ySize;i++)
        {
            u[0][i]       = boundary_u(lx,ly+i*hy,time);
            u[xSize-1][i] = boundary_u(rx,ly+i*hy,time);
        }
*/
    }

    void StepInTime_Adams()
    {
        int i,j;
        // du/dt = a \Delta u + f
        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                Curr_du_dt[i][j] = TCC*( (u[i-1][j]-2.*u[i][j]+u[i+1][j])/hx/hx  +
                                         (u[i][j-1]-2.*u[i][j]+u[i][j+1])/hy/hy) +
                                   Force[i][j];
            }
        }

        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                u[i][j] += dt*(1.5*Curr_du_dt[i][j]-0.5*Prev_du_dt[i][j]);
            }
        }

        time += dt;

//        for(i=0;i<xSize;i++)
//        {
//            u[i][0]         = boundary_u(lx+i*hx,ly,time);
//            u[i][ySize-1]   = boundary_u(lx+i*hx,ry,time);
//        }

//        for(i=0;i<ySize;i++)
//        {
//            u[0][i]       = boundary_u(lx,ly+i*hy,time);
//            u[xSize-1][i] = boundary_u(rx,ly+i*hy,time);
//        }

        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                Prev_du_dt[i][j] = Curr_du_dt[i][j];
            }
        }

    }

    void StepInTime_Crank()
    {
        int i,j,k;
        double hx_s = hx*hx,
               hy_s = hy*hy;

        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                Prev_du_dt[i][j] = TCC*( (u[i-1][j]-2.*u[i][j]+u[i+1][j])/hx/hx  +
                                         (u[i][j-1]-2.*u[i][j]+u[i][j+1])/hy/hy) +
                                   Force[i][j];

                next_u[i][j] = u[i][j] + Prev_du_dt[i][j]*dt;
            }
        }

        for(k=0;k<10;k++)
        {
            for(i = 1; i < xSize-1; i++)
            {
                for(j = 1; j < ySize-1; j++)
                {
                    next_u[i][j] = u[i][j]+ (TCC*((next_u[i-1][j]-2.*next_u[i][j]+next_u[i+1][j])/hx_s  +
                                                  (next_u[i][j-1]-2.*next_u[i][j]+next_u[i][j+1])/hy_s) +
                                                  Force[i][j] +  Prev_du_dt[i][j])*dt*0.5;
                }
            }

        }


        for(i = 1; i < xSize-1; i++)
        {
            for(j = 1; j < ySize-1; j++)
            {
                u[i][j] = next_u[i][j];
            }
        }

        time += dt;

    }
};



#endif // HEATTASK_H
