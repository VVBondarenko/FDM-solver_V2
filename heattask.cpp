#include "heattask.h"


HeatTask::HeatTask(double LX, double RX, double LY, double RY, int XSize, int YSize, double dt, double TCC)
{
    if(RX>LX)
    {
        lx = LX;
        rx = RX;
    }
    else
    {
        lx = RX;
        rx = LX;
    }

    if(RY>LY)
    {
        ly = LY;
        ry = RY;
    }
    else
    {
        ly = RY;
        ry = LY;
    }

    xSize = XSize;
    ySize = YSize;

    hx = (rx - lx)/(xSize-1);
    hy = (ry - ly)/(ySize-1);


    this->dt = dt;
    this->TCC = TCC;

    int i,j;

    u = new double* [xSize];
    for(i = 0; i < xSize; i++)
            u[i] = new double[ySize];

    Prev_du_dt = new double* [xSize];
    for(i = 0; i < xSize; i++)
        Prev_du_dt[i] = new double[ySize];

    Curr_du_dt = new double* [xSize];
    for(i = 0; i < xSize; i++)
        Curr_du_dt[i] = new double[ySize];

    for(i=0;i<xSize;i++)
    {
        for(j=0;j<ySize;j++)
        {
            u[i][j] = initial_f(lx+hx*i,ly+hy*j);
            Prev_du_dt[i][j] = 0.;
            Curr_du_dt[i][j] = 0.;
        }
    }

}

void HeatTask::StepInTime_Euler()
{
    int i,j;
    // du/dt = a \Delta u + f
    for(i = 1; i < xSize-1; i++)
    {
        for(j = 1; j < ySize-1; j++)
        {
            Prev_du_dt[i][j] = TCC*( (u[i-1][j]-2.*u[i][j]+u[i+1][j])/hx/hx  +
                                     (u[i][j-1]-2.*u[i][j]+u[i][j+1])/hy/hy) +
                                            rightpart_f(lx+i*hx,ly+j*hy,time);
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
}

void HeatTask::StepInTime_Adams()
{
    int i,j;
    // du/dt = a \Delta u + f
    for(i = 1; i < xSize-1; i++)
    {
        for(j = 1; j < ySize-1; j++)
        {
            Curr_du_dt[i][j] = TCC*( (u[i-1][j]-2.*u[i][j]+u[i+1][j])/hx/hx  +
                                     (u[i][j-1]-2.*u[i][j]+u[i][j+1])/hy/hy) +
                                            rightpart_f(lx+i*hx,ly+j*hy,time);
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

    for(i = 1; i < xSize-1; i++)
    {
        for(j = 1; j < ySize-1; j++)
        {
            Prev_du_dt[i][j] = Curr_du_dt[i][j];
        }
    }

}

void HeatTask::Output()
{
    int i, j;

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

    fclose(output);
}

void HeatTask::Plot(int Kind)
{
    Output();
    if(!Graph.isOpened())
    {
        std::vector<std::string> script;
        script.push_back("set terminal wxt");

        if(Kind == 0)
            script.push_back("splot \"result.dat\" w surface ");

        if(Kind == 1)
        {
            script.push_back("set view map");
            script.push_back("load '../parula.pal'");
            script.push_back("splot \"result.dat\" using 1:2:3 with image");
        }

        Graph.open();
        Graph.execute(script);
        PreviousPlotKind = Kind;
    }
    else
    {
        if(PreviousPlotKind == Kind)
        {
            Graph.write("replot");
            Graph.flush();
        }
        else
        {
            std::vector<std::string> script;
            script.push_back("reset");

            if(Kind == 0)
                script.push_back("splot \"result.dat\" w surface ");

            if(Kind == 1)
            {
                script.push_back("set view map");
                script.push_back("load '../parula.pal'");
                script.push_back("splot \"result.dat\" using 1:2:3 with image");
            }

//            Graph.open();
            Graph.execute(script);
            PreviousPlotKind = Kind;
        }
    }
}

void HeatTask::Animate(double FinalTime)
{
    ClosePlot();
    Output();
    std::vector<std::string> script;
    script.push_back("set terminal gif animate delay 30");
    script.push_back("set output \"heat_eq.gif\"");

    script.push_back("set xrange [0:1]");
    script.push_back("set yrange [0:1]");

    script.push_back("splot \"result.dat\"");

    Graph.open();
    Graph.execute(script);

    while(time < FinalTime)
    {
        for(int i=0; i<20; i++)
            this->StepInTime_Adams();
        Output();
//        usleep(500);
        Graph.write("splot \"result.dat\" ");
        Graph.flush();
    }

    ClosePlot();
}

void HeatTask::ClosePlot()
{
    if(Graph.isOpened())
        Graph.close();
}


