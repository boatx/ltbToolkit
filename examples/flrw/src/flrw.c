#include "flrw.h"
//blad om_k_0 ma przeciwny znak
double flrw_t0 = 0.0;
double om_m_0 = 1.0;
double om_r_0 = 0.0;
double om_k_0 = 0.0;

double ro_0 = -999;
double ro_c = -999;
double RC_0 = -999;
double M3=-999;
double f2 = -999;

void print(FILE *f, double t_min, double t_max)
{
    fprintf(f,"M3=%g\n",M3);
    fprintf(f,"f2=%g\n",f2);
    fprintf(f,"Om_m_0=%g\n",om_m_0);
    fprintf(f,"Om_r_0=%g\n",om_r_0);
    fprintf(f,"Om_k_0=%g\n",om_k_0);
    fprintf(f,"t0=%g\n",flrw_t0);
    fprintf(f,"RC_0=%g\n",RC_0);
    fprintf(f,"ro_c=%g\n",ro_c);
    fprintf(f,"ro_0=%g\n",ro_0);
    fprintf(f,"k=-1.0*f2\n");
    fprintf(f,"t_min=%g\n",t_min);
    fprintf(f,"t_max=%g\n",t_max);
    if (om_k_0 != 0)
        fprintf(f,"a_0 = %g\n",sqrt(f2/(1-om_m_0 - om_r_0))/H_0);
    else
        fprintf(f,"a_0 = %g\n",1.0);

}

void print_param_flrw(double t_min, double t_max)
{
    if (ro_0 < 0)
    {
        set_ro_0();
        set_RC_0();
        set_M3();
    }
    print(stdout,t_min,t_max);
}

void export_param_flrw(char *filename,double t_min, double t_max)
{
    FILE *f;
    f = fopen(filename,"w+");
    print(f,t_min,t_max);
    fclose(f);
}

void initialise_flrw(double om_m0, double om_r0, double t0)
{
    om_m_0 = om_m0;
    om_r_0 = om_r0;
    flrw_t0 = t0;
    set_ro_c();
    set_ro_0();
    set_f2();
    set_M3();
    if (om_k_0 != 0)
    {
        set_RC_0();
    }
}

void set_ro_0()
{
    if (ro_c < 0)
    {
        set_ro_c();
    }
    ro_0 = ro_c*om_m_0;
}

void set_ro_c()
{
    ro_c=3*H_0*H_0/(8*pi*G_grav);
}

void set_f2()
{
    //S0 (scale factor now) have to be 1.0
    //f2=-k
    if (NORM)
    {
        f2 = -1.0*(om_m_0 + om_r_0 - 1)*H_0*H_0;
    }
    else
    {
        if (om_m_0 > 1.0)
        {
            f2 = -1.0;
        }
        else if (om_m_0 < 1.0)
        {
            f2 = 1.0;
        }
        else
        {
            f2 = 0.0;
        }
    }
}

void set_M3()
{
    M3 = om_m_0/(2*H_0*G_grav);
    om_k_0 = 1 - om_m_0 - om_r_0;
    if (om_k_0 != 0)
    {
        //M3 = 4.0 * pi * ro_0 / 3.0;
        //M3=S0*S0*S0*8*pi*ro_0/6.0;
        M3 = M3*H_0*H_0*H_0;
        //M3=M3*pow(f2/om_k_0,1.5);
    }
    else
    {
        M3 = 4.0 * pi * ro_0 / 3.0;
    }
}

void set_RC_0()
{
    if (om_k_0 != 0)
    {
        RC_0 = 1.0 / (H_0*sqrt(fabs(om_m_0+om_r_0 - 1.0)));
    }
    else
    {
        RC_0 = -999.0;
    }
}


//Flat model

double M_flrw(double r)
{
    if (ro_0 < 0)
    {
        set_ro_0();
        set_M3();
    }

    if (om_k_0 != 0)
    {
        if (r >= RC_0*pi/2.0)
        {
            puts("r>RC_0*pi/2");
            r = pi*RC_0 - r;
        }
    }

    return M3*r*r*r;
}

double M_r_flrw(double r)
{
    if (om_k_0 != 0)
    {
        if (r >= RC_0*pi/2.0)
        {
            r = pi*RC_0 - r;
            return -3.0*M3*r*r;
        }
    }

    return 3.0*M3*r*r;

}

double E_flrw(double r)
{
    if (f2 == -999)
    {
        set_f2();
    }

    if (om_k_0 != 0)
    {
        if (r >= RC_0*pi/2.0)
        {
            puts("r>RC_0*pi/2");
            r = pi*RC_0 - r;
        }
    }


    return 0.5*f2*r*r;
}

double E_r_flrw(double r)
{
    if (om_k_0 != 0)
    {
        if (r >= RC_0*pi/2.0)
        {
            r = pi*RC_0 - r;
            return -0.5*f2*2.0*r;
        }
    }

    return 0.5*f2*2.0*r;
}

double tB_flrw(double __attribute__ ((unused)) r)
{

    return flrw_t0;
}

double tB_r_flrw(double __attribute__ ((unused)) r)
{

    return 0.0;
}

