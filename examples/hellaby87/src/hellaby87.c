#include "hellaby87.h"


void M0_hellaby_initialize()
{
    M0 = omm / (2*G_grav*H_0*pow(omm-1.0,1.5));
}

double M_hellaby87(double r)
{
    if (M0 < 0.0)
    {
        M0_hellaby_initialize();
    }
    return M0*(1+M1*cos(r));
}

double M_hellaby87_der(double r, void UNUSED(*params))
{
    return M_hellaby87(r);
}

double M_r_hellaby87(double r)
{
    return -M0*M1*sin(r);
}

double M_r2_hellaby87(double r)
{
    return -M0*M1*cos(r);
}

double E_hellaby87(double r)
{
    return -0.5*(1 - E1*sin(r)*sin(r));
}

double E_r_hellaby87(double r)
{
    //return E1*sin(r)*cos(r);
    return E1*0.5*sin(2.0*r);
}

double E_r2_hellaby87(double r)
{
    return E1*cos(2*r);
}

double tB_hellaby87(double r)
{
    double tB=0.0;

    if (M0 < 0.0)
    {
        M0_hellaby_initialize();
    }
    //this one i found in article http://arxiv.org/abs/1201.5845
    //return -G_grav*M_hellaby87(r)/pow(-2.0*E_hellaby87(r),1.5) + G_grav*M0*(1 - M1);
    //and thi one was used in octave code
    tB=-pi*G_grav*M_hellaby87(r)/pow(-2.0*E_hellaby87(r),1.5);
    return tB;
}

double tB_r_hellaby87(double r)
{
    double tmp=0.0;

    if (M0 < 0.0)
    {
        M0_hellaby_initialize();
    }

    double E=-2.0*E_hellaby87(r);
    tmp = M_r_hellaby87(r)*pow(E,1.5)
          + 3*M_hellaby87(r)*sqrt(E)*E_r_hellaby87(r);
    return -pi*G_grav*tmp/pow(E,3);
}

void tB_hellaby87_vector(gsl_vector *tB, gsl_vector *M, gsl_vector *E)
{
    if (M0 < 0.0)
    {
        M0_hellaby_initialize();
    }

    int i=0;
    int N=tB->size;
    double x=G_grav*M0*(1-M1);
    for (i=0;i<N;i++)
    {
        tB->data[i]=-G_grav*M->data[i]/pow(-2*E->data[i],1.5) + x;
    }
}

double tC_hellaby87(gsl_vector *M, gsl_vector *E)
{
    double tC;
    int N=E->size;
    double x=G_grav*pi*M0*(1-M1);
    int i=0;
    gsl_vector *tmp=gsl_vector_alloc(N);
    for (i=0;i<N;i++)
    {
        tmp->data[i]=G_grav*pi*pow(-2.0*E->data[i],-1.5);
    }
    gsl_vector_mul(tmp,M);
    gsl_vector_add_constant(tmp,x);
    tC=gsl_vector_max(tmp);
    gsl_vector_free(tmp);
    return tC;
}

