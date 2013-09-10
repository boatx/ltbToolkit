#include "interp_data.h"

size_t eta_n = 10000;
double dxi=5;

//initialise for etaInterpolation
Interp_data* interp_data_alloc(double max, double min, unsigned int n, int cur, double (*f)(double,int))
{
    if (n == 0)
    {
        GSL_ERROR_VAL("n is zero",GSL_EINVAL,0);
    }

    unsigned int i=0;
    double delta=(max-min)/(double)n;

    Interp_data *idata = malloc(sizeof(Interp_data));
    if (!idata)
    {
        GSL_ERROR_VAL("idata not allocated",GSL_ENOMEM,0);
    }

    //initialise fields of struct
    idata->n = n;
    idata->cur = cur;
    idata->max = max;
    idata->min = min;
    idata->x = NULL;
    idata->y = NULL;
    idata->acc = NULL;
    idata->spline = NULL;

    idata->x = malloc(sizeof(double)*n);
    if (!idata->x)
    {
        interp_data_free(idata);
        GSL_ERROR_VAL("idata not allocated",GSL_ENOMEM,0);

    }

    idata->y=malloc(sizeof(double)*n);
    if (!idata->y)
    {
        interp_data_free(idata);
        GSL_ERROR_VAL("idata not allocated",GSL_ENOMEM,0);

    }

    for (i=0; i<n; i++)
    {
        idata->x[i]=min + (double)i*delta;
        idata->y[i]=f(idata->x[i],cur);
    }

    idata->acc = gsl_interp_accel_alloc();
    if (!idata->acc)
    {
        interp_data_free(idata);
        GSL_ERROR_VAL("idata not allocated",GSL_ENOMEM,0);
    }

    idata->spline = gsl_interp_alloc(gsl_interp_cspline, n);
    if (!idata->spline)
    {
        interp_data_free(idata);
        GSL_ERROR_VAL("idata not allocated",GSL_ENOMEM,0);
    }

    if (gsl_interp_init(idata->spline, idata->y, idata->x, n))
    {
        interp_data_free(idata);
        GSL_ERROR_VAL("idata spline error",GSL_EINVAL,0);
    }
    return idata;
}

void interp_data_free(Interp_data *idata)
{
    if (!idata)
        return;

    if (idata->x)
    {
        free(idata->x);
    }
    if (idata->y)
    {
        free(idata->y);
    }
    if (idata->spline)
    {
        gsl_interp_free(idata->spline);
    }
    if (idata->acc)
    {
        gsl_interp_accel_free(idata->acc);
    }

    idata->n = 0;
    free(idata);
}

double set_eta_max(double xi_max)
{
    double eta_max = 0;

    //eta_max must be greater than xi_max
    if (xi_max > 0)
    {
        eta_max = (1 + dxi) * xi_max;
    }
    else
    {
        eta_max = (1 - dxi) * xi_max;
    }

    return eta_max;
}

double set_eta_min(double xi_min)
{
    double eta_min = 0;

    //eta_min must be lower than xi_min
    if (xi_min > 0)
    {
        eta_min = (1 - dxi) * xi_min;
    }
    else
    {
        eta_min = (1 + dxi) * xi_min;
    }

    return eta_min;
}

int check_interp_data(Interp_data *eta_interp,double xi_min, double xi_max, int cur)
{
    size_t n = eta_interp->n;

    printf("xi_min = %f xi_max = %f\n cur=%d\n",xi_min,xi_max,cur);
    printf("xi_min = %f xi_max = %f\n cur=%d\n",eta_interp->min,eta_interp->max,eta_interp->cur);
    if (eta_interp->cur != cur || eta_interp->max < xi_max ||
        eta_interp->min > xi_min)
    {
        GSL_ERROR("eta_interp wrong parameters",GSL_EINVAL);
    }

    return GSL_SUCCESS;
}
