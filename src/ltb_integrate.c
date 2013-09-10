#include "ltb_integrate.h"

int gsl_matrix_integrate_trap_range(const gsl_matrix *mat,
                                    double dr,
                                    unsigned int rstart,
                                    unsigned int rstop,
                                    unsigned int tstart,
                                    unsigned int tstop,
                                    gsl_matrix *integral
                                   )
{
    CHECK_PTR(integral,GSL_EFAULT);
    CHECK_PTR(mat,GSL_EFAULT);
    CHECK_M_DIM(mat,integral,"Mat matrix and integral matrix have diffrent dimensions");

    if (rstart>=rstop)
    {
        GSL_ERROR("rstart <= rstop",GSL_EINVAL);
    }

    if (tstart>=tstop)
    {
        GSL_ERROR("tstart <= tstop",GSL_EINVAL);
    }

    if (tstop > mat->size1)
    {
        GSL_ERROR("tstop >= mat->size1",GSL_EINVAL);
    }

    if (rstop > mat->size2)
    {
        GSL_ERROR("rstop >= mat->size2",GSL_EINVAL);
    }

    unsigned i,j;
    double sum;

    gsl_matrix_set_all(integral,0.0);
    for (i=tstart;i<tstop;++i)
    {
        sum = 0.0;

        for (j=rstart;j<rstop-1;++j)
        {
            sum += gsl_matrix_get(mat,i,j) + gsl_matrix_get(mat,i,j+1);
            gsl_matrix_set(integral,i,j+1,sum*dr*0.5);

        }
    }

    return GSL_SUCCESS;
}
