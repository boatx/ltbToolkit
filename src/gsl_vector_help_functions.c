#include "gsl_vector_help_functions.h"

//generate r values
int gsl_vector_gen(gsl_vector *vec, unsigned int n, double min, double max)
{
    CHECK_PTR(vec,GSL_EFAULT);

    unsigned int i=0;
    if (n > 1 && n <= vec->size)
    {
        double delta=(max-min)/((double)n-1);
        for (i=0;i<n;i++)
        {
            gsl_vector_set(vec,i,min + i*delta);
        }

        return GSL_SUCCESS;
    }
    else
    {
        return GSL_EINVAL;
    }
}

//generate values - f(r)
int gsl_vector_gen_fun(gsl_vector *vec, unsigned int n, gsl_vector *r, double (*f)(double))
{
    CHECK_PTR(r,GSL_EFAULT);
    CHECK_PTR(vec,GSL_EFAULT);
    CHECK_PTR(f,GSL_EFAULT);


    unsigned int i=0;
    if (vec->size == r->size)
    {
        if (n <= vec->size)
        {
            for (i=0;i<n;i++)
            {
                gsl_vector_set(vec,i,f(gsl_vector_get(r,i)));
            }

            return GSL_SUCCESS;
        }
        else
        {
            return GSL_EINVAL;
        }
    }
    else
    {
        return GSL_EBADLEN;
    }
}

int min(gsl_vector *vec, double x)
{
    CHECK_PTR(vec,GSL_EFAULT);

    size_t i=0;
    if (vec)
    {
        for (i=0; i<vec->size;i++)
        {
            if (gsl_vector_get(vec,i) > x)
            {
                gsl_vector_set(vec,i,x);
            }
        }

        return GSL_SUCCESS;
    }
    else
    {
        return GSL_EFAULT;
    }
}

int max(gsl_vector *vec, double x)
{
    CHECK_PTR(vec,GSL_EFAULT);

    size_t i=0;
    if(vec)
    {
        for (i=0; i<vec->size;i++)
        {
            if (gsl_vector_get(vec,i) < x)
            {
                gsl_vector_set(vec,i,x);
            }
        }
        return GSL_SUCCESS;
    }
    else
    {
        return GSL_EFAULT;
    }
}

int gsl_vector_substract(gsl_vector *a, gsl_vector *b)
{
    CHECK_PTR(a,GSL_EFAULT);
    CHECK_PTR(b,GSL_EFAULT);

    if (a->size != b->size)
    {
        return GSL_EBADLEN;
    }
    else
    {
        size_t i=0;
        for(i=0; i<a->size;++i)
        {
            a->data[i]-=b->data[i];
        }
    }

    return GSL_SUCCESS;
}

int gsl_vector_substract_out(gsl_vector *a, gsl_vector *b, gsl_vector *out)
{

    CHECK_PTR(a,GSL_EFAULT);
    CHECK_PTR(b,GSL_EFAULT);

    int status=GSL_SUCCESS;
    if (a->size != b->size || a->size != out->size)
    {
        return GSL_EBADLEN;
    }
    status = gsl_vector_memcpy(out,a);
    CHECK_STAT(status);
    status = gsl_vector_substract(out,b);
    return status;
}

int gsl_vector_add_constant_out(gsl_vector *dest, gsl_vector *src, double x)
{

    CHECK_PTR(dest,GSL_EFAULT);
    CHECK_PTR(src,GSL_EFAULT);

    if (dest->size != src->size)
    {
        return GSL_EBADLEN;
    }
    int status=GSL_SUCCESS;
    status = gsl_vector_memcpy(dest,src);
    CHECK_STAT(status);
    status=gsl_vector_add_constant(dest,x);
    return status;
}

int gsl_vector_load_from_file(gsl_vector *vec,char *fileName)
{
    CHECK_PTR(vec,GSL_EFAULT);
    int status=GSL_SUCCESS;
    //from file
    FILE *f=fopen(fileName,"r");
    if (f != NULL)
    {
        status = gsl_vector_fscanf(f,vec);
        fclose(f);
    }
    else
    {
        status=GSL_EFAULT;
    }

    return status;
}

int gsl_vector_set_range(gsl_vector *vec, unsigned int min, unsigned int max, double val)
{
    CHECK_PTR(vec,GSL_EFAULT);

    if (min > max || min < 0)
    {
        return GSL_EINVAL;
    }
    if (max >= vec->size)
    {
        return GSL_EBADLEN;
    }

    unsigned int i=0;

    for(i=min;i<max+1;++i)
    {
        vec->data[i]=val;
    }

    return GSL_SUCCESS;
}

//write Matrix to ascii file for gnuplot etc
int gsl_write_mat(gsl_matrix *matrix, char* fileName, char* header)
{
  CHECK_PTR(matrix,GSL_EFAULT);
  unsigned int i=0,j=0;
  FILE *f=fopen(fileName,"w");
  CHECK_PTR(f,GSL_EFAULT);

  if (header)
  {
    fprintf(f,header);
    fprintf(f,"\n");
  }

  for (i = 0; i < matrix->size1; i++)
  {
    for (j = 0; j < matrix->size2; j++)
    {
        fprintf (f,"%17.8f ",gsl_matrix_get(matrix,i,j));
    }
    fprintf(f,"\n");
  }
  fprintf(f,"\n");

  fclose(f);

  return GSL_SUCCESS;
}

//write two vector in two columns
int gsl_write_two_vec(const gsl_vector *vec, const gsl_vector *vec2, char* fileName)
{
    CHECK_PTR(vec,GSL_EFAULT);
    CHECK_PTR(vec2,GSL_EFAULT);
    int status=GSL_SUCCESS;
    gsl_matrix *mat = gsl_matrix_alloc(vec->size,2);
    CHECK_PTR(mat,GSL_ENOMEM);
    status=gsl_matrix_set_col(mat,0,vec);
    CHECK_STAT(status);
    status=gsl_matrix_set_col(mat,1,vec2);
    CHECK_STAT(status);
    status=gsl_write_mat(mat,fileName,NULL);
    CHECK_STAT(status);
    gsl_matrix_free(mat);

    return GSL_SUCCESS;
}

int gsl_write_vec(const gsl_vector *vec,char *fileName)
{
  CHECK_PTR(vec,GSL_EFAULT);
  FILE *f=fopen(fileName,"w");
  CHECK_PTR(f,GSL_EFAULT);
  unsigned int i=0;
  for (i=0;i<vec->size;++i)
  {
      fprintf(f,"%f\n",gsl_vector_get(vec,i));
  }

  fclose(f);
}

int find_indices_leq(const double *data, const double x, unsigned int n, unsigned int start)
{
    if(data)
    {
        if (start >= n)
        {
            return -1;
        }

        unsigned i=0;

        for (i=start;i<n;++i)
        {
            if(data[i]<=x)
            {
                return i;
            }
        }

        return -1;
    }
    else
    {
        return -1;
    }
}

int gsl_vector_find_indices_leq(const gsl_vector *vec, const double x, unsigned int start)
{
    if (start >= vec->size)
    {
        return GSL_EINVAL;
    }

    return find_indices_leq(vec->data,x,vec->size,start);
}

int find_indices_geq(const double *data, const double x, unsigned int n, unsigned int start)
{
    if(data)
    {
        if (start >= n)
        {
            return -1;
        }

        unsigned i=0;

        for (i=start;i<n;++i)
        {
            if(data[i]>=x)
            {
                return i;
            }
        }

        return -1;
    }
    else
    {
        return -1;
    }
}

int gsl_vector_find_indices_geq(const gsl_vector *vec, const double x, unsigned int start)
{
    if (start >= vec->size)
    {
        return GSL_EINVAL;
    }

    return find_indices_geq(vec->data,x,vec->size,start);
}

int find_indices_range(const double *data, const double xa, const double xb,
                       unsigned int n ,unsigned int start, int *ia, int *ib)
{

    if (data)
    {
        *ia=-1;
        *ib=-1;

        *ia=find_indices_geq(data,xa,n,start);
        if (*ia >= 0)
        {
            *ib=find_indices_geq(data,xb,n,(unsigned int)(*ia));
        }

        return GSL_SUCCESS;
    }
    else
    {
        return GSL_EFAULT;
    }
}

int gsl_vector_indices_range(const gsl_vector *vec, const double xa, const double xb,
                             unsigned int start, int *ia, int *ib)
{
    if (start >= vec->size)
    {
        return GSL_EINVAL;
    }

    return find_indices_range(vec->data,xa,xb,vec->size,start,ia,ib);
}

int data_to_string(double *data, char* format, unsigned int n, char *out)
{
    int status=0;
    if(data)
    {
        unsigned int i=0;
        for(i=0;i<n;++i)
        {
            status=sprintf(out,format,data[i]);
            if(status)
            {
                return status;
            }
        }

    }
    else
    {
        return -1;
    }

    return GSL_SUCCESS;
}

int gsl_vector_to_string(gsl_vector *vec, char* format, char* out)
{
    CHECK_PTR(vec,GSL_EFAULT);
    return data_to_string(vec->data,format,vec->size,out);
}

int gsl_matrix_row_normalize(gsl_matrix *mat)
{
    CHECK_PTR(mat,GSL_EFAULT);

    if(mat)
    {
        unsigned int i=0;
        double max=0;
        gsl_vector_view vv;

        for (i=0;i<mat->size1;++i)
        {
            vv = gsl_matrix_row(mat,i);
            max = gsl_vector_max(&(vv.vector));
            gsl_vector_scale(&(vv.vector),1.0/max);
        }

        return GSL_SUCCESS;

    }
    else
    {
        return GSL_EFAULT;
    }
}



