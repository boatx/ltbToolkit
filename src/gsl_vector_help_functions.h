#ifndef __GSL_VEC_MAT_HELP_FUN__
#define __GSL_VEC_MAT_HELP_FUN__

#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "ltb_error.h"

/*! \file gsl_vector_help_functions.h
 *  \brief Custom extensions of gsl_vector, gsl_matrix from GSL library.
 *
 *  Most of function from this file was used during testes.
 *  I only use functions for writing vectors and matrices and function for
 *  generating values from function (\link gsl_vector_gen_fun()\endlink and
 *  \link gsl_vector_gen() \endlink)
 */

#define M(m,i,j) (gsl_matrix_get(m,i,j))
#define MS(m,i,j,x) (gsl_matrix_set(m,i,j,x))

#define V(v,i) (gsl_vector_get(v,i))
#define VS(v,i,x) (gsl_vector_set(v,i,x))


/*! \brief Generate n values from min to max and kept them in gsl_vector vec.
 */
int gsl_vector_gen(gsl_vector *vec, unsigned int n, double min, double max);

/*! \brief Generate n values of function f(r) and kept them in gsl_vector vec.
 */
int gsl_vector_gen_fun(gsl_vector *vec, unsigned int n, gsl_vector *r, double (*f)(double));

/*! \brief If ith element of vector vec is greater than x -> set x.
 */
int min(gsl_vector *vec, double x);

/*! \brief If ith element of vector vec is lesser than x -> set x.
 */
int max(gsl_vector *vec, double x);

/*! \brief Subtraction of gsl_vector's a=a-b.
 */
int gsl_vector_substract(gsl_vector *a, gsl_vector *b);

/*! \brief Subtraction of gsl_vector's out=a-b.
 */
int gsl_vector_substract_out(gsl_vector *a, gsl_vector *b, gsl_vector *out);

/*! \brief Create new gsl_vector dest = src + x.
 */
int gsl_vector_add_constant_out(gsl_vector *dest, gsl_vector *src, double x);

/*! \brief Set value val in vec over indexes from min to max.
 */
int gsl_vector_set_range(gsl_vector *vec, unsigned int min, unsigned int max, double val);

/*! \brief Load gsl_vector vec from ascii file.
 */
int gsl_vector_load_from_file(gsl_vector *vec,char *fileName);

/*! \brief Write gsl_matrix matrix to ascii file fileName with header.
 */
int gsl_write_mat(gsl_matrix *matrix, char* fileName, char* header);

/*! \brief Write to vectors to file in two columns
 */
int gsl_write_two_vec(const gsl_vector *vec, const gsl_vector *vec2, char* fileName);

/*! \brief Find first indices of array data less equal than x.
 */
int find_indices_leq(const double *data, const double x, unsigned int n, unsigned int start);

/*! \brief Find first indices of gsl_vector vec less equal than x.
 */
int gsl_vector_find_indices_leq(const gsl_vector *vec, const double x, unsigned int start);

/*! \brief Find first indices of array data greater equal than x.
 */
int find_indices_geq(const double *data, const double x, unsigned int n, unsigned int start);

/*! \brief Find first indices of gsl_vector vec greater equal than x.
 */
int gsl_vector_find_indices_geq(const gsl_vector *vec, const double x, unsigned int start);

/*! \brief Find indices indices from range from xa to xb in array data.
 */
int find_indices_range(const double *data, const double xa, const double xb,
                       unsigned int n ,unsigned int start, int *ia, int *ib);


/*! \brief Find indices indices from range from xa to xb in gsl_vector vec.
 */
int gsl_vector_indices_range(const gsl_vector *vec, const double xa, const double xb,
                             unsigned int start, int *ia, int *ib);

/*! \brief Convert double array to string. */
int data_to_string(double *data, char* format, unsigned int n, char *out);

/*! \brief Convert gsl_vector to string. */
int gsl_vector_to_string(gsl_vector *vec, char* format, char* out);

/*! \brief Normalize rows of matrix to 1 */
int gsl_matrix_row_normalize(gsl_matrix *mat); 

#endif /*__GSL_VEC_MAT_HELP_FUN__*/
