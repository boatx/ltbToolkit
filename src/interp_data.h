#include <gsl/gsl_spline.h>
#include "ltb_error.h"

/*! \file interp_data.h
 * \brief Structure and functions needed for interpolation \a eta(xi)
 * \author Lukasz Matczynski
 * \copyright GNU Public License
*/

/*! Interp_data typdef stucture
 *
 * Used in inverting xi function to get eta
 */
typedef struct
{
    int n ;/*! Numbers of point to initialize spline procedure */
    int cur; /*! curvature of model for which we want interpolation */
    double *x; /*! x points for spline initialization */
    double *y; /*! values of function in point x for spline initialization */
    double max; /*! maximum value of points to initialization */
    double min;/*! minimum value op points to initialization */
    gsl_interp_accel *acc; /*! LUT (look up table) for accelerating of interpolation */
    gsl_interp *spline; /*! gsl interp structure needed for interpolation */
}Interp_data;

size_t eta_n; /*! Default numbers of points to initialize spline interpolation*/
/*!\brief Parameter used to calculate min and max value in Interp_data. 
 * It telling how much bigger (lower) max (min) should be than max (min) value of xi.
 * Used in set_eta_max() and set_eta_min() functions.
 *
 * \sa set_eta_max()
 * \sa set_eta_min()
 */

double dxi; 
/*! \brief Allocating Interp_data
 *
 * \param [in] max --- maximum value of points to initialize spline function
 * \param [in] min --- minimum value of points to initialize spline function
 * \param [in] n --- number of points to initialize spline function
 * \param [in] cur  curvature of model for which we want interpolate \eta(xi)
 * \param [in] f pointer to function calculating eta
 *
 * \return pointer to allocate Interp_data stuct or NULL pointer
 */
Interp_data* interp_data_alloc(double max,double min, unsigned int n, int cur,double (*f)(double,int));

/*! \brief Freeing Interp_data stuct
 *
 * \param [in] idata pointer to Interp_data structure which we want to free
 */
void interp_data_free(Interp_data *idata);

/*! \brief Calculating max value for Interp_data structure
 *
 * \param [in] xi_max maximum value of xi
 */
double set_eta_max(double xi_max);

/*! \brief Calculating min value for Interp_data structure
 *
 * \param [in] xi_mim minimum value of xi
 */
double set_eta_min(double xi_min);

/*! \brief Check Interp_data 
 *
 * Return error if \c cur != \c eta_interp->cur, \c xi_min < \c eta_interp->min
 * or \c xi_max > \c eta_interp->max
 *
 * \param [in] eta_interp --- Interp_data to check
 * \param [in] xi_min 
 * \param [in] xi_max
 * \param [in] cur curvature parameter
 *
 */
int check_interp_data(Interp_data *eta_interp,double xi_min, double xi_max, int cur);

