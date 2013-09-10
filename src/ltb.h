#ifndef LTB_MAIN
#define LTB_MAIN
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_spline.h>
#include "ltb_error.h"
#include "gsl_vector_help_functions.h"
#include "phys_constant.h"
#include "interp_data.h"
#include <omp.h>

/*! \file ltb.h
 * \brief Ltb model function
 * \details Details(TODO)
 * \author Lukasz Matczynski
 * \copyright GNU Public License
*/

extern int tB_norm;
extern unsigned int chunk_size;

/*! Ltb_model typdef struct
 * Keeping all important information about specific LTB Model.
 */
typedef struct
{
    int cur; /*! Curvature of model values=-1,0,1 or -999 (error). */
    size_t r_size; /*!< Number of r points in model, size of r_vec. */
    size_t t_size; /*!< Number of t points in model. size of t_vec  */
    gsl_vector *r_vec; /*!< values of r --- gsl_vector. */
    gsl_vector *t_vec; /*!< values of t --- gsl_vector. */
    gsl_vector *tB_vec;/*!< values of function t_{b} --- gsl_vector*/
    gsl_vector *E_vec;/*!< values of function E --- gsl_vector*/
    gsl_vector *M_vec;/*!< values of function M --- gsl_vector*/
    gsl_matrix *R_mat;/*!< Matrix of R (areal radius) --- gsl_matrix of size r_vec x r_mat*/
    double (*tB_fun)(double); /*!< tB function */
    double (*E_fun)(double); /*!< E function */
    double (*M_fun)(double); /*!< M function */
    double tB_max; /*!<maximum value of tB function, needed for normalization */
}Ltb_model;

/*! \brief Setting curvature of model
 *
 * Checking all values of E gsl_vector If all values have the same sign set cur
 * parameter and return GSL_SUCCESS. If E values have different sign set cur
 * parameter to error value -999 and return LTB_ECURV.
 * \param [in] E vector of values of \f$E(r)\f$ LTB function
 * \param [out] cur curvature parameter
 * \return GSL_SUCCESS or LTB_ECURV
 * \sa set_cur()
 */
int set_cur(const gsl_vector *E_vec,int *cur);

/*!
 * \brief Calculate curvature from E
 *
 * \param [in] E value of \f$E(r)\f$ function from LTB model
 * \return curvature parameter from set {-1,0,1}
 */
int cur_fun(double E);

/*! \brief Calculate single values of Ricci Tensor.
 * \param [in] E value of LTB \f$E(r)\f$ function
 * \param [in] E_r first derivative of E
 * \param [in] R areal radius
 * \param [in] R_r partial derivative over r of R
 * \return value of Ricci Tensor (double)
 */
double Ricci_fun(double E, double E_r, double R, double R_r);

/*!\brief Calculate matrix of Ricci Tensor values.
 *
 * All pointers should be allocated otherwise return GSL_EFAULT.
 * All matrix should have equals dimensions all vector should have size equal
 * matrices second dimension otherwise return GSL_EBADLEN.
 *
 * \param [in] ltb pointer to \c Ltb_model (using \c ltb->R_mat)
 * \param [in] E_r first derivative of E
 * \param [in] R_r partial derivative over r of R
 * \param [out] Ric pointer to Ricci matrix
 *
 * \return error code or GSL_SUCCESS;
 * \sa Ricci_fun()
 */
int ltb_model_Ricci(Ltb_model *ltb, gsl_vector *E_r, gsl_matrix *R_r,gsl_matrix *Ric);

/*! Calculate density in point (\a t,\a r).
 *
 * Calculate density using equation \cite Bonnor1985 :
 * \f[
 *  \rho = \frac{M_{,r}}{4\pi R^2R_{,r}}
 * \f]
 *
 * \param M_r first derivative of \a M(r) function in point \a r
 * \param R value of areal radius in point (\a t,\a r)
 * \param R_r partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial r}\f$) in point (\a t,\a r)
 * \return value of density
 * \sa ltb_model_rho()
 */
double rho_fun(double M_r, double R, double R_r);

/*!\brief Calculate matrix of density in LTB model.
 *
 * All pointers should be allocated otherwise return GSL_EFAULT.
 * All matrix should have equals dimensions all vector should have size equal
 * matrices second dimension otherwise return GSL_EBADLEN.
 *
 * \param [in] ltb pointer to \c Ltb_model (using \c ltb->R_mat)
 * \param [in] M_r first derivative of \a M(r) function in point \a r
 * \param [in] R_r partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial r}\f$) in point (\a t,\a r)
 * \param [out] rho_mat matrix of density
 *
 * \return error code or GSL_SUCCESS;
 * \sa rho_fun()
 */
int ltb_model_rho(Ltb_model *ltb, gsl_vector *M_r, gsl_matrix *R_r, gsl_matrix *rho_mat);

/*! \brief Calculate root square of radial part of LTB metrics.
 *
 * \param [in] R_r partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial r}\f$) in point (\a t,\a r)
 * \param [in] E value of \f$E(r)\f$ LTB function
 *
 * \return value or sqrt of radial part of LTB metrics
 */
double grr_sq_fun(double R_r, double E);

/*! \brief Calculate matrix of root square of radial part of LTB model metrics .
 *
 * \param [in] E value of \f$E(r)\f$ LTB function
 * \param [in] R_r partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial r}\f$) in point (\a t,\a r)
 * \param [out] matrix of grr_sq values
 *
 * \return error code or GSL_SUCCESS;
 * \sa grr_sq_fun()
 */
int ltb_model_grr_sq(gsl_vector *E, gsl_matrix *R_r, gsl_matrix *grr_sq);

/*!\brief Calculate expansion parameter in LTB model.
 *
 * \param [in] R_r partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial r}\f$) in point (\a t,\a r)
 * \param [in] R_t partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial t}\f$) in point (\a t,\a r)
 * \param [in] R_rt partial derivative of areal radius \a R (\f$\frac{\partial^2 R}{\partial r \partial t}\f$) in point (\a t,\a r)
 *
 *
 * \return single value of expansion parameter
 */
double exp_fun(double R, double R_r, double R_t, double R_rt);

/*! \brief Calculate matrix of expansion parameters in LTB model.
 *
 *  \param [in] ltb pointer to \c Ltb_model (using \c ltb->R_mat)
 *  \param [in] R_r partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial r}\f$) in point (\a t,\a r)
 *  \param [in] R_t partial derivative of areal radius \a R (\f$\frac{\partial R}{\partial t}\f$) in point (\a t,\a r)
 *  \param [in] R_rt partial derivative of areal radius \a R (\f$\frac{\partial^2 R}{\partial r \partial t}\f$) in point (\a t,\a r)
 *  \param [out] exp matrix of expansion parameters
 *
 * \return error code or GSL_SUCCESS;
 *
 *  \sa expansion_fun()
 */
int ltb_model_exp(Ltb_model *ltb, gsl_matrix *R_r,
                  gsl_matrix *R_t, gsl_matrix *R_rt, gsl_matrix *exp);

/*! \brief Calculate vector of chi values
 *
 * \param [in] cur curvature parameter
 * \param [in] E value of \f$E(r)\f$ LTB function
 * \param [out] chi_vec vector of chi values
 *
 * \return error code or GSL_SUCCESS;
 */
int get_chi(int cur, gsl_vector *E, gsl_vector *chi_vec);

/*! \brief Calculate values of fi in model with curvature = 0
 *  \return fi value
 *  \sa fi()
 */
double fiFlat(double eta);

/*! \brief Calculate values of fi in model with curvature > 0
 *
 * \return fi value
 * \sa fi()
 */
double fiSpherical(double eta);

/*! \brief Calculate values of fi in model with curvature < 0
 *
 * \return fi value
 * \sa fi()
 */
double fiHyperbolic(double eta);

/*! \brief Calculate values of fi
 *
 * Depending on cur parameter choosing function to calculate \a fi(eta)
 *
 * \param [in] eta value of eta
 * \param [in] cur curvature parameter
 *
 * \return fi value
 * \sa fiFlat()
 * \sa fiSpherical()
 * \sa fiHyperbolic()
 */
double fi(double eta, int cur);

/*! \brief Calculet xi value using value of tB,M,E functions
 *
 * \param [in] tB value of \a tB function
 * \param [in] M value of \a M function
 * \param [in] E value of \a E function
 * \param [in] t value of time coordinate
 * \param [in] cur curvature parameter
 *
 * \return value of xi
 *
 */
double xi_fun(double tB, double M, double E, double t, int cur);

/*! \brief Calculate xi matrix using values of \a tB,\a M, \aE functions
 *
 * If cur parameter is not in {-1,0,1}, return GSL_NAN and throw up error
 *
 * \param [in] tB values of \a tB function
 * \param [in] M values of \a M function
 * \param [in] E values of \a E function
 * \param [in] t values of time coordinates
 * \param [out] xi_mat values of xi
 *
 * \return GSL_SUCESS or error code
 * \sa xi_fun_from_eta()
 */

int get_xi(gsl_vector *tB,gsl_vector *M, gsl_vector *t, gsl_vector *E,
           int cur,gsl_matrix *xi_mat);

/*! \brief Calculate xi value in model with curvature < 0
 *  \return xi value
 *  \sa xi_fun_from_eta()
 */
double xiHyperbolic(double eta);

/*! \brief Calculate values of xi in model with curvature = 0
 *  \return xi value
 *  \sa xi_fun_from_eta()
 */
double xiFlat(double eta);

/*! \brief Calculate values of xi in model with curvature > 0
 *  \return xi value
 *  \sa xi_fun_from_eta()
 */
double xiSpherical(double eta);

/*! \brief Calculate values of xi
 *
 * Depending on cur parameter choosing function to calculate \a xi(eta)
 *
 * \param [in] eta value of eta
 * \param [in] cur curvature parameter
 *
 * \return fi value
 * \sa xiFlat()
 * \sa xiSpherical()
 * \sa xiHyperbolic()
 */
double xi_fun_from_eta(double eta, int cur);

/*! \brief Allocating interpolation data for specific Ltb_model
 *
 * \param [in] model pointer to Ltb_model
 * \param [in] n numbers of point to initialize spline
 *
 * \return pointer to Interp_data* (or null pointer)
 *
 */
Interp_data* interp_data_alloc_for_ltb(const Ltb_model *model, size_t n);

/*! \brief Finding value of max and min for Interp_data
 *
 * \param [in] model pointer to Ltb_model
 * \param [out] max ---  max value of points to initialize spline procedure
 * \pram [out] min  --- mina value of points to initialize spline procedure
 *
 * \return GSL_SUCCESS or error code
 */
int get_min_max_for_interp(const Ltb_model *model, double *max, double *min);

/*! \brief Calculate values of eta from xi
 *
 * If curvature (cur = 0) is flat --- use analytic formula for \a eta(xi).
 * Otherwise use eta_interp parameter for spline interpolation of \a eta(xi).
 * If cur != 0 && eta_interp == NULL throw up error and return GSL_NAN.
 *
 * \param [in] xi value of xi
 * \param [in] cur curvature parameter
 * \param [in] eta_interp Interp_data needed for spline interpolation
 *  Using gsl spline interpolation
 *
 * \return Single eta value.
 */
double eta_fun(double xi,int cur, Interp_data *eta_interp);

/*! \brief Calculate \a R (areal radius) value
 *
 * \param [in] M value of \a M(r) function
 * \param [in] E value of \a E(r) function
 * \param [in] tB value of \a tB(r) function
 * \param [in] eta value of eta
 * \param [in] cur value of curvature parameter
 *
 * \return value of R
 * \sa ltb_model_set_R()
 */
double R_fun(double M, double E, double eta, int cur);

/*! \brief Allocating Ltb_model structure
 *
 * \c r_size and \c t_size have to be > 0 otherwise throw up error and return
 * NULL ptr. If can't allocate any of Ltb_model pointers return NULL.
 *
 * \param [in] _r_size number of \a r points (radial coordinates)
 * \param [in] _t_size number of \a t points (time coordinates)
 * \param [in] M pointer to \a M(r) function
 * \param [in] E pointer to \a E(r) function
 * \param [in] tB pointer to \a tB(r) function
 *
 * \return pointer to Ltb_model
 *
 * \sa ltb_model_alloc_from_vec()
 * \sa ltb_model_cpy()
*/
Ltb_model* ltb_model_alloc(size_t _r_size, size_t _t_size,
                           double (*M)(double),double (*E)(double),
                           double (*tB)(double));

/*! \brief Allocating Ltb_model structure from r vector and t vector
 *
 * After allocating Ltb_model with ltb_model_alloc() function (with r parameter
 * equal to _r_vec->size and t parameter to _t_vec->size) copy _r_vec and
 * _t_vec to Ltb_model \c r_vec and \c t_vec. If copying data fail throw up
 * error and return null.
 *
 * \param [in] _r_vec vector of \a r values (radial coordinate)
 * \param [in] _t_vec vector of \a t values (time coordinate)
 * \param [in] M pointer to \a M(r) function
 * \param [in] E pointer to \a E(r) function
 * \param [in] tB pointer to \a tB(r) function
 *
 * \return pointer to Ltb_model
 *
 * \sa ltb_model_alloc()
 * \sa ltb_model_cpy()
*/
Ltb_model* ltb_model_alloc_from_vec(gsl_vector *_r_vec, gsl_vector *_t_vec,
                                    double (*M)(double),double (*E)(double),
                                    double (*tB)(double));

/*! \brief Allocating Ltb_model making copy of src
 *
 * Allocate Ltb_model with ltb_model_alloc() function, after this copying of
 * data from Ltb_model *src to newly created. If copying of any data fail throw
 * up error and return NULL.
 *
 * \param [in] src Ltb_model to copy
 * \return pointer to Ltb_model
 *
 * \sa ltb_model_alloc_from_vec()
 * \sa ltb_model_alloc()
*/
Ltb_model* ltb_model_cpy(Ltb_model *src);

/*! \brief Freeing Ltb_model structure
 *
 * \param [in] model pointer to Ltb_model
 */
void ltb_model_free(Ltb_model *model);

/*! \brief Setting values of tB, M, E function for Ltb_model model
 * 
 * Calculating values of LTB functions (\a M, \a E, \a tB) using pointers tu
 * function \c model->M_fun, \c model->E_fun, \c model->tB_fun . If one of this
 * pointers is NULL return GSL_EFAULT. Calculated values are stored in
 * Ltb_model gsl_vectors (\c M_vec, \c E_vec and \c tB_vec). For tB is
 * calculate \c model->tB_max if global parameter \link tB_norm\endlink is set values in
 * \c tB_vec will be normalized (\c gsl_vector_max(tB_vec) == 0). Depending on
 * \c E_vec curvature parameter is calculated (E < 0 - cur = 1, E > 0 -- cur =
 * -1 , E == 0 cur=0, otherwise error will be throw up and cur == -999).
 *
 *  \param [in] model pointer to Ltb_model
 *  \return GSL_SUCCESS or error code
 *  \sa ltb_model_alloc()
 *  \sa ltb_model_set_R
 */
int ltb_model_set_fun(Ltb_model *model);

/*! \brief Calculating matrix of function \a R(t,r)
 *
 * eta_interp pointer is needed to calculate eta in curved model. If you only
 * interested in flat model you can pass NULL in place of eta_interp.
 *
 * \param [in] model pointer to Ltb_model
 * \param [in] eta_interp pointer to Inter_data structure
 *
 * \return GSL_SUCCESS or error code
 * \sa ltb_model_set_fun()
 * \sa ltb_model_alloc()
 */
int ltb_model_set_R(Ltb_model *model,Interp_data *eta_interp);

/*! Write basis function E,tB,M of Ltb_model model to ascii file file_name
 *
 * Writing to file content of \c r_vec, \c M_vec, \c E_vec, \c tB_vec in four
 * columns.
 * \param [in] model pointer to Ltb_model with functions to write
 * \param [in] file_name name of file
 *
 * \return GSL_SUCCESS or error code
 * \sa ltb_model_write_fun()
 */
int ltb_model_write_fun(const Ltb_model *model, char* file_name);

/*! \brief Calculating partial derivative of \a R over \a t
 *
 * \param [in] R areal radius value
 * \param [in] M value of \a M(r) function
 * \param [in] E value of \a E(r) function
 * \param [in] sgn sign of derivative (have to be calculated separately)
 *
 * \return Partial derivative of \a R over \a t
 * \sa ltb_model_R_t()
 *
 */
double R_t_analytic_fun(double R, double M, double E, int sgn);

/*! \brief Calculating partial derivatives of \a R over \a t
 *
 * \param [in] model pointer to Ltb_model
 * \param [in] eta_interp pointer to Interp_data needed for calculating sign of
 * derivative.
 * \param [out] R_t matrix of derivatives
 *
 * \return GSL_SUCCESS or error code
 * \sa R_t_analytic_fun()
 */
int ltb_model_R_t(const Ltb_model *model, Interp_data *eta_interp, gsl_matrix *R_t);

/*! \brief Calculating partial derivative of \a R over \a r
 *
 * \param [in] R areal radius value
 * \param [in] M value of \a M(r) function
 * \param [in] E value of \a E(r) function
 * \param [in] sgn sign of derivative (have to be calculated separately)
 *
 * \return Partial derivative of \a R over \a r
 * \sa ltb_model_R_r()
 */
double R_r_analytic_fun(double R,double R_t,double E,double E_r,
                        double M,double M_r,double tB,double tB_r,double t);

/*! \brief Calculating value of partial derivatives of \a R over \a r
 *
 * \param [in] R areal radius value
 * \param [in] R_t partial derivative of R over t
 * \param [in] E \a LTB E(r) function
 * \param [in] E_r derivative of \a E(r) function
 * \param [in] M \a LTB M(r) function
 * \param [in] M_r derivative of \a M(r) function
 * \param [in] tB LTB \a tB(r) function
 * \param [in] tB_r derivative of \a tB(r) function
 * \param [out] R_r matrix of partial derivative of \a R over \a r
 *
 * \return GSL_SUCCESS or error code
 * \sa ltb_model_R_r_analytic_fun()
 */
int ltb_model_R_r(const Ltb_model *model,const gsl_matrix *R_t,
                      const gsl_vector *E_r_vec, const gsl_vector *M_r_vec,
                      const gsl_vector *tB_r_vec, gsl_matrix *R_r);

/*! \brief Calculating value of partial mixed derivative of \a R over \a r and \a t
 *
 * \param [in] R areal radius value
 * \param [in] R_r partial derivative of R over t
 * \param [in] R_t partial derivative of R over t
 * \param [in] E_r derivative of \a E(r) function
 * \param [in] M \a LTB M(r) function
 * \param [in] M_r derivative of \a M(r) function
 *
 * \return value of partial mixed derivative of \a R over \a r and \a t
 * \sa ltb_model_R_rt()
 */
double R_rt_analytic_fun(double R, double R_r, double R_t, double E_r,double M, double M_r);

/*! \brief Calculating value of partial mixed derivative of \a R over \a r and \a t
 *
 * \param [in] model pointer to Ltb_model
 * \param [in] R_r partial derivative of R over t
 * \param [in] R_t partial derivative of R over t
 * \param [in] E_r derivative of \a E(r) function
 * \param [in] M_r derivative of \a M(r) function
 * \param [out] R_rt matrix of mixed partial derivative
 *
 * \return GSL_SUCCESS or error code
 * \sa R_rt_analytic_fun()
 */
int ltb_model_R_rt(const Ltb_model *model,
                      const gsl_matrix *R_r, const gsl_matrix *R_t,
                      const gsl_vector *E_r, const gsl_vector *M_r,
                      gsl_matrix *R_rt);

#endif /*LTB_MAIN*/
