#include <gsl/gsl_matrix.h>
#include "ltb_error.h"

/*! \file ltb_integrate.h
 * \brief Functions for integrating
 * \details Details(TODO)
 * \author Lukasz Matczynski
 * \copyright GNU Public License
*/


/*! \brief Integrate using trapezoid method points kept in gsl_matrix
 *
 * This function is needed for calculating radial proper length d(t,r). To get
 * it we mast integrate points in rows of grr_sq matrix (\link
 * ltb_model_grr_sq() \endlink ). To integrate I'm using trapezoid method:
 * \f[
 *  \int_a^b f(x)\textrm{dx} \approx \frac{h}{2} \sum_{i=1}^{n}\left(f(x_i) +
 *  f(x_{i+1})\right)
 * \f]
 *
 *  Where \f$h\f$ size of integration interval.
 *  When r dimension of matrix grr_sq is n we calculate n integrals
 *  \f[
 *  \int_{m[0]}^{m[0]} f(x)\textrm{dx}\textrm{,} \int_{m[0]}^{m[1]}
 *  f(x)\textrm{dx}\textrm{,...,}\int_{m[0]}^{m[n-1]} f(x)\textrm{dx}
 *  \f]
 * Where \f$m[i]\f$ are matrix elements. First of this integrals is always 0.
 *
 * \param [in] mat matrix to integrate
 * \param [in] dr size of integration interval
 * \param [in] rstart element of row from which start to integrate
 * \param [in] rstop element of row on which stop integrating
 * \param [in] tstart row from which start to integrate
 * \param [in] tstop row on which stop integrating
 * \param [out] integral matrix where result of integration is kept
 *
 * \sa ltb_model_grr_sq()
 */
int gsl_matrix_integrate_trap_range(const gsl_matrix *mat,
                                    double dr,
                                    unsigned int rstart,
                                    unsigned int rstop,
                                    unsigned int tstart,
                                    unsigned int tstop,
                                    gsl_matrix *integral
                                   );
