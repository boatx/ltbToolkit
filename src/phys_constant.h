#include "macros.h"
#ifndef __PHYS_CONST__
#define __PHYS_CONST__

/*! \file phys_constant.h
 *  \brief Physics constant
 */

/*! pi constant */
static const double UNUSED(pi) = 3.14159265;
/*! Gravity constant G */
extern double G_grav;
/*! Hubbel parameter */
extern double H_0;
/*! c speed of light */
extern double c;

/*! \brief Default G parameter
 *c=1, Gpc / 10^10M_{sun}
 */
#define G_DEF 6.67300e-11 * 1.98892e40 / 3.08568025e25 / (299792458.0*299792458.0)

/*! \brief Default Hubble parameter
 * c=1 units Gpc*h^{-1}
 */
#define H_0_DEF 100.0 / 299792.458 * 1000

/*! \brief Default c --- speed of  light
 */
#define C_DEF 299792458.0;

/*! \brief Gyr/Gpc value
 * Parameter used to conver time values to distance values
 */
static const double UNUSED(Gyr_per_Gpc) =   3.08568025e16 / 299792458.0 / (365.25*24.0*3600.0);

/*! \brief Write physical constant to file
    \param [in] file
 */
void export_phys_const(char *file_name);

/*! \brief Set physical constant
 * \param [in] G_in Gravity Const.
 * \param [in] H_0_in actual value of Hubbel Parameter
 * \param [in] c_in speed of light
 */
void set_phys_const(double G_in, double H_0_in, double c_in);
#endif /*__PHYS_CONSTANT__*/
