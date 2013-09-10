#ifndef __HELLABY87__
#define __HELLABY87__

#include <phys_constant.h>
#include <macros.h>
#include <gsl/gsl_vector.h>
#include <math.h>

/*! \file hellaby87.h
 * Hellaby87 LTB functions
 */
static double UNUSED(M1)=5e-5;
static double UNUSED(E1)=3e-6;
static double UNUSED(M0)=-666;
/*! Cosmological parameter Omega matter */
static double omm = 1.015;

/*! LTB M function
 */
double M_hellaby87(double r);

double M_hellaby87_der(double r, void* params);

/*! LTB M function first derivative (analytic formula)
 */
double M_r_hellaby87(double r);

/*! LTB M function second derivative (analytic formula)
 */
double M_r2_hellaby87(double r);

/*! LTB E function
 */
double E_hellaby87(double r);

/*! LTB E function first derivative (analytic formula)
 */
double E_r_hellaby87(double r);

/*! LTB E function second derivative (analytic formula)
 */
double E_r2_hellaby87(double r);

/*! LTB tB function
 */
double tB_hellaby87(double r);

/*! LTB tB function first derivative (analytic formula)
 */
double tB_r_hellaby87(double r);

/*! LTB tB function generate from vectors M and E*/
void tB_hellaby87_vector(gsl_vector *tB, gsl_vector *M, gsl_vector *E);

/*! Big Crunch time
 * \return Crunch time */
double tC_hellaby87(gsl_vector *M, gsl_vector *E);

/*! Initialize M0 parameter*/
void M0_hellaby_initialize();

#endif /* __HELLABY87__ */
