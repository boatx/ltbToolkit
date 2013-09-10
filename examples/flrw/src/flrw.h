#ifndef __FLRW__
#define __FLRW__

#include <phys_constant.h>
#include <math.h>
#include <stdio.h>

/*! \file flrw.h
 * Friedman Robertson Lemaitre Walker basis function for LTB model
 */
#define NORM 1

double RC_0;

void print_param_flrw(double t_min, double t_max);

void export_param_flrw(char *filename,double t_min, double t_max);

void initialise_flrw(double om_m0, double om_r0, double t0);

void set_ro_0();

void set_ro_c();

void set_M3();

void set_f2();

void set_RC_0();

/*! M function
 */
double M_flrw(double r);

/*! E function
 */
double E_flrw(double r);

/*! tB function
 */
double tB_flrw(double __attribute__ ((unused)) r);

double M_r_flrw(double r);

double E_r_flrw(double r);

double tB_r_flrw(double __attribute__ ((unused)) r);

#endif /* __FLRW__*/
