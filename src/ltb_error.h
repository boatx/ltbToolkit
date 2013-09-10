#ifndef __LTB_ERROR__
#define __LTB_ERROR__

#include <gsl/gsl_errno.h>
#include <gsl/gsl_nan.h>


/*! \file ltb_error.h
 *  \brief Errors Handling macros and functions.
 */

/*! \brif Ltb erro codes
 *
 */
enum {
    LTB_ECURV = 100,/**<Wrong value of curvature*/
    LTB_EEFUN = 101/**<E(r)<-0.5 metric degeneration*/

};

/*! Checking status macro - if status is > 0 then return status
 */
#define CHECK_STAT(status) \
    if(status){ \
        return status; \
    }

/*! Checking pointer macro - if pointer is NULL return specific error errno
 */
#define CHECK_PTR(stptr,errno) \
    if(!stptr){ \
        gsl_error ("Bad pointer", __FILE__, __LINE__, errno) ; \
        return errno; \
    }

/*! Checking if two matrices have the same shape */
#define CHECK_M_DIM(m1,m2,reason) \
    if (m1->size1 != m2->size1 || m1->size2 != m2->size2){ \
        gsl_error(reason,__FILE__, __LINE__, GSL_EBADLEN);\
        return GSL_EBADLEN; \
    }

/*! Checking if two values differ*/
#define CHECK_SIZE(s1,s2,reason) \
    if (s1 != s2){ \
        gsl_error(reason,__FILE__, __LINE__, GSL_EBADLEN);\
        return GSL_EBADLEN; \
    }

#endif /*__LTB_ERROR__*/
