#ifndef __HELP_MACROS__
#define __HELP_MACROS__

/*! \file macros.h
 * \brief Additional helpful macros
 */

/*! \brief Prevent compilator from giving warrnig about unused variable 
 */
#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#else
# define UNUSED(x) x
#endif



#endif /*__HELP_MACROS__*/
