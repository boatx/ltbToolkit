AC_INIT([ltbToolkit], [1.0], [lukmatcz@gmail.com])
AC_CONFIG_SRCDIR(src/ltb.c)
AM_INIT_AUTOMAKE([-Wall gnu])
AC_PROG_CC([gcc cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
                 ltbtoolkit.pc
                 Makefile
                 src/Makefile
                 examples/Makefile
                 ])

AC_CONFIG_MACRO_DIR([m4])

#DX_HTML_FEATURE(ON)
#DX_CHM_FEATURE(OFF)
#DX_CHI_FEATURE(OFF)
#DX_MAN_FEATURE(OFF)
#DX_RTF_FEATURE(OFF)
#DX_XML_FEATURE(OFF)
#DX_PDF_FEATURE(ON)
#DX_PS_FEATURE(OFF)
#DX_INIT_DOXYGEN(ltbToolkit, doxygen.cfg, doc)

AM_PROG_AR
LT_INIT
AC_PROG_LIBTOOL
AM_PROG_LIBTOOL
AC_OPENMP

AC_CHECK_LIB([gslcblas],[cblas_dgemv])
AC_CHECK_LIB([gsl],[gsl_spline_init],,
        [AC_MSG_ERROR(Could not find gsl. Hint (1)  find /usr -name "*gsl*" -print   Hint (2)  install gsl*.a  e.g. from  http://www.gnu.org/software/gsl/ and put it in /usr/local/lib/ or set environment variable LDFLAGS to be -L/my/directory/lib/ with the library and CPPFLAGS to be -I/my/directory/include for the include files.)
        ])

AC_CHECK_LIB([m], [cos])
AC_CHECK_FUNCS([floor pow sqrt cbrt])

AC_FUNC_MALLOC

AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_HEADERS([string.h])
AC_TYPE_SIZE_T

AC_OUTPUT
