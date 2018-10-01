#
Perlude.
AC_INIT([VtCombiner], [0.0.1], [stephanie.georges@genetics.utah.edu])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

#Fixing file creation time
find ${srcdir} -name aclocal.m4 -exec touch {} \;
find ${srcdir} -name \*.in -exec touch {} \;
find ${srcdir} -name configure -exec touch {} \;

# Checks for programs
AC_PROG_CXX
AM_PROG_AR
LT_PREREQ([2.2])
LT_INIT()
AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_14

AC_CONFIG_HEADERS([config.h])

CFLAGS=
CPPFLAGS=
CXXFLAGS=

## START: DEBUG / RELEASE ##
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging, default: no])],
              [case "${enableval}" in
                yes) debug=true ;;
                no)  debug=false ;;
                *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
               esac],
              [debug=false])

AM_CONDITIONAL([DEBUG], [test x"$debug" = x"true"])
AS_IF([test "x$debug" = xtrue], [AC_DEFINE([LOGLVL], [LOGLV_DEBUG], [Define log level to be debug])])
## END: DEBUG / RELEASE ##

## START: PROFILING ##
AC_ARG_ENABLE([profiling],
              [AS_HELP_STRING([--enable-profiling],
                              [enable profiling, default: no])],
              [case "${enableval}" in
                yes) profiling=true ;;
                no)  profiling=false ;;
                *)   AC_MSG_ERROR([bad value ${enableval} for --enable-profiling]) ;;
               esac],
              [profiling=false])

AM_CONDITIONAL([PROFILING], [test x"$profiling" = x"true"])
## END: DEBUG / RELEASE ##

AC_ARG_WITH([htslib],
            [AS_HELP_STRING([--with-htslib=<path-to-htslib>],
                            [Specify the path to which htslib this program will link against])])

AS_IF([test "x$with_htslib" != x], [CPPFLAGS="$CPPFLAGS -I$with_htslib"])
AS_IF([test "x$with_htslib" != x], [LDFLAGS="$LDFLAGS -L$with_htslib"])

AC_ARG_WITH([gsl],
            [AS_HELP_STRING([--with-gsl=<path-to-gsl>],
                            [Specify the path to which gnu scientific library this program will link against])])

AS_IF([test "x$with_gsl" != x], [CPPFLAGS="$CPPFLAGS -I$with_gsl/include"])
AS_IF([test "x$with_gsl" != x], [LDFLAGS="$LDFLAGS -L$with_gsl/lib"])


## Checks for libraries
AC_CHECK_LIB([m], [cos])
AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], [AC_MSG_ERROR([Could not find gslcblas shared library])])
AC_CHECK_LIB([gsl], [gsl_stats_mean], [], [AC_MSG_ERROR([Could not find gsl shared library])])

AC_CHECK_LIB([hts], [bcf_read], [], [AC_MSG_ERROR([Could not find htslib shared library])])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

echo \
"----------------------------------------------
${PACKAGE_NAME} Version ${PACKAGE_VERSION}

Prefix: '${prefix}'
Compiler: '${CXX}'
CFLAGS: '${CFLAGS}'
CPPFLAGS: '${CPPFLAGS}'
CXXFLAGS: '${CXXFLAGS}'

Extra variables:
LDFLAGS='${LDFLAGS}'
LDADD='${LDADD}'
LIBS='${LIBS}'

PTHREAD_LIBS='${PTHREAD_LIBS}'
PTHREAD_CFLAGS='${PTHREAD_CFLAGS}'
PTHREAD_CC='${PTHREAD_CC}'

DATADIR=${DATADIR}
debug=${debug}
"
