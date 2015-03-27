AC_INIT([spgs], [1.0])

dnl Find the compilers and flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`dnl check for C and C++ compilers

dnl check for headers, types and type sizes
AC_CHECK_HEADERS([limits.h])
AC_CHECK_TYPES([uint64_t])
AC_CHECK_TYPES([uint_least64_t])
AC_CHECK_TYPES([unsigned long long])
AC_CHECK_SIZEOF([unsigned long])

dnl Make variable substitutions and specify output files
AC_CONFIG_FILES([src/Makevars])
AC_SUBST(DEFS)

AC_OUTPUT
