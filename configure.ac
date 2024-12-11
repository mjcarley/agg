#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

##AC_PREREQ([2.69])
AC_INIT([agg],[0.1.0],[codes@paraffinalia.co.uk])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/body.c])
AC_CONFIG_HEADERS([config.h])

AGG_MAJOR_VERSION=0
AGG_MINOR_VERSION=1
AGG_MICRO_VERSION=0
AGG_INTERFACE_AGE=0
AGG_BINARY_AGE=0
AGG_VERSION=$AGG_MAJOR_VERSION.$AGG_MINOR_VERSION.$AGG_MICRO_VERSION

AC_SUBST(AGG_MAJOR_VERSION)
AC_SUBST(AGG_MINOR_VERSION)
AC_SUBST(AGG_MICRO_VERSION)
AC_SUBST(AGG_VERSION)

# libtool versioning
LT_RELEASE=$AGG_MAJOR_VERSION.$AGG_MINOR_VERSION
LT_CURRENT=`expr $AGG_MICRO_VERSION - $AGG_INTERFACE_AGE`
LT_REVISION=$AGG_INTERFACE_AGE
LT_AGE=`expr $AGG_BINARY_AGE - $AGG_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

dnl Initialize libtool
LT_INIT

# For automake.
VERSION=$AGG_VERSION
PACKAGE=agg

# Checks for programs.
AC_PROG_CC

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -D_GNU_SOURCE"
fi
  
AC_PROG_AWK
AC_SEARCH_LIBS([strerror],[cposix])
##AC_HEADER_STDC

PKG_CHECK_MODULES([GLIB], [glib-2.0], [],
[
  echo \#################################################################
  echo \# 
  echo \# GLIB-2.0 not installed. You can install it using your 
  echo \# distribution's package manager or from:
  echo \#
  echo \# https://gitlab.gnome.org/GNOME/glib/
  echo \#
  echo \#################################################################
  exit
])

AC_CHECK_PROG(BLASWRAP, blaswrap-config, yes, no)
if test x$BLASWRAP = xyes ; then
   AC_DEFINE(HAVE_BLASWRAP, 1, [BLAS wrapper header available])
   blaswrap_cflags=`blaswrap-config --cflags`
   blaswrap_libs=`blaswrap-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# BLAS wrapper library not installed or not properly installed.
  echo \# You can get it from
  echo \#
  echo \# https://github.com/mjcarley/blaswrap
  echo \#
  echo \#################################################################
  exit
fi

# Checks for header files.
AC_CHECK_HEADERS([triangle_api.h])

# Checks for typedefs, structures, and compiler characteristics.
AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx"
fi
if test x$hasavx2 = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx2"
fi
if test x$hasfma = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mfma"
fi

AC_CHECK_LIB(m, cos)
AC_CHECK_LIB(matheval, evaluator_create, [], [
  echo \#################################################################
  echo \#
  echo \# libmatheval library not installed or not in library path
  echo \# You can get it from
  echo \#
  echo \# https://www.gnu.org/software/libmatheval/
  echo \#
  echo \# or from
  echo \#
  echo \# https://github.com/mjcarley/libmatheval-no-guile
  echo \#
  echo \#################################################################
  exit
])
AC_CHECK_LIB(triangle-api, triangle_mesh_create, [], [
  echo \#################################################################
  echo \#
  echo \# triangle-api library not installed or not in library path
  echo \# You can get it from
  echo \#
  echo \# https://github.com/wo80/Triangle/
  echo \#
  echo \#################################################################
  exit
])

AC_CONFIG_MACRO_DIR([m4])

# Checks for library functions.
AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS

CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $GLIB_CFLAGS $blaswrap_cflags"
LIBS="$LIBS -ltriangle-api $blaswrap_libs $GLIB_LIBS -lstdc++"

AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(LIBS)

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/agg-config
		 tools/Makefile])
AC_OUTPUT
