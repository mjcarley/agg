#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

##AC_PREREQ([2.69])
AC_INIT([agg],[0.1.0],[codes@paraffinalia.co.uk])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/agg.c])
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
AC_HEADER_STDC

# Checks for libraries.
AM_PATH_GLIB_2_0(2.4.0, has_glib2=yes, has_glib2=no, gthread gmodule)
if test x$has_glib2 = xyes; then
  glib_cflags=`$PKG_CONFIG glib-2.0 --cflags`
  glib_thread_cflags=`$PKG_CONFIG glib-2.0 --cflags gthread-2.0`
  glib_module_cflags=`$PKG_CONFIG glib-2.0 --cflags gmodule-2.0`
  glib_libs=`$PKG_CONFIG glib-2.0 --libs`
  glib_thread_libs=`$PKG_CONFIG glib-2.0 --libs gthread-2.0`
  glib_module_libs=`$PKG_CONFIG glib-2.0 --libs gmodule-2.0`
  GLIB_LIBS="$glib_libs"
  GLIB_DEPLIBS="$glib_libs"
else
  if test x$with_glib = xyes ; then
    AC_MSG_ERROR([
*** Directory must be specified for --with-glib])
  fi

  if test x$with_glib = x ; then 
    # Look for separately installed glib

    GLIB_REQUIRED_VERSION=1.2.8
    AM_PATH_GLIB($GLIB_REQUIRED_VERSION,,
      AC_MSG_ERROR([
*** GLIB $GLIB_REQUIRED_VERSION or better is required. The latest version of GLIB
*** is always available from ftp://ftp.gtk.org/.]),
      gmodule gthread)

    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`
    GLIB_LIBS="$glib_libs"
    GLIB_DEPLIBS="$glib_libs"
  else
    # Use uninstalled glib (assume they got the version right)

    GLIB_CONFIG=$with_glib/glib-config
    if test -x $GLIB_CONFIG ; then 
      :
    else
      AC_MSG_ERROR([GLIB directory ($with_glib) not present or not configured])
    fi

    glib_cflags=`$GLIB_CONFIG --cflags`
    glib_thread_cflags=`$GLIB_CONFIG --cflags gthread`
    glib_module_cflags=`$GLIB_CONFIG --cflags gmodule`
    glib_libs=`$GLIB_CONFIG --libs`
    glib_thread_libs=`$GLIB_CONFIG --libs gthread`
    glib_module_libs=`$GLIB_CONFIG --libs gmodule`

    glib_release=`$GLIB_CONFIG --version | sed 's%\\.[[0-9]]*$%%'`

    # canonicalize relative paths
    case $with_glib in 
      /*)
        glib_dir=$with_glib
        ;;
      *)
        glib_dir="\$(top_builddir)/$with_glib"
        ;;
    esac

    GLIB_CFLAGS="-I$glib_dir -I$glib_dir/gmodule"
    GLIB_LIBS="$glib_dir/libglib.la $glib_dir/gmodule.la"
    GLIB_DEPLIBS=

    AC_SUBST(GLIB_CFLAGS)
    AC_SUBST(GLIB_LIBS)
  fi
fi

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

AC_CONFIG_MACRO_DIR([m4])

# Checks for library functions.
AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS

CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $glib_cflags $blaswrap_cflags"
LDFLAGS="$LDFLAGS $blaswrap_libs $glib_libs $gqr_libs $mop_libs $gtv_libs -lstdc++"

AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
#AC_SUBST(DX_RULES)

LIBS="$glib_libs $LIBS $FLIBS"

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 src/agg-config
		 tools/Makefile])
AC_OUTPUT
