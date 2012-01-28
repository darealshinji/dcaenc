#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([dcaenc],[2],[patrakov@gmail.com],,[http://aepatrakov.narod.ru/dcaenc/])
AC_CONFIG_SRCDIR([dcaenc.c])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign check-news])
LT_INIT([disable-static])

# Checks for programs.
AC_PROG_CC
PKG_PROG_PKG_CONFIG 

# Checks for libraries.
AC_ARG_ENABLE([alsa],
      AS_HELP_STRING([--disable-alsa], [disable building of ALSA plugin]))

ac_save_LIBS="$LIBS"

if test "x$enable_alsa" != "xno"; then
	PKG_CHECK_MODULES([ALSA], [alsa >= 1.0.11], [HAVE_ALSA=yes], [HAVE_ALSA=no])
	AC_CHECK_LIB(asound, snd_pcm_extplug_create,,
             [HAVE_ALSA=no], -ldl)
	if test "x$enable_alsa" = "xyes" -a "x$HAVE_ALSA" = "xno"; then
		AC_MSG_ERROR([ALSA libraries not found])
	fi
fi
AM_CONDITIONAL(HAVE_ALSA, test "x$enable_alsa" != "xno" -a "x$HAVE_ALSA" = "xyes")
LIBS="$ac_save_LIBS"

# Checks for header files.
AC_HEADER_STDC

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
ac_save_LIBS="$LIBS"
LT_LIB_M
LIBS="$ac_save_LIBS"

AC_CONFIG_FILES([Makefile dcaenc.pc])
AC_OUTPUT
