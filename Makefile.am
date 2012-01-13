ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTIONS=foreign
EXTRA_DIST = dcaenc.sym dca.conf

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dcaenc.pc

bin_PROGRAMS = dcaenc
noinst_PROGRAMS = gentables
lib_LTLIBRARIES = libdcaenc.la

include_HEADERS = dcaenc.h
libdcaenc_la_SOURCES = dcaenc.c \
    dcaenc_private.h int_data.h math_tables.h softfloat.h
nodist_libdcaenc_la_SOURCES = math_tables.c
libdcaenc_la_LDFLAGS = -version-info 0:0:0 \
    -export-symbols $(top_srcdir)/dcaenc.sym

if HAVE_ALSA
plugin_LTLIBRARIES = libasound_module_pcm_dca.la
plugindir = $(libdir)/alsa-lib

libasound_module_pcm_dca_la_SOURCES = alsaplugin.c dcaenc.h
libasound_module_pcm_dca_la_LIBADD = libdcaenc.la @ALSA_LIBS@
libasound_module_pcm_dca_la_LDFLAGS = -no-undefined -avoid-version

alsaconfdir = $(datadir)/alsa/pcm
alsaconf_DATA = dca.conf
endif

dcaenc_SOURCES = main.c wavfile.c dcaenc.h wavfile.h
dcaenc_LDADD = libdcaenc.la

gentables_SOURCES = gentables.c \
    int_data.h float_data.h softfloat.h

gentables_LDADD = @LIBM@

CLEANFILES = math_tables.c
math_tables.c: Makefile gentables$(EXEEXT)
	./gentables$(EXEEXT) >$@
