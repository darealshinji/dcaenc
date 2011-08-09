ACLOCAL_AMFLAGS = -I m4

AUTOMAKE_OPTIONS=foreign
EXTRA_DIST = dcaenc.sym

bin_PROGRAMS = dcaenc
noinst_PROGRAMS = gentables
lib_LTLIBRARIES = libdcaenc.la

plugin_LTLIBRARIES = libasound_module_pcm_dca.la
plugindir = $(libdir)/alsa-lib

include_HEADERS = dcaenc.h
libdcaenc_la_SOURCES = dcaenc.c \
    dcaenc_private.h int_data.h math_tables.h softfloat.h
nodist_libdcaenc_la_SOURCES = math_tables.c
libdcaenc_la_LDFLAGS = -version-info 0:0:0 \
    -export-symbols $(top_srcdir)/dcaenc.sym

libasound_module_pcm_dca_la_SOURCES = alsaplugin.c dcaenc.h
libasound_module_pcm_dca_la_LIBADD = libdcaenc.la -lasound
libasound_module_pcm_dca_la_LDFLAGS = -no-undefined -avoid-version

dcaenc_SOURCES = main.c wavfile.c dcaenc.h wavfile.h
dcaenc_LDADD = libdcaenc.la

gentables_SOURCES = gentables.c \
    int_data.h float_data.h softfloat.h

gentables_LDADD = @LIBM@

CLEANFILES = math_tables.c
math_tables.c: Makefile gentables$(EXEEXT)
	./gentables$(EXEEXT) >$@
