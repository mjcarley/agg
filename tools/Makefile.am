## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
DEPS = $(top_builddir)/src/libagg.la
LDADD = $(top_builddir)/src/.libs/libagg.a 


bin_PROGRAMS = \
	agg-test \
	agg-parse
