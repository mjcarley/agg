## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)

DEPS = $(top_builddir)/src/libagg.la

bin_PROGRAMS = \
	agg-generate \
	agg-test

agg_generate_LDADD = \
	$(top_builddir)/src/.libs/libagg.a 

agg_test_LDADD = \
	$(top_builddir)/src/.libs/libagg.a 
