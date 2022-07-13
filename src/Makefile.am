AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"AGG\"

bin_SCRIPTS=agg-config

BUILT_SOURCES = \
	agg-config

agg-config: agg-config.in

lib_LTLIBRARIES = libagg.la

libagg_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libagg_la_SOURCES = \
	agg.c \
	binomials.h \
	binomials.c \
	naca.c \
	agg-private.h \
	transform.c \
	distribution.c \
	mesh.c \
	spacing.c \
	tinyexpr.c \
	tinyexpr.h \
	parser.c \
	body.c \
	shape.c \
	interp.c \
	functions.c

include_HEADERS = \
	agg.h

CLEANFILES = $(BUILT_SOURCES)

