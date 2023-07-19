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
	bernstein.c \
	binomials.h \
	binomials.c \
	tinyexpr.c \
	tinyexpr.h \
	section.c \
	transform.c \
	variables.c \
	surface.c \
	naca.c \
	intersection.c \
	hefsi.c \
	hefsi.h \
	patch.c \
	moller-tri.c \
	mesh.c \
	body.c \
	ico.c \
	blend.c

include_HEADERS = \
	agg.h

CLEANFILES = $(BUILT_SOURCES)

