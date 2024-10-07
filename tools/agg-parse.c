/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2023 Michael Carley
 *
 * AGG is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  AGG is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with AGG.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <agg.h>

char *progname ;

static void print_help_message(FILE *f, gint pps)

{
  fprintf(f, "Usage: %s [options] input\n\n", progname) ;

  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit;\n"
	  "  -e echo processing of file to stderr\n"
	  "  -G # GMSH .geo output file name\n"
	  "  -p # number of points per spline in parsed geometry mesh (%d)\n"
	  "  -S list available sections\n"
	  "  -s # parse and write a section to stdout\n"
	  "  -T list available transforms\n"
	  "  -t # apply a transform to a section and write to stdout\n",
	  pps) ;
  
  return ;
}  

static agg_transform_t *transform_parse(char *str)

{
  agg_transform_t *T ;
  char **tokens, *name ;
  agg_variable_t p[32] = {0} ;
  gint i, np ;

  g_assert_not_reached() ;
  
  T = agg_transform_new(8) ;

  g_strdelimit (str, "(),", ' ') ;

  tokens = g_strsplit(str, " ", 0) ;

  name = NULL ; np = 0 ;
  for ( i = 0 ; tokens[i] != NULL ; i ++ ) {
    if ( strlen(tokens[i]) != 0 ) {
      if ( name == NULL ) {
	name = g_strdup(tokens[i]) ;
      } else {
	p[np].val = atof(tokens[i]) ; np ++ ;
      }
    }
  }

  /* agg_transform_parse(T, name, p, np) ; */
  agg_transform_parse(T, p, np) ;
  
  return T ;
}

static void section_write(FILE *f, char *str, agg_transform_t *T,
			  gint npts, char *opfmt)

{
  agg_section_t *s ;
  char **tokens, *name ;
  agg_variable_t p[32] = {0} ;
  gint i, np ;
  
  s = agg_section_new(128, 128) ;

  g_strdelimit (str, "(),", ' ') ;

  tokens = g_strsplit(str, " ", 0) ;

  name = NULL ; np = 0 ;
  for ( i = 0 ; tokens[i] != NULL ; i ++ ) {
    if ( strlen(tokens[i]) != 0 ) {
      if ( name == NULL ) {
	name = g_strdup(tokens[i]) ;
      } else {
	p[np].val = atof(tokens[i]) ; np ++ ;
      }
    }
  }

  agg_section_parse(s, name, p, np) ;

  if ( opfmt == NULL ) {
    agg_section_write(f, s, T, npts) ;
    return ;
  }

  if ( strcmp(opfmt, "mpost") == 0 ) {
    agg_section_format_write(f, s, T,
			     "(%1.3lgu,%1.3lgu)--",
			     "(%1.3lgu,%1.3lgu) ;\n", npts) ;
    
    return ;
  }

  fprintf(stderr, "%s: unrecognized output format %s\n", progname, opfmt) ;
  
  return ;
}

static void body_write_points(FILE *f, agg_body_t *b,
			      agg_surface_workspace_t *w)

{
  gint i, j, k, nu, nv ;
  gdouble u, v, x[3] ;
  agg_surface_t *S ;

  nu = 32 ; nv = 128 ;
  for ( i = 0 ; i < agg_body_surface_number(b) ; i ++ ) {
    S = agg_body_surface(b, i) ;

    for ( j = 0 ; j <= nu ; j ++ ) {
      u = agg_surface_umin(S) +
	(agg_surface_umax(S) - agg_surface_umin(S))*j/nu ;
      for ( k = 0 ; k <= nv ; k ++ ) {
	v = -1 + 2.0*k/nv ;
	agg_surface_point_eval(S, u, v, x, w) ;
	fprintf(f, "%e %e %e\n", x[0], x[1], x[2]) ;
      }
    }
    
  }
	  
  return ;
}

gint main(gint argc, char **argv)

{
  agg_body_t *b ;
  agg_mesh_t *m ;
  char *file, *gfile, ch, *section, *transform, *opfmt ;
  gint pps, offp, offsp, offs, i ;
  agg_surface_workspace_t *w ;
  agg_transform_t *T ;
  FILE *output ;
  gboolean echo, write_points ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;

  echo = FALSE ; write_points = FALSE ;

  pps = 2 ;
  offp = offsp = offs = 1 ;
  gfile = NULL ;
  section = NULL ;
  transform = NULL ;
  opfmt = NULL ;
  T = NULL ;
  
  while ( (ch = getopt(argc, argv, "heG:o:Pp:Ss:Tt:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'h': print_help_message(stderr, pps) ; return 0 ; break ;
    case 'e': echo = TRUE ; break ;
    case 'G': gfile = g_strdup(optarg) ; break ;
    case 'o': opfmt = g_strdup(optarg) ; break ;
    case 'P': write_points = TRUE ; break ;
    case 'p': pps = atoi(optarg) ; break ;
    case 's': section = g_strdup(optarg) ; break ;
    case 'S':
      fprintf(stderr, "%s: available sections\n\n", progname) ;
      agg_sections_list(stderr) ;
      return 0 ;
    case 'T':
      fprintf(stderr, "%s: available transforms\n\n", progname) ;
      agg_affine_list(stderr, "  ", "\n") ;
      return 0 ;
      break ;
    case 't': transform = g_strdup(optarg) ; break ;
    }
  }

  if ( transform != NULL ) {
    T = transform_parse(transform) ;
  }
  
  if ( section != NULL ) {
    section_write(stdout, section, T, 128, opfmt) ;
		  
    return 0 ;
  }
  
  if ( optind >= argc ) {
    fprintf(stderr, "%s: input file name required\n", progname) ;
    exit(1) ;
  }

  file = g_strdup(argv[optind]) ;

  fprintf(stderr, "%s: parsing \"%s\"\n", progname, file) ;

  b = agg_body_new(32, 32) ;
  w = agg_surface_workspace_new() ;
  agg_body_read(b, file, echo) ;
  agg_body_globals_compile(b) ;
  agg_body_globals_eval(b) ;
  agg_body_globals_write(stderr, b) ;
  if ( echo ) {
    fprintf(stderr, "%d surface%s\n", agg_body_surface_number(b),
	    (agg_body_surface_number(b) > 1 ? "s" : "")) ;
    agg_body_surfaces_list(stderr, b) ;
  }

  if ( write_points ) {
    body_write_points(stdout, b, w) ;
    
    return 0 ;
  }
  
  m = agg_mesh_new(65536, 65536, 65536) ;
  /* agg_mesh_body(m, b, pps, w) ; */
  for ( i = 0 ; i < agg_body_surface_number(b) ; i ++ ) {
    agg_mesh_surface(m,i) = agg_body_surface(b,i) ;
    agg_mesh_surface_number(m) ++ ;
    agg_mesh_surface_triangulate(m, i,
				 agg_body_triangulation_settings(b,i), w) ;
  }

  if ( gfile != NULL ) {
    output = fopen(gfile, "w") ;
    if ( output == NULL ) {
      fprintf(stderr, "%s: cannot open GMSH file \"%s\"\n", progname, gfile) ;
      return 1 ;
    }

    fprintf(output, "lc = 0.1 ;\n") ;
    agg_mesh_write_gmsh(output, m, "lc", offp, offsp, offs, FALSE) ;

    fclose(output) ;
  }
  
  return 0 ;
}
