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

gint main(gint argc, gchar **argv)

{
#if 0
  gint fid = 0, i, j, k, npts, ntri ;
  agg_parser_t *p ;
  agg_body_t *b ;
  agg_crowd_t c ;
  agg_distribution_t *d ;
  agg_workspace_t *w ;
  agg_adaptive_grid_workspace_t *wg ;
  agg_mesh_t *m ;
  GScanner *scanner ;  
  gchar ch, *progname ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;
  
  while ( (ch = getopt(argc, argv, "T:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    /* case 'T': test = parse_test(optarg) ; break ; */
    }
  }
  w = agg_workspace_alloc(32) ;
  p = agg_parser_alloc() ;
  scanner = agg_scanner_alloc() ;
  g_scanner_input_file(scanner, fid) ;

  c.p = p ;

  agg_parser_crowd_read(scanner, &c) ;

  agg_crowd_list_bodies(stderr, &c) ;

  fprintf(stderr, "%s: %d bod%s in crowd\n",
	  progname,
	  agg_crowd_body_number(&c),
	  (agg_crowd_body_number(&c) == 1 ? "y" : "ies")) ;
  for ( i = 0 ; i < agg_crowd_body_number(&c) ; i ++ ) {
    b = agg_crowd_body(&c, i) ;
    fprintf(stderr, "%s: gridding body %d\n", progname, i) ;
    if ( agg_body_grid(b)->t == AGG_GRID_ADAPTIVE ) {
      for ( j = 0 ; j < agg_body_distribution_number(b) ; j ++ ) {
	d = agg_body_distribution(b,j) ;
	agg_distribution_interpolation_weights(d) ;
      }
      wg = agg_adaptive_grid_workspace_alloc(32768, 8192, 32768) ;  
      agg_adaptive_grid_make(c.b[i]->g, b, w, wg) ;
    }
  }
  
  npts = ntri = 0 ;
  for ( i = 0 ; i < agg_crowd_body_number(&c) ; i ++ ) {
    if ( c.b[i]->g != NULL ) {
      fprintf(stderr, "  body %d: %d grid points; %d triangles\n", i,
	      agg_grid_point_number_max(c.b[i]->g),
	      agg_grid_triangle_number_max(c.b[i]->g)) ;
      npts += agg_grid_point_number_max(c.b[i]->g) ;
      ntri += agg_grid_triangle_number_max(c.b[i]->g) ;
    }
  }

  fprintf(stderr, "   total: %d grid points; %d triangles\n", npts, ntri) ;

  agg_parser_expressions_evaluate(c.p) ;
  
  m = agg_mesh_alloc(npts, ntri, 1, 0, 0) ;
  for ( i = 0 ; i < agg_crowd_body_number(&c) ; i ++ ) {
    fprintf(stderr, "%s: meshing body %d\n", progname, i) ;
    b = agg_crowd_body(&c, i) ;
    for ( j = 0 ; j < agg_body_distribution_number(b) ; j ++ ) {
      d = agg_body_distribution(b,j) ;
      agg_distribution_interpolation_weights(d) ;
    }
    if ( b->g != NULL ) agg_body_mesh_grid(b, b->g, m, i, w) ;
  }

  agg_mesh_trim(m, &c, w) ;
  
  agg_mesh_points_write(stdout, m) ;
  agg_mesh_tri_write(stdout, m) ;
#endif
  
  return 0 ;
}
