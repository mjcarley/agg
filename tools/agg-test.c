/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2021, 2022 Michael Carley
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

const gchar *tests[] = {"bernstein",
			"geometry",
			"transform",
			"wing",
			"parser",
			""} ;

static gint parse_test(gchar *s)

{
  gint i ;

  for ( i = 0 ; strlen(tests[i]) != 0 ; i ++ ) {
    if ( strcmp(s, tests[i]) == 0 ) return i ;
  }
  
  return -1 ;
}

static gint bernstein_basis_test(gint n)

{
  gdouble x, tot, err, S[64] ;
  gint j ;
  
  fprintf(stderr, "Bernstein basis polynomials\n") ;
  fprintf(stderr, "===========================\n") ;
  fprintf(stderr, "polynomial basis order: %d\n", n) ;

  err = 0.0 ;
  for ( x = 0.0 ; x <= 1.0 ; x += 1.0/128 ) {
    tot = 0.0 ;
    for ( j = 0 ; j <= n ; j ++ ) {
      tot += agg_bernstein_basis_eval(n, j, x) ;
    }
    err = MAX(err, fabs(tot-1.0)) ;
  }

  fprintf(stderr, "termwise maximum deviation from unity: %lg\n", err) ;

  err = 0.0 ;
  for ( x = 0.0 ; x <= 1.0 ; x += 1.0/128 ) {
    tot = 0.0 ;
    agg_bernstein_basis(n, x, S, NULL) ;
    for ( j = 0 ; j <= n ; j ++ ) tot += S[j] ;
    err = MAX(err, fabs(tot-1.0)) ;
  }

  fprintf(stderr, "vector maximum deviation from unity: %lg\n", err) ;
  
  return 0 ;
}

static void geometry_test(void)

{
  agg_shape_t *c ;
  gdouble x, z, y, th, t, p, m, xc[1024], yc[1024], *work ;
  gint nu, nl, n ;
  
  /*NACA section*/
  th = 0.12 ;
  p = 0.1 ; m = 0.4 ; 
  /* p = 0.0 ; m = 0.0 ; */
  
  fprintf(stderr, "geometry evaluation\n") ;
  fprintf(stderr, "===================\n") ;

  n = 8 ;
  c = agg_shape_alloc(n) ;

  nu = nl = 0 ;
  for ( t = 0.5*M_PI ; t >= 0.0 ; t -= 0.5*M_PI/64 ) {
    xc[nu] = cos(t) ;
    yc[nu] = agg_naca_four(th, p, m, xc[nu]) ;
    nu ++ ;
  }

  for ( t = M_PI ; t >= 0.5*M_PI ; t -= 0.5*M_PI/64 ) {
    xc[nu+nl] = -cos(t) ;
    yc[nu+nl] = agg_naca_four(th, p, m, -xc[nu+nl]) ;
    nl ++ ;
  }

  work = (gdouble *)g_malloc0((2*(nu+nl)+2*(nu+nl+1)*(n+1))*
			      sizeof(gdouble)) ;
  
  agg_shape_fit(c, &(xc[0]), &(yc[0]), nu,
		&(xc[nu]), &(yc[nu]), nl,
		0.5, 1.0, yc[nu-1], n, TRUE, work) ;

  for ( t = M_PI ; t > 0.5*M_PI ; t -= 0.5*M_PI/256 ) {
    x = cos(t) ;
    z = agg_shape_eval(c, x, 0) ;
    y = agg_naca_four(th, p, m, x) ;
    fprintf(stdout, "%lg %lg %lg\n", -x, z, y) ;
  }

  for ( t = 0.5*M_PI ; t >= 0.0 ; t -= 0.5*M_PI/256 ) {
    x = cos(t) ;
    z = agg_shape_eval(c, x, 0) ;
    y = agg_naca_four(th, p, m, x) ;
    fprintf(stdout, "%lg %lg %lg\n", x, z, y) ;
  }

  return ;
}

static void transform_test(void)

{
  agg_shape_t *s ;
  agg_local_transform_t *tr ;
  gdouble x[3], th, t, p, m, *work, r[3] ;
  agg_parser_t *parser ;
  gint nu, nl, n ;
  
  /*NACA section*/
  th = 0.12 ;
  p = 0.1 ; m = 0.4 ; 
  p = 0.0 ; m = 0.0 ;

  fprintf(stderr, "local transform evaluation\n") ;
  fprintf(stderr, "==========================\n") ;

  n = 8 ;
  s = agg_shape_alloc(n) ;
  tr = agg_local_transform_alloc() ;
  parser = agg_parser_alloc() ;

  r[0] = 0.1 ; r[1] = 0.0 ; r[2] = 0.7 ;
  agg_local_transform_parse(tr, "shrink", r, 3) ;
  r[0] = 0.1 ; r[1] = 0.1 ; r[2] = 5.5 ;
  agg_local_transform_parse(tr, "shift", r, 3) ;

  agg_parser_variable_add(parser, "dx", 0.3) ; 
  /*check resetting of existing variables works correctly*/
  agg_parser_variable_add(parser, "dx", 0.2) ;
  agg_parser_variable_add(parser, "dy", -0.1) ;
  agg_local_transform_set_expression(tr,
				     agg_local_transform_parameter_index(tr,1,0),
				     "dx+dy", parser) ;  
  agg_local_transform_set_expression(tr,
				     agg_local_transform_parameter_index(tr,1,1),
				     "dy", parser) ;  
  
  nu = 16 ; nl = 26 ;

  work = (gdouble *)g_malloc0((4*(nu+nl)+2*(nu+nl+1)*(n+1))*
			      sizeof(gdouble)) ;

  agg_shape_fit_naca_four(s, n, th, p, m, nu, nl, work) ;

  for ( t = M_PI ; t > 0.5*M_PI ; t -= 0.5*M_PI/256 ) {
    x[0] = -cos(t) ;
    x[1] = agg_shape_eval(s, cos(t), 0) ;
    x[2] = 0.0 ;
    fprintf(stdout, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
  }

  for ( t = 0.5*M_PI ; t >= 0.0 ; t -= 0.5*M_PI/256 ) {
    x[0] = cos(t) ;
    x[1] = agg_shape_eval(s, cos(t), 0) ;
    x[2] = 0.0 ;
    fprintf(stdout, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
  }

  agg_local_transform_eval_parameters(tr) ;
  for ( t = M_PI ; t > 0.5*M_PI ; t -= 0.5*M_PI/256 ) {
    x[0] = -cos(t) ;
    x[1] = agg_shape_eval(s, cos(t), 0) ;
    x[2] = 0.0 ;
    agg_local_transform_apply(tr, x) ;
    fprintf(stdout, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
  }

  for ( t = 0.5*M_PI ; t >= 0.0 ; t -= 0.5*M_PI/256 ) {
    x[0] = cos(t) ;
    x[1] = agg_shape_eval(s, cos(t), 0) ;
    x[2] = 0.0 ;
    agg_local_transform_apply(tr, x) ;
    fprintf(stdout, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
  }
  
  return ;
}

static void wing_test(void)

{
  agg_distribution_t *wing ;
  agg_shape_t *sh ;
  agg_local_transform_t *T ;
  gdouble ttip, troot, *work, p, m, th, tw, t, span, pm[8] ;
  gdouble ctip, croot ;
  gint i, nt, n, nu, nl ;
  
  wing = agg_distribution_alloc(32) ;
  wing->smin = 0 ;
  wing->smax = 1 ;

  /*number of nodes per section on output*/
  /* nj = 128 ; */
  
  ttip = 0.05 ; troot = 0.2 ; nt = 4 ;
  croot = 0.8 ; ctip = 0.1 ;
  tw = 10*M_PI/180 ;
  span = 2.5 ;
  p = 0.0 ; m = 0.2 ;
  n = 8 ; nu = 16 ; nl = 16 ;
  
  work = (gdouble *)g_malloc0((4*(nu+nl)+2*(nu+nl+1)*(n+1))*
			      sizeof(gdouble)) ;
  for ( i = 0 ; i < nt ; i ++ ) {
    t = (gdouble)i/(nt-1) ;
    th = troot + (ttip - troot)*t ;
    sh = agg_shape_alloc(n) ;
    T = agg_local_transform_alloc() ;
    /* agg_transform_scale(T, croot + t*(ctip - croot)) ; */
    pm[0] = croot + t*(ctip - croot) ;
    agg_local_transform_parse(T, "scale", pm, 1) ;
    pm[0] = 0 ; pm[1] = 0 ; pm[2] = t*span ;
    agg_local_transform_parse(T, "shift", pm, 3) ;

    agg_shape_fit_naca_four(sh, n, th, p, m, nu, nl, work) ;
    agg_distribution_add_shape(wing, t, sh) ;
  }

  agg_distribution_interpolation_weights(wing) ;

  sh = agg_shape_alloc(n) ;
  T = agg_local_transform_alloc() ;
  /* agg_distribution_write_mesh(stdout, wing, */
  /* 			      0, 1, 15, AGG_SPACING_COSINE, */
  /* 			      -1, 1, 64, AGG_SPACING_COSINE, */
  /* 			      sh, T, NULL) ; */
  
  return ;
}

static void parser_test(void)

{
  gint fid = 0, i ;
  agg_parser_t *p ;
  agg_body_t *b ;
  agg_distribution_t *d ;
  agg_workspace_t *w ;
  agg_mesh_t *m ;
  GScanner *scanner ;
  
  fprintf(stderr, "geometry parser test\n") ;
  fprintf(stderr, "====================\n") ;

  w = agg_workspace_alloc(32) ;
  p = agg_parser_alloc() ;
  b = agg_body_alloc() ;
  scanner = agg_scanner_alloc() ;
  g_scanner_input_file(scanner, fid) ;

  if ( agg_parser_body_read(fid, scanner, p, b) != 0 ) {
    exit(1) ;
  }
  agg_parser_expressions_evaluate(p) ;

  /*generate distributions*/
  agg_body_distributions_list(stderr, b) ;

  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) {
    d = agg_body_distribution(b,i) ;
    agg_distribution_interpolation_weights(d) ;
  }

  m = agg_mesh_alloc(agg_grid_point_number_max(b->g),
		     agg_grid_triangle_number_max(b->g),
		     0, 0, 0) ;
  
  agg_body_mesh_grid(b, b->g, m, w) ;
  
  agg_mesh_points_write(stdout, m) ;
  agg_mesh_tri_write(stdout, m) ;
  
  return ; 
}

gint main(gint argc, gchar **argv)

{
  gchar ch, *progname ;
  gint test ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;
  test = -1 ;
  
  while ( (ch = getopt(argc, argv, "T:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'T': test = parse_test(optarg) ; break ;
    }
  }

  if ( test == -1 ) {
    fprintf(stderr, "%s: unrecognized test\n", progname) ;
    return 1 ;
  }

  if ( test == 0 ) {
    bernstein_basis_test(8) ;

    return 0 ;
  }

  if ( test == 1 ) {
    geometry_test() ;

    return 0 ;
  }

  if ( test == 2 ) {
    transform_test() ;

    return 0 ;
  }

  if ( test == 3 ) {
    wing_test() ;

    return 0 ;
  }

  if ( test == 4 ) {
    parser_test() ;

    return 0 ;
  }
  
  return 0 ;
}
