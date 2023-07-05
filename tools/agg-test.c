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

#include "hefsi.h"

const gchar *tests[] = {"bernstein",
			"circle",
			"aerofoil",
			"transform",
			"surface",
			"body",
			""} ;

static void surface_sphere(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			   gdouble r, agg_section_t *s)

{
  gdouble u, p[3] ;
  agg_transform_t *T ;
  gchar *expr[8] = {NULL} ;

  u = 0.0 ;
  agg_section_set_circle(s) ;
  T = agg_surface_transform(S) ;
  agg_surface_section_add(S, s, u) ;
  /*centre on (0,0): circle of radius 1/2*/
  p[0] = -0.5 ; p[1] = 0 ; p[2] = 0 ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;
  /*scale circle radius along z axis: circle of radius 1/2 at u=1/2 */
  p[0] =  0.0 ; p[1] = 0.0 ;
  expr[2] = "2*sqrt(u*(1-u))" ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, p, expr, 3) ;
  /*shift centre to z=0*/
  p[0] =  0.0 ; p[1] = 0 ; p[2] = 0 ;
  expr[2] = "-1+u" ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;
  /*scale to required radius: sphere of radius r*/
  p[0] = 2*r ;
  /*shift to required centre*/
  agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, p, expr, 1) ;
  p[0] = x ; p[1] = y ; p[2] = z ;
  expr[2] = NULL ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;

  agg_surface_umin(S) = 0.0 ; 
  agg_surface_umax(S) = 1.0 ; 

  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  agg_surface_weights_make(S) ;
  
  return ;
}

static void surface_tube(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			 gdouble r, gdouble len, agg_section_t *s)

{
  gdouble u, p[3] ;
  agg_transform_t *T ;
  gchar *expr[8] = {NULL}, etmp[1024] ;

  u = 0.0 ;
  agg_section_set_circle(s) ;
  T = agg_surface_transform(S) ;
  agg_surface_section_add(S, s, u) ;
  /*centre on (0,0): circle of radius 1/2*/
  p[0] = -0.5 ; p[1] = 0 ; p[2] = 0 ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;
  /*scale on radius*/
  p[0] = 2*r ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, p, expr, 1) ;
  /*translate in space to form cylinder*/
  p[0] = 0 ; p[1] = 0 ; p[2] = 0 ;
  sprintf(etmp, "-%lg/2 + %lg*u", len, len) ;
  expr[2] = etmp ;
  /* expr[2] = "-len/2 + len*u" ; */
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;

  /*translate to appropriate centre*/
  p[0] = x ; p[1] = y ; p[2] = z ;
  expr[2] = NULL ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;

  agg_surface_axes(S) = AGG_AXES_PZ_PX_PY ;
  
  agg_surface_umin(S) = 0.0 ; 
  agg_surface_umax(S) = 1.0 ; 

  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  agg_surface_weights_make(S) ;
  
  return ;
}

static void surface_wing(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			 gdouble root, gdouble span,
			 gdouble lm, gdouble swle, gdouble dh,
			 agg_section_t *s)

{
  gdouble u, p[3] ;
  agg_transform_t *T ;
  gchar *expr[8] = {NULL}, etmp[1024] ;

  u = 0.0 ;
  T = agg_surface_transform(S) ;
  agg_surface_section_add(S, s, u) ;
  
  expr[0] = expr[1] = expr[2] = NULL ; 
  /*scale on taper ratio, shrinking about leading edge*/
  p[0] = 0.0 ; p[1] = 0 ; p[2] = 0 ;
  sprintf(etmp, "1 + u*(%lg-1)", lm) ;
  expr[2] = g_strdup(etmp) ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, p, expr, 3) ;
  expr[0] = expr[1] = expr[2] = NULL ; 
  /*scale on root chord*/
  p[0] = root ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, p, expr, 1) ;
  expr[0] = expr[1] = expr[2] = NULL ; 

  /*translate along span and sweep*/
  p[0] = 0 ; p[1] = 0 ; p[2] = 0 ;
  sprintf(etmp, "%lg*sin(%lg)*u", span, swle) ;
  expr[0] = g_strdup(etmp) ;
  sprintf(etmp, "%lg*sin(%lg)*u", span, dh) ;
  expr[1] = g_strdup(etmp) ;
  sprintf(etmp, "%lg*u", span) ;
  expr[2] = g_strdup(etmp) ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;
  expr[0] = expr[1] = expr[2] = NULL ; 
  p[0] = x ; p[1] = y ; p[2] = z ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ;
  expr[0] = expr[1] = expr[2] = NULL ; 

  /* agg_surface_axes(S) = AGG_AXES_PZ_PY_PX ; */

  agg_surface_umin(S) = 0.0 ; 
  agg_surface_umax(S) = 1.0 ; 

  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  agg_surface_weights_make(S) ;
  
  return ;
}

static gint parse_test(gchar *s)

{
  gint i ;

  for ( i = 0 ; strlen(tests[i]) != 0 ; i ++ ) {
    if ( strcmp(s, tests[i]) == 0 ) return i ;
  }
  
  return -1 ;
}

static void bernstein_basis_test(gint n)

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
  
  return ;
}

static void circle_test(void)

{
  gint i, nx ;
  gdouble x, y, r, err ;
  agg_section_t *s ;

  s = agg_section_new(8, 8) ;
  agg_section_set_circle(s) ;
  
  nx = 128 ;
  err = 0.0 ;
  for ( i = 0 ; i <= nx ; i ++ ) {
    x = -1.0 + 2.0*i/nx ;
    y = agg_section_eval(s, x) ;
    fprintf(stdout, "%e %e\n", fabs(x), y) ;
    r = sqrt((fabs(x)-0.5)*(fabs(x)-0.5) + y*y) ;
    err = MAX(err, fabs(r-0.5)) ;
  }

  fprintf(stderr, "radius error: %lg\n", err) ;
  
  return ;
}

static void aerofoil_test(void)

{
  gint i, nx ;
  gdouble x, y ;
  agg_section_t *s ;

  s = agg_section_new(8, 8) ;
  agg_section_set_aerofoil(s, 0.5, 0.1, 0.01) ;
  
  nx = 128 ;
  for ( i = 0 ; i <= nx ; i ++ ) {
    x = -1.0 + 2.0*i/nx ;
    y = agg_section_eval(s, x) ;
    fprintf(stdout, "%e %e\n", fabs(x), y) ;
  }
  
  return ;
}

static void transform_test(void)

{
  gint i, j, nx, nu, nglobal ;
  gdouble x, y, xin[3], xout[3], p[32] ;
  agg_section_t *s ;
  agg_transform_t *T ;
  agg_variable_t global[4] ;
  gchar *args[32] ;

  s = agg_section_new(8, 8) ;
  agg_section_set_aerofoil(s, 0.5, 0.1, 0.0) ;

  T = agg_transform_new(64) ;

  /*
   * set up global variables: 
   * u is surface parameter for distribution
   * twist is total twist at wingtip
   */
  nglobal = 2 ;
  global[0].name = g_strdup("u") ;
  global[0].def =  NULL ;
  global[0].val = 0.0 ;
  global[1].name = g_strdup("twist") ;
  global[1].def =  g_strdup("15.0*pi/180.0") ;
  global[1].val = 0.4 ;

  /*non-linear scaling, shrinking about mid-chord*/
  args[0] = NULL ; args[1] = "0" ; args[2] = "0.5*u+(1-u)^2" ;
  p[0] = 0.5 ; p[1] = 0 ; p[2] = 0 ;
  /*linear twist*/
  agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, p, args, 3) ;
  args[0] = NULL ; args[1] = "0" ; args[2] = "twist*u" ;
  p[0] = 0.5 ; p[1] = 0 ; p[2] = 0 ;
  /*translate in z to generate "wing"*/
  agg_transform_operator_add(T, AGG_TRANSFORM_ROTATE, p, args, 3) ;
  args[0] = NULL ; args[1] = NULL ; args[2] = "u*3.0" ;
  p[0] = 0.0 ; p[1] = 0 ; p[2] = 0 ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, args, 3) ;

  /* agg_transform_operators_write(stderr, T) ; */
  
  agg_transform_add_global_variables(T, global, nglobal) ;
  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;
  agg_transform_variables_eval(T) ;

  nx = 65 ; nu = 17 ;
  for ( j = 0 ; j < nu ; j ++ ) {
    global[0].val = (gdouble)j/(nu-1) ;
    /*reevaluate all variables after changing the parameters*/
    agg_transform_variables_eval(T) ;
    for ( i = 0 ; i <= nx ; i ++ ) {
      x = -1.0 + 2.0*i/nx ;
      y = agg_section_eval(s, x) ;
      /* fprintf(stdout, "%e %e\n", fabs(x), y) ; */
      xin[0] = fabs(x) ; xin[1] = y ; xin[2] = xout[2] = 0 ;
      agg_transform_apply(T, xin, xout) ;
      fprintf(stdout, "%e %e %e\n", xout[0], xout[1], xout[2]) ;
    }
  }
  
  return ;
}

static void surface_test(void)

{
  gint nu, nv, i, j ;
  /* agg_transform_t *T ; */
  agg_surface_t *S ;
  agg_surface_workspace_t *w ;
  agg_section_t *s ;
  agg_patch_t *P ;
  gdouble u, v, x[3] ;
  
  S = agg_surface_new(64) ;
  P = agg_patch_new(1024) ;
  w = agg_surface_workspace_new() ;
  s = agg_section_new(32, 32) ;
  agg_section_set_circle(s) ;
  
  /* T = agg_surface_transform(S) ; */

  /*weird twisty thing*/
  /* p[0] = -0.5 ; p[1] = 0 ; p[2] = 0 ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ; */
  /* p[0] =  0.0 ; p[1] = 0.0 ; */
  /* expr[2] = "u^2*90*pi/180" ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_ROTATE, p, expr, 3) ; */
  /* p[0] =  0.0 ; p[1] = 0 ; p[2] = 0 ; */
  /* expr[2] = "3*u + abs(v)" ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, p, expr, 3) ; */

  /* for ( i = 0 ; i < 8 ; i ++ ) { */
  /*   u = (gdouble)i/7 ; */
  /*   agg_section_eta_left(s)  = 0.5 + (0.05-0.5)*u ; */
  /*   agg_section_eta_right(s) = 0.5 + (0.25-0.5)*u ; */
  /*   agg_surface_section_add(S, s, u) ; */
  /* } */

  /*sphere*/
  surface_sphere(S, 0.3, -0.1, 0.3, 1.3, s) ;
  agg_patch_mapping(P) = AGG_PATCH_SPHERICAL ;
  agg_patch_wrap_t(P) = TRUE ;
  
  nu = 16 ; nv = 33 ;
  for ( i = 0 ; i < nu ; i ++ ) {
    for ( j = 0 ; j < nv ; j ++ ) {
      agg_patch_map(P, (gdouble)i/(nu-1), (gdouble)j/(nv-1), &u, &v) ;
      agg_surface_point_eval(S, u, v, x, w) ;
      fprintf(stdout, "%e %e %e\n", x[0], x[1], x[2]) ;
    }
  }
  
  return ;
}
  
/* static void write_patch_grid(FILE *f, agg_patch_t *P, gint ns, gint nt) */

/* { */
/*   gint i, j, e1, e2 ; */
/*   gdouble s, t ; */

/*   for ( i = 0 ; i < ns ; i ++ ) { */
/*     for ( j = 0 ; j < nt ; j ++ ) { */
/*       agg_patch_sample(P, (gdouble)i/(ns-1), (gdouble)j/(nt-1), &s, &t) ; */
/*       fprintf(stdout, "%lg %lg\n", s, t) ; */
/*     } */
/*   } */

/*   return ; */
  
/*   e1 = 0 ; e2 = 2 ; */

/*   /\* for ( i = 0 ; i < ns ; i ++ ) { *\/ */
/*   /\*   patch_edge_interp(P, e1, (gdouble)i/(ns-1), &s, &t) ; *\/ */
/*   /\*   fprintf(f, "%lg %lg ", s, t) ; *\/ */
/*   /\*   patch_edge_interp(P, e2, 1.0 - (gdouble)i/(ns-1), &s, &t) ; *\/ */
/*   /\*   fprintf(f, "%lg %lg\n", s, t) ; *\/ */
/*   /\* } *\/ */

/*   e1 = 1 ; e2 = 3 ; */

/*   for ( i = 0 ; i < nt ; i ++ ) { */
/*     agg_patch_edge_interp(P, e1, (gdouble)i/(nt-1), &s, &t) ; */
/*     fprintf(f, "%lg %lg ", s, t) ; */
/*     agg_patch_edge_interp(P, e2, 1.0 - (gdouble)i/(nt-1), &s, &t) ; */
/*     fprintf(f, "%lg %lg\n", s, t) ; */
/*   } */
  
/*   return ; */
/* } */

/* static void write_patch_gmsh(FILE *f, agg_patch_t *P, agg_surface_t *S, */
/* 			     agg_surface_workspace_t *w) */

/* { */
/*   gint pps, i ; */
/*   gdouble x[3], u, v ; */

/*   for ( i = 0 ; i < agg_patch_point_number(P) ; i ++ ) { */
/*     /\* agg_patch_sample(P1, (gdouble)i/(nu-1), (gdouble)j/(nv-1), &u, &v) ; *\/ */
/*     u = agg_patch_point_s(P, i) ; */
/*     v = agg_patch_point_t(P, i) ; */
/*     agg_patch_map(P, u, v, &u, &v) ; */
/*     agg_surface_point_eval(S, u, v, x, w) ; */
/*     fprintf(f, "Point(%d) = {%e, %e, %e, lc} ;\n", i, x[0], x[1], x[2]) ; */
/*   } */

/*   for ( i = 0 ; i < agg_patch_point_number(P) - 1 ; i ++ ) { */
/*     fprintf(f, "Line(%d) = {%d, %d} ;\n", i, i, i+1) ; */
/*   }     */
  
/*   return ; */
/* } */

/* static void write_intersection_gmsh(FILE *f, agg_intersection_t *inter, */
/* 				    gint nsp, gint pps, */
/* 				    agg_surface_workspace_t *w) */

/* { */
/*   agg_patch_t *P ; */
/*   agg_surface_t *S ; */
/*   gint i, j ; */
/*   gdouble *x, u, v ; */
     
/*   for ( i = 0 ; i < agg_intersection_point_number(inter) ; i ++ ) { */
/*     x = agg_intersection_point(inter,i) ; */
/*     fprintf(f, "Point(%d) = {%e, %e, %e, lc} ;\n", i, x[0], x[1], x[2]) ; */
/*   } */

/*   for ( i = 0 ; i < nsp ; i ++ ) { */
/*     fprintf(f, "Spline(%d) = {", i) ; */
/*     for ( j = 0 ; j < pps-1 ; j ++ ) { */
/*       fprintf(f, "%d, ", i*(pps-1) + j) ; */
/*     } */
/*     fprintf(f, "%d} ;\n", i*(pps-1) + pps-1) ; */
/*   } */
  
/*   /\* for ( i = 0 ; i < agg_intersection_point_number(inter) - 1 ; i ++ ) { *\/ */
/*   /\*   fprintf(f, "Line(%d) = {%d, %d} ;\n", i, i, i+1) ; *\/ */
/*   /\* } *\/ */
  
/*   return ; */
/* } */

static void body_test(void)

{
  gint nsec, nseg, pps, offp, offsp, offs, nsurf, i, sgn ;
  gdouble len, rfuse, yw, zw, chord, span, taper, swle, dihedral ;
  gdouble yf, zf, chf, spf, tpf, swf, dhf ;
  agg_surface_t **S ; 
  agg_surface_workspace_t *w ;
  agg_wireframe_t *wf ;
  agg_intersection_t **inter, **resample ;
  agg_section_t *sc, *sw ;
  agg_patch_t **P ;
  FILE *output ;

  nsurf = 5 ;
  
  S = (agg_surface_t **)g_malloc0(nsurf*sizeof(agg_surface_t *)) ;
  P = (agg_patch_t **)g_malloc0(nsurf*sizeof(agg_patch_t *)) ;
  inter = (agg_intersection_t **)g_malloc0(50*sizeof(agg_intersection_t *)) ;
  resample = (agg_intersection_t **)g_malloc0(50*sizeof(agg_intersection_t *)) ;

  w = agg_surface_workspace_new() ;
  sc = agg_section_new(32, 32) ;
  agg_section_set_circle(sc) ;
  sw = agg_section_new(32, 32) ;
  agg_section_set_aerofoil(sw, 0.5, 0.1, 0) ;

  for ( i = 0 ; i < nsurf ; i ++ ) {
    S[i] = agg_surface_new(64) ;
    P[i] = agg_patch_new(1024) ;
    inter[i] = agg_intersection_new(8192) ;
    resample[i] = agg_intersection_new(8192) ;
  }    
  
  rfuse = 0.4 ; len = 5.0 ;
  /*wing parameters*/
  yw = -rfuse/4 ; zw = -1 ; chord = 1.5 ; span = 4.0 ; taper = 0.5 ;
  swle = 15*M_PI/180 ; dihedral = 5.0*M_PI/180.0 ;

  /*tailplane parameters*/
  yf = rfuse/8 ; zf = 0.5*len*0.7 ;  chf = 0.4 ; spf = 1.0 ;
  tpf = 0.9 ; swf = 18*M_PI/180 ; dhf = 3.0*M_PI/180.0 ;
  
  /* agg_patch_mapping(P[0]) = AGG_PATCH_SPHERICAL ; */
  agg_patch_mapping(P[0]) = AGG_PATCH_TUBULAR ;
  agg_patch_wrap_t(P[0]) = TRUE ;
  surface_tube(S[0], 0.0, 0.0, 0.0, rfuse, len, sc) ;

  sgn = 1 ;
  for ( i = 1 ; i < 3 ; i ++ ) {
    agg_patch_mapping(P[i]) = AGG_PATCH_SPHERICAL ;
    agg_patch_wrap_t(P[i]) = TRUE ;
    surface_wing(S[i], zw, yw, 0.0, chord, sgn*span, taper,
		 sgn*swle, sgn*dihedral, sw) ;
    sgn = -sgn ;
  }
  for ( i = 3 ; i < 5 ; i ++ ) {
    agg_patch_mapping(P[i]) = AGG_PATCH_SPHERICAL ;
    agg_patch_wrap_t(P[i]) = TRUE ;
    surface_wing(S[i], zf, yf, 0.0, chf, sgn*spf, tpf,
		 sgn*swf, sgn*dhf, sw) ;
    sgn = -sgn ;
  }

  wf = agg_wireframe_new(65536, 65536, 65536) ;
  nsec = 8 ; nseg = 65 ; pps = 2 ;
  offp = offsp = offs = 1 ;
  /*wing intersection with body calculated in different order for each
    wing to check handling of intersections when surfaces are added to
    body*/
  agg_surface_patch_intersection(inter[0], S[0], P[0], S[1], P[1], w) ;
  agg_surface_patch_intersection(inter[1], S[2], P[2], S[0], P[0], w) ;
  agg_surface_patch_intersection(inter[2], S[0], P[0], S[3], P[3], w) ;
  agg_surface_patch_intersection(inter[3], S[4], P[4], S[0], P[0], w) ;

  for ( i = 0 ; i < 4 ; i ++ ) {
    fprintf(stderr, "%d intersection points detected\n",
	    agg_patch_point_number(inter[i])) ;
    agg_intersection_resample(inter[i], nseg, pps, resample[i], w) ;
    agg_intersection_bbox_set(resample[i]) ;
    agg_wireframe_intersection_add(wf, resample[i], nseg, pps, w) ;
  }

  for ( i = 0 ; i < 5 ; i ++ ) {
    fprintf(stderr, "adding surface %d\n", i) ;
    agg_wireframe_surface_add(wf, S[i], P[i], nsec, nseg, pps, w) ;
  }

  output = fopen("surface1.geo", "w") ;
  fprintf(output, "lc = 0.1 ;\n") ;
  agg_wireframe_write_gmsh(output, wf, "lc", offp, offsp, offs, FALSE) ;
  fclose(output) ;

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
    circle_test() ;

    return 0 ;
  }

  if ( test == 2 ) {
    aerofoil_test() ;

    return 0 ;
  }

  if ( test == 3 ) {
    transform_test() ;

    return 0 ;
  }

  if ( test == 4 ) {
    surface_test() ;

    return 0 ;
  }

  if ( test == 5 ) {
    body_test() ;
    
    return 0 ;
  }
  
  return 0 ;
}

