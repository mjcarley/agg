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
			"parser",
			"ico",
			"fit",
			"section_fit",
			"derivatives",
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
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  /*scale circle radius along z axis: circle of radius 1/2 at u=1/2 */
  p[0] =  0.0 ; p[1] = 0.0 ;
  expr[2] = "2*sqrt(u*(1-u))" ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  /*shift centre to z=0*/
  p[0] =  0.0 ; p[1] = 0 ; p[2] = 0 ;
  expr[2] = "-1+u" ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  /*scale to required radius: sphere of radius r*/
  p[0] = 2*r ;
  /*shift to required centre*/
  agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, 0, 1, p, expr,
			     NULL, NULL, 1) ;
  p[0] = x ; p[1] = y ; p[2] = z ;
  expr[2] = NULL ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;

  agg_surface_umin(S) = 0.0 ; 
  agg_surface_umax(S) = 1.0 ; 

  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  agg_surface_weights_make(S) ;
  
  return ;
}

static void surface_tube(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			 agg_variable_t *global, gint nglobal,
			 agg_section_t *s)

{
  gdouble u, p[3] ;
  agg_transform_t *T ;
  gchar *expr[8] = {NULL}, etmp[1024] ;

  u = 0.0 ;
  agg_section_set_circle(s) ;
  T = agg_surface_transform(S) ;
  agg_transform_add_global_variables(T, global, nglobal) ;  

  agg_surface_section_add(S, s, u) ;
  /*centre on (0,0): circle of radius 1/2*/
  p[0] = -0.5 ; p[1] = 0 ; p[2] = 0 ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  /*scale on radius*/
  expr[0] = "2*radius" ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, 0, 1, p, expr,
			     NULL, NULL, 1) ;
  expr[0] = NULL ;
  /*translate in space to form cylinder*/
  p[0] = 0 ; p[1] = 0 ; p[2] = 0 ;
  /* sprintf(etmp, "-%lg/2 + %lg*u", len, len) ; */
  sprintf(etmp, "-len/2 + len*u") ;
  expr[2] = etmp ;
  /* expr[2] = "-len/2 + len*u" ; */
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;

  /*translate to appropriate centre*/
  p[0] = x ; p[1] = y ; p[2] = z ;
  expr[2] = NULL ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;

  agg_surface_axes(S) = AGG_AXES_PZ_PX_PY ;
  
  agg_surface_umin(S) = 0.0 ; 
  agg_surface_umax(S) = 1.0 ; 

  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  agg_surface_weights_make(S) ;
  
  return ;
}

static void surface_wing(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			 agg_variable_t *global, gint nglobal,
			 agg_section_t *s)

{
  gdouble u, p[3] ;
  agg_transform_t *T ;
  gchar etmp[3][1024], *expr[8] = {NULL} ;

  u = 0.0 ;
  T = agg_surface_transform(S) ;
  agg_surface_section_add(S, s, u) ;

  agg_transform_add_global_variables(T, global, nglobal) ;  
  
  expr[0] = expr[1] = expr[2] = NULL ; 
  /*scale on taper ratio, shrinking about leading edge*/
  p[0] = 0.0 ; p[1] = 0 ; p[2] = 0 ;
  sprintf(etmp[2], "1 + u*(taper-1)") ;
  expr[2] = etmp[2] ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  expr[0] = expr[1] = expr[2] = NULL ; 
  /*scale on root chord*/
  p[0] = 0 ;
  sprintf(etmp[0], "root") ;
  expr[0] = etmp[0] ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, 0, 1, p, expr,
			     NULL, NULL, 1) ;
  expr[0] = expr[1] = expr[2] = NULL ; 

  sprintf(etmp[0], "(1-u)^(1/8)") ;
  expr[0] = etmp[0] ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SCALE_Y, 0, 1, p, expr,
			     NULL, NULL, 1) ;
  expr[0] = expr[1] = expr[2] = NULL ; 
  
  /*translate along span and sweep; add del to stop wings
    meeting in the middle and confusing the intersection calculation*/
  p[0] = 0 ; p[1] = 0 ; p[2] = 0 ;
  sprintf(etmp[0], "span*sin(sweep)*u") ;
  expr[0] = etmp[0] ;
  sprintf(etmp[1], "span*sin(dihedral)*u") ;
  expr[1] = etmp[1] ;
  sprintf(etmp[2], "del + span*u") ;
  expr[2] = etmp[2] ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  expr[0] = expr[1] = expr[2] = NULL ; 
  p[0] = x ; p[1] = y ; p[2] = z ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  expr[0] = expr[1] = expr[2] = NULL ; 

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
  gdouble x, tot, err, S[64], dS[64], ee, err_d, dB ;
  gint j ;
  
  fprintf(stderr, "Bernstein basis polynomials\n") ;
  fprintf(stderr, "===========================\n") ;
  fprintf(stderr, "polynomial basis order: %d\n", n) ;

  ee = 1e-6 ;
  
  err = err_d = 0.0 ;
  for ( x = 0.0 ; x <= 1.0 ; x += 1.0/128 ) {
    tot = 0.0 ;
    for ( j = 0 ; j <= n ; j ++ ) {
      tot += agg_bernstein_basis_eval(n, j, x) ;
      dS[0] = agg_bernstein_derivative_eval(n, j, x) ;      
      dB = (agg_bernstein_basis_eval(n, j, x+ee/2) -
	    agg_bernstein_basis_eval(n, j, x-ee/2))/ee ;
      err_d = MAX(err_d, fabs(dS[0]-dB)) ;
    }
    err = MAX(err, fabs(tot-1.0)) ;
  }

  fprintf(stderr, "termwise maximum deviation from unity: %lg\n", err) ;
  fprintf(stderr, "maximum derivative error: %lg\n", err_d) ;

  err = 0.0 ;
  for ( x = 0.0 ; x <= 1.0 ; x += 1.0/128 ) {
    tot = 0.0 ;
    agg_bernstein_basis(n, x, S, dS) ;
    for ( j = 0 ; j <= n ; j ++ ) tot += S[j] ;
    err = MAX(err, fabs(tot-1.0)) ;
  }

  fprintf(stderr, "vector maximum deviation from unity: %lg\n", err) ;

  err_d = 0.0 ;
  for ( x = 0.0 ; x <= 1.0 ; x += 1.0/128 ) {
    tot = 0.0 ;
    agg_bernstein_basis(n, x, S, dS) ;
    for ( j = 0 ; j <= n ; j ++ ) {
      dB = (agg_bernstein_basis_eval(n, j, x+ee/2) -
	    agg_bernstein_basis_eval(n, j, x-ee/2))/ee ;
      err_d = MAX(err_d, fabs(dS[j]-dB)) ;
    }
  }

  fprintf(stderr, "maximum derivative error: %lg\n", err_d) ;
  
  return ;
}

static void circle_test(void)

{
  gint i, nx ;
  gdouble x, y, dy, r, err ;
  agg_section_t *s ;

  s = agg_section_new(8, 8) ;
  agg_section_set_circle(s) ;
  
  nx = 128 ;
  err = 0.0 ;
  for ( i = 0 ; i <= nx ; i ++ ) {
    x = -1.0 + 2.0*i/nx ;
    y = agg_section_eval(s, x) ;
    dy = agg_section_diff(s, x) ;
    fprintf(stdout, "%e %e %e\n", fabs(x), y, dy) ;
    r = sqrt((fabs(x)-0.5)*(fabs(x)-0.5) + y*y) ;
    err = MAX(err, fabs(r-0.5)) ;
  }

  fprintf(stderr, "radius error: %lg\n", err) ;
  
  return ;
}

static void aerofoil_test(void)

{
  gint i, nx ;
  gdouble x, y, dy ;
  agg_section_t *s ;

  s = agg_section_new(8, 8) ;
  agg_section_set_aerofoil(s, 0.5, 0.1, 0.01) ;
  
  nx = 1024 ;
  for ( i = 0 ; i <= nx ; i ++ ) {
    x = -1.0 + 2.0*i/nx ;
    y = agg_section_eval(s, x) ;
    dy = agg_section_diff(s, x) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e %1.16e\n", fabs(x), y, x, dy) ;
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
  agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, 0, 1, p, args,
			     NULL, NULL, 3) ;
  args[0] = NULL ; args[1] = "0" ; args[2] = "twist*u" ;
  p[0] = 0.5 ; p[1] = 0 ; p[2] = 0 ;
  /*translate in z to generate "wing"*/
  agg_transform_operator_add(T, AGG_TRANSFORM_ROTATE, 0, 1, p, args,
			     NULL, NULL, 3) ;
  args[0] = NULL ; args[1] = NULL ; args[2] = "u*3.0" ;
  p[0] = 0.0 ; p[1] = 0 ; p[2] = 0 ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, args,
			     NULL, NULL, 3) ;

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
  agg_surface_t *S ;
  agg_surface_workspace_t *w ;
  agg_section_t *s ;
  agg_patch_t *P ;
  gdouble u, v, x[3] ;
  
  S = agg_surface_new(64) ;
  P = agg_patch_new() ;
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

static void body_test(void)

{
  gint nsec, nsp, pps, offp, offsp, offs, nsurf, i, j ;
  gdouble rfuse, yw, zw ;
  agg_surface_t **S ; 
  agg_surface_workspace_t *w ;
  agg_mesh_t *wf ;
  agg_body_t *b ;
  agg_section_t *sc, *sw ;
  agg_patch_t **P ;
  FILE *output ;

  nsurf = 2 ;

  b = agg_body_new(32, 32) ;
  
  agg_body_global_add(b, "radius", NULL, 0.4) ;
  agg_body_global_add(b, "len", NULL, 5) ;
  agg_body_global_add(b, "root", NULL, 1.5) ;
  agg_body_global_add(b, "span", "len*0.7", 5) ;
  agg_body_global_add(b, "dihedral", "5*pi/180", 0) ;
  agg_body_global_add(b, "sweep", "25*pi/180", 0) ;
  agg_body_global_add(b, "taper", NULL, 0.5) ;
  agg_body_global_add(b, "del", "span*0", 0) ;

  agg_body_globals_compile(b) ;
  
  agg_body_globals_eval(b) ;
  agg_body_globals_write(stderr, b) ;
  
  S = (agg_surface_t **)g_malloc0(nsurf*sizeof(agg_surface_t *)) ;
  P = (agg_patch_t **)g_malloc0(nsurf*sizeof(agg_patch_t *)) ;

  w = agg_surface_workspace_new() ;
  sc = agg_section_new(32, 32) ;
  agg_section_set_circle(sc) ;
  sw = agg_section_new(32, 32) ;
  agg_section_set_aerofoil(sw, 0.5, 0.2, 0) ;

  for ( i = 0 ; i < nsurf ; i ++ ) {
    S[i] = agg_surface_new(64) ;
    P[i] = agg_patch_new() ;
  }
  
  rfuse = 0.4 ;
  /* len = 5.0 ; */
  /*wing parameters*/
  yw = -rfuse/4 ; zw = -1 ;

  /*tailplane parameters*/
  /* yf = rfuse/8 ; zf = 0.5*len*0.7 ;  chf = 0.4 ; spf = 1.0 ; */
  /* tpf = 0.9 ; swf = 18*M_PI/180 ; dhf = 3.0*M_PI/180.0 ; */
  
  for ( i = 0 ; i < nsurf ; i ++ ) {
    agg_patch_mapping(P[i]) = AGG_PATCH_SPHERICAL ;
    agg_patch_wrap_t(P[i]) = TRUE ;
    surface_wing(S[i], zw, yw, 0.0,
		 agg_body_globals(b), agg_body_global_number(b),
		 sw) ;
  }
  
  agg_surface_axes(S[1]) = AGG_AXES_PX_PY_MZ ;  
  agg_patch_invert(P[1]) = TRUE ;
  
  wf = agg_mesh_new(65536, 65536, 65536) ;
  nsec = 12 ; nsp = 8 ; pps = 4 ;
  offp = offsp = offs = 1 ;

  for ( i = 0 ; i < nsurf ; i ++ ) {
    fprintf(stderr, "adding surface %d\n", i) ;
    agg_mesh_surface(wf, i) = S[i] ; 
    agg_mesh_patch(wf, i) = P[i] ;
    agg_mesh_surface_number(wf) ++ ;
  }

  /* agg_mesh_body_regular(wf, b, nsec, nsp, pps, w) ;   */
  g_assert_not_reached() ;
  
  output = fopen("surface.geo", "w") ;
  fprintf(output, "lc = 0.1 ;\n") ;
  agg_mesh_write_gmsh(output, wf, "lc", offp, offsp, offs, FALSE) ;
  fclose(output) ;

  return ;
}

static void parser_test(gchar *file)

{
  agg_body_t *b ;
  agg_mesh_t *m ;
  gint nsec, nsp, pps, offp, offsp, offs ;
  agg_surface_workspace_t *w ;
  FILE *output ;
  
  b = agg_body_new(32, 32) ;
  w = agg_surface_workspace_new() ;

  agg_body_read(b, file, TRUE) ;

  agg_body_globals_compile(b) ;
  agg_body_globals_eval(b) ;
  agg_body_globals_write(stderr, b) ;
  fprintf(stderr, "%d surface%s\n", agg_body_surface_number(b),
	  (agg_body_surface_number(b) > 1 ? "s" : "")) ;
  agg_body_surfaces_list(stderr, b) ;

  m = agg_mesh_new(65536, 65536, 65536) ;
  nsec = 8 ; nsp = 65 ; pps = 2 ;
  offp = offsp = offs = 1 ;

  /* m->B[0].invert = TRUE ; */
  /* m->B[1].invert = TRUE ; */
  
  agg_mesh_body(m, b, pps, w) ;

  output = fopen("surface.geo", "w") ;
  fprintf(output, "lc = 0.1 ;\n") ;
  agg_mesh_write_gmsh(output, m, "lc", offp, offsp, offs, FALSE) ;
  fclose(output) ;
  
  return ;
}

static void derivative_test(gchar *file)

{
  agg_body_t *b ;
  gdouble u, v, de, x[3] ;
  gint nu, nv, i, j ;
  agg_surface_workspace_t *w ;
  FILE *output ;
  
  b = agg_body_new(32, 32) ;
  w = agg_surface_workspace_new() ;

  agg_body_read(b, file, TRUE) ;

  agg_body_globals_compile(b) ;
  agg_body_globals_eval(b) ;
  agg_body_globals_write(stderr, b) ;
  agg_body_surfaces_list(stderr, b) ;

  nu = 8 ; nv = 16 ;
  for ( i = 1 ; i < nu ; i ++ ) {
    u = (gdouble)i/nu ;
    for ( j = 1 ; j < nv ; j ++ ) {
      v = -1 + (gdouble)j*2.0/nv ;
      agg_surface_point_eval(agg_body_surface(b,0), u, v, x, w) ;
      fprintf(stdout, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
    }
  }
  
  return ;
}

#if 0
static void ico_test(void)

{
  gdouble st[1024], u, v, p[3], x[3] ;
  gint e[2000], t[1024], i, nt, np, ne, nl, offp, offsp, offs ;
  agg_mesh_t *m ;
  agg_patch_t *P ;
  agg_surface_t *S ;
  agg_transform_t *T ;
  agg_section_t *s ;
  gchar *expr[8] = {NULL} ;
  agg_surface_workspace_t *w ;
  gboolean physical_points ;
  FILE *output ;
  
  S = agg_surface_new(64) ;
  P = agg_patch_new() ;
  agg_patch_mapping(P) = AGG_PATCH_SPHERICAL ;
  agg_patch_wrap_t(P) = TRUE ;
  w = agg_surface_workspace_new() ;
  s = agg_section_new(32, 32) ;
  agg_section_set_circle(s) ;

  physical_points = TRUE ;
  
  agg_surface_section_add(S, s, 0) ;
  T = agg_surface_transform(S) ;
  /*centre on (0,0): circle of radius 1/2*/
  p[0] = -0.5 ; p[1] = 0 ; p[2] = 0 ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  /*scale circle radius along z axis: circle of radius 1/2 at u=1/2 */
  p[0] =  0.0 ; p[1] = 0.0 ;
  expr[2] = "2*sqrt(u*(1-u))" ;
  agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, 0, 1, p, expr,
			     NULL, NULL, 3) ;
  /*flatten the shape */
  /* p[0] =  0.0 ; p[1] = 0.0 ; */
  /* expr[0] = "sqrt(u*(1-u))" ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_SCALE_Y, p, expr, 1) ; */
  /*translate circle in z to form sphere*/
  p[0] =  0.0 ; p[1] = 0 ; p[2] = 0 ;
  expr[2] = "u" ;
  agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr,
			     NULL, NULL, 3) ;

  agg_surface_umin(S) = 0.0 ; 
  agg_surface_umax(S) = 1.0 ; 
  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  agg_surface_weights_make(S) ;

  m = agg_mesh_new(65536, 65536, 65536) ;

  agg_mesh_surface(m, 0) = S ;
  agg_mesh_patch(m, 0) = P ;
  agg_mesh_surface_number(m) == 1 ;

  agg_mesh_icosahedron(m, 0, w) ;
  
  agg_mesh_refine_loop(m, w) ;
  agg_mesh_refine_loop(m, w) ;
  agg_mesh_refine_loop(m, w) ;
  
  offp = 1 ; offsp = 0 ; offs = 0 ;
  output = fopen("surface.geo", "w") ;
  fprintf(output, "lc = 0.1 ;\n") ;
  agg_mesh_write_gmsh(output, m, "lc", offp, offsp, offs, FALSE) ;
  fclose(output) ;

  return ;
}
#endif

static void fit_test(void)

{
  gint i, nxu, nxl, nl, nu ;
  gdouble x, y, dy, r, err, xu[64], yu[64], xl[64], yl[64] ;
  gdouble t, p, m ;
  agg_section_t *s ;

  t = 0.2 ; p = 0.1 ; m = 0.3 ;
  s = agg_section_new(32, 32) ;

  nxu = 37 ;
  for ( i = 0 ; i < nxu ; i ++ ) {
    xu[i] = (gdouble)i/(nxu-1) ;
    yu[i] = agg_naca_four(t, p, m, xu[i]) ;
  }
  nxl = 50 ;
  for ( i = 0 ; i < nxl ; i ++ ) {
    xl[i] = (gdouble)i/(nxl-1) ;
    yl[i] = agg_naca_four(t, p, m, -xl[i]) ;
  }

  agg_section_fit(s, xu, 1, yu, 1, nxu, xl, 1, yl, 1, nxl, 0.5, 0, 10, 10) ;

  /* agg_section_trailing_edge_upper(s) = 0 ; */
  /* agg_section_trailing_edge_lower(s) = 0 ; */

  err = 0 ;
  for ( i = 0 ; i <= 128 ; i ++ ) {
    x = -1 + 2.0*(gdouble)i/128 ;
    y = agg_section_eval(s, x) ;
    fprintf(stdout, "%lg %lg ", ABS(x), y) ;
    y = agg_naca_four(t, p, m, x) ;
    fprintf(stdout, "%lg\n", y) ;
    err = MAX(err, ABS(y-agg_section_eval(s, x))) ;
  }

  fprintf(stderr, "maximum error: %lg\n", err) ;
  
  return ;
}

static void section_fit_test(void)

{
  gint nsec, nsp, pps, offp, offsp, offs, i, j ;
  gdouble rfuse, yw, zw, u, v, x[3], y[3], emax ;
  agg_surface_t *S ; 
  agg_surface_workspace_t *w ;
  agg_mesh_t *wf ;
  agg_body_t *b ;
  agg_section_t *sc, *sw ;
  agg_patch_t *P ;
  agg_chebyshev_t *C ;
  FILE *output ;

  b = agg_body_new(32, 32) ;
  
  agg_body_global_add(b, "radius", NULL, 0.4) ;
  agg_body_global_add(b, "len", NULL, 5) ;
  agg_body_global_add(b, "root", NULL, 1.5) ;
  agg_body_global_add(b, "span", "len*0.7", 5) ;
  agg_body_global_add(b, "dihedral", "5*pi/180", 0) ;
  agg_body_global_add(b, "sweep", "25*pi/180", 0) ;
  agg_body_global_add(b, "taper", NULL, 0.5) ;
  agg_body_global_add(b, "del", "span*0", 0) ;

  agg_body_globals_compile(b) ;
  
  agg_body_globals_eval(b) ;
  agg_body_globals_write(stderr, b) ;
  
  w = agg_surface_workspace_new() ;
  sc = agg_section_new(32, 32) ;
  agg_section_set_circle(sc) ;
  sw = agg_section_new(32, 32) ;
  /* agg_section_set_ellipse(sw, 0.25) ; */
  agg_section_set_aerofoil(sw, 0.5, 0.25, 0) ;

  S = agg_surface_new(64) ;
  P = agg_patch_new() ;
  
  rfuse = 0.4 ;
  yw = -rfuse/4 ; zw = -1 ;
  
  agg_patch_mapping(P) = AGG_PATCH_SPHERICAL ;
  agg_patch_wrap_t(P) = TRUE ;
  surface_wing(S, zw, yw, 0.0,
	       agg_body_globals(b), agg_body_global_number(b),
	       sw) ;

  C = agg_chebyshev_new(16384, 3) ;
  u = 0.5 ;
  /* agg_chebyshev_surface_section_adaptive(C, S, u, 8, 1e-9, 1e-1, w) ; */
  agg_chebyshev_surface_section(C, S, u, 32, w) ;

  fprintf(stderr, "%d intervals on Chebyshev interpolator\n",
	  agg_chebyshev_interval_number(C)) ;
  fprintf(stderr, "shortest interval: %lg\n",
	  agg_chebyshev_interval_shortest(C)) ;
  emax = 0 ;
  for ( i = 0 ; i <= 313 ; i ++ ) {
    v = -1.0 + 2.0*(gdouble)i/313 ;
    agg_surface_point_eval(S, u, v, x, w) ;
    agg_chebyshev_eval(C, v, y) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
	    v, x[0], x[1], x[2], y[0], y[1], y[2]) ;
    emax = MAX(emax, fabs(x[0] - y[0])) ;
    emax = MAX(emax, fabs(x[1] - y[1])) ;
    emax = MAX(emax, fabs(x[2] - y[2])) ;
  }

  fprintf(stderr, "maximum geometric error: %lg\n", emax) ;
  
  return ;
}

gint main(gint argc, gchar **argv)

{
  gchar ch, *progname, *file = "test.agg" ;
  gint test ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;
  test = -1 ;
  
  while ( (ch = getopt(argc, argv, "i:T:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'i': file = g_strdup(optarg) ; break ;
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

  if ( test == 6 ) {
    parser_test(file) ;
    
    return 0 ;
  }

  /* if ( test == 7 ) { */
  /*   ico_test() ; */
    
  /*   return 0 ; */
  /* } */

  if ( test == 8 ) {
    fit_test() ;
    
    return 0 ;
  }

  if ( test == 9 ) {
    section_fit_test() ;
    
    return 0 ;
  }
  
  if ( test == 10 ) {
    derivative_test(file) ;
    
    return 0 ;
  }

  return 0 ;
}

