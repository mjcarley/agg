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

#include "agg-private.h"

/* #include "hefsi.h" */

const char *tests[] = {"bernstein",
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
			"affine",
			"intersection",
			"section_interp",
			"surface_diff",
			""} ;

static void surface_sphere(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			   gdouble r, agg_section_t *s)

{
  gdouble u, p[3] ;
  gint i ;
  agg_transform_t *T ;
  agg_affine_t *A ;
  char *expr[8] = {NULL} ;

  u = 0.0 ;
  agg_section_set_circle(s) ;
  T = agg_surface_transform(S) ;
  agg_surface_section_add(S, s, u) ;
  /*centre on (0,0): circle of radius 1/2*/
  A = agg_affine_new(0) ;
  p[0] = -0.5 ; p[1] = 0 ; p[2] = 0 ;
  expr[0] = NULL ; expr[1] = NULL ; expr[2] = NULL ;
  agg_affine_translation(A, p, expr, 3) ;
  agg_transform_affine_add(T, A) ;
  /*scale circle radius along z axis: circle of radius 1/2 at u=1/2 */
  A = agg_affine_new(0) ;
  p[0] =  1.0 ; p[1] = 1.0 ;
  expr[0] = NULL ; expr[1] = NULL ;
  expr[0] = "2*sqrt(u*(1-u))" ;
  expr[1] = "2*sqrt(u*(1-u))" ;
  expr[2] = NULL ;
  
  agg_affine_scale(A, p, expr, 3) ;
  agg_transform_affine_add(T, A) ;
  /*shift centre to z=0*/
  A = agg_affine_new(0) ;
  p[0] =  0.0 ; p[1] = 0 ; p[2] = 0 ;
  expr[0] = NULL ; expr[1] = NULL ; expr[2] = "(-0.5+u)" ;
  agg_affine_translation(A, p, expr, 3) ;
  agg_transform_affine_add(T, A) ;
  /*scale to required radius: sphere of radius r*/
  A = agg_affine_new(0) ;
  p[0] = 2*r ; expr[0] = NULL ;
  agg_affine_scale(A, p, expr, 1) ;
  agg_transform_affine_add(T, A) ;
  /*shift to required centre*/
  A = agg_affine_new(0) ;
  p[0] = x ; p[1] = y ; p[2] = z ;
  expr[0] = NULL ; expr[1] = NULL ; expr[2] = NULL ;
  agg_affine_translation(A, p, expr, 3) ;
  agg_transform_affine_add(T, A) ;

  agg_surface_umin(S) = 0.0 ; 
  agg_surface_umax(S) = 1.0 ; 

  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
    A = agg_transform_affine(T,i) ;
    agg_affine_expressions_compile(A, T->e) ;
  }
  
  agg_surface_weights_make(S) ;
  
  return ;
}

static void surface_tube(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			 agg_variable_t *global, gint nglobal,
			 agg_section_t *s)

{
  gdouble u, p[3] ;
  agg_transform_t *T ;
  char *expr[8] = {NULL}, etmp[1024] ;

  g_assert_not_reached() ;
  
  /* u = 0.0 ; */
  /* agg_section_set_circle(s) ; */
  /* T = agg_surface_transform(S) ; */
  /* agg_transform_add_global_variables(T, global, nglobal) ;   */

  /* agg_surface_section_add(S, s, u) ; */
  /* /\*centre on (0,0): circle of radius 1/2*\/ */
  /* p[0] = -0.5 ; p[1] = 0 ; p[2] = 0 ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr, */
  /* 			     NULL, NULL, 3) ; */
  /* /\*scale on radius*\/ */
  /* expr[0] = "2*radius" ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, 0, 1, p, expr, */
  /* 			     NULL, NULL, 1) ; */
  /* expr[0] = NULL ; */
  /* /\*translate in space to form cylinder*\/ */
  /* p[0] = 0 ; p[1] = 0 ; p[2] = 0 ; */
  /* /\* sprintf(etmp, "-%lg/2 + %lg*u", len, len) ; *\/ */
  /* sprintf(etmp, "-len/2 + len*u") ; */
  /* expr[2] = etmp ; */
  /* /\* expr[2] = "-len/2 + len*u" ; *\/ */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr, */
  /* 			     NULL, NULL, 3) ; */

  /* /\*translate to appropriate centre*\/ */
  /* p[0] = x ; p[1] = y ; p[2] = z ; */
  /* expr[2] = NULL ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr, */
  /* 			     NULL, NULL, 3) ; */

  /* agg_surface_axes(S) = AGG_AXES_PZ_PX_PY ; */
  
  /* agg_surface_umin(S) = 0.0 ;  */
  /* agg_surface_umax(S) = 1.0 ;  */

  /* agg_expression_data_compile(T->e) ; */
  /* agg_transform_expressions_compile(T) ; */

  /* agg_surface_weights_make(S) ; */
  
  return ;
}

static void surface_wing(agg_surface_t *S, gdouble x, gdouble y, gdouble z,
			 agg_variable_t *global, gint nglobal,
			 agg_section_t *s)

{
  gdouble u, p[3] ;
  agg_transform_t *T ;
  char etmp[3][1024], *expr[8] = {NULL} ;

  g_assert_not_reached() ;
  
  /* u = 0.0 ; */
  /* T = agg_surface_transform(S) ; */
  /* agg_surface_section_add(S, s, u) ; */

  /* agg_transform_add_global_variables(T, global, nglobal) ;   */
  
  /* expr[0] = expr[1] = expr[2] = NULL ;  */
  /* /\*scale on taper ratio, shrinking about leading edge*\/ */
  /* p[0] = 0.0 ; p[1] = 0 ; p[2] = 0 ; */
  /* sprintf(etmp[2], "1 + u*(taper-1)") ; */
  /* expr[2] = etmp[2] ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, 0, 1, p, expr, */
  /* 			     NULL, NULL, 3) ; */
  /* expr[0] = expr[1] = expr[2] = NULL ;  */
  /* /\*scale on root chord*\/ */
  /* p[0] = 0 ; */
  /* sprintf(etmp[0], "root") ; */
  /* expr[0] = etmp[0] ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_SCALE, 0, 1, p, expr, */
  /* 			     NULL, NULL, 1) ; */
  /* expr[0] = expr[1] = expr[2] = NULL ;  */

  /* sprintf(etmp[0], "(1-u)^(1/8)") ; */
  /* expr[0] = etmp[0] ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_SCALE_Y, 0, 1, p, expr, */
  /* 			     NULL, NULL, 1) ; */
  /* expr[0] = expr[1] = expr[2] = NULL ;  */
  
  /* /\*translate along span and sweep; add del to stop wings */
  /*   meeting in the middle and confusing the intersection calculation*\/ */
  /* p[0] = 0 ; p[1] = 0 ; p[2] = 0 ; */
  /* sprintf(etmp[0], "span*sin(sweep)*u") ; */
  /* expr[0] = etmp[0] ; */
  /* sprintf(etmp[1], "span*sin(dihedral)*u") ; */
  /* expr[1] = etmp[1] ; */
  /* sprintf(etmp[2], "del + span*u") ; */
  /* expr[2] = etmp[2] ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr, */
  /* 			     NULL, NULL, 3) ; */
  /* expr[0] = expr[1] = expr[2] = NULL ;  */
  /* p[0] = x ; p[1] = y ; p[2] = z ; */
  /* agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, expr, */
  /* 			     NULL, NULL, 3) ; */
  /* expr[0] = expr[1] = expr[2] = NULL ;  */

  /* agg_surface_umin(S) = 0.0 ;  */
  /* agg_surface_umax(S) = 1.0 ;  */

  /* agg_expression_data_compile(T->e) ; */
  /* agg_transform_expressions_compile(T) ; */

  /* agg_surface_weights_make(S) ; */
  
  return ;
}

static gint parse_test(char *s)

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

/* static void transform_test(void) */

/* { */
/*   gint i, j, nx, nu, nglobal ; */
/*   gdouble x, y, xin[3], xout[3], p[32] ; */
/*   agg_section_t *s ; */
/*   agg_transform_t *T ; */
/*   agg_variable_t global[4] ; */
/*   char *args[32] ; */

/*   s = agg_section_new(8, 8) ; */
/*   agg_section_set_aerofoil(s, 0.5, 0.1, 0.0) ; */

/*   T = agg_transform_new(64) ; */

/*   /\* */
/*    * set up global variables:  */
/*    * u is surface parameter for distribution */
/*    * twist is total twist at wingtip */
/*    *\/ */
/*   nglobal = 2 ; */
/*   global[0].name = g_strdup("u") ; */
/*   global[0].def =  NULL ; */
/*   global[0].val = 0.0 ; */
/*   global[1].name = g_strdup("twist") ; */
/*   global[1].def =  g_strdup("15.0*pi/180.0") ; */
/*   global[1].val = 0.4 ; */

/*   /\*non-linear scaling, shrinking about mid-chord*\/ */
/*   args[0] = NULL ; args[1] = "0" ; args[2] = "0.5*u+(1-u)^2" ; */
/*   p[0] = 0.5 ; p[1] = 0 ; p[2] = 0 ; */
/*   /\*linear twist*\/ */
/*   agg_transform_operator_add(T, AGG_TRANSFORM_SHRINK, 0, 1, p, args, */
/* 			     NULL, NULL, 3) ; */
/*   args[0] = NULL ; args[1] = "0" ; args[2] = "twist*u" ; */
/*   p[0] = 0.5 ; p[1] = 0 ; p[2] = 0 ; */
/*   /\*translate in z to generate "wing"*\/ */
/*   agg_transform_operator_add(T, AGG_TRANSFORM_ROTATE, 0, 1, p, args, */
/* 			     NULL, NULL, 3) ; */
/*   args[0] = NULL ; args[1] = NULL ; args[2] = "u*3.0" ; */
/*   p[0] = 0.0 ; p[1] = 0 ; p[2] = 0 ; */
/*   agg_transform_operator_add(T, AGG_TRANSFORM_TRANSLATE, 0, 1, p, args, */
/* 			     NULL, NULL, 3) ; */

/*   /\* agg_transform_operators_write(stderr, T) ; *\/ */
  
/*   agg_transform_add_global_variables(T, global, nglobal) ; */
/*   agg_expression_data_compile(T->e) ; */
/*   agg_transform_expressions_compile(T) ; */
/*   agg_transform_variables_eval(T) ; */

/*   nx = 65 ; nu = 17 ; */
/*   for ( j = 0 ; j < nu ; j ++ ) { */
/*     global[0].val = (gdouble)j/(nu-1) ; */
/*     /\*reevaluate all variables after changing the parameters*\/ */
/*     agg_transform_variables_eval(T) ; */
/*     for ( i = 0 ; i <= nx ; i ++ ) { */
/*       x = -1.0 + 2.0*i/nx ; */
/*       y = agg_section_eval(s, x) ; */
/*       /\* fprintf(stdout, "%e %e\n", fabs(x), y) ; *\/ */
/*       xin[0] = fabs(x) ; xin[1] = y ; xin[2] = xout[2] = 0 ; */
/*       agg_transform_apply(T, xin, xout) ; */
/*       fprintf(stdout, "%e %e %e\n", xout[0], xout[1], xout[2]) ; */
/*     } */
/*   } */
  
/*   return ; */
/* } */

static void surface_test(void)

{
  gint nu, nv, i, j ;
  agg_surface_t *S ;
  agg_surface_workspace_t *w ;
  agg_section_t *s ;
  /* agg_patch_t *P ; */
  gdouble u, v, x[3] ;
  
  S = agg_surface_new(64) ;
  /* P = agg_patch_new() ; */
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
  /* surface_sphere(S, 0.3, -0.1, 0.3, 1.3, s) ; */
  /* agg_patch_mapping(P) = AGG_PATCH_SPHERICAL ; */
  /* agg_patch_wrap_t(P) = TRUE ; */
  
  /* nu = 16 ; nv = 33 ; */
  /* for ( i = 0 ; i < nu ; i ++ ) { */
  /*   for ( j = 0 ; j < nv ; j ++ ) { */
  /*     agg_patch_map(P, (gdouble)i/(nu-1), (gdouble)j/(nv-1), &u, &v) ; */
  /*     agg_surface_point_eval(S, u, v, x, w) ; */
  /*     fprintf(stdout, "%e %e %e\n", x[0], x[1], x[2]) ; */
  /*   } */
  /* } */
  
  return ;
}

static void body_test(void)

{
#if 0
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
#endif
  return ;
}

static void parser_test(char *file)

{
  agg_body_t *b ;
  agg_mesh_t *m ;
  gint nsec, nsp, pps, offp, offsp, offs, i ;
  agg_surface_workspace_t *w ;
  agg_triangulation_settings_t settings ;
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
  
  /* agg_mesh_body(m, b, pps, w) ; */

  agg_triangulation_type(&settings) = AGG_TRIANGULATION_GRID_REGULAR ;
  agg_triangulation_section_number(&settings) = 7 ;  
  agg_triangulation_section_point_number(&settings) = 38 ;
  agg_triangulation_points_per_spline(&settings) = 5 ;
  agg_triangulation_invert(&settings) = FALSE ;
  agg_triangulation_sampling_s(&settings) = AGG_SAMPLING_LINEAR ;
  agg_triangulation_sampling_t(&settings) = AGG_SAMPLING_COSINE ;
  
  for ( i = 0 ; i < agg_body_surface_number(b) ; i ++ ) {
    agg_mesh_surface(m,i) = agg_body_surface(b,i) ;
    agg_mesh_surface_number(m) ++ ;
    agg_mesh_surface_triangulate(m, i,
				 agg_body_triangulation_settings(b,i), w) ;
  }
  
  output = fopen("surface.geo", "w") ;
  fprintf(output, "lc = 0.1 ;\n") ;
  agg_mesh_write_gmsh(output, m, "lc", offp, offsp, offs, FALSE) ;
  fclose(output) ;
  
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
  char *expr[8] = {NULL} ;
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
#if 0
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
#endif  
  return ;
}

static void affine_test(void)

{
  gint i, j, nx, nu, nglobal ;
  gdouble x, y, xin[3], xout[3], p[32] ;
  agg_section_t *s ;
  agg_transform_t *T ;
  agg_variable_t global[4], params[32] ;
  agg_affine_t *A ;
  char *args[32] ;

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

  agg_transform_add_global_variables(T, global, nglobal) ;
  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;
  agg_transform_variables_eval(T) ;

  A = agg_affine_new(0) ;
  agg_variable_definition(&(params[0])) = "translate" ;
  agg_variable_definition(&(params[1])) = NULL ;
  agg_variable_value(&(params[1]))      = 0.0 ;
  
  agg_variable_definition(&(params[2])) = NULL ;
  agg_variable_value(&(params[2]))      = 0.5 ;
  agg_variable_definition(&(params[3])) = NULL ;
  agg_variable_value(&(params[3]))      = 0.0 ;
  agg_affine_parse(A, params, 4) ;
  
  agg_affine_expressions_compile(A, T->e) ;
  agg_transform_affine_add(T, A) ;
  A = agg_affine_new(0) ;
  args[0] = "2*pi*u" ; p[0] = 0 ;
  agg_affine_rotation_x(A, p, args, 1) ;
  agg_affine_expressions_compile(A, T->e) ;
  agg_transform_affine_add(T, A) ;

  nx = 65 ; nu = 17 ;
  for ( j = 0 ; j < nu ; j ++ ) {
    gdouble tform[16] = {
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1} ;
    global[0].val = (gdouble)j/(nu-1) ;
    /*reevaluate all variables after changing the parameters*/
    agg_transform_variables_eval(T) ;
    for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
      A = agg_transform_affine(T,i) ;
      agg_affine_matrices_evaluate(A) ;
    }
    for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
      A = agg_transform_affine(T,i) ;
      agg_mat_mat_mul_4(tform, agg_affine_matrix(A,0), tform) ;
      /* matmat_mul_4(tform, agg_affine_matrix(A,0), tform) ; */
    }
    for ( i = 0 ; i <= nx ; i ++ ) {
      x = -1.0 + 2.0*i/nx ;
      y = agg_section_eval(s, x) ;
      /* fprintf(stdout, "%e %e\n", fabs(x), y) ; */
      xin[0] = fabs(x) ; xin[1] = y ; xin[2] = 0.0 ;
      agg_affine_point_transform(xout, tform, xin) ;
      fprintf(stdout, "%e %e %e\n", xout[0], xout[1], xout[2]) ;
    }
  }
  
  return ;
}

/* static gint hefsi_func(gdouble s, gdouble t, gdouble *x, gpointer data) */

/* { */
/*   gpointer *hdata = data ; */
/*   agg_surface_t *S = hdata[0] ; */
/*   agg_surface_workspace_t *w = hdata[2] ; */
/*   gdouble u, v ; */

/*   agg_surface_point_eval(S, s, t, x, w) ; */

/*   return 0 ; */
/* } */

/* static void intersection_test(void) */

/* { */
/*   gint nu, nv, i, j ; */
/*   agg_surface_t *S1, *S2 ; */
/*   agg_surface_workspace_t *w ; */
/*   agg_section_t *s ; */
/*   gdouble u, v, x[3] ; */
/*   agg_intersection_t *inter ; */
/*   hefsi_surface_t *h1, *h2 ; */
/*   hefsi_workspace_t *wh ; */
/*   hefsi_segment_t *seg1 ; */
/*   gpointer data1[3], data2[3] ; */
/*   gint dmin, dmax ; */
/*   gdouble scale, tol ; */
/*   GSList *il ; */
  
/*   S1 = agg_surface_new(64) ; */
/*   S2 = agg_surface_new(64) ; */
/*   w = agg_surface_workspace_new() ; */
/*   s = agg_section_new(32, 32) ; */
/*   agg_section_set_circle(s) ; */

/*   data1[0] = S1 ; data1[2] = w ; */
/*   data2[0] = S2 ; data2[2] = w ; */
  
/*   /\*sphere*\/ */
/*   surface_sphere(S1, 0.3, -0.1, 0.3, 1.3, s) ; */
/*   surface_sphere(S2, 0.4,  0.7, 0.9, 1.3, s) ; */
  
/*   nu = 16 ; nv = 33 ; */
/*   for ( i = 0 ; i <= nu ; i ++ ) { */
/*     u = (gdouble)i/nu ; */
/*     for ( j = 0 ; j <= nv ; j ++ ) { */
/*       v = -1 + 2.0*(gdouble)j/nv ; */
/*       agg_surface_point_eval(S1, u, v, x, w) ; */
/*       fprintf(stdout, "%e %e %e\n", x[0], x[1], x[2]) ; */
/*     } */
/*   } */

/*   for ( i = 0 ; i <= nu ; i ++ ) { */
/*     u = (gdouble)i/nu ; */
/*     for ( j = 0 ; j <= nv ; j ++ ) { */
/*       v = -1 + 2.0*(gdouble)j/nv ; */
/*       agg_surface_point_eval(S2, u, v, x, w) ; */
/*       fprintf(stdout, "%e %e %e\n", x[0], x[1], x[2]) ; */
/*     } */
/*   } */

/*   inter = agg_intersection_new(8192) ; */

/*   wh = hefsi_workspace_new() ; */

/*   dmin = 6 ; dmax = 10 ; scale = 18/16.0 ; tol = 1e-6 ; */
/*   h1 = hefsi_surface_new(hefsi_func, data1, FALSE, 0, 1, 0, 1) ; */
/*   h2 = hefsi_surface_new(hefsi_func, data2, FALSE, 0, 1, 0, 1) ; */
/*   hefsi_surface_initialize(h1, dmin, dmax) ; */
/*   hefsi_set_bounding_boxes(h1, scale) ; */
/*   hefsi_surface_initialize(h2, dmin, dmax) ; */
/*   hefsi_set_bounding_boxes(h2, scale) ; */
/*   hefsi_surface_intersections(h1, h2, tol, wh) ; */

/*   for ( i = 0 ; i < wh->c->len ; i ++ ) { */
/*     for ( il = hefsi_workspace_curve(wh,i) ; il != NULL ; */
/* 	  il = il->next ) { */
/*       j = GPOINTER_TO_INT(il->data) ; */
/*       seg1 = hefsi_workspace_segment(wh,j) ; */
/*       fprintf(stdout, "%e %e %e\n", seg1->x1[0], seg1->x1[1], seg1->x1[2]) ; */
/*       /\* for ( j = 0 ; j < 7 ; j ++ ) { *\/ */
/*       /\* 	seg1->x1[j] = 	0.5*(seg1->x1[j] + seg1->x2[j]) ; *\/ */
/*       /\* } *\/ */
/*     } */
/*   } */
  
/*   return ; */
/* } */

static void affine_derivative_test(void)

{
  agg_affine_t *A ;
  agg_transform_t *T ;
  agg_variable_t global[4], params[32] ;
  gint nglobal, i, j ;
  char *args[32] ;
  gdouble p[32], u, du, T0[16], T1[16], dT[16] ;
  
  fprintf(stderr, "differentiation of affine matrices\n") ;
  fprintf(stderr, "==================================\n") ;

  T = agg_transform_new(64) ;

  nglobal = 2 ;
  global[0].name = g_strdup("twist") ;
  global[0].def =  g_strdup("15.0*pi/180.0") ;
  global[0].val = 0.4 ;
  global[1].name = g_strdup("span") ;
  global[1].def =  NULL ;
  global[1].val = 3.4 ;

  agg_transform_variable_add(T, "u", NULL, 0.0) ;
  
  agg_transform_add_global_variables(T, global, nglobal) ;
  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;
  agg_transform_variables_eval(T) ;

  A = agg_affine_new(3) ;

  args[0] = "twist*u" ; p[0] = 0 ;
  agg_affine_rotation_x(A, p, args, 1) ;
  agg_affine_differentiate(A, "u") ;
  agg_affine_expressions_compile(A, T->e) ;
  agg_transform_affine_add(T, A) ;

  A = agg_affine_new(3) ;

  args[0] = NULL ; p[0] = 0 ;
  args[1] = NULL ; p[1] = 0 ;
  args[2] = "span*u" ; p[2] = 0 ;
  agg_affine_translation(A, p, args, 3) ;
  agg_affine_differentiate(A, "u") ;
  agg_affine_expressions_compile(A, T->e) ;
  agg_transform_affine_add(T, A) ;
  
  u = 0.7 ; du = 1e-3 ;

  agg_transform_evaluate(T, u+0.5*du, 0, T1) ;
  agg_transform_evaluate(T, u-0.5*du, 0, T0) ;
  agg_transform_evaluate(T, u       , 1, dT) ;

  for ( i = 0 ; i < 16 ; i ++ ) T1[i] = (T1[i] - T0[i])/du ;
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    for ( j = 0 ; j < 4 ; j ++ ) {
      fprintf(stderr, "%lg (%lg) ", dT[4*i+j], T1[4*i+j]) ;
    }
    fprintf(stderr, "\n") ;
  }
  
  /* agg_affine_write(stderr, A) ; */
  
  return ;
}

static void section_interp_test(void)

{
  agg_section_t *s, *si, *ds, *s0, *s1 ;
  agg_surface_t *S ;
  gdouble u, du, err[32] = {0} ;
  gint i ;

  fprintf(stderr, "section interpolation test\n") ;
  fprintf(stderr, "==========================\n") ;
  
  S = agg_surface_new(64) ;

  ds = agg_section_new(8, 8) ;
  si = agg_section_new(8, 8) ;
  s0 = agg_section_new(8, 8) ;
  s1 = agg_section_new(8, 8) ;
  
  s = agg_section_new(8, 8) ;
  for ( u = 0 ; u <= 1.0 ; u += 1.0/16 ) {
    agg_section_set_aerofoil(s, 0.25 + 0.3*u, 0.2-0.1*u, 0.013*u) ;
    agg_surface_section_add(S, s, u) ;
  }

  agg_surface_weights_make(S) ;

  u = 5.0/16 ; du = 1e-6 ;
  agg_surface_section_diff(S, u, si, ds) ;
  agg_surface_section_interp(S, u+du/2, s1) ;
  agg_surface_section_interp(S, u-du/2, s0) ;

  for ( i = 0 ; i <= agg_section_order_upper(s1) ; i ++ ) {
    agg_section_coefficient_upper(s1,i) =
      (agg_section_coefficient_upper(s1,i) -
       agg_section_coefficient_upper(s0,i))/du ;
  }
  for ( i = 0 ; i <= agg_section_order_lower(s1) ; i ++ ) {
    agg_section_coefficient_lower(s1,i) =
      (agg_section_coefficient_lower(s1,i) -
       agg_section_coefficient_lower(s0,i))/du ;
  }
  agg_section_eta_left(s1)  = (agg_section_eta_left(s1) -
			       agg_section_eta_left(s0))/du ;
  agg_section_eta_right(s1) = (agg_section_eta_right(s1) -
			       agg_section_eta_right(s0))/du ;
  agg_section_trailing_edge_upper(s1) =
    (agg_section_trailing_edge_upper(s1) -
     agg_section_trailing_edge_upper(s0))/du ;
  agg_section_trailing_edge_lower(s1) =
    (agg_section_trailing_edge_lower(s1) -
     agg_section_trailing_edge_lower(s0))/du ;
  
  agg_surface_section_interp(S, u, s0) ;
  
  for ( i = 0 ; i <= agg_section_order_upper(si) ; i ++ ) {
    fprintf(stderr, "%lg %lg %lg %lg\n",
	    agg_section_coefficient_upper(si,i), 
	    agg_section_coefficient_upper(s0,i), 
	    agg_section_coefficient_upper(ds,i),
	    agg_section_coefficient_upper(s1,i)) ;
    err[0] = MAX(err[0],ABS(agg_section_coefficient_upper(si,i)-
			    agg_section_coefficient_upper(s0,i))) ;
    err[1] = MAX(err[1],ABS(agg_section_coefficient_upper(ds,i)-
			    agg_section_coefficient_upper(s1,i))) ;
  }
  for ( i = 0 ; i <= agg_section_order_lower(si) ; i ++ ) {
    fprintf(stderr, "%lg %lg %lg %lg\n",
	    agg_section_coefficient_lower(si,i), 
	    agg_section_coefficient_lower(s0,i), 
	    agg_section_coefficient_lower(ds,i),
	    agg_section_coefficient_lower(s1,i)) ; 
    err[2] = MAX(err[2],ABS(agg_section_coefficient_lower(si,i)-
			    agg_section_coefficient_lower(s0,i))) ;
    err[3] = MAX(err[3],ABS(agg_section_coefficient_lower(ds,i)-
			    agg_section_coefficient_lower(s1,i))) ;
  }

  fprintf(stderr, "%lg %lg %lg %lg\n",
	  agg_section_eta_left(si), agg_section_eta_left(s0), 
	  agg_section_eta_left(ds), agg_section_eta_left(s1)) ; 
  err[4] = MAX(err[4],ABS(agg_section_eta_left(si)-
			  agg_section_eta_left(s0))) ;
  err[5] = MAX(err[5],ABS(agg_section_eta_left(ds)-
			  agg_section_eta_left(s1))) ;
  fprintf(stderr, "%lg %lg %lg %lg\n",
	  agg_section_eta_right(si), agg_section_eta_right(s0), 
	  agg_section_eta_right(ds), agg_section_eta_right(s1)) ; 
  err[6] = MAX(err[6],ABS(agg_section_eta_right(si)-
			  agg_section_eta_right(s0))) ;
  err[7] = MAX(err[7],ABS(agg_section_eta_right(ds)-
			  agg_section_eta_right(s1))) ;

  fprintf(stderr, "%lg %lg %lg %lg\n",
	  agg_section_trailing_edge_upper(si), agg_section_trailing_edge_upper(s0), 
	  agg_section_trailing_edge_upper(ds), agg_section_trailing_edge_upper(s1)) ; 
  err[8] = MAX(err[8],ABS(agg_section_trailing_edge_upper(si)-
			  agg_section_trailing_edge_upper(s0))) ;
  err[9] = MAX(err[9],ABS(agg_section_trailing_edge_upper(ds)-
			  agg_section_trailing_edge_upper(s1))) ;
  fprintf(stderr, "%lg %lg %lg %lg\n",
	  agg_section_trailing_edge_lower(si), agg_section_trailing_edge_lower(s0), 
	  agg_section_trailing_edge_lower(ds), agg_section_trailing_edge_lower(s1)) ; 
  err[10] = MAX(err[10],ABS(agg_section_trailing_edge_lower(si)-
			  agg_section_trailing_edge_lower(s0))) ;
  err[11] = MAX(err[11],ABS(agg_section_trailing_edge_lower(ds)-
			  agg_section_trailing_edge_lower(s1))) ;
  
  fprintf(stderr, "errors:") ;
  for ( i = 0 ; i < 12 ; i ++ ) {
    fprintf(stderr, " %lg", err[i]) ;
  }
  fprintf(stderr, "\n") ;
  
  return ;
}

static void surface_diff_test(char *file)

{
  agg_body_t *b ;
  agg_surface_workspace_t *w ;
  agg_surface_t *S ;
  gint i ;
  gdouble u, v, del, x[3], dx[3], xu[3], xv[3], yu[3], yv[3] ;
  
  fprintf(stderr, "surface_differentiation test\n") ;
  fprintf(stderr, "============================\n") ;

  b = agg_body_new(32, 32) ;
  w = agg_surface_workspace_new() ;

  agg_body_read(b, file, FALSE) ;

  agg_body_globals_compile(b) ;
  agg_body_globals_eval(b) ;
  agg_body_globals_write(stderr, b) ;
  fprintf(stderr, "%d surface%s\n", agg_body_surface_number(b),
	  (agg_body_surface_number(b) > 1 ? "s" : "")) ;
  agg_body_surfaces_list(stderr, b) ;

  S = agg_body_surface(b, 1) ;

  u = 0.3 ; v = 0.0 ;

  agg_surface_point_diff(S, u, v, x, xu, xv, w) ;

  fprintf(stderr, "analytical\n") ;
  fprintf(stderr, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
  fprintf(stderr, "%lg %lg %lg\n", xu[0], xu[1], xu[2]) ;
  fprintf(stderr, "%lg %lg %lg\n", xv[0], xv[1], xv[2]) ;

  agg_surface_point_diff_numerical(S, u, v, x, xu, xv, w) ;

  fprintf(stderr, "numerical\n") ;
  fprintf(stderr, "%lg %lg %lg\n", x[0], x[1], x[2]) ;
  fprintf(stderr, "%lg %lg %lg\n", xu[0], xu[1], xu[2]) ;
  fprintf(stderr, "%lg %lg %lg\n", xv[0], xv[1], xv[2]) ;
  
  return ;
}

gint main(gint argc, char **argv)

{
  char ch, *progname, *file = "test.agg" ;
  gint test ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;
  test = -1 ;
  
  while ( (ch = getopt(argc, argv, "i:t:")) != EOF ) {
    switch (ch) {
    default: g_assert_not_reached() ; break ;
    case 'i': file = g_strdup(optarg) ; break ;
    case 't': test = parse_test(optarg) ; break ;
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

  /* if ( test == 3 ) { */
  /*   transform_test() ; */

  /*   return 0 ; */
  /* } */

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
  
  if ( test == 11 ) {
    affine_test() ;
    
    return 0 ;
  }

  /* if ( test == 12 ) { */
  /*   intersection_test() ; */
    
  /*   return 0 ; */
  /* } */

  if ( test == 10 ) {
    affine_derivative_test() ;

    return 0 ;
  }

  if ( test == 13 ) {
    section_interp_test() ;

    return 0 ;
  }

  if ( test == 14 ) {
    surface_diff_test(file) ;
    
    return 0 ;
  }
  
  return 0 ;
}

