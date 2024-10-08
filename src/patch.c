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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <glib.h>

#include <agg.h>

#include "agg-private.h"

static const struct {
  char *name ;
  agg_patch_mapping_t map ;
} mapping_list[] =
  {
    {"bilinear",      AGG_PATCH_BILINEAR},
    {"spherical",     AGG_PATCH_SPHERICAL},
    {"hemispherical", AGG_PATCH_HEMISPHERICAL},
    {"tubular",       AGG_PATCH_TUBULAR},
    {NULL,        -1}
  } ;

static void shapefunc(gdouble s, gdouble t, gdouble L[])

{
  L[0] = 1.0 - s - t ;
  L[1] =       s     ;
  L[2] =           t ;
    
  return ;
}

/** 
 * @{ 
 *
 * @ingroup patches
 */

/** 
 * Allocate a new patch for mapping of surface parameters
 * 
 * @return newly allocated ::agg_patch_t
 */

agg_patch_t *agg_patch_new(void)

{
  agg_patch_t *P ;

  P = (agg_patch_t *)g_malloc0(sizeof(agg_patch_t)) ;

  agg_patch_mapping(P) = AGG_PATCH_BILINEAR ;
  agg_patch_wrap_s(P) = FALSE ;
  agg_patch_wrap_t(P) = FALSE ;
  agg_patch_invert(P) = FALSE ;

  agg_patch_hole_number(P) = 2 ;
  agg_curve_type(agg_patch_curve_smin(P))  = AGG_CURVE_POLYNOMIAL ;
  agg_curve_order(agg_patch_curve_smin(P)) = 0 ;
  (agg_curve_data(agg_patch_curve_smin(P)))[0] = 0 ;
  agg_curve_type(agg_patch_curve_smax(P))  = AGG_CURVE_POLYNOMIAL ;
  agg_curve_order(agg_patch_curve_smax(P)) = 0 ;
  (agg_curve_data(agg_patch_curve_smax(P)))[0] = 1 ;
  
  return P ;
}

static gint agg_bilinear_map(agg_surface_t *S,
			     agg_patch_t *P, gdouble s, gdouble t,
			     gdouble *u, gdouble *v)
{
  /* gdouble w0, w1, w2, w3 ; */

  /* /\*basic bilinear mapping for now */
  /*   https://en.wikipedia.org/wiki/Bilinear_interpolation */
  /* *\/ */
  /* w0 = (1.0-s)*(1.0-t) ; */
  /* w1 = (    s)*(1.0-t) ; */
  /* w2 = (    s)*(    t) ; */
  /* w3 = (1.0-s)*(    t) ; */
  
  /* *u = */
  /*   agg_patch_corner_u(P,0)*w0 + */
  /*   agg_patch_corner_u(P,1)*w1 + */
  /*   agg_patch_corner_u(P,2)*w2 + */
  /*   agg_patch_corner_u(P,3)*w3 ; */
  /* *v = */
  /*   agg_patch_corner_v(P,0)*w0 + */
  /*   agg_patch_corner_v(P,1)*w1 + */
  /*   agg_patch_corner_v(P,2)*w2 + */
  /*   agg_patch_corner_v(P,3)*w3 ; */

  g_assert_not_reached() ; /*untested code*/
  
  *u = agg_surface_umin(S) + s*(agg_surface_umax(S) - agg_surface_umin(S)) ;
  *v = -1.0 + 2.0*t ;
  
  return 0 ;
}

static gint agg_spherical_map(agg_surface_t *S,
			      agg_patch_t *P, gdouble s, gdouble t,
			      gdouble *u, gdouble *v)

{

  if ( t > 1.0 ) t -= 1.0 ;
  if ( t < 0.0 ) t += 1.0 ;
  g_assert(0.0 <= t && t <= 1.0) ;

  *u = 0.5*(1.0-cos(M_PI*s)) ;
  
  if ( t < 0.5 ) {
    *v = -0.5*(1.0 + cos(2.0*M_PI*t)) ;
  } else {
    *v =  0.5*(1.0 + cos(2.0*M_PI*t)) ;
  }
  
  return 0 ;
}

static gint agg_hemispherical_map(agg_surface_t *S,
				  agg_patch_t *P, gdouble s, gdouble t,
				  gdouble *u, gdouble *v)

{
  agg_curve_t *cmin, *cmax ;
  gdouble smin, smax, tmp, ut ;

  if ( t > 1.0 ) t -= 1.0 ;
  if ( t < 0.0 ) t += 1.0 ;
  g_assert(0.0 <= t && t <= 1.0) ;

  cmin = agg_patch_curve_smin(P) ;
  cmax = agg_patch_curve_smax(P) ;
  agg_curve_eval(cmin, t, &smin, &tmp) ;
  agg_curve_eval(cmax, t, &smax, &tmp) ;

  /* smin = 0.0 ; smax = 1.0 ; */
  
  ut = sin(0.5*M_PI*(smin + s*(smax-smin))) ;

  *u = agg_surface_umin(S) + ut*(agg_surface_umax(S) - agg_surface_umin(S)) ;

  
  if ( t < 0.5 ) {
    *v = -0.5*(1.0 + cos(2.0*M_PI*t)) ;
  } else {
    *v =  0.5*(1.0 + cos(2.0*M_PI*t)) ;
  }
  
  return 0 ;
}

static gint agg_tubular_map(agg_surface_t *S,
			    agg_patch_t *P, gdouble s, gdouble t,
			    gdouble *u, gdouble *v)

{
  if ( t > 1.0 ) t -= 1.0 ;
  if ( t < 0.0 ) t += 1.0 ;
  g_assert(0.0 <= t && t <= 1.0) ;

  *u = agg_surface_umin(S) + s*(agg_surface_umax(S) - agg_surface_umin(S)) ;
  /* *u = s ; */

  if ( t < 0.5 ) {
    *v = -0.5*(1.0 + cos(2.0*M_PI*t)) ;
  } else {
    *v =  0.5*(1.0 + cos(2.0*M_PI*t)) ;
  }
  
  return 0 ;
}

/** 
 * Apply mapping in a patch to find parametric variables on surface.
 * 
 * @param S surface to which P is applied;
 * @param P patch whose mapping is to be applied;
 * @param s variable on \a P, \f$0\leq s \leq 1\f$;
 * @param t variable on \a P, \f$0\leq t \leq 1\f$;
 * @param u surface parametric variable on \a S;
 * @param v surface parametric variable on \a S.
 * 
 * @return 0 on success.
 */

gint agg_patch_map(agg_surface_t *S,
		   agg_patch_t *P,
		   gdouble s, gdouble t,
		   gdouble *u, gdouble *v)

{
  if ( s < 0 && agg_patch_wrap_s(P) ) { s = 1.0 + s ; }
  if ( s > 1 && agg_patch_wrap_s(P) ) { s = s - 1.0; }
  if ( t < 0 && agg_patch_wrap_t(P) ) { t = 1.0 + t ; }
  if ( t > 1 && agg_patch_wrap_t(P) ) { t = t - 1.0; }

  if ( s < 0 || s > 1 )
    g_error("%s: s (%lg) out of range", __FUNCTION__, s) ;
  if ( t < 0 || t > 1 )
    g_error("%s: t (%lg) out of range", __FUNCTION__, t) ;

  if ( agg_patch_mapping(P) == AGG_PATCH_BILINEAR ) {
    return agg_bilinear_map(S, P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_SPHERICAL ) {
    return agg_spherical_map(S, P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_HEMISPHERICAL ) {
    return agg_hemispherical_map(S, P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_TUBULAR ) {
    return agg_tubular_map(S, P, s, t, u, v) ;
  }
  
  g_error("%s: unrecognized mapping %d", __FUNCTION__, agg_patch_mapping(P)) ;
  
  return 0 ;
}

static gint wrapping_parse(char *str, gboolean *wrap)

{
  if ( strcmp("wrap", str) == 0 ) { *wrap = TRUE ; return 0 ; }
  if ( strcmp("nowrap", str) == 0 ) { *wrap = FALSE ; return 0 ; }

  return -1 ;
}

/** 
 * Parse string data for an ::agg_patch_t.
 * 
 * @param P on exit, is set according to input variables;
 * @param p array of variables for parameters setting \a P;
 * @param np number of variables in \a p.
 * 
 * @return 0 on success.
 */

gint agg_patch_parse(agg_patch_t *P, agg_variable_t *p, gint np)

{
  gboolean wrap_s, wrap_t, invert ;
  gint i ;

  if ( np < 3 ) {
    g_error("%s: at least three parameters required for mapping",
	    __FUNCTION__) ;
  }

  invert = FALSE ;
  if ( wrapping_parse(agg_variable_definition(&(p[1])), &wrap_s) != 0 ) {
    g_error("%s: unrecognized mapping \"%s\"",
	    __FUNCTION__, agg_variable_definition(&(p[1]))) ;
  }

  if ( wrapping_parse(agg_variable_definition(&(p[2])), &wrap_t) != 0 ) {
    g_error("%s: unrecognized mapping \"%s\"",
	    __FUNCTION__, agg_variable_definition(&(p[2]))) ;
  }

  if ( np > 3 ) {
    if ( strcmp(agg_variable_definition(&(p[3])), "invert") == 0 ) {
      invert = TRUE ;
    } else {
      g_error("%s: unrecognized inversion parameter \"%s\"",
	      __FUNCTION__, agg_variable_definition(&(p[3]))) ;
    }
  }
  
  for ( i = 0 ; mapping_list[i].name != NULL ; i ++ ) {
    if ( strcmp(mapping_list[i].name, agg_variable_definition(&(p[0]))) == 0 ) {
      agg_patch_mapping(P) = mapping_list[i].map ;
      agg_patch_wrap_s(P)  = wrap_s ;
      agg_patch_wrap_t(P)  = wrap_t ;
      agg_patch_invert(P)  = invert ;
      return 0 ;
    }
  }

  g_error("%s: unrecognized mapping \"%s\"", __FUNCTION__,
	  agg_variable_definition(&(p[0]))) ;
  
  return 0 ;
}


/** 
 * Convert derivatives with respect to surface coordinates \f$(u,v)\f$
 * to derivatives with respect to patch variables \f$(s,t)\f$.
 * 
 * @param P surface patch;
 * @param s parameter on \a P;
 * @param t parameter on \a P;
 * @param xu \f$\partial\mathbf{x}/\partial u\f$ at \f$(u(s,t),v(s,t))\f$;
 * @param xv \f$\partial\mathbf{x}/\partial v\f$ at \f$(u(s,t),v(s,t))\f$;
 * @param xs on exit contains \f$\partial\mathbf{x}/\partial s\f$;
 * @param xt on exit contains \f$\partial\mathbf{x}/\partial t\f$.
 * 
 * @return 0 on success.
 */

gint agg_patch_surface_diff(agg_patch_t *P,
			    gdouble s, gdouble t,
			    gdouble *xu, gdouble *xv,
			    gdouble *xs, gdouble *xt)

{
  gdouble xtu[]={xu[0], xu[1], xu[2]}, xtv[]={xv[0], xv[1], xv[2]} ;
  gdouble us, ut, vs, vt ;

  us = ut = vs = vt = 0 ;
  switch ( agg_patch_mapping(P) ) {
  default: g_assert_not_reached() ; break ;
  case AGG_PATCH_BILINEAR:
    us = 1.0 ; vt = 2.0 ;
    break ;
  case AGG_PATCH_SPHERICAL:
    us = 0.5*M_PI*sin(M_PI*s) ;
    if ( t < 0.5 ) {
      vt = 0.5*2.0*M_PI*(sin(2.0*M_PI*t)) ;
    } else {
      vt =  -0.5*2.0*M_PI*(sin(2.0*M_PI*t)) ;
    }
    break ;
  case AGG_PATCH_TUBULAR:
    us = 1.0 ;
    if ( t < 0.5 ) {
      vt = 0.5*2.0*M_PI*(sin(2.0*M_PI*t)) ;
    } else {
      vt =  -0.5*2.0*M_PI*(sin(2.0*M_PI*t)) ;
    }
    break ;
  }

  xs[0] = xtu[0]*us + xtv[0]*vs ;
  xs[1] = xtu[1]*us + xtv[1]*vs ;
  xs[2] = xtu[2]*us + xtv[2]*vs ;

  xt[0] = xtu[0]*ut + xtv[0]*vt ;
  xt[1] = xtu[1]*ut + xtv[1]*vt ;
  xt[2] = xtu[2]*ut + xtv[2]*vt ;
  
  return 0 ;
}

/** 
 * Estimate the derivatives on a surface with respect to parametric
 * coordinates \f$(s,t)\f$.
 * 
 * @param S surface to differentiate;
 * @param P mapping patch;
 * @param s coordinate on \a P;
 * @param t coordinate on \a P;
 * @param x on exit, contains point \f$\mathbf{x}(u(s,t),v(s,t))\f$;
 * @param xs on exit, contains 
 * \f$\partial\mathbf{x}(u(s,t),v(s,t))/\partial s\f$;
 * @param xt on exit, contains 
 * \f$\partial\mathbf{x}(u(s,t),v(s,t))/\partial t\f$;
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */

gint agg_patch_point_diff(agg_surface_t *S, agg_patch_t *P,
			  gdouble s, gdouble t,
			  gdouble *x, gdouble *xs, gdouble *xt,
			  agg_surface_workspace_t *w)

{
  gdouble ee, u, v ;

  ee = 1e-3 ;

  agg_patch_map(S, P, s, t, &u, &v) ;
  agg_surface_point_eval(S, u, v, x, w) ;

  if ( s > ee ) {
    agg_patch_map(S, P, s-ee, t, &u, &v) ;
    agg_surface_point_eval(S, u, v, xs, w) ;
    agg_vector_diff(xs, x, xs) ;
  } else {
    agg_patch_map(S, P, s+ee, t, &u, &v) ;
    agg_surface_point_eval(S, u, v, xs, w) ;
    agg_vector_diff(xs, xs, x) ;
  }

  xs[0] /= ee ; xs[1] /= ee ; xs[2] /= ee ; 

  if ( t > ee ) {
    agg_patch_map(S, P, s, t-ee, &u, &v) ;
    agg_surface_point_eval(S, u, v, xt, w) ;
    agg_vector_diff(xt, x, xt) ;
  } else {
    agg_patch_map(S, P, s, t+ee, &u, &v) ;
    agg_surface_point_eval(S, u, v, xt, w) ;
    agg_vector_diff(xt, xt, x) ;
  }

  xt[0] /= ee ; xt[1] /= ee ; xt[2] /= ee ; 
  
  return 0 ;
}

static void spherical_st_to_x(gdouble s, gdouble t, gdouble *x)

{
  gdouble th, ph ;

  if ( t > 1 ) t -= 1 ;
  if ( t < 0 ) t += 1 ;
  
  th = (1.0 - t)*2.0*M_PI ; ph = (1.0 - s)*M_PI ;
  x[0] = cos(th)*sin(ph) ; x[1] = sin(th)*sin(ph) ; x[2] = cos(ph) ;

  return ;
}

static void spherical_x_to_st(gdouble *x, gdouble *s, gdouble *t)

{
  gdouble r ;

  r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) ;
  *t = atan2(x[1], x[0]) ;
  if ( (*t) < 0 ) (*t) += 2.0*M_PI ;
  *s = 1.0 - acos(x[2]/r)/M_PI ;
  *t = 1.0 - (*t)*0.5/M_PI ;

  return ;
}

static void hemispherical_st_to_x(gdouble s, gdouble t, gdouble *x)

{
  gdouble th, ph ;

  if ( t > 1 ) t -= 1 ;
  if ( t < 0 ) t += 1 ;
  
  th = (1.0 - t)*2.0*M_PI ; ph = 0.5*(1.0 - s)*M_PI ;
  x[0] = cos(th)*sin(ph) ; x[1] = sin(th)*sin(ph) ; x[2] = cos(ph) ;

  return ;
}

static void hemispherical_x_to_st(gdouble *x, gdouble *s, gdouble *t)

{
  gdouble r ;

  r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) ;
  *t = atan2(x[1], x[0]) ;
  if ( (*t) < 0 ) (*t) += 2.0*M_PI ;
  *s = 1.0 - 2.0*acos(x[2]/r)/M_PI ;
  *t = 1.0 - (*t)*0.5/M_PI ;

  return ;
}

/** 
 * Split an edge, the line segment joining two points in \f$(s,t)\f$
 * parametric plane, remapping as required to respect the underlying
 * mapping of the patch. For a linear patch
 * \f$s=s_{0}+a(s_{1}-s_{0})\f$, \f$t=t_{0}+a(t_{1}-t_{0})\f$; for
 * other mappings, such as spherical or hemispherical, the
 * interpolation is different.
 * 
 * @param P surface mapping patch;
 * @param s0 coordinate of first point on segment;
 * @param t0 coordinate of first point on segment;
 * @param s1 coordinate of second point on segment;
 * @param t1 coordinate of second point on segment;
 * @param a parametric distance on segment;
 * @param s on output, coordinate of interpolated point on segment;
 * @param t on output, coordinate of interpolated point on segment.
 * 
 * @return 0 on success.
 */

gint agg_patch_edge_split(agg_patch_t *P,
			  gdouble s0, gdouble t0, gdouble s1, gdouble t1, 
			  gdouble a, gdouble *s, gdouble *t)

{
  gdouble x0[3], x1[3], x[3] ;

  switch ( agg_patch_mapping(P) ) {
  default:
    g_assert_not_reached() ;
    break ;
  case AGG_PATCH_SPHERICAL:
    spherical_st_to_x(s0, t0, x0) ;
    spherical_st_to_x(s1, t1, x1) ;

    x[0] = x0[0] + a*(x1[0] - x0[0]) ;
    x[1] = x0[1] + a*(x1[1] - x0[1]) ;
    x[2] = x0[2] + a*(x1[2] - x0[2]) ;

    spherical_x_to_st(x, s, t) ;

    if ( t0 > 1 || t1 > 1 ) (*t) += 1 ;
    if ( t0 < 0 || t1 < 0 ) (*t) -= 1 ;
    
    break ;
  case AGG_PATCH_HEMISPHERICAL:
    hemispherical_st_to_x(s0, t0, x0) ;
    hemispherical_st_to_x(s1, t1, x1) ;

    x[0] = x0[0] + a*(x1[0] - x0[0]) ;
    x[1] = x0[1] + a*(x1[1] - x0[1]) ;
    x[2] = x0[2] + a*(x1[2] - x0[2]) ;

    hemispherical_x_to_st(x, s, t) ;

    if ( t0 > 1 || t1 > 1 ) (*t) += 1 ;
    if ( t0 < 0 || t1 < 0 ) (*t) -= 1 ;
    
    break ;
  }
  
  return 0 ;
}

/** 
 * Split a triangle in the parametric plane, remapping as required to
 * respect the underlying mapping of the patch.
 * 
 * @param P surface mapping patch;
 * @param s0 coordinate of first point on triangle;
 * @param t0 coordinate of first point on triangle;
 * @param s1 coordinate of second point on triangle;
 * @param t1 coordinate of second point on triangle;
 * @param s2 coordinate of third point on triangle;
 * @param t2 coordinate of third point on triangle;
 * @param u coordinate on triangle;
 * @param v coordinate on triangle;
 * @param s on output, coordinate of interpolated point on triangle;
 * @param t on output, coordinate of interpolated point on triangle.
 * 
 * @return 0 on success.
 */

gint agg_patch_triangle_interp(agg_patch_t *P,
			       gdouble s0, gdouble t0,
			       gdouble s1, gdouble t1, 
			       gdouble s2, gdouble t2, 
			       gdouble u,  gdouble v,
			       gdouble *s, gdouble *t)

{
  gdouble x0[3], x1[3], x2[3], x[3], L[3] ;

  switch ( agg_patch_mapping(P) ) {
  default: g_assert_not_reached() ; break ;
  case AGG_PATCH_SPHERICAL:
    spherical_st_to_x(s0, t0, x0) ;
    spherical_st_to_x(s1, t1, x1) ;
    spherical_st_to_x(s2, t2, x2) ;

    shapefunc(u, v, L) ;
    x[0] = L[0]*x0[0] + L[1]*x1[0] + L[2]*x2[0] ; 
    x[1] = L[0]*x0[1] + L[1]*x1[1] + L[2]*x2[1] ; 
    x[2] = L[0]*x0[2] + L[1]*x1[2] + L[2]*x2[2] ; 

    spherical_x_to_st(x, s, t) ;

    break ;
  case AGG_PATCH_HEMISPHERICAL:
    hemispherical_st_to_x(s0, t0, x0) ;
    hemispherical_st_to_x(s1, t1, x1) ;
    hemispherical_st_to_x(s2, t2, x2) ;

    shapefunc(u, v, L) ;
    x[0] = L[0]*x0[0] + L[1]*x1[0] + L[2]*x2[0] ; 
    x[1] = L[0]*x0[1] + L[1]*x1[1] + L[2]*x2[1] ; 
    x[2] = L[0]*x0[2] + L[1]*x1[2] + L[2]*x2[2] ; 

    hemispherical_x_to_st(x, s, t) ;

    break ;
  }
  
  return 0 ;
}

/**
 * @}
 */
