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
#include <math.h>
#include <string.h>

#include <glib.h>

#include <agg.h>

#include "agg-private.h"

static const struct {
  gchar *name ;
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

agg_patch_t *agg_patch_new(gint nstmax)

{
  agg_patch_t *P ;

  P = (agg_patch_t *)g_malloc0(sizeof(agg_patch_t)) ;

  agg_patch_mapping(P) = AGG_PATCH_BILINEAR ;
  agg_patch_wrap_s(P) = FALSE ;
  agg_patch_wrap_t(P) = FALSE ;
  agg_patch_invert(P) = FALSE ;

  agg_patch_hole_number(P) = 0 ;
  agg_curve_type(agg_patch_curve_smin(P))  = AGG_CURVE_POLYNOMIAL ;
  agg_curve_order(agg_patch_curve_smin(P)) = 0 ;
  (agg_curve_data(agg_patch_curve_smin(P)))[0] = 0 ;
  agg_curve_type(agg_patch_curve_smax(P))  = AGG_CURVE_POLYNOMIAL ;
  agg_curve_order(agg_patch_curve_smax(P)) = 0 ;
  (agg_curve_data(agg_patch_curve_smax(P)))[0] = 1 ;
  
  return P ;
}

static gint agg_bilinear_map(agg_patch_t *P, gdouble s, gdouble t,
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
  
  *u = s ;
  *v = -1.0 + 2.0*t ;
  
  return 0 ;
}

static gint agg_spherical_map(agg_patch_t *P, gdouble s, gdouble t,
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


static gint agg_hemispherical_map(agg_patch_t *P, gdouble s, gdouble t,
				  gdouble *u, gdouble *v)

{
  agg_curve_t *cmin, *cmax ;
  gdouble smin, smax, tmp ;

  if ( t > 1.0 ) t -= 1.0 ;
  if ( t < 0.0 ) t += 1.0 ;
  g_assert(0.0 <= t && t <= 1.0) ;

  cmin = agg_patch_curve_smin(P) ;
  cmax = agg_patch_curve_smax(P) ;
  agg_curve_eval(cmin, t, &smin, &tmp) ;
  agg_curve_eval(cmax, t, &smax, &tmp) ;

  *u = sin(0.5*M_PI*(smin + s*(smax-smin))) ;

  if ( t < 0.5 ) {
    *v = -0.5*(1.0 + cos(2.0*M_PI*t)) ;
  } else {
    *v =  0.5*(1.0 + cos(2.0*M_PI*t)) ;
  }
  
  return 0 ;
}

static gint agg_tubular_map(agg_patch_t *P, gdouble s, gdouble t,
			    gdouble *u, gdouble *v)

{
  if ( t > 1.0 ) t -= 1.0 ;
  if ( t < 0.0 ) t += 1.0 ;
  g_assert(0.0 <= t && t <= 1.0) ;

  *u = s ;

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
 * @param P patch whose mapping is to be applied;
 * @param s variable on \a P, \f$0\leq s \leq 1\f$;
 * @param t variable on \a P, \f$0\leq t \leq 1\f$;
 * @param u surface parametric variable;
 * @param v surface parametric variable.
 * 
 * @return 0 on success.
 */

gint agg_patch_map(agg_patch_t *P, gdouble s, gdouble t,
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
    return agg_bilinear_map(P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_SPHERICAL ) {
    return agg_spherical_map(P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_HEMISPHERICAL ) {
    return agg_hemispherical_map(P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_TUBULAR ) {
    return agg_tubular_map(P, s, t, u, v) ;
  }
  
  g_error("%s: unrecognized mapping %d", __FUNCTION__, agg_patch_mapping(P)) ;
  
  return 0 ;
}

static gint wrapping_parse(gchar *str, gboolean *wrap)

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

#if 0
/** 
 * Evaluate patch coordinates on a clipping curve
 * 
 * @param c contains patch clipping;
 * @param u parameter on \a c;
 * @param s on exit contains \f$s(u)\f$ on patch;
 * @param t on exit contains \f$t(u)\f$ on patch.
 * 
 * @return 0 on success.
 */

gint agg_patch_clip_eval(agg_patch_clipping_t *c, gdouble u,
			 gdouble *s, gdouble *t)

{
  g_assert(agg_patch_clipping_orientation(c) ==  1 ||
	   agg_patch_clipping_orientation(c) == -1) ;
  
  if ( agg_patch_clipping_type(c) == AGG_CLIP_CONSTANT_S ) {
    *t = u ; *s = agg_patch_clipping_data(c,0) ;

    return 0 ;
  }

  if ( agg_patch_clipping_type(c) == AGG_CLIP_CONSTANT_T ) {
    *t = agg_patch_clipping_data(c,0) ; *s = u ;

    return 0 ;
  }

  if ( agg_patch_clipping_type(c) == AGG_CLIP_ELLIPSE ) {
    gdouble th = agg_patch_clipping_orientation(c)*2.0*M_PI*u ;

    *s = agg_patch_clipping_data(c,0) + agg_patch_clipping_data(c,2)*cos(th) ;
    *t = agg_patch_clipping_data(c,1) + agg_patch_clipping_data(c,3)*sin(th) ;
      
    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}
#endif

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

#if 0
/** 
 * Set orientation of two clipping curves so that they are traversed
 * in the same sense in physical space. The test is performed by
 * evaluating approximate normals to the curve planes (this assumes
 * that the curves are not too far from planar) and checking the
 * normals are approximately parallel, using the sign of the scalar
 * product. If necessary, the orientation of one curve is switched,
 * with priority given to changing the orientation of a closed curve,
 * corresponding to a hole in a surface, over a cut on constant
 * \f$s\f$ or \f$t\f$. If curves are not mutually oriented and are
 * both constant parameter cuts, execution stops with an error (until
 * I can decide how best to handle this case).
 * 
 * @param c1 first patch clipping;
 * @param P1 first parametric patch;
 * @param S1 first physical surface;
 * @param c2 second patch clipping;  
 * @param P2 second parametric patch;
 * @param S2 second physical surface;
 * @param w workspace for surface point evaluation.
 * 
 * @return 0 on success.
 */

gint agg_clipping_orient(agg_patch_clipping_t *c1, agg_patch_t *P1,
			 agg_surface_t *S1,
			 agg_patch_clipping_t *c2, agg_patch_t *P2,
			 agg_surface_t *S2,
			 agg_surface_workspace_t *w)

{
  gdouble u, v, s, t, t0[] = {0.25, 0.5, 0.75}, x1[9], x2[9], N1[3], N2[3] ;
  gint i ;

  for ( i = 0 ; i < 3 ; i ++ ) {
    agg_patch_clip_eval(c1, t0[i], &s, &t) ;
    agg_patch_map(P1, s, t, &u, &v) ;
    agg_surface_point_eval(S1, u, v, &(x1[3*i]), w) ;

    agg_patch_clip_eval(c2, t0[i], &s, &t) ;
    agg_patch_map(P2, s, t, &u, &v) ;
    agg_surface_point_eval(S2, u, v, &(x2[3*i]), w) ;
  }

  agg_vector_diff(&(x1[3*0]), &(x1[3*0]), &(x1[3*2])) ;
  agg_vector_diff(&(x1[3*1]), &(x1[3*1]), &(x1[3*2])) ;
  agg_vector_diff(&(x2[3*0]), &(x2[3*0]), &(x2[3*2])) ;
  agg_vector_diff(&(x2[3*1]), &(x2[3*1]), &(x2[3*2])) ;

  agg_vector_cross(N1, &(x1[3*0]), &(x1[3*1])) ;
  agg_vector_cross(N2, &(x2[3*0]), &(x2[3*1])) ;

  /*normals are (about) the same direction*/
  if ( agg_vector_scalar(N1, N2) > 0 ) return 0 ;
  
  if ( agg_patch_clipping_type(c1) == AGG_CLIP_ELLIPSE ) {
    c1->ornt = -(c1->ornt) ;
    return 0 ;
  }

  if ( agg_patch_clipping_type(c2) == AGG_CLIP_ELLIPSE ) {
    c2->ornt = -(c2->ornt) ;
    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}
#endif

gint agg_patch_point_diff(agg_surface_t *S, agg_patch_t *P,
			  gdouble s, gdouble t,
			  gdouble *x, gdouble *xs, gdouble *xt,
			  agg_surface_workspace_t *w)

{
  gdouble ee, u, v ;

  ee = 1e-6 ;

  agg_patch_map(P, s, t, &u, &v) ;
  agg_surface_point_eval(S, u, v, x, w) ;

  if ( s > ee ) {
    agg_patch_map(P, s-ee, t, &u, &v) ;
    agg_surface_point_eval(S, u, v, xs, w) ;
    agg_vector_diff(xs, x, xs) ;
  } else {
    agg_patch_map(P, s+ee, t, &u, &v) ;
    agg_surface_point_eval(S, u, v, xs, w) ;
    agg_vector_diff(xs, xs, x) ;
  }

  xs[0] /= ee ; xs[1] /= ee ; xs[2] /= ee ; 

  if ( t > ee ) {
    agg_patch_map(P, s, t-ee, &u, &v) ;
    agg_surface_point_eval(S, u, v, xt, w) ;
    agg_vector_diff(xt, x, xt) ;
  } else {
    agg_patch_map(P, s, t+ee, &u, &v) ;
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

gint agg_patch_triangle_interp(agg_patch_t *P,
			       gdouble s0, gdouble t0,
			       gdouble s1, gdouble t1, 
			       gdouble s2, gdouble t2, 
			       gdouble u,  gdouble v,
			       gdouble *s, gdouble *t)

{
  switch ( agg_patch_mapping(P) ) {
  default: g_assert_not_reached() ; break ;
  case AGG_PATCH_SPHERICAL:
    gdouble x0[3], x1[3], x2[3], x[3], L[3] ;

    spherical_st_to_x(s0, t0, x0) ;
    spherical_st_to_x(s1, t1, x1) ;
    spherical_st_to_x(s2, t2, x2) ;

    shapefunc(u, v, L) ;
    x[0] = L[0]*x0[0] + L[1]*x1[0] + L[2]*x2[0] ; 
    x[1] = L[0]*x0[1] + L[1]*x1[1] + L[2]*x2[1] ; 
    x[2] = L[0]*x0[2] + L[1]*x1[2] + L[2]*x2[2] ; 

    spherical_x_to_st(x, s, t) ;

    break ;
  }
  
  return 0 ;
}


/**
 * @}
 */
