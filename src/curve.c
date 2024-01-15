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

#include <blaswrap.h>

#include "agg-private.h"

/** 
 * @{ 
 *
 * @ingroup curves
 */

static gint fourier_eval(gdouble *C, gint nc, gdouble x,
			 gdouble *s, gdouble *t)

{
  gint i ;

  *s = C[0] ;
  for ( i = 1 ; i <= nc ; i ++ ) {
    (*s) += C[2*i-1]*cos(-2.0*M_PI*i*x) + C[2*i+0]*sin(-2.0*M_PI*i*x) ;
  }
  *t = x ;
  
  return 0 ;
}

static gint polynomial_eval(gdouble *C, gint nc, gdouble x,
			    gdouble *s, gdouble *t)

{
  gint i ;

  *s = C[nc] ;
  
  for ( i = nc-1 ; i >= 0 ; i -- ) {
    g_assert_not_reached() ;
    *s = (*s)*x + C[i] ;
  }
  *t = x ;
  
  return 0 ;
}

static gint ellipse_eval(gdouble *C, gdouble x,
			 gdouble *s, gdouble *t)

{
  gdouble a, b, r, th ;

  a = C[2] ; b = C[3] ;
  th = 2.0*M_PI*x ;
  
  r = a*b/sqrt(b*b*cos(th)*cos(th) + a*a*sin(th)*sin(th)) ;

  *s = C[0] + r*cos(th) ;
  *t = C[1] + r*sin(th) ;
    
  return 0 ;
}

/** 
 * Evaluate coordinates on a parameterized curve
 * 
 * @param c ::agg_curve_t to be evaluated;
 * @param x parameter on \a c;
 * @param s on exit, set to \f$s\f$ in parametric plane;
 * @param t on exit, set to \f$t\f$ in parametric plane.
 * 
 * @return 0 on success.
 */

gint agg_curve_eval(agg_curve_t *c, gdouble x, gdouble *s, gdouble *t)

{
  if ( agg_curve_type(c) == AGG_CURVE_FOURIER )
    return fourier_eval(c->data, c->n, x, s, t) ;

  if ( agg_curve_type(c) == AGG_CURVE_POLYNOMIAL )
    return polynomial_eval(c->data, c->n, x, s, t) ;
  
  if ( agg_curve_type(c) == AGG_CURVE_ELLIPSE )
    return ellipse_eval(c->data, x, s, t) ;

  g_assert_not_reached() ;
  
  return 0 ;
}

static gboolean point_in_ellipse(agg_curve_t *c, gdouble del,
				 gdouble s, gdouble t)

{
  gdouble smax, smin, tmax, tmin ;

  smax = c->data[0] + c->data[2] ; 
  smin = c->data[0] - c->data[2] ; 
  tmax = c->data[1] + c->data[3] ; 
  tmin = c->data[1] - c->data[3] ; 

  if ( s >= smax + del ) return FALSE ;
  if ( s <= smin - del ) return FALSE ;
  if ( t >= tmax + del ) return FALSE ;
  if ( t <= tmin - del ) return FALSE ;
  
  return TRUE ;
}

gboolean agg_curve_point_orientation(agg_curve_t *c, gdouble del,
				     gdouble s, gdouble t)

{
  if ( agg_curve_type(c) == AGG_CURVE_ELLIPSE )
    return point_in_ellipse(c, del, s, t) ;
  
  g_assert_not_reached() ;
  
  return FALSE ;
}

/** 
 * Estimate a normal to the plane approximately containing a curve on
 * a surface. This is a rough estimate intended to be used in checking
 * curve orientations. 
 * 
 * @param c an ::agg_curve_t;
 * @param S surface on which curve is defined;
 * @param P patch for mapping \a c onto \a S;
 * @param n on exit contains estimated normal;
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */

gint agg_curve_plane_normal(agg_curve_t *c, agg_surface_t *S, agg_patch_t *P,
			    gdouble *n, agg_surface_workspace_t *w)

{
  gdouble x[9], s, t, u, v ;
  gint i ;
  const gdouble tn[] = {0.25, 0.5, 0.75} ;
  const gint nt = 3 ;

  for ( i = 0 ; i < nt ; i ++ ) {
    agg_curve_eval(c, tn[i], &s, &t) ;
    agg_patch_map(P, s, t, &u, &v) ;
    agg_surface_point_eval(S, u, v, &(x[3*i]), w) ;
  }

  agg_vector_diff(&(x[3*0]), &(x[3*0]), &(x[3*2])) ;
  agg_vector_diff(&(x[3*1]), &(x[3*1]), &(x[3*2])) ;
  agg_vector_cross(n, &(x[3*0]), &(x[3*1])) ;

  u = agg_vector_length(n) ;
  n[0] /= u ; n[1] /= u ; n[2] /= u ; 
  
  return 0 ;
}

/**
 * @}
 */
