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

#include <blaswrap.h>

#include "agg-private.h"

static gint fourier_eval(gdouble *C, gint nc, gdouble x,
			 gdouble *s, gdouble *t)

{
  gint i ;

  *s = C[0] ;
  for ( i = 1 ; i <= nc ; i ++ ) {
    (*s) += C[2*i-1]*cos(2.0*M_PI*i*x) + C[2*i+0]*sin(2.0*M_PI*i*x) ;
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
