/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2022 Michael Carley
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

#include <blaswrap.h>

#include <agg.h>

#include "agg-private.h"

agg_shape_interpolant_t *agg_shape_interpolant_alloc(gint nmax)

{
  agg_shape_interpolant_t *a ;

  a = (agg_shape_interpolant_t *)g_malloc0(sizeof(agg_shape_interpolant_t)) ;

  a->nmax = nmax ;
  a->i = (gdouble *)g_malloc0(4*nmax*sizeof(gdouble)) ;
  
  return a ;
}

gint agg_shape_interpolant_make(agg_shape_interpolant_t *a,
				agg_shape_t *sh, gint n)

{
  gint i ;
  gdouble s ;
  
  g_assert(n <= a->nmax) ;

  i = 0 ;
  s = agg_spacing_eval(-1.0, 1.0, n, AGG_SPACING_LINEAR, i) ;
  a->i[4*i+0] = s ;
  a->i[4*i+1] = fabs(s) ;
  a->i[4*i+2] = agg_shape_eval(sh, s, 0) ;
  a->i[4*i+3] = 0.0 ;
  
  for ( i = 1 ; i < n ; i ++ ) {
    s = agg_spacing_eval(-1.0, 1.0, n, AGG_SPACING_LINEAR, i) ;
    a->i[4*i+0] = s ;
    a->i[4*i+1] = fabs(s) ;
    a->i[4*i+2] = agg_shape_eval(sh, s, 0) ;
    a->i[4*i+3] = a->i[4*(i-1)+3] +
      sqrt((a->i[4*i+1] - a->i[4*(i-1)+1])*(a->i[4*i+1] - a->i[4*(i-1)+1]) +
	   (a->i[4*i+2] - a->i[4*(i-1)+2])*(a->i[4*i+2] - a->i[4*(i-1)+2])) ;
  }

  a->n = n ;
  
  return 0 ;
}

static gint bsearch_array(gdouble *x, gint str, gint n, gdouble y)

{
  gint i, j, m ;

  i = 0 ; j = n ; /* - 1 ; */
  while ( i < j ) {
    m = (gint)floor(0.5*(i+j)) ;
    if ( x[m*str] <= y && x[(m+1)*str] > y) return m ;
    if ( x[m*str] < y ) i = m + 1 ;
    else j = m ;
  }

  return i ;
}

static gdouble interp_linear(gdouble x0, gdouble y0,
			     gdouble x1, gdouble y1,
			     gdouble x)

{
  gdouble y ;

  y = (x - x0)/(x1 - x0) ;
  y = y0 + y*(y1-y0) ;
  
  return y ;
}

gint agg_shape_interpolant_step(agg_shape_interpolant_t *a,
				gdouble *s, gdouble dl)

/*
 * estimate value of s coordinate to step a distance dl along curve
 */

{
  gdouble len ;
  gint i ;

  if ( *s >= a->i[4*(a->n-1)] ) return 1 ;
  
  /*locate s and interpolate l*/
  i = bsearch_array(a->i, 4, a->n, *s) ;

  g_assert(a->i[4*i] <= *s && a->i[4*(i+1)] >= *s) ;

  len = interp_linear(a->i[4* i   +0], a->i[4* i   +3],
		      a->i[4*(i+1)+0], a->i[4*(i+1)+3],
		      *s) ;
  len += dl ;
  
  /* len = a->i[4*i+3] + dl ; */
  while ( a->i[4*i+3] < len && i < a->n ) i ++ ;

  if ( i >= a->n-1 ) return 1.1 ;
  
  *s = interp_linear(a->i[4* i   +3], a->i[4* i   +0],
		     a->i[4*(i+1)+3], a->i[4*(i+1)+0],
		     len) ;
  
  return 0 ;
}
