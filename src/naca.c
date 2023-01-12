/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2021 Michael Carley
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

/*
 * Based on Ladson, Charles L., Brooks, Cuyler W., Jr, Hill, Acquilla
 * S., and Sproles, Darrell W., Computer Program to Obtain Ordinates
 * for NACA Airfoils, NASA Technical Memorandum 4741, 1996.
 */

gdouble agg_naca_four(gdouble t, gdouble p, gdouble m, gdouble x)
  
{
  gdouble y, sgn = 1 ;
  gdouble a[] = {0.2969, -0.1260, -0.3516, 0.2843, -0.1015} ;

  g_assert(x >= -1.0 && x <= 1.0) ;
  
  if ( x < 0 ) { x = -x ; sgn = -1 ; }

  y = a[0]*sqrt(x) ;
  y += x*(a[1] + x*(a[2] + x*(a[3] + x*a[4]))) ;
  
  y *= sgn*t/0.2 ;

  if ( p == 0.0 || m == 0.0 ) return y ;
  /*camber distribution*/
  if ( x <= m ) {
    y += p/m/m*(2*m*x - x*x) ;
  } else {
    y += p/(1.0-m)/(1.0-m)*((1.0-2*m) + 2*m*x - x*x) ;
  }
  
  return y ;
}

gint agg_shape_fit_naca_four(agg_shape_t *s, gint n,
			     gdouble th, gdouble p, gdouble m,
			     gboolean sharp,
			     gint nu, gint nl, gdouble *work)


{
  gdouble *xu, *yu, *xl, *yl, t, dy, y, xe ;
  gint i ;
  
  xu = work ; yu = &(xu[nu]) ;
  xl = &(yu[nu]) ; yl = &(xl[nl]) ;

  /*trailing edge constants from NASA TM4741, for extended trailing edge*/
  y = 0.002 ; dy = 0.234 ; xe = 1.0 + y/dy ;
  
  for ( i = 0 ; i < nu ; i ++ ) {
    t = 0.5*M_PI - 0.5*M_PI*i/(nu-1) ;
    xu[i] =  cos(t) ; yu[i] = agg_naca_four(th, p, m,  xu[i]) ;
  }
  if ( sharp ) {
    xu[nu-1] = xe ; yu[nu-1] = 0.0 ;
    for ( i = 0 ; i < nu ; i ++ ) {
      xu[i] /= xe ;
    }
  }
  
  for ( i = 0 ; i < nl ; i ++ ) {
    t = M_PI - 0.5*M_PI*i/(nl-1) ;
    xl[i] = -cos(t) ; yl[i] = agg_naca_four(th, p, m, -xl[i]) ;
  }
  xl[nl-1] = 0 ; yl[i] = agg_naca_four(th, p, m, -xl[i]) ;
  if ( sharp ) {
    xl[0] = xe ; yl[0] = 0.0 ;
    for ( i = 0 ; i < nl ; i ++ ) {
      xl[i] /= xe ;
    }
  }

  agg_shape_fit(s, xu, yu, nu, xl, yl, nl,
		0.5, 1.0, yu[nu-1], n, TRUE, &(yl[nl])) ;

  return 0 ;
}
