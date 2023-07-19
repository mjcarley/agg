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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

static void add_point(gdouble *x, gdouble *th, gdouble *ph)

{
  *ph = acos(x[2]) ;
  *th = atan2(x[1], x[0]) ;
  return ;
}

gint agg_sphere_ico_base(gdouble *th, gint tstr,
			 gdouble *ph, gint pstr,
			 gint *e, gint estr, gboolean convert)

/*
 * nodes and edges of regular icosahedron on a unit sphere
 * 
 * 12 nodes, 30 edges
 * 
 * data from https://en.wikipedia.org/wiki/Regular_icosahedron
 */
  
{
  gint i ;
  gdouble x[3], sq5 ;

  sq5 = 1.0/sqrt(5.0) ;

  i = 0 ;
  x[0] =  1.0 ; x[1] = 0.0 ; x[2] = 0.0 ;
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  x[0] = -x[0] ; x[1] = -x[1] ; x[2] = -x[2] ; 
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  
  x[0] = sq5 ; x[1] = 2.0*sq5 ; x[2] = 0.0 ; 
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  x[0] = -x[0] ; x[1] = -x[1] ; x[2] = -x[2] ; 
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;

  x[0] = sq5 ; x[1] = (1.0-sq5)/2.0 ; x[2] = sqrt((1.0+sq5)/2.0) ;
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  x[0] = -x[0] ; x[1] = -x[1] ; x[2] = -x[2] ; 
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;

  x[0] = sq5 ; x[1] = (1.0-sq5)/2.0 ; x[2] = -sqrt((1.0+sq5)/2.0) ;
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  x[0] = -x[0] ; x[1] = -x[1] ; x[2] = -x[2] ; 
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;

  x[0] = sq5 ; x[1] = (-1.0-sq5)/2.0 ; x[2] = sqrt((1.0-sq5)/2.0) ;
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  x[0] = -x[0] ; x[1] = -x[1] ; x[2] = -x[2] ; 
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  
  x[0] = sq5 ; x[1] = (-1.0-sq5)/2.0 ; x[2] = -sqrt((1.0-sq5)/2.0) ;
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;
  x[0] = -x[0] ; x[1] = -x[1] ; x[2] = -x[2] ; 
  add_point(x, &(th[i*tstr]), &(ph[i*pstr])) ; i ++ ;  

  if ( !convert ) return 0 ;

  for ( i = 0 ; i < 12 ; i ++ ) {
    ph[i*pstr] = (1.0 + cos(ph[i*pstr]))/2.0 ;
    th[i*tstr] = (1.0 + th[i*tstr]/M_PI)/2.0 ;
  }
  
  return 0 ;
}
