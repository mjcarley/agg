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

static gdouble linear_to_angle_cosine(gdouble t)

{
  gdouble th ;

  if ( t == 0 ) return M_PI ;
  
  if ( t > 0 ) {
    th = 2.0*M_PI - acos(2.0*t - 1.0) ;
    return th ;
  }

  th = acos(-2.0*t - 1.0) ;

  return th ;
}

static gdouble angle_to_linear_cosine(gdouble th)

{
  gdouble t ;
  
  if ( th > M_PI ) {
    t = 0.5*(1.0 + cos(th)) ;
    return t ;
  }
  
  t = -0.5*(1.0 + cos(th)) ;
  return t ;
}

/** 
 * Evaluate point in discretization of an interval
 *
 * Find the \f$i\f$th point of discretization between \f$t_{\min}\f$
 * and \f$\_{\max}\f$ for a given discretization spacing. 
 * 
 * @param tmin lower limit \f$t_{\min}\f$;
 * @param tmax lower limit \f$t_{\max}\f$;
 * @param nt number of discretization intervals;
 * @param s discretization;
 * @param i point of discretization.
 * 
 * @return \f$t_{i}\f$ of discretization of interval \f$(t_{\min},t_{\max})\f$. 
 */

gdouble agg_spacing_eval(gdouble tmin, gdouble tmax, gint nt,
			 agg_spacing_t s, gint i)

{
  gdouble t, th0, th1, th ;

  if ( s == AGG_SPACING_LINEAR ) {
    return (tmin + (tmax-tmin)*i/(nt-1)) ;
  }

  if ( s == AGG_SPACING_COSINE || AGG_SPACING_CIRCULAR ) {
    th0 = linear_to_angle_cosine(tmin) ;
    th1 = linear_to_angle_cosine(tmax) ;
    th = th0 + (th1-th0)*i/(nt-1) ;
    t = angle_to_linear_cosine(th) ;
    return t ;
  }
  
  t = tmin + (tmax-tmin)*i/(nt-1) ;
  switch ( s ) {
  default: g_error("%s: unknown spacing type %u", __FUNCTION__, s) ;
  case AGG_SPACING_COSINE:  return SIGN(t)*0.5*(1.0-cos(M_PI*t)) ; break ;
  case AGG_SPACING_HALFCOS: return SIGN(t)*(1.0-cos(0.5*M_PI*t)) ; break ;
  case AGG_SPACING_HALFSIN: return sin(0.5*M_PI*t) ; break ;
  }

  return 0.0 ;
}


