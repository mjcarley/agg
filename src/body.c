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

agg_body_t *agg_body_alloc(void)

{
  agg_body_t *b ;

  b = (agg_body_t *)g_malloc0(sizeof(agg_body_t)) ;

  agg_body_distribution_number(b) = 0 ;

  agg_body_parser(b) = NULL ;
  agg_body_grid(b) = NULL ;
  
  return b ;
}

gint agg_body_distribution_add(agg_body_t *b, agg_distribution_t *d,
			       gchar *name)

{
  g_assert(agg_body_distribution_number(b) <
	   AGG_BODY_DISTRIBUTION_NUMBER_MAX) ;

  b->names[agg_body_distribution_number(b)] = g_strdup(name) ;
  b->d[agg_body_distribution_number(b)] = d ;
  
  agg_body_distribution_number(b) ++ ;

  return 0 ;
}

gint agg_body_distributions_list(FILE *f, agg_body_t *b)

{
  gint i ;

  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) {
    fprintf(f, "%s\n", b->names[i]) ;
  }
  
  return 0 ;
}

gint agg_body_distribution_locate_u(agg_body_t *b, gdouble u)

{
  gint i ;
  agg_distribution_t *d ;
  
  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) {
    d = agg_body_distribution(b, i) ;
    if ( (agg_distribution_parameter_min(d) <= u &&
	  u <= agg_distribution_parameter_max(d)) ||
	 (agg_distribution_parameter_min(d) >= u &&
	  u >= agg_distribution_parameter_max(d))
	 ) return i ;
  }    
  
  return -1 ;
}

gint agg_body_point_eval(agg_body_t *b, gdouble u, gdouble v, gdouble *x,
			 agg_workspace_t *w)

{
  gint i, *axes ;
  agg_distribution_t *d ;
  gdouble y[3] ;
  agg_parser_t *p ;
  agg_shape_t *sh = agg_workspace_shape(w) ;
  agg_local_transform_t *T = agg_workspace_local_transform(w) ;
  
  g_assert( (p = agg_body_parser(b)) != NULL) ;
  i = agg_body_distribution_locate_u(b, u) ;
  d = agg_body_distribution(b, i) ;
  axes = d->axes ;
  agg_distribution_interpolate_shape(d, u, sh) ;

  p->values[AGG_PARSER_PARAMETER_RESERVED_S] = u ;
  agg_parser_expressions_evaluate(p) ;
  agg_local_transform_set_parameters(T, d->t) ;
  
  y[0] = fabs(v) ;
  y[1] = agg_shape_eval(sh, v, -1) ;
  y[2] = 0.0 ;
  agg_local_transform_apply(T, y) ;
  x[0] = y[SIGN(axes[0])*axes[0]-1]*SIGN(axes[0]) ;
  x[1] = y[SIGN(axes[1])*axes[1]-1]*SIGN(axes[1]) ;
  x[2] = y[SIGN(axes[2])*axes[2]-1]*SIGN(axes[2]) ;

  return 0 ;
}

gint agg_body_point_normal_eval(agg_body_t *b, gdouble u, gdouble v,
				gdouble *x, gdouble *n, agg_workspace_t *w)

{
  gdouble ee, xu[3], xv[3], t ;

  ee = 1e-6 ;

  agg_body_point_eval(b, u, v, x, w) ;

  t = (u > agg_body_smin(b) + ee ? u-ee : u+ee) ; 
  agg_body_point_eval(b, t, v, xu, w) ;
  xu[0] = (xu[0] - x[0])/(t - u) ;
  xu[1] = (xu[1] - x[1])/(t - u) ;
  xu[2] = (xu[2] - x[2])/(t - u) ;

  t = (v > agg_body_tmin(b) + ee ? v-ee : v+ee) ; 
  agg_body_point_eval(b, u, t, xv, w) ;
  xv[0] = (xv[0] - x[0])/(t - v) ;
  xv[1] = (xv[1] - x[1])/(t - v) ;
  xv[2] = (xv[2] - x[2])/(t - v) ;

  agg_vector_cross(n, xu, xv) ;
  t = agg_vector_length(n) ;
  n[0] /= t ; n[1] /= t ; n[2] /= t ;
  
  return 0 ;
}

gint agg_body_parameter_limits(agg_body_t *b, gdouble *umin, gdouble *umax)

{
  gint i ;
  agg_distribution_t *d ;
  
  *umin = G_MAXDOUBLE ; *umax = -G_MAXDOUBLE ;

  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) {
    d = agg_body_distribution(b,i) ;
    *umin = MIN(*umin, agg_distribution_parameter_min(d)) ;
    *umax = MAX(*umax, agg_distribution_parameter_max(d)) ;
  }
  
  return 0 ;
}
