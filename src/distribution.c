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
#include <math.h>
#include <string.h>

#include <glib.h>

#include <blaswrap.h>

#include <agg.h>

#include "agg-private.h"

static const struct {
  gchar *s ;
  gint a[3] ;
} agg_axes_list[] =
  {
   {"xyz",  {1, 2, 3}},
   {"yzx",  {2, 3, 1}},
   {"zxy",  {3, 1, 2}},
   {"xzy",  {1, 3, 2}},
   {"zyx",  {3, 2, 1}},
   {"yxz",  {2, 1, 3}},
   {NULL,   {1, 2, 3}}
  } ;

agg_distribution_t *agg_distribution_alloc(gint nsmax)

{
  agg_distribution_t *d ;
  
  d = (agg_distribution_t *)g_malloc0(sizeof(agg_distribution_t)) ;

  d->nsmax = nsmax ;
  d->ns = 0 ;
  
  d->ts = (gdouble *)g_malloc(2*nsmax*sizeof(gdouble)) ;
  d->w  = &(d->ts[nsmax]) ;
  d->sh = (agg_shape_t **)g_malloc(nsmax*sizeof(agg_shape_t *)) ;
  d->t = (agg_local_transform_t *)g_malloc0(sizeof(agg_local_transform_t)) ;
  
  d->axes[0] = 1 ; d->axes[1] = 2 ; d->axes[2] = 3 ;
  agg_distribution_invert(d) = FALSE ;
  
  return d ;
}

gint agg_distribution_add_shape(agg_distribution_t *d,
				gdouble t, agg_shape_t *s)

{
  if ( d->ns >= d->nsmax )
    g_error("%s: too many shapes (%d) in distribution (nsmax=%d)",
	    __FUNCTION__, d->ns, d->nsmax) ;
  
  d->sh[d->ns] = s ; 
  d->ts[d->ns] = t ; 

  d->ns ++ ;
  
  return 0 ;
}

gint agg_distribution_interpolation_weights(agg_distribution_t *d)

{
  gint j, k ;

  for ( j = 0 ; j < d-> ns ; j ++ ) {
    d->w[j] = 1.0 ;
    for ( k = 0 ; k < j ; k ++ ) {
      d->w[j] *= d->ts[j] - d->ts[k] ;
    }
    for ( k = j+1 ; k < d->ns ; k ++ ) {
      d->w[j] *= d->ts[j] - d->ts[k] ;
    }
    d->w[j] = 1.0/d->w[j] ;
  }
  
  return 0 ;
}

static void interpolate_aerofoil(agg_distribution_t *d, gdouble t,
				 agg_shape_t *s)

{
  gint i, j ;
  gdouble w, wtot ;
  
  g_assert(s->nsmax >= d->sh[0]->nsmax) ;
  /* g_assert(s->nb == 1) ; */
  
  s->closed = TRUE ;
  s->type = AGG_SHAPE_AEROFOIL ;

  s->i[0] = d->sh[0]->i[0] ; 
  s->i[1] = d->sh[0]->i[1] ; 
  s->i[2] = d->sh[0]->i[2] ; 
  s->i[3] = d->sh[0]->i[3] ;
  s->nb = d->sh[0]->nb ;
  
  memset(s->n1, 0, 8*sizeof(gdouble)) ;
  memset(s->n2, 0, 8*sizeof(gdouble)) ;

  memset(s->s, 0, (s->nsmax)*sizeof(gdouble)) ;

  wtot = 0.0 ;
  for ( i = 0 ; i < agg_distribution_station_number(d) ; i ++ ) {
    if ( fabs(t - d->ts[i]) < 1e-12 ) {
      s->n1[0] = d->sh[i]->n1[0] ;
      s->n2[0] = d->sh[i]->n2[0] ;
      memcpy(s->s, d->sh[i]->s, (d->sh[i]->i[3])*sizeof(gdouble)) ;
      return ;
    }
    g_assert(d->sh[i]->i[3] == d->sh[0]->i[3]) ;
    w = d->w[i]/(t - d->ts[i]) ;
    wtot += w ;
    s->n1[0] += (d->sh[i]->n1[0])*w ;
    s->n2[0] += (d->sh[i]->n2[0])*w ;
    for ( j = 0 ; j < d->sh[i]->i[3] ; j ++ ) {
      s->s[j] += w*(d->sh[i]->s[j]) ;
    }
  }
  
  s->n1[0] /= wtot ; s->n2[0] /= wtot ;
  for ( j = 0 ; j < s->i[3] ; j ++ ) s->s[j] /= wtot ;

  return ;
}

static void interpolate_ellipse(agg_distribution_t *d, gdouble t,
				agg_shape_t *s)

{
  gint i, j, nb ;
  gdouble w, wtot ;
  
  g_assert(s->nsmax >= d->sh[0]->nsmax) ;
  
  s->closed = TRUE ;
  s->type = AGG_SHAPE_ELLIPSE ;

  s->i[0] = d->sh[0]->i[0] ; 
  s->i[1] = d->sh[0]->i[1] ; 
  s->i[2] = d->sh[0]->i[2] ; 

  memset(s->n1, 0, 8*sizeof(gdouble)) ;
  memset(s->n2, 0, 8*sizeof(gdouble)) ;

  memset(s->s, 0, (s->nsmax)*sizeof(gdouble)) ;
  nb = d->sh[0]->nb ;
  memset(s->b, 0, (nb+2)*sizeof(gdouble)) ;
  s->nb = nb ;
  
  wtot = 0.0 ;
  for ( i = 0 ; i < agg_distribution_station_number(d) ; i ++ ) {
    if ( fabs(t - d->ts[i]) < 1e-12 ) {
      agg_shape_copy(s, d->sh[i]) ;
      return ;
    }
    g_assert(d->sh[i]->i[2] == d->sh[0]->i[2]) ;
    g_assert(d->sh[i]->nb == nb) ;
    w = d->w[i]/(t - d->ts[i]) ;
    wtot += w ;
    for ( j = 0 ; j < d->sh[i]->i[2] ; j ++ ) {
      s->s[j] += w*(d->sh[i]->s[j]) ;
    }
    for ( j = 0 ; j < nb+2 ; j ++ ) {
      s->b[j] += w*(d->sh[i]->b[j]) ;
    }
    for ( j = 0 ; j < nb+1 ; j ++ ) {
      s->n1[j] += (d->sh[i]->n1[j])*w ;
      s->n2[j] += (d->sh[i]->n2[j])*w ;
    }
  }
  
  for ( j = 0 ; j < nb+1 ; j ++ ) {
    s->n1[j] /= wtot ;
    s->n2[j] /= wtot ;
  }
  for ( j = 0 ; j < s->i[2] ; j ++ ) s->s[j] /= wtot ;
  for ( j = 0 ; j < nb+2 ; j ++ )    s->b[j] /= wtot ;

  return ;
}


gint agg_distribution_interpolate_shape(agg_distribution_t *d, gdouble s,
					agg_shape_t *sh)

{
  /*might want to check limits on t (but allow some extrapolation?)*/

  if ( sh != NULL ) {
    switch ( d->sh[0]->type ) {
    default: g_assert_not_reached() ; break ;
    case AGG_SHAPE_AEROFOIL: interpolate_aerofoil(d, s, sh) ; break ;
    case AGG_SHAPE_ELLIPSE : interpolate_ellipse(d, s, sh) ; break ;
    }
  }
  
  return 0 ;
}

static void split_string_sign(gchar *str, gchar s[], gint sgn[])

{
  gint i, n ;

  n = 0 ; i = 0 ;
  s[3] = '\0' ;
  while ( n < 3 && i < strlen(str) ) {
    sgn[n] = 1 ;
    if ( str[i] == '-' ) { sgn[n] = -1 ; i ++ ; } 
    if ( str[i] == '+' ) { sgn[n] = +1 ; i ++ ; } 
    
    s[n] = str[i] ; n ++ ; i ++ ;
  }
  
  return ;
}

gint agg_distribution_axes_parse(gchar *str, gint *a)

{
  gint i, sgn[3] ;
  gchar s[4] ;

  split_string_sign(str, s, sgn) ;
  
  for ( i = 0 ; agg_axes_list[i].s != NULL ; i ++ ) {
    if ( !strcmp(agg_axes_list[i].s, s) ) {
      a[0] = agg_axes_list[i].a[0]*sgn[0] ;
      a[1] = agg_axes_list[i].a[1]*sgn[1] ;
      a[2] = agg_axes_list[i].a[2]*sgn[2] ;

      return 0 ;
    }
  }
  
  return -1 ;
}

gint agg_distribution_point_eval(agg_distribution_t *d,
				 gdouble u, gdouble v,
				 agg_parser_t *p,
				 agg_shape_t *sh,
				 gdouble *x)

{
  gdouble y[3] ;
  gint *axes ;
  
  agg_distribution_interpolate_shape(d, u, sh) ;

  p->values[AGG_PARSER_PARAMETER_RESERVED_S] = u ;
  agg_parser_expressions_evaluate(p) ;
  agg_local_transform_eval_parameters(d->t) ;

  y[0] = fabs(v) ;
  y[1] = agg_shape_eval(sh, v, -1) ;
  y[2] = 0.0 ;
  agg_local_transform_apply(d->t, y) ;
  axes = d->axes ;

  x[0] = y[SIGN(axes[0])*axes[0]-1]*SIGN(axes[0]) ;
  x[1] = y[SIGN(axes[1])*axes[1]-1]*SIGN(axes[1]) ;
  x[2] = y[SIGN(axes[2])*axes[2]-1]*SIGN(axes[2]) ;  
  
  return 0 ;
}

gint agg_distribution_point_normal_eval(agg_distribution_t *d,
					gdouble u, gdouble v,
					agg_parser_t *p,
					agg_shape_t *sh,
					gdouble *x, gdouble *n, gdouble *J)

{
  gdouble y[3] ;
  gint *axes ;
  
  agg_distribution_interpolate_shape(d, u, sh) ;

  p->values[AGG_PARSER_PARAMETER_RESERVED_S] = u ;
  agg_parser_expressions_evaluate(p) ;
  agg_local_transform_eval_parameters(d->t) ;

  y[0] = fabs(v) ;
  y[1] = agg_shape_eval(sh, v, -1) ;
  y[2] = 0.0 ;
  agg_local_transform_apply(d->t, y) ;
  axes = d->axes ;

  x[0] = y[SIGN(axes[0])*axes[0]-1]*SIGN(axes[0]) ;
  x[1] = y[SIGN(axes[1])*axes[1]-1]*SIGN(axes[1]) ;
  x[2] = y[SIGN(axes[2])*axes[2]-1]*SIGN(axes[2]) ;  
  
  return 0 ;
}
