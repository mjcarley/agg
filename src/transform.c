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
  gchar *transform ;
  agg_local_transform_func_t func ;
  gint np ;
} local_transform_parse_list[] =
  {
   {"rotate", agg_local_transform_rotate, 3},
   {"shift",  agg_local_transform_shift,  3},
   {"shrink", agg_local_transform_shrink, 3},
   {"scale",  agg_local_transform_scale, 1},
   {NULL, NULL, 0}
  } ;

gint agg_local_transform_init(agg_local_transform_t *t)

{
  gint i ;
  t->nt = 0 ;
  t->p1[0] = 0 ;
  for ( i = 0 ; i < AGG_TRANSFORM_PARAMETER_NUMBER ; i ++ ) {
    t->isexpr[i] = FALSE ;
  }

  return 0 ;
}

agg_local_transform_t *agg_local_transform_alloc(void)

{
  agg_local_transform_t *t ;
  
  t = (agg_local_transform_t *)g_malloc0(sizeof(agg_local_transform_t)) ;

  agg_local_transform_init(t) ;
  
  return t ;
}

gint agg_local_transform_rotate(gdouble *x, gdouble *p, gdouble *y)

{
  gdouble xt, yt ;
  
  xt = cos(p[2])*(x[0] - p[0]) - sin(p[2])*(x[1] - p[1]) + p[0] ;
  yt = sin(p[2])*(x[0] - p[0]) + cos(p[2])*(x[1] - p[1]) + p[1] ;
  
  y[0] = xt ; y[1] = yt ;

  return 0 ;
}

gint agg_local_transform_shift(gdouble *x, gdouble *p, gdouble *y)

{
  y[0] = x[0] + p[0] ; y[1] = x[1] + p[1] ; y[2] = x[2] + p[2] ;

  return 0 ;
}

gint agg_local_transform_scale(gdouble *x, gdouble *p, gdouble *y)

{
  y[0] = x[0]*p[0] ; y[1] = x[1]*p[0] ; y[2] = x[2]*p[0] ;

  return 0 ;
}

gint agg_local_transform_shrink(gdouble *x, gdouble *p, gdouble *y)

{
  y[0] = p[0] + p[2]*(x[0] - p[0]) ;
  y[1] = p[1] + p[2]*(x[1] - p[1]) ;

  return 0 ;
}

gint agg_local_transform_parse(agg_local_transform_t *T, gchar *type,
			       gdouble *p, gint np)

/*
 * type: string identifier
 * p:    parameters (all doubles)
 * np:   number of parameters (checked for transform)
 */
  
{
  gint i, j ;

  for ( i = 0 ; local_transform_parse_list[i].transform != NULL ; i ++ ) {
    if ( !strcmp(local_transform_parse_list[i].transform, type) ) break ;
  }
       
  if ( local_transform_parse_list[i].transform == NULL ) return -1 ;

  /*check number of parameters*/
  if ( np != local_transform_parse_list[i].np ) {
    g_error("%s: transform \"%s\" should have %d parameter%s, not %d",
  	    __FUNCTION__, type,
	    local_transform_parse_list[i].np,
	    plural_character(local_transform_parse_list[i].np),
	    np) ;
  }

  T->func[T->nt] = local_transform_parse_list[i].func ;
  T->p1[T->nt+1] = T->p1[T->nt] + np ;
  for ( j = 0 ; j < np ; j ++ ) T->p[T->p1[T->nt]+j] = p[j] ;

  T->nt ++ ;
  
  return 0 ;
}

gint agg_local_transform_apply(agg_local_transform_t *T, gdouble *y)

{
  gint i, s ;

  for ( i = 0 ; i < T->nt ; i ++ ) {
    s = (T->func[i])(y, &(T->p[T->p1[i]]), y) ;
    if ( s != 0 ) {
      g_error("%s: return error (%d) in local transformation %d",
	      __FUNCTION__, s, i) ;
    }
  }
  
  return 0 ;
}

