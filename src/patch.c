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

#include <blaswrap.h>

#include <agg.h>

#include "agg-private.h"

agg_patch_t *agg_patch_new(gint nstmax)

{
  agg_patch_t *P ;

  P = (agg_patch_t *)g_malloc0(sizeof(agg_patch_t)) ;

  P->st = (gdouble *)g_malloc0(AGG_PATCH_POINT_SIZE*nstmax*sizeof(gdouble)) ;
  agg_patch_point_number(P) = 0 ;
  agg_patch_point_number_max(P) = nstmax ;

  gint i ;
  for ( i = 0 ; i  < agg_patch_point_number_max(P) ; i ++ ) {
    agg_patch_point_s(P,i) =
      agg_patch_point_t(P,i) =
      agg_patch_point_len(P,i) = 37.0 ;
  }

  agg_patch_corner_number(P) = 5 ;
  agg_patch_point_number(P) = 4 ;
  /*reasonable defaults for most surfaces*/
  agg_patch_corner_index(P,0) = 0 ;
  agg_patch_corner_index(P,1) = 1 ;
  agg_patch_corner_index(P,2) = 2 ;
  agg_patch_corner_index(P,3) = 3 ;
  agg_patch_corner_index(P,4) = 4 ;
  
  agg_patch_corner_u(P,0) =  0.0 ;
  agg_patch_corner_v(P,0) = -1.0 ;
  agg_patch_corner_u(P,1) =  1.0 ;
  agg_patch_corner_v(P,1) = -1.0 ;
  agg_patch_corner_u(P,2) =  1.0 ;
  agg_patch_corner_v(P,2) =  1.0 ;
  agg_patch_corner_u(P,3) =  0.0 ;
  agg_patch_corner_v(P,3) =  1.0 ;

  agg_patch_corner_s(P,0) = 0.0 ;
  agg_patch_corner_t(P,0) = 0.0 ;
  agg_patch_corner_s(P,1) = 1.0 ;
  agg_patch_corner_t(P,1) = 0.0 ;
  agg_patch_corner_s(P,2) = 1.0 ;
  agg_patch_corner_t(P,2) = 1.0 ;
  agg_patch_corner_s(P,3) = 0.0 ;
  agg_patch_corner_t(P,3) = 1.0 ;
  agg_patch_corner_s(P,4) = 0.0 ;
  agg_patch_corner_t(P,4) = 0.0 ;

  agg_patch_mapping(P) = AGG_PATCH_BILINEAR ;
  agg_patch_wrap_s(P) = FALSE ;
  agg_patch_wrap_t(P) = FALSE ;
  
  return P ;
}

static gint agg_bilinear_map(agg_patch_t *P, gdouble s, gdouble t,
			     gdouble *u, gdouble *v)
{
  gdouble w0, w1, w2, w3 ;

  /*basic bilinear mapping for now
    https://en.wikipedia.org/wiki/Bilinear_interpolation
  */
  w0 = (1.0-s)*(1.0-t) ;
  w1 = (    s)*(1.0-t) ;
  w2 = (    s)*(    t) ;
  w3 = (1.0-s)*(    t) ;
  
  *u =
    agg_patch_corner_u(P,0)*w0 +
    agg_patch_corner_u(P,1)*w1 +
    agg_patch_corner_u(P,2)*w2 +
    agg_patch_corner_u(P,3)*w3 ;
  *v =
    agg_patch_corner_v(P,0)*w0 +
    agg_patch_corner_v(P,1)*w1 +
    agg_patch_corner_v(P,2)*w2 +
    agg_patch_corner_v(P,3)*w3 ;

  return 0 ;
}

static gint agg_spherical_map(agg_patch_t *P, gdouble s, gdouble t,
			      gdouble *u, gdouble *v)

{
  *u = 0.5*(1.0-cos(M_PI*s)) ;

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
  *u = s ;
  /* 0.5*(1.0-cos(M_PI*s)) ; */

  if ( t < 0.5 ) {
    *v = -0.5*(1.0 + cos(2.0*M_PI*t)) ;
  } else {
    *v =  0.5*(1.0 + cos(2.0*M_PI*t)) ;
  }
  
  return 0 ;
}

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

  if ( agg_patch_mapping(P) == AGG_PATCH_TUBULAR ) {
    return agg_tubular_map(P, s, t, u, v) ;
  }
  
  g_error("%s: unrecognized mapping %d", __FUNCTION__, agg_patch_mapping(P)) ;
  
  return 0 ;
}

gint agg_patch_boundary_write(FILE *f, agg_patch_t *P, gdouble ds,
			      agg_surface_workspace_t *w)

{
  gint i ;

  for ( i = 0 ; i < agg_patch_point_number(P) ; i ++ ) {
    fprintf(f, "%lg %lg %lg\n",
	    agg_patch_point_s(P,i), agg_patch_point_t(P,i),
	    agg_patch_point_len(P,i)) ;
  }
  
  return 0 ;
}

static gdouble distance(gdouble *st1, gdouble *st2)

{
  return sqrt((st1[0] - st2[0])*(st1[0] - st2[0]) +
	      (st1[1] - st2[1])*(st1[1] - st2[1])) ;
}

gint agg_patch_edges_parameterize(agg_patch_t *P)

{
  gint i, e ;
  gdouble len ;
  
  agg_patch_point_len(P,0) = 0 ;
  for ( e = 0 ; e < agg_patch_corner_number(P) ; e ++ ) {
    len = 0 ;
    for ( i = agg_patch_corner_index(P,e+0) ;
	  i < agg_patch_corner_index(P,e+1) ; i ++ ) {
      len += distance(&(agg_patch_point_s(P,i+0)),
		      &(agg_patch_point_s(P,i+1))) ;
      agg_patch_point_len(P,i+1) = len ;
    }
    for ( i = agg_patch_corner_index(P,e+0) ;
	  i < agg_patch_corner_index(P,e+1) ; i ++ ) {
      agg_patch_point_len(P,i+1) /= len ;
    }
  }

  /*last edge needs special treatment*/
  len = 0 ;
  for ( i = agg_patch_corner_index(P,agg_patch_corner_number(P)-1) ;
	i < agg_patch_point_number(P) ; i ++ ) {
    len += distance(&(agg_patch_point_s(P,i+0)),
		    &(agg_patch_point_s(P,i+1))) ;
    agg_patch_point_len(P,i+1) = len ;
  }
  for ( i = agg_patch_corner_index(P,agg_patch_corner_number(P)-1) ;
	i < agg_patch_point_number(P) ; i ++ ) {
    agg_patch_point_len(P,i+1) /= len ;
  }
  
  return 0 ;
}  

gint agg_patch_edge_interp(agg_patch_t *P, gint e, gdouble p,
			   gdouble *s, gdouble *t)

{
  gint i ;
  gdouble dp, p0, p1 ;

  if ( p == 1 ) {
    i = agg_patch_corner_index(P,e+1) ;
    *s = agg_patch_point_s(P,i) ;
    *t = agg_patch_point_t(P,i) ;
    return 0 ;
  }
  p0 = 0 ;
  for ( i = agg_patch_corner_index(P,e+0) ;
	i < agg_patch_corner_index(P,e+1) ;
	i ++ ) {
    p1 = agg_patch_point_len(P,i+1) ;
    if (  p0 <= p && p < p1 ) break ;
    p0 = p1 ;
  }

  dp = (p - p0)/(p1 - p0) ;
  *s = dp*agg_patch_point_s(P,i+1) + (1.0-dp)*agg_patch_point_s(P,i) ;
  *t = dp*agg_patch_point_t(P,i+1) + (1.0-dp)*agg_patch_point_t(P,i) ;
  
  return 0 ;
}


gint agg_patch_sample(agg_patch_t *P, gdouble ss, gdouble ts,
		      gdouble *s, gdouble *t)

{
  gint i, ec, edges[] = {0, 1, 2, 3, 0, 1, 2} ;
  gdouble s0, t0, s1, t1 ;
  
  /*find the cut edge*/
  ec = -1 ;
  for ( i = 0 ; i < agg_patch_corner_number(P) ; i ++ ) {
    if ( agg_patch_edge_is_cut(P,i) ) ec = i ;
  }

  if ( ec == -1 ) {
    *s = ss ; *t = ts ;
    return 0 ;
  }
  
  agg_patch_edge_interp(P, ec, ss, &s0, &t0) ;
  agg_patch_edge_interp(P, edges[ec+2], 1.0-ss, &s1, &t1) ;

  *s = ts*s0 + (1.0 - ts)*s1 ; 
  *t = ts*t0 + (1.0 - ts)*t1 ; 
  
  return 0 ;
}
