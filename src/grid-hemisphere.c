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

static gint search_split(gint *split, gint ns, gint idx)

{
  gint i ;

  for ( i = 0 ; i < ns ; i ++ )
    if ( split[2*i+0] == idx ) return split[2*i+1] ;

  return -1 ;
  /* g_assert_not_reached() ; */
}

static void hemisphere_point(gdouble u, gdouble v, gdouble *x)

{
  gdouble sgn, r ;
  sgn = SIGN(v) ;

  v = fabs(v) ;

  r = sqrt(1.0 - u*u) ;
  x[0] = (v-0.5)*r ;
  x[1] = sgn*sqrt(v)*sqrt(1.0-v)*r ;
  x[2] = 0.5*u ;
  
  return ;
}

static void hemisphere_uv(gdouble *xi, gdouble *u, gdouble *v)

{
  gdouble r, x[3] ;
  
  /*scale radius to shift back to sphere surface*/
  r = sqrt((xi[0] - 0.0)*(xi[0] - 0.0) +
	   (xi[1] - 0.0)*(xi[1] - 0.0) +
	   (xi[2] - 0.0)*(xi[2] - 0.0)) ;

  x[0] = 0.0 + (xi[0] - 0.0)*0.5/r ;
  x[1] = 0.0 + (xi[1] - 0.0)*0.5/r ;
  x[2] = 0.0 + (xi[2] - 0.0)*0.5/r ;

  /*find (u,v) for x*/
  *u = 2.0*x[2] ;
  r = sqrt(1.0 - (*u)*(*u)) ;

  if ( r == 0 ) {
    *v = 0 ; return ;
  }
  
  *v = x[0]/r + 0.5 ;
  if ( *v > 1 ) *v = 1 ;
  if ( x[1] < 0.0 ) *v = -(*v) ;

  g_assert(!isnan(*u)) ;
  g_assert((*u) >=  0.0) ;
  g_assert((*u) <=  1.0) ;

  if ( isnan(*v) ) {
    g_error("%s: NaN, v=%lg; xi = (%lg, %lg, %lg)",
	    __FUNCTION__, *v, xi[0], xi[1], xi[2]) ;
  }

  if ( (*v) < -1.0 || (*v) > 1.0 ) {
    g_error("%s: v out of range (1+%lg);"
	    "xi = (%lg, %lg, %lg); "
	    "x = (%lg, %lg, %lg)",
	    __FUNCTION__, fabs(*v) - 1.0,
	    xi[0], xi[1], xi[2],
	    x[0], x[1], x[2]) ;
  }
  
  return ;
}

static void bisect_hemisphere_edge(gdouble *uv0, gdouble *uv1,
				   gdouble *u, gdouble *v)

{
  gdouble x0[3], x1[3] ;

  /*find points on sphere and bisect*/
  hemisphere_point(uv0[0], uv0[1], x0) ;
  hemisphere_point(uv1[0], uv1[1], x1) ;

  x0[0] = 0.5*(x0[0] + x1[0]) ;
  x0[1] = 0.5*(x0[1] + x1[1]) ;
  x0[2] = 0.5*(x0[2] + x1[2]) ;

  hemisphere_uv(x0, u, v) ;
    
  return ;
}

static gint split_index(gint i, gint j, gint n)

{
  if ( i > j ) return i*n + j ;

  return j*n + i ;
}

static void split_edge(agg_grid_t *g, gint np, gint i, gint j,
		       gint *split, gint *ns)

{
  gint idx, tmp ;
  gdouble u, v ;
  
  g_assert(i != j) ;

  if ( i < j ) {
    tmp = i ; i = j ; j = tmp ;
  }

  idx = split_index(i, j, np) ;

  if ( search_split(split, *ns, idx) != -1 ) return ;
  
  split[2*(*ns)+0] = idx ;
  split[2*(*ns)+1] = agg_grid_point_number(g) ;
  (*ns) ++ ;
  
  bisect_hemisphere_edge(&(agg_grid_point_u(g, i)),
			 &(agg_grid_point_u(g, j)), &u, &v) ;
  idx = agg_grid_point_number(g) ;
  
  agg_grid_point_u(g, idx) = u ;
  agg_grid_point_v(g, idx) = v ;
  agg_grid_point_number(g) ++ ;
  
  return ;
}

static void split_triangle(agg_grid_t *g, gint i, gint np, gint *split, gint ns)

{
  gint *tri, idx, new[3], *tnew ;

  tri = agg_grid_triangle(g, i) ;
  idx = split_index(tri[0], tri[1], np) ;
  new[0] = search_split(split, ns, idx) ;
  idx = split_index(tri[1], tri[2], np) ;
  new[1] = search_split(split, ns, idx) ;
  idx = split_index(tri[2], tri[0], np) ;
  new[2] = search_split(split, ns, idx) ;

  tnew = agg_grid_triangle(g, agg_grid_triangle_number(g)) ;
  tnew[0] = new[0] ; tnew[1] = tri[1] ; tnew[2] = new[1] ;
  agg_grid_triangle_number(g) ++ ;

  tnew = agg_grid_triangle(g, agg_grid_triangle_number(g)) ;
  tnew[0] = new[1] ; tnew[1] = tri[2] ; tnew[2] = new[2] ;
  agg_grid_triangle_number(g) ++ ;

  tnew = agg_grid_triangle(g, agg_grid_triangle_number(g)) ;
  tnew[0] = new[2] ; tnew[1] = tri[0] ; tnew[2] = new[0] ;
  agg_grid_triangle_number(g) ++ ;

  tri[0] = new[0] ; tri[1] = new[1] ; tri[2] = new[2] ;
  
  return ;
}

gint agg_grid_hemispherical(agg_grid_t *g, gint refine)

{
  gint i, j, *tri, nt, split[65536], ns, np, i1, i2 ;
  gint faces[] =
    {0, 2, 5,
     0, 5, 1,
     0, 1, 7,
     0, 7, 3,
     0, 3, 2,
     1, 5, 9,
     5, 2, 4,
     3, 7, 6,
     7, 1, 8,
     4, 9, 5,
     8, 6, 7,
     9, 8, 1} ;
  gdouble uv[] =
    {0.850650808352040, -0.500000000000000,
     0.850650808352040,  0.500000000000000,
     0.000000000000000, -0.762865556059567,
     0.000000000000000, -0.237134443940433,
     0.000000000000000,  1.000000000000000,
     0.525731112119134,  1.000000000000000,
     0.000000000000000,  0.000000000000000,
     0.525731112119134,  0.000000000000000,
     0.000000000000000,  0.237134443940433,
     0.000000000000000,  0.762865556059567} ;

  agg_grid_topology(g) = AGG_GRID_HEMISPHERICAL ;
  agg_grid_point_number(g) = 0 ;
  for ( i = 0 ; i < 10 ; i ++ ) {
    agg_grid_point_u(g, i) = uv[2*i+0] ;
    agg_grid_point_v(g, i) = uv[2*i+1] ;
    if ( agg_grid_point_number(g) >= agg_grid_point_number_max(g) ) {
      fprintf(stderr,
	      "%s: not enough points allocated in grid (%d)\n",
	      __FUNCTION__,
	      agg_grid_point_number_max(g)) ;
      return 1 ;
    }
    agg_grid_point_number(g) ++ ;
  }

  i1 = 1 ; i2 = 2 ;
  if ( agg_grid_invert(g) ) { i2 = 1 ; i1 = 2 ; }
  agg_grid_triangle_number(g) = 0 ;
  for ( i = 0 ; i < 12 ; i ++ ) {
    tri = agg_grid_triangle(g, agg_grid_triangle_number(g)) ;
    tri[ 0] = faces[3*i+0] ;
    tri[i1] = faces[3*i+1] ;
    tri[i2] = faces[3*i+2] ;
    if ( agg_grid_triangle_number(g) >= agg_grid_triangle_number_max(g) ) {
      fprintf(stderr,
	      "%s: not enough triangles allocated in grid (%d)\n",
	      __FUNCTION__,
	      agg_grid_triangle_number_max(g)) ;
      return 1 ;
    }
    agg_grid_triangle_number(g) ++ ;
  }

  for ( i = 0 ; i < refine ; i ++ ) {
    nt = agg_grid_triangle_number(g) ;
    np = agg_grid_point_number(g) ;
    ns = 0 ;
    for ( j = 0 ; j < nt ; j ++ ) {
      /*split edges and index new points*/
      if ( agg_grid_point_number(g)+3 >= agg_grid_point_number_max(g) ) {
	fprintf(stderr,
		"%s: not enough points allocated in grid (%d)\n",
		__FUNCTION__,
		agg_grid_point_number_max(g)) ;
	return 1 ;
      }
      tri = agg_grid_triangle(g, j) ;
      split_edge(g, np, tri[0], tri[1], split, &ns) ;
      split_edge(g, np, tri[1], tri[2], split, &ns) ;
      split_edge(g, np, tri[2], tri[0], split, &ns) ;
    }
    /*split the triangles and insert into g*/
    for ( j = 0 ; j < nt ; j ++ ) {
      if ( agg_grid_triangle_number(g)+3 >=
	   agg_grid_triangle_number_max(g) ) {
	fprintf(stderr,
		"%s: not enough triangles allocated in grid (%d)\n",
		__FUNCTION__,
		agg_grid_triangle_number_max(g)) ;
	return 1 ;
      }
      split_triangle(g, j, np, split, ns) ;
    }
  }
  
  g->interp_area = agg_grid_interp_area_hemispherical ;

  return 0 ;
}  

gint agg_grid_interp_area_hemispherical(agg_grid_t *g, gdouble *uv,
					gdouble s, gdouble t,
					gdouble *u, gdouble *v)

{
  gint *tri ;
  gdouble x[9], L[3] ;
  
  hemisphere_point(uv[0], uv[1], &(x[3*0])) ;
  hemisphere_point(uv[2], uv[3], &(x[3*1])) ;
  hemisphere_point(uv[4], uv[5], &(x[3*2])) ;

  L[0] = 1.0 - s - t ; L[1] = s ; L[2] = t ;

  x[0] = L[0]*x[0] + L[1]*x[3] + L[2]*x[6] ;
  x[1] = L[0]*x[1] + L[1]*x[4] + L[2]*x[7] ;
  x[2] = L[0]*x[2] + L[1]*x[5] + L[2]*x[8] ;

  hemisphere_uv(x, u, v) ;
  
  return 0 ;
}
