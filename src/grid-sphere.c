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

static void sphere_point(gdouble u, gdouble v, gdouble *x)

{
  gdouble sgn ;
  sgn = SIGN(v) ;

  v = fabs(v) ;
  x[0] = 2.0*sqrt(u)*sqrt(1.0 - u)*(v-0.5) ;
  x[1] = sgn*2.0*sqrt(u)*sqrt(1.0-u)*sqrt(v)*sqrt(1.0-v) ;
  x[2] = u ;
  
  return ;
}

static void sphere_uv(gdouble *xi, gdouble *u, gdouble *v)

{
  gdouble r, x[3] ;
  
  /*scale radius to shift back to sphere surface*/
  r = sqrt((xi[0] - 0.0)*(xi[0] - 0.0) +
	   (xi[1] - 0.0)*(xi[1] - 0.0) +
	   (xi[2] - 0.5)*(xi[2] - 0.5)) ;

  x[0] = 0.0 + (xi[0] - 0.0)*0.5/r ;
  x[1] = 0.0 + (xi[1] - 0.0)*0.5/r ;
  x[2] = 0.5 + (xi[2] - 0.5)*0.5/r ;

  /*find (u,v) for x*/
  *u = x[2] ;
  if ( (*u) == 0.0 || (*u) == 1.0 ) {
    *v = 0.0 ; return ;
  }
  *v = x[0]/2.0/sqrt(*u)/sqrt(1.0-*u) + 0.5 ;
  if ( (*v) > 1.0 && (*v) < 1.0 + 1e-10 ) (*v) = 1.0 ;
  
  if ( x[1] < 0 ) *v = -(*v) ;

  g_assert(!isnan(*u)) ;
  g_assert((*u) >=  0.0) ;
  g_assert((*u) <=  1.0) ;
  g_assert(!isnan(*v)) ;
  g_assert((*v) >= -1.0) ;
  g_assert((*v) <=  1.0) ;
  
  return ;
}

static void bisect_sphere_edge(gdouble *uv0, gdouble *uv1,
			       gdouble *u, gdouble *v)

{
  gdouble x0[3], x1[3] ;

  /*find points on sphere and bisect*/
  sphere_point(uv0[0], uv0[1], x0) ;
  sphere_point(uv1[0], uv1[1], x1) ;

  x0[0] = 0.5*(x0[0] + x1[0]) ;
  x0[1] = 0.5*(x0[1] + x1[1]) ;
  x0[2] = 0.5*(x0[2] + x1[2]) ;

  sphere_uv(x0, u, v) ;
    
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
  gint idx ;
  gdouble u, v ;
  
  g_assert(i != j) ;

  if ( i < j ) return ;

  idx = split_index(i, j, np) ;

  split[2*(*ns)+0] = idx ;
  split[2*(*ns)+1] = agg_grid_point_number(g) ;
  (*ns) ++ ;
  
  bisect_sphere_edge(&(agg_grid_point_u(g, i)),
  		     &(agg_grid_point_u(g, j)), &u, &v) ;
  idx = agg_grid_point_number(g) ;
  
  agg_grid_point_u(g, idx) = u ;
  agg_grid_point_v(g, idx) = v ;
  agg_grid_point_number(g) ++ ;
  
  return ;
}

/* static gint compare_int(gconstpointer a, gconstpointer b) */

/* { */
/*   gint i, j ; */

/*   i = *((gint *)a) ; j = *((gint *)b) ; */

/*   if ( i < j ) return -1 ; */
/*   if ( i > j ) return  1 ; */

/*   /\* g_assert_not_reached() ; *\/ */
  
/*   return 0 ; */
/* } */

static gint search_split(gint *split, gint ns, gint idx)

{
  gint i ;

  for ( i = 0 ; i < ns ; i ++ )
    if ( split[2*i+0] == idx ) return split[2*i+1] ;
  
  g_assert_not_reached() ;
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

#if 0
static void check_split(gint *split, gint ns)

{
  gint i, j ;

  for ( i = 0 ; i < ns ; i ++ ) {
    for ( j = i+1 ; j < ns ; j ++ ) {
      g_assert(split[2*i+0] != split[2*j+0]) ;
    }
  }
  
  return ;
}
#endif

gint agg_grid_spherical(agg_grid_t *g, gint refine)

/* 
 * Details of vertices and connectivity in a bunch of places, including here:
 * 
 * http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
 *
 */
  
{
  gint i, j, *tri, nt, split[65536], ns, np, i1, i2 ;
  gint faces[] = {0,  11,  5,
		  0 ,  5,  1,
		  0 ,  1,  7,
		  0 ,  7, 10,
		  0 , 10, 11,
		  1 ,  5,  9,
		  5 , 11,  4,
		  11, 10,  2,
		  10,  7,  6,
		  7 ,  1,  8,
		  3 ,  9,  4,
		  3 ,  4,  2,
		  3 ,  2,  6,
		  3 ,  6,  8,
		  3 ,  8,  9,
		  4 ,  9,  5,
		  2 ,  4, 11,
		  6 ,  2, 10,
		  8 ,  6,  7,
		  9 ,  8,  1} ;
  gdouble uv[] = {5.0000000000000000e-01,  2.3713444394043320e-01,
		  5.0000000000000000e-01,  7.6286555605956674e-01,
		  5.0000000000000000e-01, -2.3713444394043320e-01,
		  5.0000000000000000e-01, -7.6286555605956674e-01,
		  9.2532540417601994e-01, -5.0000000000000000e-01,
		  9.2532540417601994e-01,  5.0000000000000000e-01,
		  7.4674595823980006e-02, -5.0000000000000000e-01,
		  7.4674595823980006e-02,  5.0000000000000000e-01,
		  2.3713444394043320e-01,  1.0000000000000000e+00,
		  7.6286555605956674e-01,  1.0000000000000000e+00,
		  2.3713444394043320e-01,  0.0000000000000000e+00,
		  7.6286555605956674e-01,  0.0000000000000000e+00} ;

  agg_grid_topology(g) = AGG_GRID_SPHERICAL ;
  agg_grid_point_number(g) = 0 ;
  for ( i = 0 ; i < 12 ; i ++ ) {
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
  for ( i = 0 ; i < 20 ; i ++ ) {
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
    /*sort the indices for searching*/
    /* qsort(split, ns, 2*sizeof(gint), compare_int) ; */
    /* check_split(split, ns) ; */
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
  
  g->interp_area = agg_grid_interp_area_spherical ;
  g->interp_line = agg_grid_interp_line_spherical ;
  
  return 0 ;
}

gint agg_grid_interp_area_spherical(agg_grid_t *g, gdouble *uv,
				    gdouble s, gdouble t,
				    gdouble *u, gdouble *v)

{
  /* gint *tri ; */
  gdouble x[9], L[3] ;
  
  /* tri = agg_grid_triangle(g, i) ; */

  /* sphere_point(agg_grid_point_u(g, tri[0]), */
  /* 	       agg_grid_point_v(g, tri[0]), */
  /* 	       &(x[0])) ; */
  /* sphere_point(agg_grid_point_u(g, tri[1]), */
  /* 	       agg_grid_point_v(g, tri[1]), */
  /* 	       &(x[3])) ; */
  /* sphere_point(agg_grid_point_u(g, tri[2]), */
  /* 	       agg_grid_point_v(g, tri[2]), */
  /* 	       &(x[6])) ; */
  sphere_point(uv[0], uv[1], &(x[3*0])) ;
  sphere_point(uv[2], uv[3], &(x[3*1])) ;
  sphere_point(uv[4], uv[5], &(x[3*2])) ;

  L[0] = 1.0 - s - t ; L[1] = s ; L[2] = t ;

  x[0] = L[0]*x[0] + L[1]*x[3] + L[2]*x[6] ;
  x[1] = L[0]*x[1] + L[1]*x[4] + L[2]*x[7] ;
  x[2] = L[0]*x[2] + L[1]*x[5] + L[2]*x[8] ;

  sphere_uv(x, u, v) ;
  
  return 0 ;
}

gint agg_grid_interp_line_spherical(agg_grid_t *g, gdouble *uv,
				    gdouble s,
				    gdouble *u, gdouble *v)

{
  gdouble x[6] ;

  sphere_point(uv[0], uv[1], &(x[3*0])) ;
  sphere_point(uv[2], uv[3], &(x[3*1])) ;

  x[0] = x[0] + s*(x[3] - x[0]) ;
  x[1] = x[1] + s*(x[4] - x[1]) ;
  x[2] = x[2] + s*(x[5] - x[2]) ;

  sphere_uv(x, u, v) ;
  
  return 0 ;
}
