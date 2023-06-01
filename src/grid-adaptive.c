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

#define REAL double
#define VOID int

#include "triangle.h"

static gdouble tri_area(gdouble *x1, gdouble *x2, gdouble *x3)

{
  gdouble A, r12[3], r23[3], r[3] ;

  agg_vector_diff(r12, x2, x1) ;
  agg_vector_diff(r23, x3, x2) ;
  agg_vector_cross(r, r12, r23) ;

  A = 0.5*agg_vector_length(r) ;
  
  return A ;
}

gint triunsuitable(gdouble *triorg, gdouble *tridest, gdouble *triapex,
		   gdouble area, gpointer user_data)

{
  gpointer *data = user_data ;
  agg_body_t *b = data[0] ;
  agg_adaptive_grid_data_t *quality = data[1] ;
  agg_workspace_t *w = data[2] ;
  gdouble A, len, x1[3], x2[3], x3[3], n1[3], n2[3], n3[3] ;

  agg_body_point_normal_eval(b, triorg [0], triorg [1], x1, n1, w) ;
  agg_body_point_normal_eval(b, tridest[0], tridest[1], x2, n2, w) ;
  agg_body_point_normal_eval(b, triapex[0], triapex[1], x3, n3, w) ;

  /*check element curvature*/
  if ( agg_vector_scalar(n1,n2) < agg_adaptive_grid_ncos(quality) ) return 1 ;
  if ( agg_vector_scalar(n2,n3) < agg_adaptive_grid_ncos(quality) ) return 1 ;
  if ( agg_vector_scalar(n3,n1) < agg_adaptive_grid_ncos(quality) ) return 1 ;
  
  len = MAX(agg_vector_distance2(x1,x2),
	    MAX(agg_vector_distance2(x2,x3),
		agg_vector_distance2(x3,x1))) ;
  len = sqrt(len) ;

  if ( len > quality->lmax ) return 1 ;
  if ( len < quality->lmin ) return 0 ;

  A = tri_area(x1, x2, x3) ;

  if ( A > quality->amax ) return 1 ;
  if ( A < quality->amin ) return 0 ;
  
  return 0 ;
}

/*
 * Triangle options:
 * 
 * p to triangulate inside a defined boundary
 * q for quality mesh
 * u for user-defined triangle constraint
 * 
 * triangulateio inputs:
 * 
 * pointlist, numberofpoints, numberofpointattributes (boundary of
 * triangulation)
 * 
 * pointmarkerlist (to label intersection curves)
 * 
 * segmentlist, numberofsegments, segmentmarkerlist (to mark boundary)
 *
 * numberofholes, numberofregions (probably both zero, at least to
 * start with)
 *
 * triangulateio outputs:
 *
 * pointlist = NULL
 */

gint agg_grid_adaptive(agg_grid_t *g, agg_adaptive_grid_data_t *d)
  
{
  g->data = g_malloc0(sizeof(agg_adaptive_grid_data_t)) ;

  memcpy(g->data, d, sizeof(agg_adaptive_grid_data_t)) ;
  
  agg_grid_topology(g) = AGG_GRID_ADAPTIVE ;
  agg_grid_point_number(g) = 0 ;
  agg_grid_triangle_number(g) = 0 ;

  g->interp_area = agg_grid_interp_area_adaptive ;
  g->interp_line = agg_grid_interp_line_adaptive ;
  
  return 0 ;
}

gint agg_grid_interp_area_adaptive(agg_grid_t *g, gdouble *uv,
				   gdouble s, gdouble t,
				   gdouble *u, gdouble *v)

{
  
  return 0 ;
}

gint agg_grid_interp_line_adaptive(agg_grid_t *g, gdouble *uv,
				   gdouble s,
				   gdouble *u, gdouble *v)

{
  
  return 0 ;
}

agg_adaptive_grid_workspace_t *agg_adaptive_grid_workspace_alloc(gint np,
								 gint ns,
								 gint nt)

{
  agg_adaptive_grid_workspace_t *w ;  
  struct triangulateio *in, *out ;
  
  w = (agg_adaptive_grid_workspace_t *)
    g_malloc0(sizeof(agg_adaptive_grid_workspace_t)) ;

  in = w->in = g_malloc0(sizeof(struct triangulateio)) ;
  out = w->out = g_malloc0(sizeof(struct triangulateio)) ;

  w->npmax = np ; w->nsmax = ns ; w->ntmax = nt ;
  w->nt = 0 ; w->ns = 0 ; w->np = 0 ;

  in->pointlist     = (gdouble *)g_malloc0(2*np*sizeof(gdouble)) ;
  in->segmentlist   = (gint *)g_malloc0(2*ns*sizeof(gint)) ;
  out->pointlist    = (gdouble *)g_malloc0(2*np*sizeof(gdouble)) ;
  out->trianglelist = (gint *)g_malloc0(3*nt*sizeof(gint)) ;
  out->segmentlist   = (gint *)g_malloc0(2*ns*sizeof(gint)) ;
  
  return w ;
}

gint agg_adaptive_grid_make(agg_grid_t *g, agg_body_t *b,
			    agg_workspace_t *w,
			    agg_adaptive_grid_workspace_t *wg)

{
  struct triangulateio *in, *out ;
  gpointer data[4] ;
  gdouble *x ;
  gint i, *tri, *t ;
  
  agg_grid_point_number(g) = agg_grid_triangle_number(g) = 0 ;

  /*really basic grid to start with*/
  in = wg->in ; out = wg->out ;
  x = in->pointlist ;

  x[2*0+0] = agg_body_smin(b) ; x[2*0+1] = agg_body_tmin(b) ; 
  x[2*1+0] = agg_body_smax(b) ; x[2*1+1] = agg_body_tmin(b) ; 
  x[2*2+0] = agg_body_smax(b) ; x[2*2+1] = agg_body_tmax(b) ; 
  x[2*3+0] = agg_body_smin(b) ; x[2*3+1] = agg_body_tmax(b) ; 
  in->numberofpoints = 4 ;
  in->numberofpointattributes = 0 ;

  t = in->segmentlist ;
  in->numberofsegments = 4 ;
  in->segmentlist[2*0+0] = 0 ; in->segmentlist[2*0+1] = 1 ; 
  in->segmentlist[2*1+0] = 1 ; in->segmentlist[2*1+1] = 2 ; 
  in->segmentlist[2*2+0] = 2 ; in->segmentlist[2*2+1] = 3 ; 
  in->segmentlist[2*3+0] = 3 ; in->segmentlist[2*3+1] = 0 ; 
  
  /*Quiet, polygon boundary, quality mesh, user-defined splitting,
    zero-based indexing*/
  data[0] = b ; data[1] = g->data ; data[2] = w ;
  triangulate("Qpquz", in, out, NULL, data) ;

  /*transfer the workspace data to the agg_grid*/
  x = out->pointlist ;
  agg_grid_point_number(g) = out->numberofpoints ;
  if ( agg_grid_point_number(g) > agg_grid_point_number_max(g) ) {
    g_error("%s: too many points (%d) for maximum number allowed (%d)",
	    __FUNCTION__, agg_grid_point_number(g),
	    agg_grid_point_number_max(g)) ;
  }    
  for ( i = 0 ; i < out->numberofpoints ; i ++ ) {
    agg_grid_point_u(g,i) = x[2*i+0] ;
    agg_grid_point_v(g,i) = x[2*i+1] ;
  }

  t = out->trianglelist ;
  agg_grid_triangle_number(g) = out->numberoftriangles ;
  if ( agg_grid_triangle_number(g) > agg_grid_triangle_number_max(g) ) {
    g_error("%s: too many triangles (%d) for maximum number allowed (%d)",
	    __FUNCTION__, agg_grid_triangle_number(g),
	    agg_grid_triangle_number_max(g)) ;
  }    
  for ( i = 0 ; i < out->numberoftriangles ; i ++ ) {
    tri = agg_grid_triangle(g, i) ;
    tri[0] = t[3*i+0] ; 
    tri[1] = t[3*i+1] ; 
    tri[2] = t[3*i+2] ; 
  }
  
  return 0 ;
}
