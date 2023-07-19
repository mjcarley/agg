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

#ifdef HAVE_TRIANGLE_API_H
#include <triangle.h>
#include <triangle_api.h>
void traversalinit(struct memorypool *pool);
void *traverse(struct memorypool *pool);
vertex vertextraverse(mesh *m);
triangle *triangletraverse(mesh *m);
int plus1mod3[3] = {1, 2, 0};
int minus1mod3[3] = {2, 0, 1};
#endif /*HAVE_TRIANGLE_API_H*/

/** 
 * @{ 
 *
 * @ingroup mesh
 */

/** 
 * Allocate a new ::agg_mesh_t
 * 
 * @param npmax maximum number of points in mesh;
 * @param nspmax maximum number of splines in mesh;
 * @param nemax maximum number of elements in mesh.
 * 
 * @return pointer to newly allocated ::agg_mesh_t.
 */

agg_mesh_t *agg_mesh_new(gint npmax, gint nspmax, gint nemax)

{
  agg_mesh_t *m ;

  m = (agg_mesh_t *)g_malloc0(sizeof(agg_mesh_t)) ;

  m->p   = (gdouble *)g_malloc0(AGG_MESH_POINT_SIZE*npmax*
				sizeof(gdouble)) ;
  m->sp  = (gint *)g_malloc0(nspmax*8*sizeof(gint)) ;
  m->isp = (gint *)g_malloc0(nspmax*sizeof(gint)) ;
  m->e   = (gint *)g_malloc0(AGG_MESH_ELEMENT_SIZE*nemax*sizeof(gint)) ;
  m->ptags = (gint *)g_malloc0(npmax*sizeof(gint)) ;

  m->np = 0 ; m->npmax = npmax ;
  m->nsp = 1 ; m->nspmax = nspmax ;
  m->isp[1] = 0 ;
  m->ne = 0 ; m->nemax = nemax ;

  agg_mesh_intersection_number(m) = 0 ;  
  agg_mesh_surface_number(m) = 0 ;

  m->nb = 0 ;
  
  return m ;
}

static void write_element_gmsh(FILE *f, gint *e, gint i, gint offsp, gint offs)

{
  gint j, idx, ns ;
  
  if ( e[0] == 0 ) return ;

  ns = agg_mesh_element_size(e) ;
  fprintf(f, "Curve Loop(%d) = {", offs + i + 1) ;
  for ( j = 0 ; j < ns-1 ; j  ++ ) {
    if ( e[j] > 0 ) idx = offsp + e[j] ;
    else idx = -(offsp - e[j]) ;
    fprintf(f, "%d, ", idx) ;
  }
  j = ns-1 ;
  if ( e[j] > 0 ) idx = offsp + e[j] ;
  else idx = -(offsp - e[j]) ;
  fprintf(f, "%d} ;\n", idx) ;

  fprintf(f, "Surface(%d) = {%d} ;\n", offs + i+1, offs + i+1) ;

  return ;
}

static void write_element_gmsh_occ(FILE *f, gint *e, gint i,
				   gint offsp, gint offs)

{
  gint j, idx, ns ;
  
  if ( e[0] == 0 ) return ;

  ns = agg_mesh_element_size(e) ;
  fprintf(f, "Curve Loop(%d) = {", offs + 2*i + 1) ;
  for ( j = 0 ; j < ns-1 ; j  ++ ) {
    idx = ABS(e[j]) + offsp ;
    fprintf(f, "%d, ", idx) ;
  }
  j = ns-1 ;
  idx = ABS(e[j]) + offsp ;
  fprintf(f, "%d} ;\n", idx) ;

  fprintf(f, "Surface(%d) = {%d} ;\n", offs + 2*i + 2, offs + i + 1) ;

  return ;
}

/** 
 * Write a ::agg_mesh_t to file in GMSH .geo format.
 * 
 * @param f output file stream;
 * @param m mesh to write to file;
 * @param len string to include point definitions for local length scale;
 * @param offp point index offset (added to point indices);
 * @param offsp spline index offset (added to spline indices);
 * @param offs surface index offset (added to GMSH surface indices);
 * @param opencascade if TRUE, assume OpenCascade engine is being used.
 * 
 * @return 0 on success.
 */

gint agg_mesh_write_gmsh(FILE *f, agg_mesh_t *m,
			      gchar *len, gint offp, gint offsp, gint offs,
			      gboolean opencascade)

{
  gint i, j, n, *sp ;
  gdouble *p ;
  
  /*write points*/
  for ( i = 0 ; i < agg_mesh_point_number(m) ; i ++ ) {
    p = agg_mesh_point(m, i) ;
    fprintf(f, "Point(%d) = {%lg, %lg, %lg, %s} ;\n",
	    offp + i, p[0], p[1], p[2], len) ;
  }

  /*write splines*/
  for ( i = 1 ; i < agg_mesh_spline_number(m) ; i ++ ) {
    sp = agg_mesh_spline(m, i) ;
    n = agg_mesh_spline_length(m, i) ;
    if ( n == 2 ) {
      fprintf(f, "Line(%d) = {%d, %d} ;\n",
	      offsp+i, offp + sp[0], offp + sp[1]) ;
    } else {
      fprintf(f, "Spline(%d) = {", offsp+i) ;
      for ( j = 0 ; j < n-1 ; j ++ ) {
	fprintf(f, "%d, ", offp + sp[j]) ;
      }
      fprintf(f, "%d} ;\n", offp + sp[n-1]) ;
    }
  }

  /*write elements*/
  if ( opencascade ) {
    for ( i = 0 ; i < agg_mesh_element_number(m) ; i ++ ) {
      write_element_gmsh_occ(f, agg_mesh_element(m,i), i, offsp, offs) ;
    }
  } else {
    for ( i = 0 ; i < agg_mesh_element_number(m) ; i ++ ) {
      write_element_gmsh(f, agg_mesh_element(m,i), i, offsp, offs) ;
    }
  }
  
  return 0 ;
}

/** 
 * Find the indices of the end points of a spline in a mesh
 * 
 * @param m mesh containing the splines;
 * @param s index of spline;
 * @param p0 on exit, index of one end point of spline \a s;
 * @param p1 on exit, index of other end point of spline \a s.
 * 
 * @return 0 on success
 */

gint agg_mesh_spline_ends(agg_mesh_t *m, gint s, gint *p0, gint *p1)

{
  gint *sp, n ;

  g_assert(s < agg_mesh_spline_number(m)) ;
  g_assert(s >= 0) ;
  sp = agg_mesh_spline(m, s) ;
  n =  agg_mesh_spline_length(m, s) ;

  *p0 = sp[0] ; *p1 = sp[n-1] ;
  
  return 0 ;
}

static void mesh_spline_surface_interp(agg_mesh_t *m, gint s,
				       gint p0, gint p1, gint pps,
				       agg_surface_workspace_t *w)
{
  gint i, *sp, nsp, np ;
  gdouble s1, t1, s0, t0 ;

  nsp = agg_mesh_spline_number(m) ;
  sp = agg_mesh_spline(m, nsp) ;

  sp[0] = p0 ; sp[pps-1] = p1 ;

  s0 = agg_mesh_point_s(m, p0) ;
  t0 = agg_mesh_point_t(m, p0) ;
  s1 = agg_mesh_point_s(m, p1) ;
  t1 = agg_mesh_point_t(m, p1) ;
  
  for ( i = 1 ; i < pps-1 ; i ++ ) {
    np = agg_mesh_surface_point_add(m, s,
				    s0 + (s1 - s0)*i/(pps-1),
				    t0 + (t1 - t0)*i/(pps-1), w) ;
    sp[i] = np ;
    agg_mesh_point_tag(m,np) = s ;
  }

  agg_mesh_spline_number(m) ++ ;
  m->isp[nsp+1] = m->isp[nsp] + pps ;

  return ;
}

static void mesh_spline_surface_blend_interp(agg_mesh_t *m, gint b,
					     gint p0, gint p1, gint pps,
					     agg_surface_workspace_t *w)
{
  gint i, *sp, nsp, np ;
  gdouble s1, t1, s0, t0, s, t, *p ;
  agg_surface_blend_t *B = agg_mesh_surface_blend(m,b) ;
  
  nsp = agg_mesh_spline_number(m) ;
  sp = agg_mesh_spline(m, nsp) ;

  sp[0] = p0 ; sp[pps-1] = p1 ;

  s0 = agg_mesh_point_s(m, p0) ;
  t0 = agg_mesh_point_t(m, p0) ;
  s1 = agg_mesh_point_s(m, p1) ;
  t1 = agg_mesh_point_t(m, p1) ;
  
  for ( i = 1 ; i < pps-1 ; i ++ ) {
    s = s0 + (s1 - s0)*i/(pps-1) ;
    t = t0 + (t1 - t0)*i/(pps-1) ;
    np = agg_mesh_point_number(m) ;
    p  = agg_mesh_point(m,np) ;
    agg_surface_blend_evaluate(B, s, t, p, w) ;
    sp[i] = np ;
    agg_mesh_point_tag(m,np) = -1-b ;
    agg_mesh_point_number(m) ++ ;      
  }

  agg_mesh_spline_number(m) ++ ;
  m->isp[nsp+1] = m->isp[nsp] + pps ;

  return ;
}

/** 
 * Generate a spline by interpolating mapping between two existing points. 
 * 
 * @param m an ::agg_mesh_t;
 * @param s index in \a m of surface to use for interpolation;
 * @param p0 index of first point of new spline;
 * @param p1 index of end point of new spline;
 * @param pps total number of points on new spline;
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */

gint agg_mesh_spline_interp_points(agg_mesh_t *m, gint s,
				   gint p0, gint p1, gint pps,
				   agg_surface_workspace_t *w)
  
{
  gint b ;

  if ( s >= 0 ) {
    g_assert(s < agg_mesh_surface_number(m)) ;
    g_assert(agg_mesh_patch(m,s) != NULL) ;
    g_assert(agg_mesh_surface(m,s) != NULL) ;

    mesh_spline_surface_interp(m, s, p0, p1, pps, w) ;

    return 0 ;
  }

  /*negative s means we are looking at a blended surface*/
  b = -s - 1 ;
  g_assert(b < agg_mesh_surface_blend_number(m)) ;
  
  mesh_spline_surface_blend_interp(m, s, p0, p1, pps, w) ;
  
  return 0 ;
  
  /* nsp = agg_mesh_spline_number(m) ; */
  /* sp = agg_mesh_spline(m, nsp) ; */

  /* sp[0] = p0 ; sp[pps-1] = p1 ; */

  /* s0 = agg_mesh_point_s(m, p0) ; */
  /* t0 = agg_mesh_point_t(m, p0) ; */
  /* s1 = agg_mesh_point_s(m, p1) ; */
  /* t1 = agg_mesh_point_t(m, p1) ; */
  
  /* for ( i = 1 ; i < pps-1 ; i ++ ) { */
  /*   np = agg_mesh_surface_point_add(m, s, */
  /* 				    s0 + (s1 - s0)*i/(pps-1), */
  /* 				    t0 + (t1 - t0)*i/(pps-1), w) ; */
  /*   sp[i] = np ; */
  /*   agg_mesh_point_tag(m,np) = s ; */
  /* } */

  /* agg_mesh_spline_number(m) ++ ; */
  /* m->isp[nsp+1] = m->isp[nsp] + pps ; */
  
  /* return 0 ; */
}

/** 
 * Check if two splines in an ::agg_mesh_t are connected.
 * 
 * @param m an ::agg_mesh_t;
 * @param i0 index of first spline to check;
 * @param i1 index of second spline to check;
 * @param p on exit, set to index of common point of splines \a i0 and \a i1. 
 * 
 * @return TRUE if splines have an end point in common. 
 */

gboolean agg_mesh_splines_connected(agg_mesh_t *m,
				    gint i0, gint i1, gint *p)

{
  gint i0p0, i0p1, i1p0, i1p1 ;

  agg_mesh_spline_ends(m, i0, &i0p0, &i0p1) ;  
  agg_mesh_spline_ends(m, i1, &i1p0, &i1p1) ;  

  if ( i0p0 == i1p0 ) { *p = i0p0 ; return TRUE ; }
  if ( i0p0 == i1p1 ) { *p = i0p0 ; return TRUE ; }
  if ( i0p1 == i1p0 ) { *p = i0p1 ; return TRUE ; }
  if ( i0p1 == i1p1 ) { *p = i0p1 ; return TRUE ; }
  
  return FALSE ;
}


/** 
 * Add an element to a mesh.
 * 
 * @param m an ::agg_mesh_t to have an element added.
 * @param s0 index of first spline of element;
 * @param s1 index of second spline of element;
 * @param s2 index of third spline of element;
 * @param s3 index of fourth spline of element (if zero, a 
 * three-sided element is added)
 * 
 * @return 0 on success.
 */

gint agg_mesh_element_add(agg_mesh_t *m,
			  gint s0, gint s1, gint s2, gint s3)

{
  gint ne, *el ;

  ne = agg_mesh_element_number(m) ;
  el = agg_mesh_element(m, ne) ;

  el[0] = s0 ; el[1] = s1 ; el[2] = s2 ; el[3] = s3 ; el[4] = 0 ;

  /*check connections*/
  if ( !agg_mesh_splines_connected(m, ABS(s0), ABS(s1), &ne) ) {
    g_error("%s: splines s0 (%d) and s1 (%d) not connected",
	    __FUNCTION__, s0, s1) ;
  }
  if ( !agg_mesh_splines_connected(m, ABS(s1), ABS(s2), &ne) ) {
    g_error("%s: splines s1 (%d) and s2 (%d) not connected",
	    __FUNCTION__, s1, s2) ;
  }
  if ( s3 != 0 ) {
    if ( !agg_mesh_splines_connected(m, ABS(s2), ABS(s3), &ne) ) {
      g_error("%s: splines s2 (%d) and s3 (%d) not connected",
	      __FUNCTION__, s2, s3) ;
    }
    if ( !agg_mesh_splines_connected(m, ABS(s3), ABS(s0), &ne) ) {
      g_error("%s: splines s3 (%d) and s0 (%d) not connected",
	      __FUNCTION__, s3, s0) ;
    }
  } else {
    if ( !agg_mesh_splines_connected(m, ABS(s2), ABS(s0), &ne) ) {
      g_error("%s: splines s2 (%d) and s0 (%d) not connected",
	      __FUNCTION__, s2, s0) ;
    }
  }
  
  agg_mesh_element_number(m) ++ ;
  
  return 0 ;
}

static void reset_triangleio(triangleio *io)

{
  io->pointlist = (REAL *) NULL;
  io->pointattributelist = (REAL *) NULL;
  io->pointmarkerlist = (int *) NULL;
  io->numberofpoints = 0;
  io->numberofpointattributes = 0;
  
  io->trianglelist = (int *) NULL;
  io->triangleattributelist = (REAL *) NULL;
  io->trianglearealist = (REAL *) NULL;
  io->neighborlist = (int *) NULL;
  io->numberoftriangles = 0;
  io->numberofcorners = 0;
  io->numberoftriangleattributes = 0;
  
  io->segmentlist = (int *) NULL;
  io->segmentmarkerlist = (int *) NULL;
  io->numberofsegments = 0;
  
  io->holelist = (REAL *) NULL;
  io->numberofholes = 0;
  io->regionlist = (REAL *) NULL;
  io->numberofregions = 0;
  
  io->edgelist = (int *) NULL;
  io->edgemarkerlist = (int *) NULL;
  io->numberofedges = 0;

  return ;
}

static void add_point(gdouble *st, gint *nst, gdouble s, gdouble t)

{
  st[2*(*nst)+0] = s ; st[2*(*nst)+1] = t ;

  (*nst) ++ ;
  
  return ;
}

static void set_corners(gdouble *st, gint *nst,
			gint *segments, gint *nseg,
			gdouble smin, gdouble tmin,
			gdouble smax, gdouble tmax)

{
  gint nst0 = (*nst) ;
  
  add_point(st, nst, smin, tmin) ; add_point(st, nst, smax, tmin) ;
  add_point(st, nst, smax, tmax) ; add_point(st, nst, smin, tmax) ;
  segments[2*(*nseg)+0] = nst0+0 ; segments[2*(*nseg)+1] = nst0+1 ;
  (*nseg)++ ;
  segments[2*(*nseg)+0] = nst0+1 ; segments[2*(*nseg)+1] = nst0+2 ; 
  (*nseg)++ ;
  segments[2*(*nseg)+0] = nst0+2 ; segments[2*(*nseg)+1] = nst0+3 ; 
  (*nseg)++ ;
  segments[2*(*nseg)+0] = nst0+3 ; segments[2*(*nseg)+1] = nst0+0 ; 
  (*nseg)++ ;

  return ;
}

static void set_boundary(agg_mesh_t *m, gint isurf,
			 gdouble *st, gint *nst,
			 gint *segments, gint *nseg,
			 gint *segmentmarkers,
			 gdouble *holes, gint *nholes,
			 gdouble del)

/*
 * m:        mesh
 * isurf:    index of surface in mesh surface list
 * st:       list of points in patch coordinate system (s,t)
 * nst:      number of points
 * segments: segment list, two ints per segment (indices of points)
 * nseg:     number of segments
 * idxlist:  indices of mesh points lying on intersection curve
 * nidx:     number of points in idxlist
 * del:      offset from patch ends (to leave clear space for singular points)
 */
  
{
  gint i, j, nst0 ;
  agg_patch_t *P ;
  gboolean corners_set ;
  agg_patch_clipping_t *clip ;
  
  corners_set = FALSE ;
  
  /*check for clipped patch*/
  P = agg_mesh_patch(m, isurf) ;
  if ( agg_patch_clipping_number(P) == 0 ) {
    set_corners(st, nst, segments, nseg, 0, 0, 1, 1) ;
    return ;
  }

  for ( i = 0 ; i < agg_patch_clipping_number(P) ; i ++ ) {
    clip = agg_patch_clipping(P,i) ;
    if ( agg_patch_clipping_type(clip) == AGG_CLIP_CONSTANT_T ) {
      set_corners(st, nst, segments, nseg,
		  0, agg_patch_clipping_data(clip,0),
		  1, 1) ;
      corners_set = TRUE ;
    }
    if ( agg_patch_clipping_type(clip) == AGG_CLIP_CONSTANT_S ) {
      set_corners(st, nst, segments, nseg,
		  agg_patch_clipping_data(clip,0), 0,
		  1, 1) ;
      corners_set = TRUE ;
    }
    if ( agg_patch_clipping_type(clip) == AGG_CLIP_ELLIPSE ) {
      gint nc = 33 ;
      gdouble s, t ;
      nst0 = (*nst) ;
      for ( j = 0 ; j < nc ; j ++ ) {
	agg_patch_clip_eval(agg_patch_clipping(P,i), (gdouble)j/nc, &s, &t) ;
	add_point(st, nst, s, t) ;
      }
      
      for ( j = 0 ; j < nc - 1 ; j ++ ) {
	segments[2*(*nseg)+0] = nst0 + j + 0 ; 
	segments[2*(*nseg)+1] = nst0 + j + 1 ; 
	(*nseg) ++ ;      
      }    
      segments[2*(*nseg)+0] = nst0 + nc - 1 ;
      segments[2*(*nseg)+1] = nst0 ; 
      (*nseg) ++ ;

      /*centre of hole in mesh is centre of ellipse*/
      holes[(*nholes)*2+0] =
	agg_patch_clipping_data(agg_patch_clipping(P,i),0) ;
      holes[(*nholes)*2+1] = 
	agg_patch_clipping_data(agg_patch_clipping(P,i),1) ;
      (*nholes) ++ ;
    }
  }

  if ( !corners_set ) {
    set_corners(st, nst, segments, nseg, 0, 0, 1, 1) ;
  }

  return ;
}

static gboolean boundary_point(gint *b, gint nb, gint p)

{
  gint i ;

  for ( i = 0 ; i < nb ; i ++ ) { if ( b[i] == p ) return TRUE ; }
  
  return FALSE ;
}

/** 
 * Add a parametric surface to a mesh, working round clipped parts of
 * patch, using Christian Woltering's library version of Jonathan
 * Shewchuk's Triangle code:
 * 
 * https://github.com/wo80/Triangle/
 * https://www.cs.cmu.edu/~quake/triangle.html
 * 
 * @param msh a mesh to have a surface added;
 * @param isurf index of surface in mesh surface list;
 * @param args argument string to pass to triangle;
 * @param pps points per spline;
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */


gint agg_mesh_surface_add_triangle(agg_mesh_t *msh, gint isurf,
				   gchar *args, gint pps,
				   agg_surface_workspace_t *w)

{
  context *ctx ;
  triangleio in ;
  /*it has to be called this because of how the macros in triangle.h
    are defined*/
  mesh *m ;
  behavior *b ;
  statistics s ;
  gint np, np0, i0, i1, i2, j0, j1, j2, nbpts ;
  gint *sp, *e, i, vertexnumber ;
  gint boundary[8192], nbound, tag ;
  struct otri triangleloop, trisym;
  struct osub checkmark;
  vertex p1, p2, p3, vertexloop ;
  glong edgenumber, elementnumber ;
  triangle ptr; 
  subseg sptr;  
  gdouble del = 0.0625 ;
  gdouble s0, t0, s1, t1, s2, t2 ;
  agg_patch_t *P ;
  del = 0 ;
  
  P = agg_mesh_patch(msh,isurf) ;

  /*maximum number of boundary points*/
  nbpts = 1024 ;
  for ( i = 0 ; i < agg_mesh_intersection_number(msh) ; i ++ ) {
    nbpts += agg_mesh_intersection_points_end(msh, i) -
      agg_mesh_intersection_points_start(msh, i) ;
  }
  nbpts = MAX(nbpts, 2048) ;
  
  ctx = triangle_context_create() ;
  triangle_context_options(ctx, args) ;
  /* triangle_context_options(ctx, "czqa0.001"); */
  reset_triangleio(&in);

  in.pointlist = (gdouble *)g_malloc0(nbpts*2*sizeof(gdouble)) ;
  in.pointmarkerlist = (gint *)g_malloc0(nbpts*sizeof(gint)) ;
  in.segmentlist = (gint *)g_malloc0(nbpts*2*sizeof(gint)) ;
  in.segmentmarkerlist = (gint *)g_malloc0(nbpts*sizeof(gint)) ;
  in.holelist = (gdouble *)g_malloc0(16*2*sizeof(gdouble)) ;

  in.numberofpoints = in.numberofsegments = 0 ;
  in.numberofholes = 0 ;

  set_boundary(msh, isurf,
	       in.pointlist, &(in.numberofpoints),
	       in.segmentlist, &(in.numberofsegments),
	       in.segmentmarkerlist,
	       in.holelist, &(in.numberofholes),
	       del) ;
  
  np0 = in.numberofpoints ;

  triangle_mesh_create(ctx, &in) ;
  m = ctx->m ;
  b = ctx->b ;
  
  triangle_mesh_statistics(ctx, &s);

  /*
   * code adapted from Christian Woltering's Triangle API library:
   * https://github.com/wo80/Triangle/
   */
  traversalinit(&m->vertices);
  vertexnumber = b->firstnumber;
  g_assert(vertexnumber == 0) ;
  vertexloop = vertextraverse(m);
  np0 = agg_mesh_point_number(msh) ;
  nbound = 0 ;
  while (vertexloop != (vertex) NULL) {
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
      s0 = vertexloop[0] ; t0 = vertexloop[1] ;
      np = agg_mesh_surface_point_add(msh, isurf, s0, t0, w) ;
      agg_mesh_point_tag(msh,np) = isurf ;
      tag = vertexmark(vertexloop) ;
      setvertexmark(vertexloop, vertexnumber) ;
      if ( tag == 1 ) { boundary[nbound] = np ; nbound ++ ; }
      vertexnumber++ ;
    }
    vertexloop = vertextraverse(m);
  }

  /*generate edge splines*/
  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  edgenumber = b->firstnumber;
  /* To loop over the set of edges, loop over all triangles, and look at   */
  /*   the three edges of each triangle.  If there isn't another triangle  */
  /*   adjacent to the edge, operate on the edge.  If there is another     */
  /*   adjacent triangle, operate on the edge only if the current triangle */
  /*   has a smaller pointer than its neighbor.  This way, each edge is    */
  /*   considered only once.                                               */
  while (triangleloop.tri != (triangle *) NULL) {
    for (triangleloop.orient = 0; triangleloop.orient < 3;
	 triangleloop.orient++) {
      sym(triangleloop, trisym);
      if ((triangleloop.tri < trisym.tri) || (trisym.tri == m->dummytri)) {
	tspivot(triangleloop, checkmark);
	org(triangleloop, p1);
	dest(triangleloop, p2);
	i0 = np0 + vertexmark(p1) ; i1 = np0 + vertexmark(p2) ;
	agg_mesh_spline_interp_points(msh, isurf, i0, i1, pps, w) ;
	sp = agg_mesh_spline(msh,agg_mesh_spline_number(msh)-1) ;
	sp[0] = i0 ; sp[pps-1] = i1 ;
	edgenumber++;
      }
    }
    triangleloop.tri = triangletraverse(m);
  }

  traversalinit(&m->triangles);
  triangleloop.tri = triangletraverse(m);
  triangleloop.orient = 0;
  elementnumber = b->firstnumber;
  while (triangleloop.tri != (triangle *) NULL) {
    org(triangleloop, p1);
    dest(triangleloop, p2);
    apex(triangleloop, p3);
    /* Triangle number, indices for three vertices. */
    i0 = vertexmark(p1) ;
    i1 = vertexmark(p2) ;
    i2 = vertexmark(p3) ;

    j0 = i0 + np0 ; j1 = i1 + np0 ; j2 = i2 + np0 ; 
    
    i0 = agg_mesh_spline_from_endpoints(msh, j0, j1) ;
    if ( i0 != 0 ) {
      i1 = agg_mesh_spline_from_endpoints(msh, j1, j2) ;
      if ( i1 != 0 ) {
	i2 = agg_mesh_spline_from_endpoints(msh, j2, j0) ;
	if ( i2 != 0 ) {
	  e = agg_mesh_element(msh, agg_mesh_element_number(msh)) ;
	  if ( agg_patch_invert(P) ) {
	    e[0] = -i2 ; e[1] = -i1 ; e[2] = -i0 ; e[3] = 0 ;
	  } else {
	    e[0] = i0 ; e[1] = i1 ; e[2] = i2 ; e[3] = 0 ;
	  }
	  agg_mesh_element_number(msh) ++ ;
	}
      }
    }
    
    triangleloop.tri = triangletraverse(m);
    elementnumber++;
  }
  
  triangle_context_destroy(ctx) ;
  /*** adapted code ends here ****/
  
  if ( agg_patch_mapping(P) == AGG_PATCH_TUBULAR ||
       agg_patch_mapping(P) == AGG_PATCH_BILINEAR ) return 0 ;

  if ( del == 0 ) return 0 ;
  
  /*add the end cap to deal with singularity on surfaces*/
  for ( i = 1 ; i < agg_mesh_spline_number(msh) ; i ++ ) {
    agg_mesh_spline_ends(msh, i, &i0, &i1) ;
    
    if ( boundary_point(boundary, nbound, i0) &&
	 boundary_point(boundary, nbound, i1) ) {
      s0 = agg_mesh_point_s(msh, i0) ;
      t0 = agg_mesh_point_t(msh, i0) ;
      s1 = agg_mesh_point_s(msh, i1) ;
      t1 = agg_mesh_point_t(msh, i1) ;
      if ( s0 == 1.0-del && s1 == 1.0-del ) {
	s2 = 1 ; t2 = 0.5*(t0+t1) ;
	np = agg_mesh_surface_point_add(msh, isurf, s2, t2, w) ;
	agg_mesh_spline_interp_points(msh, isurf, i0 , np, pps, w) ;
	i0 = agg_mesh_spline_number(msh) - 1 ;
	agg_mesh_spline_interp_points(msh, isurf, i1, np, pps, w) ;
	i1 = agg_mesh_spline_number(msh) - 1 ;
	agg_mesh_element_add(msh, -i, i0, -i1, 0) ;
      }
      if ( s0 == del && s1 == del ) {
	s2 = 0 ; t2 = 0.5*(t0+t1) ;
	np = agg_mesh_surface_point_add(msh, isurf, s2, t2, w) ;
	agg_mesh_spline_interp_points(msh, isurf, i0 , np, pps, w) ;
	i0 = agg_mesh_spline_number(msh) - 1 ;
	agg_mesh_spline_interp_points(msh, isurf, i1, np, pps, w) ;
	i1 = agg_mesh_spline_number(msh) - 1 ;
	agg_mesh_element_add(msh, -i, i0, -i1, 0) ;
      }
      
    }
  }
  
  return 0 ;
}

/** 
 * Add elements from a surface blend to a mesh
 * 
 * @param m mesh to add surface elements to;
 * @param iB index of surface blend in \a m;
 * @param nsec number of sections on blend;
 * @param pps points per spline on element edges;
 * @param w workspace for point evaluation
 * 
 * @return 0 on success.
 */

gint agg_mesh_surface_blend_add(agg_mesh_t *m, gint iB, gint nsec, gint pps,
				agg_surface_workspace_t *w)

{
  gint i, j, np, np0, tag, nsp0, s0, s1, s2, s3, nsp ;
  gdouble s, t, *p ;
  agg_surface_blend_t *B ;

  g_assert(iB >= 0 && iB < agg_mesh_surface_blend_number(m)) ;
  
  B = agg_mesh_surface_blend(m,iB) ;
  nsp = agg_surface_blend_spline_number(B) ;
  agg_surface_blend_invert(B) = FALSE ;
  i = B->ic[0] ; j = B->ic[1] ; 
  if ( agg_patch_clipping_orientation(agg_patch_clipping(B->P[0], i)) !=
       agg_patch_clipping_orientation(agg_patch_clipping(B->P[1], j)) ) {
    agg_surface_blend_invert(B) = TRUE ;
  }
  
  tag = -1-iB ;
  np0 = agg_mesh_point_number(m) ;
  nsp0 = agg_mesh_spline_number(m) ;
  for ( i = 0 ; i <= nsec ; i ++ ) {
    for ( j = 0 ; j <= nsp ; j ++ ) {
      s = (gdouble)i/nsec ; t = (gdouble)j/nsp ;
      np = agg_mesh_point_number(m) ;
      p  = agg_mesh_point(m,np) ;
      agg_surface_blend_evaluate(B, s, t, p, w) ;
      agg_mesh_point_s(m, np) = s ;
      agg_mesh_point_t(m, np) = t ;
      agg_mesh_point_number(m) ++ ;      
    }
  }

  /*splines joining element corners*/
  for ( i = 0 ; i < nsec ; i ++ ) {
    for ( j = 0 ; j < nsp ; j ++ ) {
      agg_mesh_spline_interp_points(m, tag,
				    np0 + (i+0)*(nsp+1) + j + 0,
				    np0 + (i+1)*(nsp+1) + j + 0, pps, w) ;
      agg_mesh_spline_interp_points(m, tag,
				    np0 + (i+0)*(nsp+1) + j + 0,
				    np0 + (i+0)*(nsp+1) + j + 1, pps, w) ;
    }
  }
  i = nsec ;
  for ( j = 0 ; j < nsp ; j ++ ) {
    agg_mesh_spline_interp_points(m, tag,
				  np0 + (i+0)*(nsp+1) + j + 0,
				  np0 + (i+0)*(nsp+1) + j + 1, pps, w) ;
  }
  j = nsp ;
  for ( i = 0 ; i < nsec ; i ++ ) {
    agg_mesh_spline_interp_points(m, tag,
				  np0 + (i+0)*(nsp+1) + j + 0,
				  np0 + (i+1)*(nsp+1) + j + 0, pps, w) ;
  }    
  
  for ( i = np0 ; i < agg_mesh_point_number(m) ; i ++ ) {
    agg_mesh_point_tag(m, i) = tag ;
  }

    /*elements*/
  for ( i = 0 ; i < nsec - 1 ; i ++ ) {
    for ( j = 0 ; j < nsp - 1 ; j ++ ) {
      s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
      s1 = nsp0 + 2*(i+1)*nsp + 2*(j+0) + 1 ;
      s2 = nsp0 + 2*(i+0)*nsp + 2*(j+1) + 0 ;
      s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;
      
      if ( agg_surface_blend_invert(B) ) {
	agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
      } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
      }
    }
  }

  i = nsec - 1 ;
  for ( j = 0 ; j < nsp - 1 ; j ++ ) {
    s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nsp +   (j+0) + 0 ;
    s2 = nsp0 + 2*(i+0)*nsp + 2*(j+1) + 0 ;
    s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;
    if ( agg_surface_blend_invert(B) ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  j = nsp - 1 ;
  for ( i = 0 ; i < nsec - 1 ; i ++ ) {
    s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nsp + 2*(j+0) + 1 ;
    s2 = nsp0 + 2*nsec*nsp + j + 1 + i ;
    s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;
    
    if ( agg_surface_blend_invert(B) ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  i = nsec - 1 ; j = nsp - 1 ;
  s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
  s1 = nsp0 + 2*(i+1)*nsp + (j+0) + 0 ;
  s2 = nsp0 + 2*nsec*nsp + j + 1 + i ;
  s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;    
  if ( agg_surface_blend_invert(B) ) {
    agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
  } else {
    agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
  }

  return 0 ;
}

/** 
 * Find a spline from its endpoints
 * 
 * @param m a mesh containing a set of splines;
 * @param p0 index of a point of \a m;
 * @param p1 index of a point of \a m.
 * 
 * @return index of a spline which has endpoints \a p0 and \a p1, or 0
 * if there is no such spline.
 */

gint agg_mesh_spline_from_endpoints(agg_mesh_t *m, gint p0, gint p1)

{
  gint i, *sp, len ;

  for ( i = 1 ; i < agg_mesh_spline_number(m) ; i ++ ) {
    sp = agg_mesh_spline(m, i) ;
    len = agg_mesh_spline_length(m, i) ;
    if ( sp[0] == p0 && sp[len-1] == p1 ) return  i ;
    if ( sp[0] == p1 && sp[len-1] == p0 ) return -i ;
  }
  
  return 0 ;
}

/** 
 * Add a point to a mesh
 * 
 * @param m mesh to have point added;
 * @param surf index in \a m of surface on which to evaluate point;
 * @param s coordinate on patch for surface \a surf;
 * @param t coordinate on patch for surface \a surf;
 * @param w workspace for surface evaluation.
 * 
 * @return index of newly added point.
 */

gint agg_mesh_surface_point_add(agg_mesh_t *m, gint surf,
				gdouble s, gdouble t,
				agg_surface_workspace_t *w)

{
  gint np ;
  gdouble *p, u, v ;
  agg_surface_t *S ;
  agg_patch_t *P ;
  
  np = agg_mesh_point_number(m) ;

  if ( np >= agg_mesh_point_number_max(m) ) {
    g_error("%s: not enough space allocated (%d points)",
	    __FUNCTION__, agg_mesh_point_number_max(m)) ;
  }
  
  p = agg_mesh_point(m, np) ;
  P = agg_mesh_patch(m, surf) ;
  g_assert((S = agg_mesh_surface(m, surf)) != NULL) ;
  
  if ( P != NULL ) {
    agg_patch_map(P, s, t, &u, &v) ;    
  } else {
    u = s ; v = t ;
  }
  agg_surface_point_eval(S, u, v, p, w) ;    
  p[3] = s ; p[4] = t ;
  agg_mesh_point_tag(m,np) = surf ;

  agg_mesh_point_number(m) ++ ;
  
  return np ;
}

/** 
 * Add surface elements to a mesh, on a square grid in the parametric space 
 * 
 * @param m mesh to add elements to;
 * @param iS index of surface in surface list of \a m;
 * @param nsec number of sections on surface;
 * @param nsp number of splines per section;
 * @param pps points per spline on sections;
 * @param w workspace for point evaluation.
 * 
 * @return 0 on success.
 */

gint agg_mesh_surface_add_grid(agg_mesh_t *m, gint iS,
			       gint nsec, gint nsp, gint pps,
			       agg_surface_workspace_t *w)

{
  gint i, j, np0, nsp0, s0, s1, s2, s3 ;
  gdouble s, t, tmin, tmax, smin, smax ;
  agg_patch_t *P ;
  
  np0 = agg_mesh_point_number(m) ;
  nsp0 = agg_mesh_spline_number(m) ;

  smin = tmin = 0.0 ; smax = tmax = 1.0 ;
  P = agg_mesh_patch(m, iS) ;
  if ( agg_patch_clipping_number(P) > 0 ) {
    if ( agg_patch_clipping_type(agg_patch_clipping(P,0)) ==
	 AGG_CLIP_CONSTANT_S ) {
      smin = agg_patch_clipping_data(agg_patch_clipping(P,0),0) ;
    }
    if ( agg_patch_clipping_type(agg_patch_clipping(P,0)) ==
	 AGG_CLIP_CONSTANT_T ) {
      tmin = agg_patch_clipping_data(agg_patch_clipping(P,0),0) ;
    }
  }
  
  for ( i = 0 ; i <= nsec ; i ++ ) {
    for ( j = 0 ; j <= nsp ; j ++ ) {
      s = smin + (smax - smin)*(gdouble)i/nsec ;
      t = tmin + (tmax - tmin)*(gdouble)j/nsp ;
      agg_mesh_surface_point_add(m, iS, s, t, w) ;
    }
  }

  /*splines joining element corners*/
  for ( i = 0 ; i < nsec ; i ++ ) {
    for ( j = 0 ; j < nsp ; j ++ ) {
      agg_mesh_spline_interp_points(m, iS,
				    np0 + (i+0)*(nsp+1) + j + 0,
				    np0 + (i+1)*(nsp+1) + j + 0, pps, w) ;
      agg_mesh_spline_interp_points(m, iS,
				    np0 + (i+0)*(nsp+1) + j + 0,
				    np0 + (i+0)*(nsp+1) + j + 1, pps, w) ;
    }
  }
  i = nsec ;
  for ( j = 0 ; j < nsp ; j ++ ) {
    agg_mesh_spline_interp_points(m, iS,
				  np0 + (i+0)*(nsp+1) + j + 0,
				  np0 + (i+0)*(nsp+1) + j + 1, pps, w) ;
  }
  j = nsp ;
  for ( i = 0 ; i < nsec ; i ++ ) {
    agg_mesh_spline_interp_points(m, iS,
				  np0 + (i+0)*(nsp+1) + j + 0,
				  np0 + (i+1)*(nsp+1) + j + 0, pps, w) ;
  }    

  for ( i = np0 ; i < agg_mesh_point_number(m) ; i ++ ) {
    agg_mesh_point_tag(m, i) = iS ;
  }

  /*elements*/
  for ( i = 0 ; i < nsec - 1 ; i ++ ) {
    for ( j = 0 ; j < nsp - 1 ; j ++ ) {
      s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
      s1 = nsp0 + 2*(i+1)*nsp + 2*(j+0) + 1 ;
      s2 = nsp0 + 2*(i+0)*nsp + 2*(j+1) + 0 ;
      s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;
      
      if ( agg_patch_invert(agg_mesh_patch(m, iS)) ) {
	agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
      } else {
	agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
      }
    }
  }

  i = nsec - 1 ;
  for ( j = 0 ; j < nsp - 1 ; j ++ ) {
    s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nsp +   (j+0) + 0 ;
    s2 = nsp0 + 2*(i+0)*nsp + 2*(j+1) + 0 ;
    s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;
    if ( agg_patch_invert(agg_mesh_patch(m, iS)) ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  j = nsp - 1 ;
  for ( i = 0 ; i < nsec - 1 ; i ++ ) {
    s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nsp + 2*(j+0) + 1 ;
    s2 = nsp0 + 2*nsec*nsp + j + 1 + i ;
    s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;
    
    if ( agg_patch_invert(agg_mesh_patch(m, iS)) ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  i = nsec - 1 ; j = nsp - 1 ;
  s0 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 0 ;
  s1 = nsp0 + 2*(i+1)*nsp + (j+0) + 0 ;
  s2 = nsp0 + 2*nsec*nsp + j + 1 + i ;
  s3 = nsp0 + 2*(i+0)*nsp + 2*(j+0) + 1 ;    
  if ( agg_patch_invert(agg_mesh_patch(m, iS)) ) {
    agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
  } else {
    agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
  }
  
  return 0 ;
}

/** 
 * Generate a mesh for surfaces and surface blends on a body
 * 
 * @param m mesh to generate;
 * @param b body containing surface data;
 * @param pps points per spline on element edges;
 * @param w workspace for point evaluation.
 * 
 * @return 0 on success.
 */

gint agg_mesh_body(agg_mesh_t *m, agg_body_t *b, gint pps,
		   agg_surface_workspace_t *w)

{
  gint i, j, nsec, nsp ;
  gdouble area ;
  gchar args[64] ;
  agg_intersection_t *inter ;

  nsp = 0 ;
  agg_mesh_surface_number(m) = 0 ;
  for ( i = 0 ; i < agg_body_surface_number(b); i ++ ) {
    agg_mesh_surface(m,i) = agg_body_surface(b,i) ;
    agg_mesh_patch(m,i) = agg_body_patch(b,i) ;
    agg_mesh_surface_number(m) ++ ;
  }
  
  inter = agg_intersection_new(8192) ;
  for ( i = 0 ; i < agg_body_surface_number(b); i ++ ) {
    for ( j = i+1 ; j < agg_body_surface_number(b); j ++ ) {
      agg_surface_patch_intersection(inter,
				     agg_body_surface(b,i),
				     agg_body_patch(b,i),
				     agg_body_surface(b,j),
				     agg_body_patch(b,j),
				     &(m->B[m->nb]), w) ;
      if ( agg_intersection_point_number(inter) != 0 ) {
	inter = agg_intersection_new(8192) ;
	m->nb ++ ;
      }
    }
  }

  for ( i = 0 ; i < agg_mesh_surface_number(m) ; i ++ ) {
    g_assert(agg_surface_grid(agg_body_surface(b,i)) != AGG_GRID_UNDEFINED) ;
    if ( agg_surface_grid(agg_body_surface(b,i)) == AGG_GRID_REGULAR ) {
      nsec = agg_surface_grid_section_number(agg_body_surface(b,i)) ;
      nsp = agg_surface_grid_spline_number(agg_body_surface(b,i)) ;      
      agg_mesh_surface_add_grid(m, i, nsec, nsp, pps, w) ;
    }
    if ( agg_surface_grid(agg_body_surface(b,i)) == AGG_GRID_TRIANGLE ) {
      area = agg_surface_grid_element_area(agg_body_surface(b,i)) ;
      sprintf(args, "pzqa%lg", area) ;
      agg_mesh_surface_add_triangle(m, i, args, pps, w) ;
    }
  }

  for ( i = 0 ; i < agg_mesh_surface_blend_number(m) ; i ++ ) {
    agg_mesh_surface_blend_add(m, i, 4, pps, w) ;
  }

  return 0 ;
}

/** 
 * Find the oriented nodes of an element in an ::agg_mesh_t
 * 
 * @param m mesh containing elements to be queried;
 * @param e element index;
 * @param nodes on exit contains nodes of element \e in cyclic order;
 * @param nnodes on exit contains number of nodes in \a nodes;
 * @param s index in \a m of surface containing nodes.
 * 
 * @return 0 on success.
 */

gint agg_mesh_element_nodes(agg_mesh_t *m, gint e,
			    gint *nodes, gint *nnodes, gint *s)

{
  gint i, *el, *sp, len, tag ;

  el = agg_mesh_element(m, e) ;
  *nnodes = agg_mesh_element_size(el) ;

  for ( i = 0 ; i < (*nnodes) ; i ++ ) {
    sp = agg_mesh_spline(m, ABS(el[i])) ;
    len = agg_mesh_spline_length(m, ABS(el[i])) ;
    if ( el[i] > 0 ) {
      nodes[i] = sp[0] ;
    } else {
      nodes[i] = sp[len-1] ;
    }
    tag = agg_mesh_point_tag(m, nodes[i]) ;
    if ( i > 0 ) {
      if ( (*s) != tag )
	g_error("%s: element nodes not all from same surface",
		__FUNCTION__) ;
    }
    *s = tag ;
  }

  for ( i = 0 ; i < (*nnodes)-1 ; i ++ ) {
    gint sc ;
    sc = agg_mesh_spline_from_endpoints(m, nodes[i], nodes[i+1]) ;
    if ( ABS(sc) != ABS(el[i]) ) {
      g_error("%s: nodes do not match spline endpoints", __FUNCTION__) ;
    }
  }
  
  return 0 ;
}

/**
 * @}
 */
