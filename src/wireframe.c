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

agg_wireframe_t *agg_wireframe_new(gint npmax, gint nspmax, gint nemax)

{
  agg_wireframe_t *w ;

  w = (agg_wireframe_t *)g_malloc0(sizeof(agg_wireframe_t)) ;

  w->p   = (gdouble *)g_malloc0(AGG_WIREFRAME_POINT_SIZE*npmax*
				sizeof(gdouble)) ;
  w->sp  = (gint *)g_malloc0(nspmax*8*sizeof(gint)) ;
  w->isp = (gint *)g_malloc0(nspmax*sizeof(gint)) ;
  w->e   = (gint *)g_malloc0(AGG_WIREFRAME_ELEMENT_SIZE*nemax*sizeof(gint)) ;
  w->ptags = (gint *)g_malloc0(npmax*sizeof(gint)) ;

  w->np = 0 ; w->npmax = npmax ;
  w->nsp = 1 ; w->nspmax = nspmax ;
  w->isp[1] = 0 ;
  w->ne = 0 ; w->nemax = nemax ;

  agg_wireframe_intersection_number(w) = 0 ;  
  agg_wireframe_surface_number(w) = 0 ;
  
  return w ;
}

static void write_element_gmsh(FILE *f, gint *e, gint i, gint offsp, gint offs)

{
  gint j, idx, ns ;
  
  if ( e[0] == 0 ) return ;

  ns = agg_wireframe_element_size(e) ;
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

  ns = agg_wireframe_element_size(e) ;
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

gint agg_wireframe_write_gmsh(FILE *f, agg_wireframe_t *w,
			      gchar *len, gint offp, gint offsp, gint offs,
			      gboolean opencascade)

{
  gint i, j, n, *sp ;
  gdouble *p ;
  
  /*write points*/
  for ( i = 0 ; i < agg_wireframe_point_number(w) ; i ++ ) {
    p = agg_wireframe_point(w, i) ;
    fprintf(f, "Point(%d) = {%lg, %lg, %lg, %s} ;\n",
	    offp + i, p[0], p[1], p[2], len) ;
  }

  /*write splines*/
  for ( i = 1 ; i < agg_wireframe_spline_number(w) ; i ++ ) {
    sp = agg_wireframe_spline(w, i) ;
    n = agg_wireframe_spline_length(w, i) ;
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
    for ( i = 0 ; i < agg_wireframe_element_number(w) ; i ++ ) {
      write_element_gmsh_occ(f, agg_wireframe_element(w,i), i, offsp, offs) ;
    }
  } else {
    for ( i = 0 ; i < agg_wireframe_element_number(w) ; i ++ ) {
      write_element_gmsh(f, agg_wireframe_element(w,i), i, offsp, offs) ;
    }
  }
  
  return 0 ;
}

gint agg_wireframe_spline_ends(agg_wireframe_t *w, gint s, gint *p0, gint *p1)

{
  gint *sp, n ;

  g_assert(s < agg_wireframe_spline_number(w)) ;
  g_assert(s >= 0) ;
  sp = agg_wireframe_spline(w, s) ;
  n =  agg_wireframe_spline_length(w, s) ;

  *p0 = sp[0] ; *p1 = sp[n-1] ;
  
  return 0 ;
}
  
gint agg_wireframe_spline_interp_points(agg_wireframe_t *w, gint s,
					gint p0, gint p1, gint pps,
					agg_surface_workspace_t *ws)

{
  gint i, *sp, nsp, np ;
  gdouble s1, t1, s0, t0 ;

  g_assert(s < agg_wireframe_surface_number(w)) ;
  g_assert(agg_wireframe_patch(w,s) != NULL) ;
  g_assert(agg_wireframe_surface(w,s) != NULL) ;

  nsp = agg_wireframe_spline_number(w) ;
  sp = agg_wireframe_spline(w, nsp) ;

  sp[0] = p0 ; sp[pps-1] = p1 ;

  s0 = agg_wireframe_point_s(w, p0) ;
  t0 = agg_wireframe_point_t(w, p0) ;
  s1 = agg_wireframe_point_s(w, p1) ;
  t1 = agg_wireframe_point_t(w, p1) ;
  
  for ( i = 1 ; i < pps-1 ; i ++ ) {
    np = agg_wireframe_surface_point_add(w, s,
					 s0 + (s1 - s0)*i/(pps-1),
					 t0 + (t1 - t0)*i/(pps-1), ws) ;
    sp[i] = np ;
    agg_wireframe_point_tag(w,np) = s ;
  }

  agg_wireframe_spline_number(w) ++ ;
  w->isp[nsp+1] = w->isp[nsp] + pps ;
  
  return 0 ;
}

gboolean agg_wireframe_splines_connected(agg_wireframe_t *w,
					 gint i0, gint i1, gint *p)

{
  gint i0p0, i0p1, i1p0, i1p1 ;

  agg_wireframe_spline_ends(w, i0, &i0p0, &i0p1) ;  
  agg_wireframe_spline_ends(w, i1, &i1p0, &i1p1) ;  

  if ( i0p0 == i1p0 ) { *p = i0p0 ; return TRUE ; }
  if ( i0p0 == i1p1 ) { *p = i0p0 ; return TRUE ; }
  if ( i0p1 == i1p0 ) { *p = i0p1 ; return TRUE ; }
  if ( i0p1 == i1p1 ) { *p = i0p1 ; return TRUE ; }
  
  return FALSE ;
}

gint agg_wireframe_element_add(agg_wireframe_t *w,
			       gint s0, gint s1, gint s2, gint s3)

{
  gint ne, *el ;

  ne = agg_wireframe_element_number(w) ;
  el = agg_wireframe_element(w, ne) ;

  el[0] = s0 ; el[1] = s1 ; el[2] = s2 ; el[3] = s3 ; el[4] = 0 ;

  /*check connections*/
  if ( !agg_wireframe_splines_connected(w, ABS(s0), ABS(s1), &ne) ) {
    g_error("%s: splines s0 (%d) and s1 (%d) not connected",
	    __FUNCTION__, s0, s1) ;
  }
  if ( !agg_wireframe_splines_connected(w, ABS(s1), ABS(s2), &ne) ) {
    g_error("%s: splines s1 (%d) and s2 (%d) not connected",
	    __FUNCTION__, s1, s2) ;
  }
  if ( s3 != 0 ) {
    if ( !agg_wireframe_splines_connected(w, ABS(s2), ABS(s3), &ne) ) {
      g_error("%s: splines s2 (%d) and s3 (%d) not connected",
	      __FUNCTION__, s2, s3) ;
    }
    if ( !agg_wireframe_splines_connected(w, ABS(s3), ABS(s0), &ne) ) {
      g_error("%s: splines s3 (%d) and s0 (%d) not connected",
	      __FUNCTION__, s3, s0) ;
    }
  } else {
    if ( !agg_wireframe_splines_connected(w, ABS(s2), ABS(s0), &ne) ) {
      g_error("%s: splines s2 (%d) and s0 (%d) not connected",
	      __FUNCTION__, s2, s0) ;
    }
  }
  
  agg_wireframe_element_number(w) ++ ;
  
  return 0 ;
}

gboolean agg_wireframe_spline_degenerate(agg_wireframe_t *w, gint s)

{
  gint p0, p1 ;
  gdouble *pt0, *pt1 ;

  agg_wireframe_spline_ends(w, s, &p0, &p1) ;  
  pt0 = agg_wireframe_point(w, p0) ;
  pt1 = agg_wireframe_point(w, p1) ;

  if ( pt0[0] == pt1[2] && pt0[2] == pt1[2] && pt0[2] == pt1[2] ) return TRUE ;
  
  return FALSE ;
}

static void add_spline(agg_wireframe_t *w, gint i0, gint i1)

{
  gint i, nsp, *sp ;

  nsp = agg_wireframe_spline_number(w) ;

  sp = agg_wireframe_spline(w, nsp) ;
  for ( i = 0 ; i <= i1-i0 ; i ++ ) {
    sp[i] = i0 + i ;
  }

  w->isp[nsp+1] = w->isp[nsp] + i1-i0+1 ;

  agg_wireframe_spline_number(w) ++ ;  

  return ;
}

gint agg_wireframe_intersection_add(agg_wireframe_t *w,
				    agg_intersection_t *inter,
				    gint nsp, gint pps,
				    agg_surface_workspace_t *ws)

{
  gint ni, i, np, np0 ;
  gdouble *p ;

  if ( agg_intersection_point_number(inter) == 0 ) return 0 ;

  /*find limits on intersection*/
  agg_intersection_bbox_set(inter) ;

  ni = agg_wireframe_intersection_number(w) ;
  agg_wireframe_intersection(w, ni) = inter ;

  g_assert(agg_intersection_point_number(inter) == nsp*(pps-1) + 1) ;

  np0 = agg_wireframe_point_number(w) ;
  agg_wireframe_intersection_points_start(w, ni) = np0 ;
    
  for ( i = 0 ; i < agg_intersection_point_number(inter) ; i ++ ) {
    np = agg_wireframe_point_number(w) ;
    p = agg_wireframe_point(w, np) ;
    memcpy(p, agg_intersection_point(inter,i), 3*sizeof(gdouble)) ;
    agg_wireframe_point_tag(w,np) = -1 - ni ;
    agg_wireframe_point_number(w) ++ ;    
  }
  agg_wireframe_intersection_points_end(w, ni) =
    agg_wireframe_point_number(w) ;

  for ( i = 0 ; i < nsp ; i ++ ) {
    add_spline(w, np0 + i*(pps-1), np0 + (i+1)*(pps-1)) ;
  }

  agg_wireframe_intersection_points_per_spline(w,ni) = pps ;
  agg_wireframe_intersection_spline_number(w,ni) = nsp ;
  
  agg_wireframe_intersection_number(w) ++ ;

  return 0 ;
}

static void surface_intersections(agg_wireframe_t *w, agg_surface_t *S,
				  gint *inter, gint *ninter)

{
  gint i ;

  *ninter = 0 ;
  for ( i = 0 ; i < agg_wireframe_intersection_number(w) ; i ++ ) {
    if ( agg_intersection_surface1(w->inter[i]) == S ||
	 agg_intersection_surface2(w->inter[i]) == S ) {
      inter[(*ninter)] = i ; (*ninter) ++ ;
    }
  }
  
  return ;
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

static void set_boundary(agg_wireframe_t *w, gint isurf,
			 gint *inter, gint ninter,
			 gdouble *st, gint *nst,
			 gint *segments, gint *nseg,
			 gint *segmentmarkers,
			 gint *idxlist, gint *nidx)
			 
{
  agg_surface_t *S ;
  agg_intersection_t *iS ;
  gdouble *sti, tmin, tmax, smin, smax ;
  gint i, pps, nst0 ;
  
  if ( ninter == 0 ) {
    segments[2*(*nseg)+0] = (*nst)+0 ; segments[2*(*nseg)+1] = (*nst)+1 ;
    (*nseg) ++ ;
    segments[2*(*nseg)+0] = (*nst)+1 ; segments[2*(*nseg)+1] = (*nst)+2 ;
    (*nseg) ++ ;
    segments[2*(*nseg)+0] = (*nst)+2 ; segments[2*(*nseg)+1] = (*nst)+3 ;
    (*nseg) ++ ;
    segments[2*(*nseg)+0] = (*nst)+3 ; segments[2*(*nseg)+1] = (*nst)+0 ;
    (*nseg) ++ ;
    
    st[2*(*nst)+0] = 0 ; st[2*(*nst)+1] = 0 ; (*nst) ++ ;
    st[2*(*nst)+0] = 1 ; st[2*(*nst)+1] = 0 ; (*nst) ++ ;
    st[2*(*nst)+0] = 1 ; st[2*(*nst)+1] = 1 ; (*nst) ++ ;
    st[2*(*nst)+0] = 0 ; st[2*(*nst)+1] = 1 ; (*nst) ++ ;

    return ;
  }

  /*look for correct list of boundary points in parametric space of surface*/
  S = agg_wireframe_surface(w, isurf) ;
  iS = agg_wireframe_intersection(w, inter[0]) ;

  if ( agg_intersection_surface1(iS) == S) {
    sti = &(iS->st[0]) ;
    smin = agg_intersection_bbox_s1_min(iS) ;
    smax = agg_intersection_bbox_s1_max(iS) ;
    tmin = agg_intersection_bbox_t1_min(iS) ;
    tmax = agg_intersection_bbox_t1_max(iS) ;
  } else {
    sti = &(iS->st[2]) ;
    smin = agg_intersection_bbox_s2_min(iS) ;
    smax = agg_intersection_bbox_s2_max(iS) ;
    tmin = agg_intersection_bbox_t2_min(iS) ;
    tmax = agg_intersection_bbox_t2_max(iS) ;
  }

  /*check limits to see if this is a hole or a boundary*/
  if ( smin > 0 && smax < 1 && tmin > 0 && tmax < 1 ) {
    set_boundary(w, isurf, inter, 0, st, nst, segments, nseg, segmentmarkers,
		 idxlist, nidx) ;
    return ;
  }
  
  /*copy spline end points into boundary*/
  pps = agg_wireframe_intersection_points_per_spline(w,inter[0]) ;
  nst0 = *nst ;
  for ( i = 0 ; i <= agg_wireframe_intersection_spline_number(w, inter[0]) ;
	i ++ ) {
    st[2*(*nst) + 0] = sti[i*(pps-1)*AGG_INTERSECTION_DATA_SIZE+0] ;
    st[2*(*nst) + 1] = sti[i*(pps-1)*AGG_INTERSECTION_DATA_SIZE+1] ;

    segments[2*(*nseg)+0] = (*nst)+0 ; segments[2*(*nseg)+1] = (*nst)+1 ; 
    
    (*nst) ++ ; (*nseg) ++ ;
    idxlist[(*nidx)] =
      agg_wireframe_intersection_points_start(w, inter[0]) + i*(pps-1) ;
    (*nidx) ++ ;
  }

  segments[2*(*nseg)+0] = (*nst)+0 ; segments[2*(*nseg)+1] = (*nst)+1 ;
  (*nseg) ++ ;
  st[2*(*nst) + 0] = 1 ; st[2*(*nst) + 1] = 1 ; (*nst) ++ ;
  segments[2*(*nseg)+0] = (*nst)+0 ; segments[2*(*nseg)+1] = (*nst)+1 ;
  (*nseg) ++ ;
  st[2*(*nst) + 0] = 1 ; st[2*(*nst) + 1] = 0 ; (*nst) ++ ;
  segments[2*(*nseg)+0] = (*nst)-1 ; segments[2*(*nseg)+1] = nst0 ;
  (*nseg) ++ ;
  
  return ;
}

static void add_hole(agg_wireframe_t *w, gint isurf,
		     gint inter, gint ninter,
		     gdouble *st, gint *nst,
		     gint *segments, gint *nseg,
		     gint *segmentmarkers, 
		     gdouble *hole, gint *nhole,
		     gint *idxlist, gint *nidx)
{
  agg_surface_t *S ;
  agg_intersection_t *iS ;
  gdouble *sti, tmin, tmax, smin, smax ;
  gint i, pps, p0 ;
  
  if ( ninter == 0 ) return ;

  /*look for correct list of boundary points in parametric space of surface*/
  S = agg_wireframe_surface(w, isurf) ;
  iS = agg_wireframe_intersection(w, inter) ;

  if ( agg_intersection_surface1(iS) == S) {
    sti = &(iS->st[0]) ;
    smin = agg_intersection_bbox_s1_min(iS) ;
    smax = agg_intersection_bbox_s1_max(iS) ;
    tmin = agg_intersection_bbox_t1_min(iS) ;
    tmax = agg_intersection_bbox_t1_max(iS) ;
  } else {
    sti = &(iS->st[2]) ;
    smin = agg_intersection_bbox_s2_min(iS) ;
    smax = agg_intersection_bbox_s2_max(iS) ;
    tmin = agg_intersection_bbox_t2_min(iS) ;
    tmax = agg_intersection_bbox_t2_max(iS) ;
  }

  /*check limits to see if this is a hole or a boundary*/
  if ( !(smin > 0 && smax < 1 && tmin > 0 && tmax < 1) ) return ;

  /*the intersection is a hole in this surface*/
  
  /*copy spline end points into boundary*/
  pps = agg_wireframe_intersection_points_per_spline(w,inter) ;
  p0 = (*nst) ;
  hole[2*(*nhole)+0] = hole[2*(*nhole)+1] = 0 ;
  for ( i = 0 ; i <= agg_wireframe_intersection_spline_number(w, inter) ;
	i ++ ) {
    st[2*(*nst) + 0] = sti[i*(pps-1)*AGG_INTERSECTION_DATA_SIZE+0] ;
    st[2*(*nst) + 1] = sti[i*(pps-1)*AGG_INTERSECTION_DATA_SIZE+1] ;

    hole[2*(*nhole)+0] += st[2*(*nst) + 0] ;
    hole[2*(*nhole)+1] += st[2*(*nst) + 1] ;

    (*nst) ++ ;

    idxlist[(*nidx)] =
      agg_wireframe_intersection_points_start(w, inter) + i*(pps-1) ;
    (*nidx) ++ ;    
  }

  hole[2*(*nhole)+0] /= i ;
  hole[2*(*nhole)+1] /= i ;
  
  (*nhole) ++ ;

  for ( i = p0 ; i < (*nst) ; i ++ ) {
    segments[2*(*nseg)+0] = i + 0 ; segments[2*(*nseg)+1] = i + 1 ;
    segmentmarkers[(*nseg)] = 3 ;
    (*nseg) ++ ;
  }

  segments[2*(*nseg)+0] = (*nst) - 1 ; 
  segments[2*(*nseg)+1] = p0 ;
  segmentmarkers[(*nseg)] = 3 ;
  (*nseg) ++ ;
  
  return ;
}


gint agg_wireframe_surface_add(agg_wireframe_t *w,
			       agg_surface_t *S, agg_patch_t *P,
			       gint nsec, gint nseg, gint pps,
			       agg_surface_workspace_t *ws)

{
  context *ctx ;
  triangleio in ;
  mesh *m ;
  behavior *b ;
  statistics s ;
  gint np, np0, i0, i1, i2, j0, j1, j2, isurf, inter[32], ninter, nbpts ;
  gint npre, idxpre[1024], *sp, *e, i, vertexnumber ;
  struct otri triangleloop, trisym;
  struct osub checkmark;
  vertex p1, p2, p3, vertexloop ;
  glong edgenumber, elementnumber ;
  triangle ptr; 
  subseg sptr;  

  isurf = agg_wireframe_surface_number(w) ;
  agg_wireframe_patch(w,isurf) = P ;
  agg_wireframe_surface(w,isurf) = S ;
  agg_wireframe_surface_number(w) ++ ;

  ninter = 0 ;
  surface_intersections(w, S, inter, &ninter) ;
  /* g_assert(ninter < 2) ; */

  /*maximum number of boundary points*/
  nbpts = 4 ;
  if ( ninter > 0 ) {
    nbpts += agg_wireframe_intersection_spline_number(w, inter[0]) + 4096 ;
  }

  nbpts = MAX(nbpts, 2048) ;
  
  ctx = triangle_context_create() ;
  triangle_context_options(ctx, "pzqa0.001");
  reset_triangleio(&in);

  in.pointlist = (gdouble *)g_malloc0(nbpts*2*sizeof(gdouble)) ;
  in.pointmarkerlist = (gint *)g_malloc0(nbpts*sizeof(gint)) ;
  in.segmentlist = (gint *)g_malloc0(nbpts*2*sizeof(gint)) ;
  in.segmentmarkerlist = (gint *)g_malloc0(nbpts*sizeof(gint)) ;
  in.holelist = (gdouble *)g_malloc0(16*2*sizeof(gdouble)) ;

  npre = 0 ;
  in.numberofpoints = in.numberofsegments = 0 ;
  in.numberofholes = 0 ;
  
  for ( i = 0 ; i < ninter ; i ++ ) {
    add_hole(w, isurf, inter[i], ninter,
	     in.pointlist, &(in.numberofpoints),
	     in.segmentlist, &(in.numberofsegments),
	     in.segmentmarkerlist,
	     in.holelist, &(in.numberofholes),
	     idxpre, &npre) ;
  }
 
  set_boundary(w, isurf, inter, ninter,
	       in.pointlist, &(in.numberofpoints), 
	       in.segmentlist, &(in.numberofsegments),
	       in.segmentmarkerlist,
	       idxpre, &npre) ;
  
  np0 = in.numberofpoints ;

  /* for ( i = 0 ; i < 20 ; i ++ ) { */
  /*   for ( j = 0 ; j < 20 ; j ++ ) { */
  /*     in.pointlist[(in.numberofpoints)*2+0] = 0.5 + 0.5*(gdouble)i/20 ; */
  /*     in.pointlist[(in.numberofpoints)*2+1] = (gdouble)j/20 ; */
  /*     in.numberofpoints ++ ; */
  /*   } */
  /* } */

 /* for ( i = 0 ; i < 19 ; i ++ ) { */
 /*   for ( j = 0 ; j < 19 ; j ++ ) { */
 /*     in.segmentlist[(in.numberofsegments)*2+0] = np0 + 20*i+j ; */
 /*     in.segmentlist[(in.numberofsegments)*2+1] = np0 + 20*i+j+1 ; */
 /*     in.numberofsegments ++ ; */
 /*     in.segmentlist[(in.numberofsegments)*2+0] = np0 + 20*i+j ; */
 /*     in.segmentlist[(in.numberofsegments)*2+1] = np0 + 20*(i+1)+j ; */
 /*     in.numberofsegments ++ ; */
 /*   } */
 /* } */
 
  triangle_mesh_create(ctx, &in) ;
  m = ctx->m ;
  b = ctx->b ;
  
  triangle_mesh_statistics(ctx, &s);

  traversalinit(&m->vertices);
  vertexnumber = b->firstnumber;
  vertexloop = vertextraverse(m);
  np0 = agg_wireframe_point_number(w) ;
  while (vertexloop != (vertex) NULL) {
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
      np = agg_wireframe_surface_point_add(w, isurf,
					   vertexloop[0], vertexloop[1], ws) ;
      agg_wireframe_point_tag(w,np) = isurf ;
      setvertexmark(vertexloop, vertexnumber);
      vertexnumber++;
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
	i0 = vertexmark(p1) ; i1 = vertexmark(p2) ;
	if ( i0 >= npre ) { j0 = i0 + np0 ; } else { j0 = idxpre[i0] ; }
	if ( i1 >= npre ) { j1 = i1 + np0 ; } else { j1 = idxpre[i1] ; }
	if ( i0 >= npre || i1 >= npre ) {
	  i0 += np0 ; i1 += np0 ;
	  agg_wireframe_spline_interp_points(w, isurf, i0, i1, pps, ws) ;
	  sp = agg_wireframe_spline(w,agg_wireframe_spline_number(w)-1) ;
	  sp[0] = j0 ; sp[pps-1] = j1 ;
	  edgenumber++;
	}
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

    if ( i0 >= npre ) { j0 = i0 + np0 ; } else { j0 = idxpre[i0] ; }
    if ( i1 >= npre ) { j1 = i1 + np0 ; } else { j1 = idxpre[i1] ; }
    if ( i2 >= npre ) { j2 = i2 + np0 ; } else { j2 = idxpre[i2] ; }
    
    i0 = agg_wireframe_spline_from_endpoints(w, j0, j1) ;
    if ( i0 != 0 ) {
      i1 = agg_wireframe_spline_from_endpoints(w, j1, j2) ;
      if ( i1 != 0 ) {
	i2 = agg_wireframe_spline_from_endpoints(w, j2, j0) ;
	if ( i2 != 0 ) {
	  e = agg_wireframe_element(w, agg_wireframe_element_number(w)) ;
	  e[0] = i0 ; e[1] = i1 ; e[2] = i2 ; e[3] = 0 ;
	  agg_wireframe_element_number(w) ++ ;
	}
      }
    }
    
    triangleloop.tri = triangletraverse(m);
    elementnumber++;
  }
  
  return 0 ;
}

gint agg_wireframe_spline_from_endpoints(agg_wireframe_t *w,
					 gint p0, gint p1)

{
  gint i, *sp, len ;

  for ( i = 1 ; i < agg_wireframe_spline_number(w) ; i ++ ) {
    sp = agg_wireframe_spline(w, i) ;
    len = agg_wireframe_spline_length(w, i) ;
    if ( sp[0] == p0 && sp[len-1] == p1 ) return  i ;
    if ( sp[0] == p1 && sp[len-1] == p0 ) return -i ;
  }
  
  return 0 ;
}

gint agg_wireframe_surface_point_add(agg_wireframe_t *w, gint surf,
				     gdouble s, gdouble t,
				     agg_surface_workspace_t *ws)

{
  gint np ;
  gdouble *p, u, v ;
  agg_surface_t *S ;
  agg_patch_t *P ;
  
  np = agg_wireframe_point_number(w) ;

  if ( np >= agg_wireframe_point_number_max(w) ) {
    g_error("%s: not enough space allocated (%d points)",
	    __FUNCTION__, agg_wireframe_point_number_max(w)) ;
  }
  
  p = agg_wireframe_point(w, np) ;
  P = agg_wireframe_patch(w, surf) ;
  g_assert((S = agg_wireframe_surface(w, surf)) != NULL) ;
  
  if ( P != NULL ) {
    agg_patch_map(P, s, t, &u, &v) ;    
  } else {
    u = s ; v = t ;
  }
  agg_surface_point_eval(S, u, v, p, ws) ;    
  p[3] = s ; p[4] = t ;

  agg_wireframe_point_number(w) ++ ;
  
  return np ;
}
