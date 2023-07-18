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
  agg_mesh_t *w ;

  w = (agg_mesh_t *)g_malloc0(sizeof(agg_mesh_t)) ;

  w->p   = (gdouble *)g_malloc0(AGG_MESH_POINT_SIZE*npmax*
				sizeof(gdouble)) ;
  w->sp  = (gint *)g_malloc0(nspmax*8*sizeof(gint)) ;
  w->isp = (gint *)g_malloc0(nspmax*sizeof(gint)) ;
  w->e   = (gint *)g_malloc0(AGG_MESH_ELEMENT_SIZE*nemax*sizeof(gint)) ;
  w->ptags = (gint *)g_malloc0(npmax*sizeof(gint)) ;

  w->np = 0 ; w->npmax = npmax ;
  w->nsp = 1 ; w->nspmax = nspmax ;
  w->isp[1] = 0 ;
  w->ne = 0 ; w->nemax = nemax ;

  agg_mesh_intersection_number(w) = 0 ;  
  agg_mesh_surface_number(w) = 0 ;
  
  return w ;
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
  gint i, *sp, nsp, np ;
  gdouble s1, t1, s0, t0 ;

  g_assert(s < agg_mesh_surface_number(m)) ;
  g_assert(agg_mesh_patch(m,s) != NULL) ;
  g_assert(agg_mesh_surface(m,s) != NULL) ;

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
  
  return 0 ;
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

static void add_spline(agg_mesh_t *m, gint i0, gint i1)

{
  gint i, nsp, *sp ;

  nsp = agg_mesh_spline_number(m) ;

  sp = agg_mesh_spline(m, nsp) ;
  for ( i = 0 ; i <= i1-i0 ; i ++ ) {
    sp[i] = i0 + i ;
  }

  m->isp[nsp+1] = m->isp[nsp] + i1-i0+1 ;

  agg_mesh_spline_number(m) ++ ;  

  return ;
}

/** 
 * Add an intersection to a mesh. This should be done before adding
 * surfaces, so that they can be cut if necessary. 
 * 
 * @param m an ::agg_mesh_t;
 * @param inter a surface-surface intersection;
 * @param nsp number of splines in \a inter;
 * @param pps number of points on each spline in \a inter;
 * 
 * @return 0 on success.
 */

gint agg_mesh_intersection_add(agg_mesh_t *m,
			       agg_intersection_t *inter,
			       gint nsp, gint pps)

{
  gint ni, i, np, np0 ;
  gdouble *p ;

  if ( agg_intersection_point_number(inter) == 0 ) return 0 ;

  /*find limits on intersection*/
  agg_intersection_bbox_set(inter) ;

  ni = agg_mesh_intersection_number(m) ;
  agg_mesh_intersection(m, ni) = inter ;
  
  agg_mesh_intersection_number(m) ++ ;

  return 0 ;
  
  g_assert(agg_intersection_point_number(inter) == nsp*(pps-1) + 1) ;

  np0 = agg_mesh_point_number(m) ;
  agg_mesh_intersection_points_start(m, ni) = np0 ;
    
  for ( i = 0 ; i < agg_intersection_point_number(inter) ; i ++ ) {
    np = agg_mesh_point_number(m) ;
    p = agg_mesh_point(m, np) ;
    memcpy(p, agg_intersection_point(inter,i), 3*sizeof(gdouble)) ;
    agg_mesh_point_tag(m,np) = -1 - ni ;
    agg_mesh_point_number(m) ++ ;    
  }
  agg_mesh_intersection_points_end(m, ni) =
    agg_mesh_point_number(m) ;

  for ( i = 0 ; i < nsp ; i ++ ) {
    add_spline(m, np0 + i*(pps-1), np0 + (i+1)*(pps-1)) ;
  }

  agg_mesh_intersection_points_per_spline(m,ni) = pps ;
  agg_mesh_intersection_spline_number(m,ni) = nsp ;
  
  agg_mesh_intersection_number(m) ++ ;

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

#define CONTACT_DATA_SIZE 5

static void edge_contacts(agg_mesh_t *m, gint isurf, gint inter,
			  gint *contacts, gint *nc,
			  gint *holes, gint *nh)

/*
 * if surface isurf on mesh m is intersected on curve inter, add to
 * contacts array:
 * 
 * [(inter) (index of first contact point) (index of second contact point)
 *  (0 for hole, 1 for boundary) (curve number on intersection, 0 or 1)]
 */
  
{
  gint i, j, c ;
  gdouble s, t ;
  agg_surface_t *S ;
  agg_intersection_t *iS ;
  
  g_assert((*nc) < 8) ;
  
  S = agg_mesh_surface(m, isurf) ;
  iS = agg_mesh_intersection(m, inter) ;
  
  if ( agg_intersection_surface1(iS) != S &&
       agg_intersection_surface2(iS) != S ) return ;

  /*which curve are we using on the intersection (first or second surface)*/
  c = ( agg_intersection_surface1(iS) == S ? 0 : 1 ) ;
  /*assume a hole to start with*/
  contacts[(*nc)*CONTACT_DATA_SIZE+3] = 0 ;
  contacts[(*nc)*CONTACT_DATA_SIZE+4] = c ;

  /*look for points lying on edges of mesh*/
  for ( i = j = 0 ; i < agg_intersection_point_number(iS) ; i ++ ) {
    s = agg_intersection_point_s(iS, i, c) ;
    t = agg_intersection_point_t(iS, i, c) ;
    if ( s == 0 || s == 1 || t == 0 || t == 1 ) {
      contacts[(*nc)*CONTACT_DATA_SIZE+  0] = inter ;
      contacts[(*nc)*CONTACT_DATA_SIZE+1+j] = i ;
      j ++ ;
    }
  }

  g_assert(j < 3) ;

  if ( j == 2 ) {
    /*two contacts, boundary curve*/
    contacts[(*nc)*CONTACT_DATA_SIZE+3] = 1 ;
    (*nc) ++ ;
  } else {
    /*no contacts, hole*/
    holes[(*nh)*2+0] = inter ;
    holes[(*nh)*2+1] = c ;
    (*nh) ++ ;
  }

  /*this is a slit in the patch and we have no way to deal with it*/
  g_assert(j != 1 ) ;
  
  return ;
}

static gint insert_edge_point(gdouble *epoints, gint ne,
			      gdouble s, gdouble t,
			      gint inter, gint idx, gint c)
{
  epoints[ne*CONTACT_DATA_SIZE+0] = s ;
  epoints[ne*CONTACT_DATA_SIZE+1] = t ;
  epoints[ne*CONTACT_DATA_SIZE+2] = inter ;
  epoints[ne*CONTACT_DATA_SIZE+3] = idx ;
  epoints[ne*CONTACT_DATA_SIZE+4] = c ;
  
  return ne + 1 ;
}

static void edge_points(agg_mesh_t *m, gint isurf,
			gint *contacts, gint nc,
			gdouble es, gdouble et, gdouble e,
			gdouble ws, gdouble wt,
			gboolean *cut,
			gdouble *epoints, gint *ne)

/*
 * on exit, epoints contains:
 * 
 * [s t (intersection index) (index of intersection point) 0]
 */
  
{
  gint i, c, idx, ne0, inter ;
  agg_intersection_t *iS ;
  gdouble s, t ;
  
  ne0 = (*ne) ; *cut = FALSE ;
  for ( i = 0 ; i < nc ; i ++ ) {
    if ( contacts[i*CONTACT_DATA_SIZE+3] == 1 ) {
      /*this is a boundary and not a hole*/
      inter = contacts[i*CONTACT_DATA_SIZE+0] ;
      iS = agg_mesh_intersection(m, inter) ;
      c = contacts[i*CONTACT_DATA_SIZE+4] ;
      /*check contact points to see if they are on this edge*/
      idx = contacts[i*CONTACT_DATA_SIZE+1] ;
      s = agg_intersection_point_s(iS, idx, c) ;
      t = agg_intersection_point_t(iS, idx, c) ;
      if ( es*s + et*t == e ) {
	(*ne) = insert_edge_point(epoints, *ne, s, t, inter, idx, c) ;
	*cut = TRUE ;
      }

      idx = contacts[i*CONTACT_DATA_SIZE+2] ;
      s = agg_intersection_point_s(iS, idx, c) ;
      t = agg_intersection_point_t(iS, idx, c) ;
      if ( es*s + et*t == e ) {
	(*ne) = insert_edge_point(epoints, *ne, s, t, inter, idx, c) ;
	*cut = TRUE ;
      }      
    }
  }

  /*only handle single intersection surfaces for now*/
  g_assert((*ne)-ne0 < 2) ;

  /*will need to do a sort on newly added points here when multiple
    intersections are handled*/
  
  return ;
}

static gdouble boundary_length(gdouble s, gdouble t)

{
  if ( t == 0 ) return s ;
  if ( s == 1 ) return 1.0 + t ; 
  if ( t == 1 ) return 2.0 + (1.0-s) ; 
  if ( s == 0 ) return 3.0 + (1.0-t) ;

  g_assert_not_reached() ;
  
  return 5 ;
}
  
static gint compare_boundary(gconstpointer p1, gconstpointer p2)

{
  const gdouble *st1 = p1, *st2 = p2 ;
  gdouble len1, len2 ;

  len1 = boundary_length(st1[0], st1[1]) ;
  len2 = boundary_length(st2[0], st2[1]) ;

  if ( len1 < len2 ) return -1 ;
  if ( len1 > len2 ) return  1 ;
  
  return 0 ;
}

static void set_boundary(agg_mesh_t *m, gint isurf,
			 gdouble *st, gint *nst,
			 gint *segments, gint *nseg,
			 gint *segmentmarkers,
			 gdouble *holes, gint *nholes,
			 /* gint *idxlist, gint *nidx, */
			 gint ppe, gdouble del)

/*
 * m:        mesh
 * isurf:    index of surface in mesh surface list
 * st:       list of points in patch coordinate system (s,t)
 * nst:      number of points
 * segments: segment list, two ints per segment (indices of points)
 * nseg:     number of segments
 * idxlist:  indices of mesh points lying on intersection curve
 * nidx:     number of points in idxlist
 * ppe:      points per patch edge
 * del:      offset from patch ends (to leave clear space for singular points)
 */
  
{
  agg_intersection_t *iS ;
  gdouble epoints[CONTACT_DATA_SIZE*8], shole, thole ;
  gint i, j, j0, j1, inter, c, nst0, contacts[CONTACT_DATA_SIZE*8],
    nc, ne ;
  gint iholes[2*8], nh ;
  gboolean edges_cut[4] = {FALSE} ;
  
  nst0 = (*nst) ;
  
  /*find contact edges and sense of intersections with patch edges*/
  for ( i = nc = nh = 0 ; i < agg_mesh_intersection_number(m) ; i ++ ) {
    edge_contacts(m, isurf, i, contacts, &nc, iholes, &nh) ;
  }
  
  /*single intersections only for now*/
  g_assert(nc == 1 || nc == 0) ;
  /*generate ordered list around the boundary of edge contact points*/
  for ( i = ne = 0 ; i < nc ; i ++ ) {
    edge_points(m, isurf, contacts, nc, 0, 1, 0,  1,  0,
		&(edges_cut[0]), epoints, &ne) ;
    edge_points(m, isurf, contacts, nc, 1, 0, 1,  0,  1,
		&(edges_cut[1]), epoints, &ne) ;
    edge_points(m, isurf, contacts, nc, 0, 1, 1, -1,  0,
		&(edges_cut[2]), epoints, &ne) ;
    edge_points(m, isurf, contacts, nc, 1, 0, 0,  0, -1,
		&(edges_cut[3]), epoints, &ne) ;
  }

  /*epoints now contains an ordered list of intersection points round
    the boundary: insert the corners where required*/
  if ( edges_cut[0] ) {
    if ( edges_cut[1] ) {
      ne = insert_edge_point(epoints, ne, 0, 0, -1, -1, -1) ;
      ne = insert_edge_point(epoints, ne, 1, 1, -1, -1, -1) ;
      ne = insert_edge_point(epoints, ne, 0, 1, -1, -1, -1) ;      
    } else {
      if ( edges_cut[2] ) {
	ne = insert_edge_point(epoints, ne, 1, 0, -1, -1, -1) ;
	ne = insert_edge_point(epoints, ne, 1, 1, -1, -1, -1) ;
      } else {
	if ( edges_cut[3] ) {
	  ne = insert_edge_point(epoints, ne, 1, 0, -1, -1, -1) ;
	  ne = insert_edge_point(epoints, ne, 1, 1, -1, -1, -1) ;
	  ne = insert_edge_point(epoints, ne, 0, 1, -1, -1, -1) ;
	} else {
	  g_assert_not_reached() ;
	}
      }      
    }
  } else {
    if ( edges_cut[1] ) {
      g_assert(!edges_cut[2]) ;
      g_assert(!edges_cut[3]) ;
    } else {
      if ( edges_cut[2] ) {
	g_assert(!edges_cut[3]) ;
      } else {
	/*no cut edges, put all the corners in*/
	ne = insert_edge_point(epoints, ne,   del, 0, -1, -1, -1) ;
	ne = insert_edge_point(epoints, ne, 1-del, 0, -1, -1, -1) ;
	ne = insert_edge_point(epoints, ne, 1-del, 1, -1, -1, -1) ;
	ne = insert_edge_point(epoints, ne,   del, 1, -1, -1, -1) ;
      }
    }
  }
  /*sort the boundary points*/
  qsort(epoints, ne, CONTACT_DATA_SIZE*sizeof(gdouble), compare_boundary) ;

  /*insert segments and points from epoints*/
  for ( i = 0 ; i < ne ; i ++ ) {
    if ( epoints[(i+0)*CONTACT_DATA_SIZE+2] != -1 &&
	 epoints[(i+1)*CONTACT_DATA_SIZE+2] != -1 ) {
      /*insert the intersection points*/
      /*identify the intersection and which points to use*/
      inter = epoints[i*CONTACT_DATA_SIZE+2] ;
      iS = agg_mesh_intersection(m, inter) ;
      /* pps = agg_mesh_intersection_points_per_spline(m,inter) ; */
      c     = epoints[i*CONTACT_DATA_SIZE+4] ;
      g_assert(iS->S[c] == m->S[isurf]) ;
      j0 = epoints[(i+0)*CONTACT_DATA_SIZE+3] ;
      j1 = epoints[(i+1)*CONTACT_DATA_SIZE+3] ;
      if ( j0 < j1 ) {
	for ( j = j0 ; j <= j1 ; j ++ )
	  add_point(st, nst,
		    agg_intersection_point_s(iS, j, c),
		    agg_intersection_point_t(iS, j, c)) ;
	(*nst) -- ;
      } else {
	for ( j = j0 ; j >= j1 ; j -- )
	  add_point(st, nst,
		    agg_intersection_point_s(iS, j, c),
		    agg_intersection_point_t(iS, j, c)) ;
	(*nst) -- ;
      }
    } else {
      /*corner point*/
      add_point(st, nst,
		epoints[i*CONTACT_DATA_SIZE+0],
		epoints[i*CONTACT_DATA_SIZE+1]) ;
    }
  }

  for ( i = nst0 ; i < (*nst)-1 ; i ++ ) {
    segments[2*(*nseg)+0] = i   ; 
    segments[2*(*nseg)+1] = i+1 ;
    (*nseg) ++ ;
  }

  segments[2*(*nseg)+0] = (*nst)-1   ; 
  segments[2*(*nseg)+1] = nst0 ;
  (*nseg) ++ ;

  /*add holes*/
  for ( i = 0 ; i < nh ; i ++ ) {
    inter = iholes[2*i+0] ; c = iholes[2*i+1] ;
    iS = agg_mesh_intersection(m, inter) ;
    nst0 = (*nst) ;
    shole = thole = 0.0 ;
    for ( j = 0 ; j < agg_intersection_point_number(iS)-1 ; j ++ ) {
      add_point(st, nst,
		agg_intersection_point_s(iS, j, c),
		agg_intersection_point_t(iS, j, c)) ;
      shole += agg_intersection_point_s(iS, j, c) ;
      thole += agg_intersection_point_t(iS, j, c) ;
    }
    shole /= (gdouble)agg_intersection_point_number(iS) ;
    thole /= (gdouble)agg_intersection_point_number(iS) ;

    for ( j = 0 ; j < agg_intersection_point_number(iS)-1 ; j ++ ) {
      segments[2*(*nseg)+0] = nst0 + j + 0 ; 
      segments[2*(*nseg)+1] = nst0 + j + 1 ; 
      (*nseg) ++ ;      
    }    
    segments[2*(*nseg)+0] = nst0 + agg_intersection_point_number(iS) - 2 ;
    segments[2*(*nseg)+1] = nst0 ; 
    (*nseg) ++ ;

    holes[(*nholes)*2+0] = shole ;
    holes[(*nholes)*2+1] = thole ;
    (*nholes) ++ ;
  }
  
  return ;
}

static gboolean boundary_point(gint *b, gint nb, gint p)

{
  gint i ;

  for ( i = 0 ; i < nb ; i ++ ) {
    if ( b[i] == p ) return TRUE ;
  }
  
  return FALSE ;
}

/** 
 * Add a parametric surface to a mesh, cutting on intersections if
 * necessary. This requires that any intersections be added in
 * advance, using ::agg_mesh_intersection_add. 
 * 
 * @param msh a mesh to have a surface added;
 * @param S the surface to add;
 * @param P the patch mapping \a S;
 * @param nsec number of sections to compute;
 * @param nseg number of segments on each section; 
 * @param pps number of points on each segment;
 * @param args argument string to pass to triangle;
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */

gint agg_mesh_surface_add(agg_mesh_t *msh,
			  agg_surface_t *S, agg_patch_t *P,
			  gint nsec, gint nseg, gint pps,
			  gchar *args,
			  agg_surface_workspace_t *w)

{
  context *ctx ;
  triangleio in ;
  /*it has to be called this because of how the macros in triangle.h
    are defined*/
  mesh *m ;
  behavior *b ;
  statistics s ;
  gint np, np0, i0, i1, i2, j0, j1, j2, isurf, nbpts ;
  gint *sp, *e, i, vertexnumber ;
  gint boundary[8192], nbound, tag ;
  struct otri triangleloop, trisym;
  struct osub checkmark;
  vertex p1, p2, p3, vertexloop ;
  glong edgenumber, elementnumber ;
  triangle ptr; 
  subseg sptr;  
  gint ncap = 4 ;
  gdouble del = 0.0625 ;
  gdouble s0, t0, s1, t1, s2, t2 ;
  del = 0 ;
  
  isurf = agg_mesh_surface_number(msh) ;
  agg_mesh_patch(msh,isurf) = P ;
  agg_mesh_surface(msh,isurf) = S ;
  agg_mesh_surface_number(msh) ++ ;

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
	       ncap, del) ;
  
  np0 = in.numberofpoints ;

  triangle_mesh_create(ctx, &in) ;
  m = ctx->m ;
  b = ctx->b ;
  
  triangle_mesh_statistics(ctx, &s);

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
	i0 = vertexmark(p1) ; i1 = vertexmark(p2) ;
	j0 = i0 + np0 ; j1 = i1 + np0 ; 
	/* i0 += np0 ; i1 += np0 ; */
	agg_mesh_spline_interp_points(msh, isurf, j0, j1, pps, w) ;
	sp = agg_mesh_spline(msh,agg_mesh_spline_number(msh)-1) ;
	sp[0] = j0 ; sp[pps-1] = j1 ;
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
 * Find a spline from its endpoints
 * 
 * @param m a mesh containing a set of splines;
 * @param p0 index of a point of \a m;
 * @param p1 index of another point of \a m.
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
i * @param t coordinate on patch for surface \a surf;
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
 * Mesh a ::agg_body_t containing a collection of
 * (possibly-intersecting) surfaces, using Christian Woltering's
 * library version of Jonathan Shewchuk's Triangle code:
 * 
 * https://github.com/wo80/Triangle/
 * https://www.cs.cmu.edu/~quake/triangle.html
 * 
 * @param m an ::agg_mesh_t to contain mesh on output;
 * @param b body to be meshed;
 * @param nsec number of sections to generate on each surface of \a b;
 * @param nsp number of splines on each section;
 * @param pps number of points per spline;
 * @param args argument string to pass to Triangle (see Triangle link
 * for details);
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */

gint agg_mesh_body_triangle(agg_mesh_t *m, agg_body_t *b,
			    gint nsec, gint nsp, gint pps,
			    gchar *args,
			    agg_surface_workspace_t *w)

{
  gint i, j ;
  agg_intersection_t *inter, *resample ;

  inter    = agg_intersection_new(8192) ;
  resample = agg_intersection_new(8192) ;
  for ( i = 0 ; i < agg_body_surface_number(b); i ++ ) {
    for ( j = i+1 ; j < agg_body_surface_number(b); j ++ ) {
      agg_surface_patch_intersection(inter,
				     agg_body_surface(b,i),
				     agg_body_patch(b,i),
				     agg_body_surface(b,j),
				     agg_body_patch(b,j), w) ;
      if ( agg_intersection_point_number(inter) != 0 ) {
	agg_intersection_resample(inter, nsp, pps, resample, w) ;
	agg_intersection_bbox_set(resample) ;
	agg_mesh_intersection_add(m, resample, nsp, pps) ;
	inter    = agg_intersection_new(8192) ;
	resample = agg_intersection_new(8192) ;
      }
    }
  }
  
  for ( i = 0 ; i < agg_body_surface_number(b); i ++ ) {
    agg_mesh_surface_add(m, agg_body_surface(b,i), agg_body_patch(b,i),
			 nsec, nsp, pps, args, w) ;
  }

  return 0 ;
}

gint agg_mesh_body_regular(agg_mesh_t *m, agg_body_t *b,
			    gint nsec, gint nsp, gint pps,
			    agg_surface_workspace_t *w)

{
  gint i, j, iS, np0, nsp0, s0, s1, s2, s3 ;
  gdouble s, t ;
  
  if ( agg_mesh_intersection_number(m) > 0 ) {
    g_error("%s: cannot generate regular grids on intersected surfaces",
	    __FUNCTION__) ;
  }

  for ( iS = 0 ; iS < agg_mesh_surface_number(m) ; iS ++ ) {
    np0 = agg_mesh_point_number(m) ;
    nsp0 = agg_mesh_spline_number(m) ;
    for ( i = 0 ; i <= nsec ; i ++ ) {
      for ( j = 0 ; j <= nsp ; j ++ ) {
	s = (gdouble)i/nsec ; t = (gdouble)j/nsp ;
	agg_mesh_surface_point_add(m, iS, s, t, w) ;
      }
    }

    /*splines joining element corners*/
    for ( i = 0 ; i < nsec ; i ++ ) {
      for ( j = 0 ; j < nsp ; j ++ ) {
	agg_mesh_spline_interp_points(m, iS,
				      np0 + (i+0)*(nsp+1) + j + 0,
				      np0 + (i+1)*(nsp+1) + j + 0, pps, w) ;
	/* fprintf(stderr, "%d %d\n", */
	/* 	agg_mesh_spline_number(m)-1, (nsp0+2*(i+0)*nsp+2*(j+0)+0)) ; */
	agg_mesh_spline_interp_points(m, iS,
				      np0 + (i+0)*(nsp+1) + j + 0,
				      np0 + (i+0)*(nsp+1) + j + 1, pps, w) ;
	/* fprintf(stderr, "%d %d\n", */
	/* 	agg_mesh_spline_number(m)-1, nsp0+2*(i+0)*nsp+2*(j+0)+1) ; */
      }
    }
    i = nsec ;
    for ( j = 0 ; j < nsp ; j ++ ) {
      agg_mesh_spline_interp_points(m, iS,
				    np0 + (i+0)*(nsp+1) + j + 0,
				    np0 + (i+0)*(nsp+1) + j + 1, pps, w) ;
      /* fprintf(stderr, "%d %d\n", */
      /* 	      agg_mesh_spline_number(m)-1,  */
      /* 		      (nsp0+2*(i+0)*nsp+(j+0)+0)) ; */
    }
    j = nsp ;
    for ( i = 0 ; i < nsec ; i ++ ) {
      agg_mesh_spline_interp_points(m, iS,
				    np0 + (i+0)*(nsp+1) + j + 0,
				    np0 + (i+1)*(nsp+1) + j + 0, pps, w) ;
      /* fprintf(stderr, "%d %d\n", */
      /* 	      agg_mesh_spline_number(m)-1, */
      /* 	      nsp0 + 2*nsec*nsp + j + i) ; */
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
      /* agg_mesh_element_add(m, s0, s1, -s2, -s3) ; */
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
          
      /* agg_mesh_element_add(m, s0, s1, -s2, -s3) ; */
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
    /* agg_mesh_element_add(m, s0, s1, -s2, -s3) ; */
    if ( agg_patch_invert(agg_mesh_patch(m, iS)) ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
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
