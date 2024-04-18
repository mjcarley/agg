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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#include <triangle_api.h>
#pragma GCC diagnostic pop

void traversalinit(struct memorypool *pool);
void *traverse(struct memorypool *pool);
vertex vertextraverse(mesh *m);
triangle *triangletraverse(mesh *m);
int plus1mod3[3] = {1, 2, 0};
int minus1mod3[3] = {2, 0, 1};
#endif /*HAVE_TRIANGLE_API_H*/

static gint compare_double(gconstpointer a, gconstpointer b)

{
  gdouble x, y ;

  x = *((gdouble *)a) ; y = *((gdouble  *)b) ;

  if ( x < y ) return -1 ;
  if ( x > y ) return  1 ;
  
  return 0 ;
}

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
  
  mesh_spline_surface_blend_interp(m, b, p0, p1, pps, w) ;
  
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

static void insert_grid_segments(gdouble *st, gint *nst,
				 gint *segments, gint *nseg,
				 gint ns, gint nt)
{
  gint i, j, nst0 ;

  nst0 = *nst ;
  for ( i = 0 ; i <= ns ; i ++ ) {
    for ( j = 0 ; j <= nt ; j ++ ) {
      st[2*(*nst) + 0] = (gdouble)i/ns ;
      st[2*(*nst) + 1] = (gdouble)j/nt ;
      (*nst) ++ ;
    }
  }

  for ( i = 0 ; i <= ns ; i ++ ) {
    for ( j = 0 ; j < nt ; j ++ ) {      
      segments[2*(*nseg)+0] = nst0 + i*(nt+1) + j + 0 ; 
      segments[2*(*nseg)+1] = nst0 + i*(nt+1) + j + 1 ;
      (*nseg) ++ ;
    }
  }

  for ( j = 0 ; j <= nt ; j ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      segments[2*(*nseg)+0] = nst0 + (i+0)*(nt+1) + j ;
      segments[2*(*nseg)+1] = nst0 + (i+1)*(nt+1) + j ;
      (*nseg) ++ ;
    }
  }
  
  return ;
}

static gboolean edges_connected(gint *e, gint i, gint j,
				gint *p0, gint *p1)

{
  if ( e[2*i+0] == e[2*j+0] ) {
    *p0 = e[2*i+1] ; *p1 = e[2*j+1] ;
    
    return TRUE ;
  }

  if ( e[2*i+0] == e[2*j+1] ) {
    *p0 = e[2*i+1] ; *p1 = e[2*j+0] ;
    
    return TRUE ;
  }

  if ( e[2*i+1] == e[2*j+0] ) {
    *p0 = e[2*i+0] ; *p1 = e[2*j+1] ;
    
    return TRUE ;
  }

  if ( e[2*i+1] == e[2*j+1] ) {
    *p0 = e[2*i+0] ; *p1 = e[2*j+0] ;
    
    return TRUE ;
  }
  
  return FALSE ;
}

static void split_triangle(gint *e, gint *ne, gint ne0,
			   gint *t, gint *nt, gint t0)

{
  gint i0, i1, j0, j1, k0, k1, tri[12], ntri, p[6], np ;

  /*triangle split edges*/
  i0 = t[3*t0+0] ; i1 = i0 + ne0 ; 
  j0 = t[3*t0+1] ; j1 = j0 + ne0 ; 
  k0 = t[3*t0+2] ; k1 = k0 + ne0 ; 

  ntri = 0 ; np = 0 ;

  if ( edges_connected(e, i0, j0, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i0 ; tri[3*ntri+1] = j0 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, i0, j1, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i0 ; tri[3*ntri+1] = j1 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, i0, k0, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i0 ; tri[3*ntri+1] = k0 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, i0, k1, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i0 ; tri[3*ntri+1] = k1 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, i1, j0, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i1 ; tri[3*ntri+1] = j0 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, i1, j1, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i1 ; tri[3*ntri+1] = j1 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, i1, k0, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i1 ; tri[3*ntri+1] = k0 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, i1, k1, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = i1 ; tri[3*ntri+1] = k1 ; ntri ++ ; np ++ ;
  }

  if ( edges_connected(e, j0, k0, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = j0 ; tri[3*ntri+1] = k0 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, j0, k1, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = j0 ; tri[3*ntri+1] = k1 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, j1, k0, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = j1 ; tri[3*ntri+1] = k0 ; ntri ++ ; np ++ ;
  }
  if ( edges_connected(e, j1, k1, &(p[2*np+0]), &(p[2*np+1])) ) {
    tri[3*ntri+0] = j1 ; tri[3*ntri+1] = k1 ; ntri ++ ; np ++ ;
  }

  if ( ntri != 3 )
    g_error("%s: invalid triangle number %d", __FUNCTION__, t0) ;

  g_assert(np == ntri) ;

  /*three new edges from the open jaws*/
  e[2*((*ne)+0)+0] = p[0] ; e[2*((*ne)+0)+1] = p[1] ; 
  e[2*((*ne)+1)+0] = p[2] ; e[2*((*ne)+1)+1] = p[3] ; 
  e[2*((*ne)+2)+0] = p[4] ; e[2*((*ne)+2)+1] = p[5] ; 

  t[3*(*nt)+0] = tri[3*0+0] ; t[3*(*nt)+1] = tri[3*0+1] ;
  t[3*(*nt)+2] = (*ne)+0 ; (*nt) ++ ;
  t[3*(*nt)+0] = tri[3*1+0] ; t[3*(*nt)+1] = tri[3*1+1] ;
  t[3*(*nt)+2] = (*ne)+1 ; (*nt) ++ ;
  t[3*(*nt)+0] = tri[3*2+0] ; t[3*(*nt)+1] = tri[3*2+1] ;
  t[3*(*nt)+2] = (*ne)+2 ; (*nt) ++ ;

  t[3*t0+0] = (*ne) + 0 ; 
  t[3*t0+1] = (*ne) + 1 ; 
  t[3*t0+2] = (*ne) + 2 ; 
  
  (*ne) += 3 ;
  
  return ;
}

static void subdivide_loop(agg_patch_t *P,
			   gdouble *st, gint *nst,
			   gint *e, gint *ne,
			   gint *tri, gint *nt)

{
  gint i, j, k, ne0, nt0 ;

  ne0 = *ne ; nt0 = *nt ;
  /*split triangle edges with new points*/
  for ( i = 0 ; i < ne0 ; i ++ ) {
    j = e[2*i+0] ; k = e[2*i+1] ;

    agg_patch_edge_split(P, st[2*j+0], st[2*j+1], st[2*k+0], st[2*k+1],
			 0.5, &(st[2*(*nst)+0]), &(st[2*(*nst)+1])) ;
    
    e[2*i+1] = (*nst) ;
    e[2*(*ne)+0] = (*nst) ; e[2*(*ne)+1] = k ;
    g_assert((*ne) == i+ne0) ;
    (*nst) ++ ; (*ne) ++ ;
  }

  /*split triangles and add new internal edges*/
  for ( i = 0 ; i < nt0 ; i ++ ) {
    split_triangle(e, ne, ne0, tri, nt, i) ;
  }
  
  return ;
}

static void cartesian_to_unit_sphere(gdouble *x, gdouble *ph, gdouble *th)

{
  gdouble r ;

  r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) ;
  *ph = acos(x[2]/r) ;
  *th = atan2(x[1], x[0]) ;

  if ( *th < 0 ) (*th) += 2.0*M_PI ;

  return ;
}

/*
 * icosahedron data taken from:
 * 
 * https://danielsieger.com/blog/2021/01/03/generating-platonic-solids.html
 * 
 */

static void sphere_ico_base(gdouble *st, gint *ed, gint *tr,
			    gint offp, gint offe,
			    gint *np, gint *ne, gint *nt)

{
  gdouble x[3], phi, a, b, r ;
  gint i ;
  const gint edges[] = {
    0, 1,   0,  2,    0, 18,    0,  7,    0,  9,   0, 19,    1,  2,
    1, 3,   1,  7,    1,  8,    2,  3,    2,  5,   2,  9,    3, 16,
    3, 5,   3,  8,    3, 17,    4,  5,    4, 10,   4, 11,    4, 13,
    5, 9,   5, 11,    6,  9,    6, 10,    6, 11,   6, 12,    7,  8,
    7, 15,  8, 14,    9, 11,   10, 11,   10, 12,  10, 13,   12, 13,
    4, 16,  5, 16,   14, 17,    8, 17,    6, 18,   9, 18,   15, 19,
    7, 19
  } ;
  const gint tri[] = {
    41, 42, 28,     5,  3, 42,     4,  2, 40,    39, 23, 40,
    8,   9, 27,     0,  8,  3,     1,  6,  0,    12,  1,  4,
    21, 11, 12,    22, 21, 30,    23, 25, 30,    24, 31, 25,
    32, 24, 26,    34, 33, 32,     7, 15,  9,    10,  7,  6,
    11, 14, 10,    19, 17, 22,    18, 19, 31,    33, 20, 18,
    38, 37, 29,    15, 16, 38,    36, 13, 14,    17, 35, 36
  } ;
  phi = 0.5*(1.0 + sqrt(5.0)) ;
  a = 1.0 ; b = 1.0/phi ;
  r = sqrt(a*a + b*b) ;

  a /= r ; b /= r ;

  (*np) = 0 ;
  x[0] =  0 ; x[1] =  b ; x[2] = -a ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  b ; x[1] =  a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -b ; x[1] =  a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  0 ; x[1] =  b ; x[2] =  a ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  0 ; x[1] = -b ; x[2] =  a ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -a ; x[1] =  0 ; x[2] =  b ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  0 ; x[1] = -b ; x[2] = -a ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  a ; x[1] =  0 ; x[2] = -b ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  a ; x[1] =  0 ; x[2] =  b ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -a ; x[1] =  0 ; x[2] = -b ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  b ; x[1] = -a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -b ; x[1] = -a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;

  /*12 points in regular icosahedron at this stage*/  
  for ( i = 0 ; i < (*np) ; i ++ ) {
    st[2*i+0] = 1.0 - st[2*i+0]/M_PI ;
    st[2*i+1] = 1.0 - st[2*i+1]*0.5/M_PI ;
  }

  /*extra points to wrap around patch boundaries in parametric space*/
  st[2*(*np)+0] = st[2*7+0] ; st[2*(*np)+1] = 0 ; (*np) ++ ;
  st[2*(*np)+0] = st[2*8+0] ; st[2*(*np)+1] = 0 ; (*np) ++ ;
  st[2*(*np)+0] = st[2*4+0] ; st[2*(*np)+1] = st[2*4+1] + 1.0 ; (*np) ++ ;
  st[2*(*np)+0] = st[2*6+0] ; st[2*(*np)+1] = st[2*6+1] + 1.0 ; (*np) ++ ;

  st[2*(*np)+0] = 1.0 ; st[2*(*np)+1] = 0.5 ; (*np) ++ ;
  st[2*(*np)+0] = 1.0 ; st[2*(*np)+1] = 1.0 ; (*np) ++ ;
  st[2*(*np)+0] = 0.0 ; st[2*(*np)+1] = 0.5 ; (*np) ++ ;
  st[2*(*np)+0] = 0.0 ; st[2*(*np)+1] = 1.0 ; (*np) ++ ;
  
  *ne = 43 ;
  for ( i = 0 ; i < (*ne) ; i ++ ) {
    ed[2*i+0] = edges[2*i+0] + offp ;
    ed[2*i+1] = edges[2*i+1] + offp ;
  }

  *nt = 24 ;
  for ( i = 0 ; i < (*nt) ; i ++ ) {
    tr[3*i+0] = tri[3*i+0] + offe ;
    tr[3*i+1] = tri[3*i+1] + offe ;
    tr[3*i+2] = tri[3*i+2] + offe ;
  }
  
  return ;
}

static void insert_sphere_ico_segments(agg_patch_t *P,
				       gdouble *st, gint *nst,
				       gint *segments, gint *nseg, gint sub)

{
  gint i, nst0, nseg0, ne, np, nt, tri[8192*3] ;

  nst0 = *nst ; nseg0 = *nseg ;
  sphere_ico_base(&(st[2*nst0]), &(segments[2*nseg0]), tri, nst0, nseg0,
		  &np, &ne, &nt) ;
  (*nst) += np ; (*nseg) += ne ;

  for ( i = 0 ; i < sub ; i ++ ) {
    subdivide_loop(P, &(st[2*nst0]), &np, &(segments[2*nseg0]), &ne, tri, &nt) ;
  }

  (*nst) = nst0 + np ;
  (*nseg) = nseg0 + ne ;
    
  return ;
}

/*
 * Regular gridding of spherical surface:
 * 
 * https://danielsieger.com/blog/2021/03/27/generating-spheres.html
 * 
 */

static void insert_sphere_uv_segments(agg_patch_t *P,
				      gdouble *st, gint *nst,
				      gint *segments, gint *nseg,
				      gint ns, gint nt)

{
  gint i, j, nseg0 ;

  nseg0 = *nseg ;
  for ( i = 1 ; i < ns ; i ++ ) {
    for ( j = 0 ; j <= nt ; j ++ ) {
      st[2*(*nst)+0] = (gdouble)i/ns ;
      st[2*(*nst)+1] = (gdouble)j/nt ;
      (*nst) ++ ;
    }
  }

  for ( i = 1 ; i < ns-1 ; i ++ ) {
    for ( j = 0 ; j <= nt ; j ++ ) {
      segments[2*(*nseg)+0] = nseg0 + (i-1+0)*(nt+1) + j ;
      segments[2*(*nseg)+1] = nseg0 + (i-1+1)*(nt+1) + j ;
      (*nseg) ++ ;
    }
  }

  for ( j = 0 ; j < nt ; j ++ ) {
    for ( i = 1 ; i < ns ; i ++ ) {
      segments[2*(*nseg)+0] = nseg0 + (i-1)*(nt+1) + j + 0 ;
      segments[2*(*nseg)+1] = nseg0 + (i-1)*(nt+1) + j + 1 ;
      (*nseg) ++ ;
    }
  }

  i = 1 ;
  for ( j = 0 ; j < nt ; j ++ ) {
    st[2*(*nst)+0] = 0.0 ;
    st[2*(*nst)+1] = (gdouble)(j+0.5)/nt ;
    (*nst) ++ ;
    segments[2*(*nseg)+0] = nseg0 + (i-1)*(nt+1) + j + 0 ;
    segments[2*(*nseg)+1] = (*nst) - 1 ;
    (*nseg) ++ ;   
    segments[2*(*nseg)+0] = (*nst) - 1 ;
    segments[2*(*nseg)+1] = nseg0 + (i-1)*(nt+1) + j + 1 ;
    (*nseg) ++ ;   
  }

  i = ns - 1 ;
  for ( j = 0 ; j < nt ; j ++ ) {
    st[2*(*nst)+0] = 1.0 ;
    st[2*(*nst)+1] = (gdouble)(j+0.5)/nt ;
    (*nst) ++ ;
    segments[2*(*nseg)+0] = nseg0 + (i-1)*(nt+1) + j + 0 ;
    segments[2*(*nseg)+1] = (*nst) - 1 ;
    (*nseg) ++ ;   
    segments[2*(*nseg)+0] = (*nst) - 1 ;
    segments[2*(*nseg)+1] = nseg0 + (i-1)*(nt+1) + j + 1 ;
    (*nseg) ++ ;   
  }
  
  return ;
}

static void hemisphere_ico_base(gdouble *st, gint *ed, gint *tr,
				gint offp, gint offe,
				gint *np, gint *ne, gint *nt)

{
  gdouble x[3], phi, a, b, r ;
  gint i ;
  const gint edges[] = {
    0,  1,    0,  2,    0,  5,    0,  6,    1,  2,    1,  4,
    1,  7,    2,  4,    2,  6,    2, 13,    2, 14,    3,  4,
    3,  8,    3,  9,    3, 10,    3, 13,    4,  7,    4,  9,
    4, 13,    5,  6,    6, 12,    6, 14,    7,  9,    8,  9,
    8, 10,    8, 11,   10, 11,   12, 14
  } ;
  gint tri[] = {
     2,  3, 19,     1,  8,  3,     0,  4,  1,    5,  7,  4,
    16,  5,  6,    22, 17, 16,    13, 11, 17,   12, 13, 23,
    24, 14, 12,    26, 24, 25,    15, 18, 11,   18,  9,  7,
    10, 21,  8,    21, 27, 20
  } ;
  phi = 0.5*(1.0 + sqrt(5.0)) ;
  a = 1.0 ; b = 1.0/phi ;
  r = sqrt(a*a + b*b) ;

  a /= r ; b /= r ;
  
  (*np) = 0 ;
  x[0] =  b ; x[1] =  a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -b ; x[1] =  a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  0 ; x[1] =  b ; x[2] =  a ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  0 ; x[1] = -b ; x[2] =  a ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -a ; x[1] =  0 ; x[2] =  b ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  a ; x[1] =  0 ; x[2] = -b*0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  a ; x[1] =  0 ; x[2] =  b ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -a ; x[1] =  0 ; x[2] = -b*0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] =  b ; x[1] = -a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;
  x[0] = -b ; x[1] = -a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &(st[2*(*np)+0]), &(st[2*(*np)+1])) ; (*np) ++ ;

  for ( i = 0 ; i < (*np) ; i ++ ) {
    st[2*i+0] = 1.0 - st[2*i+0]/M_PI*2 ;
    st[2*i+1] = 1.0 - st[2*i+1]*0.5/M_PI ;
  }

  /*extra points to wrap around patch boundaries in parametric space*/  
  st[2*(*np)+0] = st[2*6+0] ; st[2*(*np)+1] = 0 ; (*np) ++ ;
  st[2*(*np)+0] = st[2*7+0] ; st[2*(*np)+1] = 0 ; (*np) ++ ;
  st[2*(*np)+0] = st[2*3+0] ; st[2*(*np)+1] = st[2*3+1] + 1.0 ; (*np) ++ ;
  
  st[2*(*np)+0] = 1.0 ; st[2*(*np)+1] = 0.5 ; (*np) ++ ;
  st[2*(*np)+0] = 1.0 ; st[2*(*np)+1] = 1.0 ;
  (*np) ++ ;

  *ne = 28 ;
  for ( i = 0 ; i < (*ne) ; i ++ ) {
    ed[2*i+0] = edges[2*i+0] + offp ;
    ed[2*i+1] = edges[2*i+1] + offp ;
  }

  *nt = 14 ;
  for ( i = 0 ; i < (*nt) ; i ++ ) {
    tr[3*i+0] = tri[3*i+0] + offe ;
    tr[3*i+1] = tri[3*i+1] + offe ;
    tr[3*i+2] = tri[3*i+2] + offe ;
  }
  
  return ;
}

static void insert_hemisphere_ico_segments(agg_patch_t *P,
					   gdouble *st, gint *nst,
					   gint *segments, gint *nseg, gint sub)

{
  gint i, nst0, nseg0, ne, np, nt, tri[8192*3] ;

  nst0 = *nst ; nseg0 = *nseg ;

  hemisphere_ico_base(&(st[2*nst0]), &(segments[2*nseg0]), tri, nst0, nseg0,
		      &np, &ne, &nt) ;
  (*nst) += np ; (*nseg) += ne ;

  for ( i = 0 ; i < sub ; i ++ ) {
    subdivide_loop(P, &(st[2*nst0]), &np, &(segments[2*nseg0]), &ne, tri, &nt) ;
  }

  (*nst) = nst0 + np ;
  (*nseg) = nseg0 + ne ;
  
  return ;
}

/*
 * Regular gridding of hemispherical surface:
 * 
 * https://danielsieger.com/blog/2021/03/27/generating-spheres.html
 * 
 */

static void insert_hemisphere_uv_segments(agg_patch_t *P,
					  gdouble *st, gint *nst,
					  gint *segments, gint *nseg,
					  gint ns, gint nt)

{
  gint i, j, nseg0 ;

  nseg0 = *nseg ;
  for ( i = 0 ; i < ns ; i ++ ) {
    for ( j = 0 ; j <= nt ; j ++ ) {
      st[2*(*nst)+0] = (gdouble)i/ns ;
      st[2*(*nst)+1] = (gdouble)j/nt ;
      (*nst) ++ ;
    }
  }

  for ( i = 0 ; i < ns-1 ; i ++ ) {
    for ( j = 0 ; j <= nt ; j ++ ) {
      segments[2*(*nseg)+0] = nseg0 + (i+0)*(nt+1) + j ;
      segments[2*(*nseg)+1] = nseg0 + (i+1)*(nt+1) + j ;
      (*nseg) ++ ;
    }
  }

  for ( j = 0 ; j < nt ; j ++ ) {
    for ( i = 0 ; i < ns ; i ++ ) {
      segments[2*(*nseg)+0] = nseg0 + (i-1)*(nt+1) + j + 0 ;
      segments[2*(*nseg)+1] = nseg0 + (i-1)*(nt+1) + j + 1 ;
      (*nseg) ++ ;
    }
  }

  i = ns - 1 ;
  for ( j = 0 ; j < nt ; j ++ ) {
    st[2*(*nst)+0] = 1.0 ;
    st[2*(*nst)+1] = (gdouble)(j+0.5)/nt ;
    (*nst) ++ ;
    segments[2*(*nseg)+0] = nseg0 + (i-1)*(nt+1) + j + 0 ;
    segments[2*(*nseg)+1] = (*nst) - 1 ;
    (*nseg) ++ ;   
    segments[2*(*nseg)+0] = (*nst) - 1 ;
    segments[2*(*nseg)+1] = nseg0 + (i-1)*(nt+1) + j + 1 ;
    (*nseg) ++ ;   
  }
  
  return ;
}

static gboolean in_patch_hole(agg_patch_t *P, gdouble *st, gdouble del)

{
  gint i ;
  agg_curve_t *c ;
  
  for ( i = 2 ; i < agg_patch_hole_number(P) ; i ++ ) {
    c = agg_patch_hole(P, i) ;
    if ( agg_curve_point_orientation(c, del, st[0], st[1]) == TRUE )
      return TRUE ;
  }
  
  return FALSE ;
}

static void trim_to_holes(gdouble *st, gint *nst,
			  gint *segments, gint *nseg,
			  gdouble del, agg_patch_t *P)

{
  gint i, idx0, idx1 ;

  for ( i = 0 ; i < (*nseg) ; i ++ ) {
    idx0 = segments[2*i+0] ; idx1 = segments[2*i+1] ; 
    if ( in_patch_hole(P, &(st[2*idx0]), del) ||
	 in_patch_hole(P, &(st[2*idx1]), del) ) {
      /*replace the segment with one from the end of the list*/
      segments[2*i+0] = segments[2*(*nseg)+0] ;
      segments[2*i+1] = segments[2*(*nseg)+1] ;
      (*nseg) -- ; i -- ;
    }
  }
  
  return ;
}

static gint hemisphere_boundary_points(gdouble *t, gint sub)

{
  gdouble x[3], phi, a, b, r, s ;
  gint i, j, n = 0 ;
  
  phi = 0.5*(1.0 + sqrt(5.0)) ;
  a = 1.0 ; b = 1.0/phi ;
  r = sqrt(a*a + b*b) ;

  a /= r ; b /= r ;
  
  n = 0 ;
  x[0] =  b ; x[1] =  a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &s, &(t[n])) ; n ++ ;
  x[0] = -b ; x[1] =  a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &s, &(t[n])) ; n ++ ;
  x[0] =  a ; x[1] =  0 ; x[2] = 0 ;
  cartesian_to_unit_sphere(x, &s, &(t[n])) ; n ++ ;
  x[0] = -a ; x[1] =  0 ; x[2] = 0 ;
  cartesian_to_unit_sphere(x, &s, &(t[n])) ; n ++ ;
  x[0] =  b ; x[1] = -a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &s, &(t[n])) ; n ++ ;
  x[0] = -b ; x[1] = -a ; x[2] =  0 ;
  cartesian_to_unit_sphere(x, &s, &(t[n])) ; n ++ ;

  for ( i = 0 ; i < n ; i ++ ) {
    t[i] = 1.0 - t[i]*0.5/M_PI ;
    /* if ( t[i] > 1 ) t[i] -= 1 ; */
    /* if ( t[i] < 0 ) t[i] += 1 ; */
  }
  t[n] = 0 ; n ++ ;
  
  qsort(t, n, sizeof(gdouble), compare_double) ;

  for ( i = 0 ; i < sub ; i ++ ) {
    for ( j = 0 ; j < n-1 ; j ++ ) {
      t[n+j] = 0.5*(t[j]+t[j+1]) ;
    }
    n += n-1 ;
    qsort(t, n, sizeof(gdouble), compare_double) ;
  }

  for ( i = 0 ; i < n ; i ++ ) {
    if ( t[i] > 1 ) t[i] -= 1 ;
    if ( t[i] < 0 ) t[i] += 1 ;
  }
  qsort(t, n, sizeof(gdouble), compare_double) ;
  
  return n ;
}

static gint rail_curve_points(agg_surface_blend_t *B, gint n, gdouble *t)

{
  gint j ;
  agg_surface_t *S ;
  
  S = agg_surface_blend_surface(B,0) ;
  if ( agg_surface_grid(S) == AGG_GRID_HEMISPHERE_ICO ) {
    j = agg_surface_blend_hole(B,0) ;
    if ( j == 0 ) {
      /*rail curve is the end of the patch*/
      return hemisphere_boundary_points(t, agg_surface_grid_subdivision(S)) ;
    }
  }
  if ( agg_surface_grid(S) == AGG_GRID_HEMISPHERE_UV ) {
    n = agg_surface_grid_spline_number(S) ;
  }
  S = agg_surface_blend_surface(B,1) ;
  if ( agg_surface_grid(S) == AGG_GRID_HEMISPHERE_ICO ) {
    j = agg_surface_blend_hole(B,1) ;
    if ( j == 0 ) {
      /*rail curve is the end of the patch*/
      return hemisphere_boundary_points(t, agg_surface_grid_subdivision(S)) ;
    }
  }
  if ( agg_surface_grid(S) == AGG_GRID_HEMISPHERE_UV ) {
    n = agg_surface_grid_spline_number(S) ;
  }

  /*default is even spacing on the curve*/
  for ( j = 0 ; j <= n ; j ++ ) {
    t[j] = (gdouble)j/n ;
  }

  return n+1 ;
}

/** 
 * Add a parametric surface to a mesh, working round holes in the
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
  gint np, np0, ns0, i0, i1, i2, j0, j1, j2, nbpts ;
  gint *sp, *e, i, j, vertexnumber ;
  struct otri triangleloop, trisym;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
  struct osub checkmark;
#pragma GCC diagnostic pop
  vertex p1, p2, p3, vertexloop ;
  glong edgenumber, elementnumber ;
  triangle ptr; 
  subseg sptr;
  gdouble s0, t0, t[1024] ;
  gint nseg = 64 ;
  agg_curve_t *c ;
  agg_patch_t *P ;
  agg_surface_t *S ;
  agg_surface_blend_t *B ;
  
  P = agg_mesh_patch(msh, isurf) ;
  S = agg_mesh_surface(msh, isurf) ;
  
  /*maximum number of boundary points*/
  nbpts = 1024 ;
  for ( i = 0 ; i < agg_mesh_intersection_number(msh) ; i ++ ) {
    nbpts += agg_mesh_intersection_points_end(msh, i) -
      agg_mesh_intersection_points_start(msh, i) ;
  }
  nbpts = MAX(nbpts, 16384) ;
  
  ctx = triangle_context_create() ;
  triangle_context_options(ctx, args) ;
  reset_triangleio(&in);

  in.pointlist = (gdouble *)g_malloc0(nbpts*2*sizeof(gdouble)) ;
  in.pointmarkerlist = (gint *)g_malloc0(nbpts*sizeof(gint)) ;
  in.segmentlist = (gint *)g_malloc0(nbpts*2*sizeof(gint)) ;
  in.segmentmarkerlist = (gint *)g_malloc0(nbpts*sizeof(gint)) ;
  in.holelist = (gdouble *)g_malloc0(16*2*sizeof(gdouble)) ;

  in.numberofpoints = in.numberofsegments = 0 ;
  in.numberofholes = 0 ;

  if ( agg_surface_grid(S) == AGG_GRID_REGULAR ) {
    insert_grid_segments(in.pointlist, &(in.numberofpoints),
			 in.segmentlist, &(in.numberofsegments),
			 agg_surface_grid_section_number(S),
			 agg_surface_grid_spline_number(S)) ;
  }
  
  if ( agg_surface_grid(S) == AGG_GRID_SPHERE_ICO ) {
    insert_sphere_ico_segments(P, in.pointlist, &(in.numberofpoints),
			       in.segmentlist, &(in.numberofsegments),
			       agg_surface_grid_subdivision(S)) ;
  }

  if ( agg_surface_grid(S) == AGG_GRID_SPHERE_UV ) {
    insert_sphere_uv_segments(P, in.pointlist, &(in.numberofpoints),
			      in.segmentlist, &(in.numberofsegments),
			      agg_surface_grid_section_number(S),
			      agg_surface_grid_spline_number(S)) ;
  }

  if ( agg_surface_grid(S) == AGG_GRID_HEMISPHERE_ICO ) {
    insert_hemisphere_ico_segments(P, in.pointlist, &(in.numberofpoints),
				   in.segmentlist, &(in.numberofsegments),
				   agg_surface_grid_subdivision(S)) ;
  }
  
  if ( agg_surface_grid(S) == AGG_GRID_HEMISPHERE_UV ) {
    insert_hemisphere_uv_segments(P, in.pointlist, &(in.numberofpoints),
				  in.segmentlist, &(in.numberofsegments),
				  agg_surface_grid_section_number(S),
				  agg_surface_grid_spline_number(S)) ;
  }
  
  trim_to_holes(in.pointlist, &(in.numberofpoints),
		in.segmentlist, &(in.numberofsegments),
		0.01, P) ;
  
  for ( i = 0 ; i < agg_patch_hole_number(P)-2 ; i ++ ) {
    np0 = in.numberofpoints ;
    ns0 = in.numberofsegments ;

    c = agg_patch_hole(P,i+2) ;
    B = agg_patch_blend(P,i+2) ;
    nseg = rail_curve_points(B, nseg, t) ;
    for ( j = 0 ; j < nseg-1 ; j ++ ) {
      agg_curve_eval(c, t[j], &s0, &t0) ;
      
      in.pointlist[2*np0+2*j+0] = s0 ;
      in.pointlist[2*np0+2*j+1] = t0 ;
      in.segmentlist[2*ns0+2*j+0] = np0+j+0 ; 
      in.segmentlist[2*ns0+2*j+1] = np0+j+1 ; 
    }
    in.segmentlist[2*ns0+2*(nseg-2)+1] = np0 ;    
    in.holelist[2*i+0] = c->data[0] ;
    in.holelist[2*i+1] = c->data[1] ; 
    in.numberofholes ++ ;
    in.numberofpoints += nseg ;
    in.numberofsegments += nseg ;
  }
  
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
  while (vertexloop != (vertex) NULL) {
    if (!b->jettison || (vertextype(vertexloop) != UNDEADVERTEX)) {
      s0 = vertexloop[0] ; t0 = vertexloop[1] ;
      np = agg_mesh_surface_point_add(msh, isurf, s0, t0, w) ;
      agg_mesh_point_tag(msh,np) = isurf ;
      setvertexmark(vertexloop, vertexnumber) ;
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
  gint i, j, np, np0, tag, nsp0, s0, s1, s2, s3, nsp, nt ;
  gdouble s, t[2048], *p ;
  agg_surface_blend_t *B ;

  g_assert(iB >= 0 && iB < agg_mesh_surface_blend_number(m)) ;
  
  B = agg_mesh_surface_blend(m,iB) ;
  nsp = 16 ;
  agg_surface_blend_invert(B) = FALSE ;

  nt = nsp ;
  nt = rail_curve_points(B, nsp, t) ;
  nt -- ;
  
  tag = -1-iB ;
  np0 = agg_mesh_point_number(m) ;
  nsp0 = agg_mesh_spline_number(m) ;
  for ( i = 0 ; i <= nsec ; i ++ ) {
    for ( j = 0 ; j <= nt ; j ++ ) {
      s = (gdouble)i/nsec ;
      np = agg_mesh_point_number(m) ;
      p  = agg_mesh_point(m,np) ;
      agg_surface_blend_evaluate(B, s, t[j], p, w) ;
      agg_mesh_point_s(m, np) = s ;
      agg_mesh_point_t(m, np) = t[j] ;
      agg_mesh_point_number(m) ++ ;      
    }
  }

  /*splines joining element corners*/
  for ( i = 0 ; i < nsec ; i ++ ) {
    for ( j = 0 ; j < nt ; j ++ ) {
      agg_mesh_spline_interp_points(m, tag,
				    np0 + (i+0)*(nt+1) + j + 0,
				    np0 + (i+1)*(nt+1) + j + 0, pps, w) ;
      agg_mesh_spline_interp_points(m, tag,
				    np0 + (i+0)*(nt+1) + j + 0,
				    np0 + (i+0)*(nt+1) + j + 1, pps, w) ;
    }
  }
  i = nsec ;
  for ( j = 0 ; j < nt ; j ++ ) {
    agg_mesh_spline_interp_points(m, tag,
				  np0 + (i+0)*(nt+1) + j + 0,
				  np0 + (i+0)*(nt+1) + j + 1, pps, w) ;
  }
  j = nt ;
  for ( i = 0 ; i < nsec ; i ++ ) {
    agg_mesh_spline_interp_points(m, tag,
				  np0 + (i+0)*(nt+1) + j + 0,
				  np0 + (i+1)*(nt+1) + j + 0, pps, w) ;
  }    
  
  for ( i = np0 ; i < agg_mesh_point_number(m) ; i ++ ) {
    agg_mesh_point_tag(m, i) = tag ;
  }

  /*elements*/
  for ( i = 0 ; i < nsec - 1 ; i ++ ) {
    for ( j = 0 ; j < nt - 1 ; j ++ ) {
      s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
      s1 = nsp0 + 2*(i+1)*nt + 2*(j+0) + 1 ;
      s2 = nsp0 + 2*(i+0)*nt + 2*(j+1) + 0 ;
      s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;
      
      if ( agg_surface_blend_invert(B) ) {
	agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
      } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
      }
    }
  }

  i = nsec - 1 ;
  for ( j = 0 ; j < nt - 1 ; j ++ ) {
    s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nt +   (j+0) + 0 ;
    s2 = nsp0 + 2*(i+0)*nt + 2*(j+1) + 0 ;
    s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;
    if ( agg_surface_blend_invert(B) ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  j = nt - 1 ;
  for ( i = 0 ; i < nsec - 1 ; i ++ ) {
    s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nt + 2*(j+0) + 1 ;
    s2 = nsp0 + 2*nsec*nt + j + 1 + i ;
    s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;
    
    if ( agg_surface_blend_invert(B) ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  i = nsec - 1 ; j = nt - 1 ;
  s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
  s1 = nsp0 + 2*(i+1)*nt + (j+0) + 0 ;
  s2 = nsp0 + 2*nsec*nt + j + 1 + i ;
  s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;    
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
  
  np0 = agg_mesh_point_number(m) ;
  nsp0 = agg_mesh_spline_number(m) ;

  smin = tmin = 0.0 ; smax = tmax = 1.0 ;
  
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
  gint i, j ;
  gdouble area ;
  gchar args[64] ;

  agg_mesh_surface_number(m) = 0 ;
  for ( i = 0 ; i < agg_body_surface_number(b); i ++ ) {
    agg_mesh_surface(m,i) = agg_body_surface(b,i) ;
    agg_mesh_patch(m,i) = agg_body_patch(b,i) ;
    agg_mesh_surface_number(m) ++ ;
  }

  fprintf(stderr, "%s: added %d surfaces\n", __FUNCTION__,
	  agg_mesh_surface_number(m)) ;
  
  for ( i = 0 ; i < agg_body_surface_number(b); i ++ ) {
    fprintf(stderr, "%s: adding surface %d\n", __FUNCTION__, i) ;
    for ( j = i+1 ; j < agg_body_surface_number(b); j ++ ) {
      if ( agg_surface_patch_trim(agg_body_surface(b,i), agg_body_patch(b,i),
				  0.05,
				  agg_body_surface(b,j), agg_body_patch(b,j),
				  0.05, &(m->B[m->nb]), w) ) {
	m->nb ++ ;
      }
    }
  }
  
  for ( i = 0 ; i < agg_mesh_surface_number(m) ; i ++ ) {
    g_assert(agg_surface_grid(agg_body_surface(b,i)) != AGG_GRID_UNDEFINED) ;
    area = agg_surface_grid_element_area(agg_body_surface(b,i)) ;
    area = 0.01 ;
    sprintf(args, "pza%lg", area) ;
    fprintf(stderr, "%s: triangulating surface %d\n", __FUNCTION__, i) ;    
    agg_mesh_surface_add_triangle(m, i, args, pps, w) ;
  }

  for ( i = 0 ; i < agg_mesh_surface_blend_number(m) ; i ++ ) {
    fprintf(stderr, "%s: building blend %d\n", __FUNCTION__, i) ;
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
