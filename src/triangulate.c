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

/* void traversalinit(struct memorypool *pool); */
/* void *traverse(struct memorypool *pool); */
/* vertex vertextraverse(mesh *m); */
/* triangle *triangletraverse(mesh *m); */
/* int plus1mod3[3] = {1, 2, 0}; */
/* int minus1mod3[3] = {2, 0, 1};
 */

extern void traversalinit(struct memorypool *pool);
extern void *traverse(struct memorypool *pool);
extern vertex vertextraverse(mesh *m);
extern triangle *triangletraverse(mesh *m);
extern int plus1mod3[] ;
extern int minus1mod3[] ;
#endif /*HAVE_TRIANGLE_API_H*/

static gdouble point_sample(agg_sampling_t sample,
			    gdouble tmin, gdouble tmax, gint nt, gint i)

{
  gdouble t, tc ;
  
  if ( sample == AGG_SAMPLING_LINEAR ) {
    t = tmin + (tmax - tmin)*i/nt ;
    return t ;
  }
  
  if ( sample == AGG_SAMPLING_COSINE ) {
    t = tmax - (tmax - tmin)*(1+cos(M_PI*i/nt))*0.5 ;
    /* t = tmin + 0.5*(tmax - tmin)*(1 - cos(M_PI*i/nt)) ; */
    /* t = cos(M_PI*((gdouble)i/nt-1)) ; */

    return t ;
  }

  if ( sample == AGG_SAMPLING_COSINE_DOUBLE ) {
    tc = (gdouble)i/nt ;
    if ( tc < 0.5 ) {
      t = -0.5*(1.0 + cos(2.0*M_PI*tc)) ;
    } else {
      t =  0.5*(1.0 + cos(2.0*M_PI*tc)) ;
    }

    t = tmin + (tmax - tmin)*(1+t)*0.5 ;

    return t ;
  }

  g_assert_not_reached() ;
  
  return 0.0 ;
}

static void tri_spline_surface_interp(agg_mesh_t *m, gint s,
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

static void triangulate_regular(agg_mesh_t *m,
				gint iS,
				agg_triangulation_settings_t *settings,
				agg_surface_workspace_t *w)

{
  gint i0, i1, np, pps, *sp ;
  gdouble u, v, s, t, smin, smax, tmin, tmax ;
  agg_surface_t *S ;
  gint ns, nt ;
  gint i, j, np0, nsp0, nspc, nspt, s0, s1, s2, s3 ;
  gboolean invert ;
  
  S = agg_mesh_surface(m,iS) ; 

  ns = agg_triangulation_section_number(settings) ;  
  nt = agg_triangulation_section_point_number(settings) ;
  pps = agg_triangulation_points_per_spline(settings) ;

  invert = agg_triangulation_invert(settings) ;
  
  tmin = -1 ; tmax = 1 ;
  smin = agg_surface_umin(agg_mesh_surface(m,iS)) ;
  smax = agg_surface_umax(agg_mesh_surface(m,iS)) ;

  np0 = agg_mesh_point_number(m) ;
  nsp0 = agg_mesh_spline_number(m) ;
  for ( i = 0 ; i <= ns ; i ++ ) {
    s = point_sample(agg_triangulation_sampling_s(settings), smin, smax,
		     ns, i) ;
    for ( j = 0 ; j <= nt ; j ++ ) {
      t = point_sample(agg_triangulation_sampling_t(settings), tmin, tmax,
		       nt, j) ;
      agg_mesh_surface_point_add(m, iS, s, t, w) ;
    }
  }

  /*splines joining element corners*/
  for ( i = 0 ; i < ns ; i ++ ) {
    for ( j = 0 ; j < nt ; j ++ ) {
      tri_spline_surface_interp(m, iS,
				np0 + (i+0)*(nt+1) + j + 0,
				np0 + (i+1)*(nt+1) + j + 0, pps, w) ;
      tri_spline_surface_interp(m, iS,
				np0 + (i+0)*(nt+1) + j + 0,
				np0 + (i+0)*(nt+1) + j + 1, pps, w) ;
    }
  }
  i = ns ;
  for ( j = 0 ; j < nt ; j ++ ) {
    tri_spline_surface_interp(m, iS,
			      np0 + (i+0)*(nt+1) + j + 0,
			      np0 + (i+0)*(nt+1) + j + 1, pps, w) ;
  }
  j = nt ;
  nspt = agg_mesh_spline_number(m) ;
  for ( i = 0 ; i < ns ; i ++ ) {
    tri_spline_surface_interp(m, iS,
			      np0 + (i+0)*(nt+1) + j + 0,
			      np0 + (i+1)*(nt+1) + j + 0, pps, w) ;
  }    

  for ( i = np0 ; i < agg_mesh_point_number(m) ; i ++ ) {
    agg_mesh_point_tag(m, i) = iS ;
  }

  /*elements*/
  for ( i = 0 ; i < ns - 1 ; i ++ ) {
    for ( j = 0 ; j < nt - 1 ; j ++ ) {
      s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
      s1 = nsp0 + 2*(i+1)*nt + 2*(j+0) + 1 ;
      s2 = nsp0 + 2*(i+0)*nt + 2*(j+1) + 0 ;
      s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;
      
      if ( invert ) {
	agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
      } else {
	agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
      }
    }
  }

  i = ns - 1 ;
  for ( j = 0 ; j < nt - 1 ; j ++ ) {
    s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nt +   (j+0) + 0 ;
    s2 = nsp0 + 2*(i+0)*nt + 2*(j+1) + 0 ;
    s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;
    if ( invert ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  j = nt - 1 ;
  for ( i = 0 ; i < ns - 1 ; i ++ ) {
    s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
    s1 = nsp0 + 2*(i+1)*nt + 2*(j+0) + 1 ;
    s2 = nsp0 + 2*ns*nt + j + 1 + i ;
    s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;
    
    if ( invert ) {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    } else {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    }
  }
  i = ns - 1 ; j = nt - 1 ;
  s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
  s1 = nsp0 + 2*(i+1)*nt + (j+0) + 0 ;
  s2 = nsp0 + 2*ns*nt + j + 1 + i ;
  s3 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 1 ;    
  if ( invert ) {
    agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
  } else {
    agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
  }

  if ( !agg_surface_close(S) ) return ;

  nspc = agg_mesh_spline_number(m) ;
  
  j = 0 ;
  for ( i = 0 ; i <= ns ; i ++ ) {
    tri_spline_surface_interp(m, iS,
			      np0 + (i+0)*(nt+1) + j + 0,
			      np0 + (i+0)*(nt+1) + j + nt, 2, w) ;
  }

  j = 0 ;
  for ( i = 0 ; i < ns ; i ++ ) {
    s0 = nsp0 + 2*(i+0)*nt + 2*(j+0) + 0 ;
    s1 = nspc + i + 1 ;
    s2 = nspt + i ;
    s3 = nspc + i + 0 ;
    if ( invert ) {
      agg_mesh_element_add(m, s0, s1, -s2, -s3) ;
    } else {
      agg_mesh_element_add(m, s3, s2, -s1, -s0) ;
    }
  }
  
  return ;
}
  
gint agg_mesh_surface_triangulate(agg_mesh_t *m, gint isurf,
				  agg_triangulation_settings_t *settings,
				  agg_surface_workspace_t *w)

{

  switch ( agg_triangulation_type(settings) ) {
  default:
    g_error("%s: unrecognised triangulation type in settings", __FUNCTION__) ;
    break ;
  case AGG_TRIANGULATION_UNDEFINED:
    g_error("%s: triangulation type not set", __FUNCTION__) ;
    break ;
  case AGG_TRIANGULATION_GRID_REGULAR:
    triangulate_regular(m, isurf, settings, w) ;
    break ;
  }

  return 0 ;
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

/** 
 * Add a point to a mesh
 * 
 * @param m mesh to have point added;
 * @param surf index in \a m of surface on which to evaluate point;
 * @param s coordinate for surface \a surf;
 * @param t coordinate for surface \a surf;
 * @param w workspace for surface evaluation.
 * 
 * @return index of newly added point.
 */

gint agg_mesh_surface_point_add(agg_mesh_t *m, gint surf,
				gdouble s, gdouble t,
				agg_surface_workspace_t *w)

{
  gint np ;
  gdouble *p ;
  agg_surface_t *S ;
  
  np = agg_mesh_point_number(m) ;

  if ( np >= agg_mesh_point_number_max(m) ) {
    g_error("%s: not enough space allocated (%d points)",
	    __FUNCTION__, agg_mesh_point_number_max(m)) ;
  }
  
  p = agg_mesh_point(m, np) ;
  g_assert((S = agg_mesh_surface(m, surf)) != NULL) ;

  agg_surface_point_eval(S, s, t, p, w) ;    
  p[3] = s ; p[4] = t ;
  agg_mesh_point_tag(m,np) = surf ;

  agg_mesh_point_number(m) ++ ;
  
  return np ;
}

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
			      char *len, gint offp, gint offsp, gint offs,
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

gboolean agg_mesh_spline_zero_length(agg_mesh_t *m, gint s)

{
  gint i0, i1 ;
  gdouble *x0, *x1 ;
  
  agg_mesh_spline_ends(m, s, &i0, &i1) ;

  x0 = agg_mesh_point(m, i0) ;
  x1 = agg_mesh_point(m, i1) ;

  if ( x0[0] == x1[0] && x0[1] == x1[1] && x0[2] == x1[2] ) return TRUE ;
  
  return FALSE ;
}
