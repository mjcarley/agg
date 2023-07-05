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

agg_wireframe_t *agg_wireframe_new(gint npmax, gint nspmax, gint nemax)

{
  agg_wireframe_t *w ;

  w = (agg_wireframe_t *)g_malloc0(sizeof(agg_wireframe_t)) ;

  w->p   = (gdouble *)g_malloc0(AGG_WIREFRAME_POINT_SIZE*npmax*
				sizeof(gdouble)) ;
  w->sp  = (gint *)g_malloc0(nspmax*8*sizeof(gint)) ;
  w->isp = (gint *)g_malloc0(nspmax*sizeof(gint)) ;
  w->e   = (gint *)g_malloc0(4*nemax*sizeof(gint)) ;
  w->ptags = (gint *)g_malloc0(npmax*sizeof(gint)) ;

  w->np = 0 ; w->npmax = npmax ;
  w->nsp = 0 ; w->nspmax = nspmax ;
  w->ne = 0 ; w->nemax = nemax ;
  
  return w ;
}

static void interp_spline_points(agg_patch_t *P, agg_surface_t *S,
				 agg_surface_workspace_t *w,
				 gdouble *p, gint *np, gint *sp, gint pps)

{
  gint i, i0, i1 ;
  gdouble s, t, s1, t1, s0, t0, u, v ;

  i0 = sp[0] ; i1 = sp[pps-1] ;
  s0 = p[i0*AGG_WIREFRAME_POINT_SIZE+3] ;
  t0 = p[i0*AGG_WIREFRAME_POINT_SIZE+4] ;
  s1 = p[i1*AGG_WIREFRAME_POINT_SIZE+3] ;
  t1 = p[i1*AGG_WIREFRAME_POINT_SIZE+4] ;
  for ( i = 1 ; i < pps-1 ; i ++ ) {
    s = s0 + (s1 - s0)*i/(pps-1) ;
    t = t0 + (t1 - t0)*i/(pps-1) ;
    agg_patch_map(P, s, t, &u, &v) ;
    agg_surface_point_eval(S, u, v, &(p[(*np)*5]), w) ;
    p[(*np)*AGG_WIREFRAME_POINT_SIZE+3] = s ;
    p[(*np)*AGG_WIREFRAME_POINT_SIZE+4] = t ;
    sp[i] = (*np) ;
    (*np) ++ ;
  }
  
  return ;
}

static gboolean segments_intersect(gdouble *x1, gdouble *x2,
				   gdouble *y1, gdouble *y2,
				   gdouble *s, gdouble *t)

{
  gdouble A[4], b[2], det ;

  A[0] = x2[0] - x1[0] ; A[1] = -(y2[0] - y1[0]) ; 
  A[2] = x2[1] - x1[1] ; A[3] = -(y2[1] - y1[1]) ; 
  agg_invert2x2(A, A, &det) ;

  if ( fabs(det) < 1e-12 ) return FALSE ;
  b[0] = y1[0] - x1[0] ; 
  b[1] = y1[1] - x1[1] ; 

  *s = A[0]*b[0] + A[1]*b[1] ; 
  if ( *s < 0 || *s > 1 ) return FALSE ;
  *t = A[2]*b[0] + A[3]*b[1] ; 
  if ( *t < 0 || *t > 1 ) return FALSE ;
  
  return TRUE ;
}

static gboolean add_intersection(gdouble *sti0, gint i0,
				 gdouble *sti1, gint i1,
				 gdouble *stj0, gint j0,
				 gdouble *stj1, gint j1,
				 gdouble *p, gint *np,
				 gint *ie, gint *ni,
				 agg_patch_t *P, agg_surface_t *S,
				 agg_surface_workspace_t *w)

/*
 * stXX are points in the (s,t) parameter space
 * 
 * check for intersection of segment stiX with segment stjX and add to
 * point and intersection list if there is one
 */
  
{
  gdouble d, a, b, u, v ;

  /*wrapping tolerance check*/
  /* wtol = 0.1 ; */
  
  d = (stj1[0] - stj0[0])*(stj1[0] - stj0[0]) +
      (stj1[1] - stj0[1])*(stj1[1] - stj0[1]) ;
  
  if ( d < 3e-1 && segments_intersect(sti0, sti1, stj0, stj1, &a, &b) ) {
    /*add a new point at the intersection*/
    u = stj0[0] + b*(stj1[0] - stj0[0]) ;
    v = stj0[1] + b*(stj1[1] - stj0[1]) ;
    p[AGG_WIREFRAME_POINT_SIZE*(*np)+3] = u ;
    p[AGG_WIREFRAME_POINT_SIZE*(*np)+4] = v ;
    agg_patch_st_correct(P, &(p[AGG_WIREFRAME_POINT_SIZE*(*np)+3])) ;
    agg_patch_map(P,
		  p[AGG_WIREFRAME_POINT_SIZE*(*np)+3],
		  p[AGG_WIREFRAME_POINT_SIZE*(*np)+4],
		  &u, &v) ;
    agg_surface_point_eval(S, u, v, &(p[AGG_WIREFRAME_POINT_SIZE*(*np)]), w) ;
    
    ie[5*(*ni)+0] = i0 ; 
    ie[5*(*ni)+1] = i1 ; 
    ie[5*(*ni)+2] = (*np) ;
    ie[5*(*ni)+3] = j0 ;
    ie[5*(*ni)+4] = j1 ;
    
    (*np) ++ ; (*ni) ++ ;
    return TRUE ;
  }

  return FALSE ;
}

static void intersection_point(agg_surface_t *S, agg_patch_t *P,
			       gdouble *st, gint nst,
			       gdouble *p, gint *np, gint nsec, gint nseg,
			       gint i0, gint i1, gint *ie, gint *ni,
			       agg_surface_workspace_t *w)

{
  gdouble sti0[2], sti1[2], stj0[2], stj1[2] ;
  gint j ;
  
  sti0[0] = p[AGG_WIREFRAME_POINT_SIZE*i0+3] ;
  sti0[1] = p[AGG_WIREFRAME_POINT_SIZE*i0+4] ;
  sti1[0] = p[AGG_WIREFRAME_POINT_SIZE*i1+3] ;
  sti1[1] = p[AGG_WIREFRAME_POINT_SIZE*i1+4] ;

  for ( j = 0 ; j < nst-1 ; j ++ ) {
    stj0[0] = st[AGG_INTERSECTION_DATA_SIZE*(j+0)+0] ;
    stj0[1] = st[AGG_INTERSECTION_DATA_SIZE*(j+0)+1] ;
    stj1[0] = st[AGG_INTERSECTION_DATA_SIZE*(j+1)+0] ;
    stj1[1] = st[AGG_INTERSECTION_DATA_SIZE*(j+1)+1] ;
    if ( add_intersection(sti0, i0, sti1, i1, stj0, j, stj1, j+1, p, np, ie, ni,
			  P, S, w) )
      break ;
  }

  j = nst - 1 ;
  stj0[0] = st[AGG_INTERSECTION_DATA_SIZE*(j+0)+0] ;
  stj0[1] = st[AGG_INTERSECTION_DATA_SIZE*(j+0)+1] ;
  stj1[0] = st[AGG_INTERSECTION_DATA_SIZE*(  0)+0] ;
  stj1[1] = st[AGG_INTERSECTION_DATA_SIZE*(  0)+1] ;
  add_intersection(sti0, i0, sti1, i1, stj0, j, stj1, 0, p, np, ie, ni,
		   P, S, w) ;  

  return ;
}


static gboolean edge_intersection(gint e0, gint e1,
				  gint *ie, gint ne, gint *iint)

{
  gint i ;

  for ( i = 0 ; i < ne ; i ++ ) {
    if ( (ie[5*i+0] == e0 && ie[5*i+1] == e1) ||
	 (ie[5*i+1] == e0 && ie[5*i+0] == e1) ) {
      *iint = ie[5*i+2] ;
      return TRUE ;
    }
  }
  
  return FALSE ;
}

static gboolean element_intersected(agg_wireframe_t *w, gint i,
				    gint *ie, gint ne, gint ic[])
				    /* gint *i0, gint *i1) */

{
  gint ni, e0, e1, *e, j ;

  e = &(w->e[4*i]) ;

  for ( j = ni = 0 ; j < 4 ; j ++ ) {
    agg_wireframe_spline_ends(w, e[j], &e0, &e1) ;
    if ( edge_intersection(e0, e1, ie, ne, &(ic[2*ni+1])) ) {
      ic[2*ni+0] = j ;
      ni ++ ;
    }
  }

  /* agg_wireframe_spline_ends(w, e[0], &e0, &e1) ; */
  /* if ( edge_intersection(e0, e1, ie, ne, &(ic[ni])) ) ni ++ ; */
  /* agg_wireframe_spline_ends(w, e[1], &e0, &e1) ; */
  /* if ( edge_intersection(e0, e1, ie, ne, &(ic[ni])) ) ni ++ ; */
  /* agg_wireframe_spline_ends(w, e[2], &e0, &e1) ; */
  /* if ( edge_intersection(e0, e1, ie, ne, &(ic[ni])) ) ni ++ ; */
  /* agg_wireframe_spline_ends(w, e[3], &e0, &e1) ; */
  /* if ( edge_intersection(e0, e1, ie, ne, &(ic[ni])) ) ni ++ ; */

  if ( ni < 2 ) return FALSE ;

  return TRUE ;
}

static void cut_element(agg_wireframe_t *w, gint i, gint *ne,
			gint *ie, gint ni,
			agg_surface_t *S, agg_patch_t *P,
			agg_surface_workspace_t *ws,
			gint pps)

{
  gint *el, nsp, *sp, ic[8] ;

  if ( !element_intersected(w, i, ie, ni, ic) ) return ;

  el = &(w->e[4*i]) ;
  /*add the new spline*/
  nsp = agg_wireframe_spline_number(w) ;
  sp = agg_wireframe_spline(w, nsp) ;
  sp[0] = ic[2*0+1] ; sp[pps-1] = ic[2*1+1] ;
  w->isp[nsp+1] = w->isp[nsp] + pps ;
  interp_spline_points(P, S, ws, w->p, &(w->np), sp, pps) ;
  agg_wireframe_spline_number(w) ++ ;

  /*rotate the element indices */
  
  return ;
}

gint agg_wireframe_surface(agg_wireframe_t *w,
			   agg_surface_t *S, agg_patch_t *P,
			   agg_intersection_t *inter,
			   gint nsec, gint nseg, gint pps,
			   agg_surface_workspace_t *ws)

{
  gdouble u, v, s, t, *p, *st ;
  gint i, j, np, nsp, ne, i0, i1, ie[8192], ni, *sp, nst ;

  /*points on grid nodes*/
  for ( i = np = 0 ; i <= nsec ; i ++ ) {
    s = (gdouble)i/nsec ;
    for ( j = 0 ; j <= nseg ; j ++ ) {
      t = (gdouble)j/nseg ;
      agg_patch_map(P, s, t, &u, &v) ;
      p = agg_wireframe_point(w,np) ;
      agg_surface_point_eval(S, u, v, p, ws) ;
      p[3] = s ; p[4] = t ;

      g_assert(np == i*(nseg+1) + j) ;
      
      np ++ ;
    }
  }

  /*make the splines*/
  p = agg_wireframe_point(w,0) ;  
  sp = w->sp ;
  w->isp[0] = 0 ;
  for ( i = nsp = 0 ; i <= nsec ; i ++ ) {
    for ( j = 0 ; j < nseg ; j ++ ) {
      i0 = (i+0)*(nseg+1) + j + 0 ; 
      i1 = (i+0)*(nseg+1) + j + 1 ;
      sp[nsp*pps+    0] = i0 ;
      sp[nsp*pps+pps-1] = i1 ;

      interp_spline_points(P, S, ws, p, &np, &(sp[nsp*pps]), pps) ;

      g_assert(nsp == i*(nseg) + j) ;

      w->isp[nsp+1] = w->isp[nsp] + pps ;    
      nsp ++ ;
    }
  }
  g_assert(nsp == (nsec+1)*nseg) ;
  for ( j = 0 ; j <= nseg ; j ++ ) {
    for ( i = 0 ; i < nsec ; i ++ ) {
      i0 = (i+0)*(nseg+1) + j + 0 ; 
      i1 = (i+1)*(nseg+1) + j + 0 ;
      sp[nsp*pps+    0] = i0 ;
      sp[nsp*pps+pps-1] = i1 ;

      interp_spline_points(P, S, ws, p, &np, &(sp[nsp*pps]), pps) ;
      
      g_assert(nsp == (nsec+1)*(nseg) + j*nsec + i) ;

      w->isp[nsp+1] = w->isp[nsp] + pps ;    
      nsp ++ ;
    }
  }
  agg_wireframe_point_number(w) = np ;
  agg_wireframe_spline_number(w) = nsp ;
  
  /*make the elements (ordered lists of splines)*/
  for ( i = ne = 0 ; i < nsec ; i ++ ) {
    for ( j = 0 ; j < nseg ; j ++ ) {
      w->e[4*ne + 0] = (i+0)*nseg + j ;
      w->e[4*ne + 1] = (j+1)*nsec + i + (nsec+1)*nseg ;
      w->e[4*ne + 2] = (i+1)*nseg + j ;
      w->e[4*ne + 3] = (j+0)*nsec + i + (nsec+1)*nseg ;
      ne ++ ;
    }
  }

  agg_wireframe_element_number(w) = ne ;

  if ( inter == NULL ) return 0 ;

  nst = agg_intersection_point_number(inter) ;
  if ( agg_intersection_surface1(inter) == S ) {
    g_assert(agg_intersection_patch1(inter) == P) ;
    st = &(inter->st[0]) ;
  } else {
    g_assert(agg_intersection_surface2(inter) == S) ;
    g_assert(agg_intersection_patch2(inter) == P) ;
    st = &(inter->st[2]) ;
  }
  
  /*intersection of element edges with surface intersection curve*/
  p = agg_wireframe_point(w,0) ;
  np = agg_wireframe_point_number(w) ;
  for ( i = ni = 0 ; i < nsec ; i ++ ) {
    for ( j = 0 ; j < nseg ; j ++ ) {
      i0 = (i+0)*(nseg+1) + j + 0 ; i1 = (i+0)*(nseg+1) + j + 1 ;
      intersection_point(S, P, st, nst,
			 p, &np, nsec, nseg, i0, i1, ie, &ni, ws) ;
      i0 = (i+0)*(nseg+1) + j + 0 ; i1 = (i+1)*(nseg+1) + j + 0 ;
      intersection_point(S, P, st, nst,
			 p, &np, nsec, nseg, i0, i1, ie, &ni, ws) ;
    }
    j = nseg ;
    i0 = (i+0)*(nseg+1) + j + 0 ; i1 = (i+1)*(nseg+1) + j + 0 ;
    intersection_point(S, P, st, nst,
		       p, &np, nsec, nseg, i0, i1, ie, &ni, ws) ;
  }

  i = nsec ;
  for ( j = 0 ; j < nseg ; j ++ ) {
    i0 = (i+0)*(nseg+1) + j + 0 ; i1 = (i+0)*(nseg+1) + j + 1 ;
    intersection_point(S, P, st, nst,
		       p, &np, nsec, nseg, i0, i1, ie, &ni, ws) ;
  }
  
  agg_wireframe_point_number(w) = np ;

  /*split any elements cut by the intersection curve*/
  ne = agg_wireframe_element_number(w) ;
  for ( i = 0 ; i < agg_wireframe_element_number(w) ; i ++ ) {
    cut_element(w, i, &ne, ie, ni, S, P, ws, pps) ;
  }  
  
  return 0 ;
}

static void write_element_gmsh(FILE *f, gint *e, gint i, gint offsp, gint offs)

{
  if ( e[3] < 0 ) {
    /*three-sided element*/
    fprintf(f, "Curve Loop(%d) = {%d, %d, -%d} ;\n", offs + i,
	    offsp + e[0], offsp + e[1], offsp + e[2]) ;
    fprintf(f, "Surface(%d) = {%d} ;\n", offs + i, offs + i) ;
    return ;
  }
  fprintf(f, "Curve Loop(%d) = {%d, %d, -%d, -%d} ;\n", offs + i,
	  offsp + e[0], offsp + e[1], offsp + e[2], offsp + e[3]) ;
  fprintf(f, "Surface(%d) = {%d} ;\n", offs + i, offs + i) ;

  return ;
}

gint agg_wireframe_write_gmsh(FILE *f, agg_wireframe_t *w,
			      gchar *len, gint offp, gint offsp, gint offs)

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
  for ( i = 0 ; i < agg_wireframe_spline_number(w) ; i ++ ) {
    sp = agg_wireframe_spline(w, i) ;
    n = agg_wireframe_spline_length(w, i) ;
    fprintf(f, "Spline(%d) = {", offsp+i) ;
    for ( j = 0 ; j < n-1 ; j ++ ) {
      fprintf(f, "%d, ", offp + sp[j]) ;
    }
    fprintf(f, "%d} ;\n", offp + sp[n-1]) ;
  }

  /*write elements*/
  for ( i = 0 ; i < agg_wireframe_element_number(w) ; i ++ ) {
    write_element_gmsh(f, &(w->e[4*i]), i, offsp, offs) ;
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
  
