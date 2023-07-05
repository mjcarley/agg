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

agg_grid_t *agg_grid_new(gint nuvmax, gint nijmax)

{
  agg_grid_t *g ;

  g = (agg_grid_t *)g_malloc0(sizeof(agg_grid_t)) ;

  g->uv = (gdouble *)g_malloc0(nuvmax*AGG_GRID_POINT_SIZE*sizeof(gdouble)) ;
  g->ij = (gint *)g_malloc0(nuvmax*AGG_GRID_ELEMENT_SIZE*sizeof(gint)) ;

  agg_grid_point_number_max(g) = nuvmax ;
  agg_grid_element_number_max(g) = nijmax ;
  
  return g ;
}

gint agg_grid_basic(agg_grid_t *g,
		    gdouble umin, gdouble umax, gint nu,
		    gdouble vmin, gdouble vmax, gint nv)

{
  gint i, j, n, *e ;
  gdouble u, v ;

  if ( agg_grid_point_number_max(g) < (nu+1)*(nv+1) )
    g_error("%s: not enough space for requested number of points",
	    __FUNCTION__) ;
  if ( agg_grid_element_number_max(g) < nu*nv) 
    g_error("%s: not enough space for requested number of elements",
	    __FUNCTION__) ;
  
  agg_grid_point_number(g) = 0 ;

  for ( i = n = 0 ; i <= nu ; i ++ ) {
    u = umin + (umax - umin)*i/nu ;
    for ( j = 0 ; j <= nv ; j ++ ) {
      v = vmin + (vmax - vmin)*j/nv ;
      agg_grid_point_u(g,n) = u ;
      agg_grid_point_v(g,n) = v ;
      n ++ ;
    }
  }
  agg_grid_point_number(g) = n ;

  for ( i = n = 0 ; i < nu ; i ++ ) {
    for ( j = 0 ; j < nv ; j ++ ) {
      e = agg_grid_element(g, n) ;
      e[0] = (i+0)*(nv+1) + j + 0 ;
      e[1] = (i+1)*(nv+1) + j + 0 ;
      e[2] = (i+1)*(nv+1) + j + 1 ;
      e[3] = (i+0)*(nv+1) + j + 1 ;
      n ++ ;
    }
  }

  agg_grid_element_number(g) = n ;
  
  return 0 ;
}

static void write_spline(FILE *f, agg_grid_t *g, gint p1, gint p2, gint pps,
			 agg_surface_t *S, agg_surface_workspace_t *w)

{
  gint i ;
  gdouble u, v, x[3] ;

  for ( i = 0 ; i <= pps ; i ++ ) {
    u = agg_grid_point_u(g, p1) +
      (agg_grid_point_u(g, p2) - agg_grid_point_u(g, p1))*i/pps ;
    v = agg_grid_point_v(g, p1) +
      (agg_grid_point_v(g, p2) - agg_grid_point_v(g, p1))*i/pps ;
    agg_surface_point_eval(S, u, v, x, w) ;
    fprintf(f, "%e %e %e %e %e\n", x[0], x[1], x[2], u, v) ;
  }
  
  return ;
}

gint agg_grid_surface_write(FILE *f, agg_grid_t *g, agg_surface_t *S,
			    agg_surface_workspace_t *w)

{
  gint pps, i, j, *e ;
  gdouble u, v ;

  pps = 3 ;

  for ( i = 0 ; i < agg_grid_element_number(g) ; i ++ ) {
    e = agg_grid_element(g, i) ;
    write_spline(f, g, e[0], e[1], pps, S, w) ;
  }    
  
  return 0 ;
}

static gint spline_lookup(gint *sp, gint pps, gint nsp, gint p1, gint p2)

{
  gint i ;

  for ( i = 0 ; i < nsp ; i ++ ) {
    if (sp[pps*i+0] == p1 && sp[pps*i+pps-1] == p2) return   i+1  ;
    if (sp[pps*i+0] == p2 && sp[pps*i+pps-1] == p1) return -(i+1) ;
  }
  
  return 0 ;
}

static gint add_spline(gint *sp, gint pps, gint nsp, agg_grid_t *g,
		       gint p1, gint p2, agg_surface_t *S,
		       agg_surface_workspace_t *w, gint *np, FILE *f)

{
  gint i ;
  gdouble u, v, x[3] ;

  sp[pps*nsp + 0] = p1 ; sp[pps*nsp + pps-1] = p2 ;
  for ( i = 1 ; i < pps-1 ; i ++ ) {
    u = agg_grid_point_u(g, p1) +
      (agg_grid_point_u(g, p2) - agg_grid_point_u(g, p1))*i/(pps-1) ;
    v = agg_grid_point_v(g, p1) +
      (agg_grid_point_v(g, p2) - agg_grid_point_v(g, p1))*i/(pps-1) ;
    agg_surface_point_eval(S, u, v, x, w) ;
    fprintf(f, "Point(%d) = {%e, %e, %e, lc} ;\n",
	    *np, x[0], x[1], x[2]) ;
    sp[pps*nsp + i] = (*np) ;
    (*np) ++ ;
  }

  return nsp + 1 ;
}

gint agg_grid_surface_gmsh_write(FILE *f, agg_grid_t *g, gint pps,
				 agg_surface_t *S,
				 agg_surface_workspace_t *w)

{
  gint splines[16384], i, *e, np, nsp, j, ie, npt ;
  gdouble x[3], u, v ;
  gchar *curve = "Spline" ;

  np = 0 ;
  for ( i = 0 ; i < agg_grid_point_number(g) ; i ++ ) {
    u = agg_grid_point_u(g, i) ;
    v = agg_grid_point_v(g, i) ;
    agg_surface_point_eval(S, u, v, x, w) ;
    fprintf(f, "Point(%d) = {%e, %e, %e, lc} ;\n",
	    np, x[0], x[1], x[2]) ;
    np ++ ;
  }

  for ( i = 0 ; i < agg_grid_element_number(g) ; i ++ ) {
    e = agg_grid_element(g, i) ;
    for ( j = 0 ; j < 3 ; j ++ ) {
      ie = spline_lookup(splines, pps, nsp, e[j], e[j+1]) ;
      if ( ie == 0 )
	nsp = add_spline(splines, pps, nsp, g, e[j], e[j+1], S, w, &np, f) ;
    }
    ie = spline_lookup(splines, pps, nsp, e[3], e[0]) ;
    if ( ie == 0 )
      nsp = add_spline(splines, pps, nsp, g, e[3], e[0], S, w, &np, f) ;
  }

  /*write the splines*/
  for ( i = 0 ; i < nsp ; i ++ ) {
    fprintf(f, "%s(%d) = {", curve, i) ;
    for ( j = 0 ; j < pps-1 ; j ++ ) {
      fprintf(f, "%d, ", splines[i*pps + j]) ;
    }
    fprintf(f, "%d} ;\n", splines[i*pps + pps-1]) ;
  }

  /*write the surfaces*/
  for ( i = 0 ; i < agg_grid_element_number(g) ; i ++ ) {
    e = agg_grid_element(g, i) ;
    fprintf(f, "Curve Loop(%d) = {", i) ;
    for ( j = 0 ; j < 3 ; j ++ ) {
      ie = spline_lookup(splines, pps, nsp, e[j], e[j+1]) ;
      if ( ie > 0 ) {
	fprintf(f, "%d,", ie-1) ;
      } else {
	ie = -ie ;
	fprintf(f, "-%d,", ie-1) ;
      }
    }
    ie = spline_lookup(splines, pps, nsp, e[3], e[0]) ;
    if ( ie > 0 ) {
      fprintf(f, "%d} ;\n", ie-1) ;
    } else {
      ie = -ie ;
      fprintf(f, "-%d} ;\n", ie-1) ;
    }
    fprintf(f, "Surface(%d) = {%d} ;\n", i, i) ;
  }
	
  return 0 ;
}
