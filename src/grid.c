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

static gboolean check_constant(gchar *str, gpointer expr)

{
  if ( expr != NULL ) {
    fprintf(stderr, "%s must be a numerical constant\n", str) ;    
    return FALSE ;
  }

  return TRUE ;
}

static gboolean check_non_negative(gchar *str, gdouble p, gdouble *v)

{
  if ( p < 0 ) {
    fprintf(stderr, "%s must be non-negative (%lg)\n", str, p) ;
    return FALSE ;
  }
  
  *v = p ;
  return TRUE ;
}

static gboolean check_positive(gchar *str, gdouble p, gdouble *v)

{
  if ( p <= 0 ) {
    fprintf(stderr, "%s must be greater than zero (%lg)\n", str, p) ;
    return FALSE ;
  }
  if ( v != NULL ) *v = p ;

  return TRUE ;
}

static gboolean check_int(gchar *str, gdouble p, gint *v)

{
  if ( p != round(p) ) {
    fprintf(stderr, "%s must be an integer (%lg)\n", str, p) ;
    return FALSE ;
  }
  
  *v = (gint)p ;
  
  return TRUE ;
}

static gdouble linear_to_angle(gdouble x)

{
  gdouble t ;
  
  t = acos(2.0*fabs(x) - 1.0) ;
  if ( x > 0 ) t = 2.0*M_PI - t ;

  return t ;
}

agg_grid_t *agg_grid_alloc(gint np, gint nt)

{
  agg_grid_t *g ;

  g = (agg_grid_t *)g_malloc0(sizeof(agg_grid_t)) ;

  g->tri = (gint *)g_malloc0(3*nt*sizeof(gint)) ;
  g->uv  = (gdouble *)g_malloc0(2*np*sizeof(gdouble)) ;

  agg_grid_point_number_max(g) = np ;
  agg_grid_triangle_number_max(g) = nt ;
  agg_grid_invert(g) = FALSE ;
  agg_grid_topology(g) = AGG_GRID_NONE ;
  agg_grid_data(g) = NULL ;

  return g ;
}

gint agg_grid_init(agg_grid_t *g)

{
  agg_grid_point_number(g) = 0 ;
  agg_grid_triangle_number(g) = 0 ;
  
  return 0 ;
}

gint agg_grid_square(agg_grid_t *g,
		     gdouble umin, gdouble umax, gint nu, agg_spacing_t su,
		     gdouble vmin, gdouble vmax, gint nv, agg_spacing_t sv)

{
  gint i, j, *tri ;
  gdouble u, v ;
  
  g_assert(agg_grid_point_number_max(g) >= nu*nv) ;
  g_assert(agg_grid_triangle_number_max(g) >= 2*nu*nv) ;
  
  agg_grid_point_number(g) = 0 ;

  for ( i = 0 ; i < nu ; i ++ ) {
    u = agg_spacing_eval(umin, umax, nu, su, i) ;
    for ( j = 0 ; j < nv ; j ++ ) {
      v = agg_spacing_eval(vmin, vmax, nv, sv, j) ;
      agg_grid_point_u(g,i*nv+j) = u ;
      agg_grid_point_v(g,i*nv+j) = v ;
      agg_grid_point_number(g) ++ ;
    }
  }

  agg_grid_triangle_number(g) = 0 ;

  for ( i = 0 ; i < nu-1 ; i ++ ) {
    for ( j = 0 ; j < nv-1 ; j ++ ) {
      tri = agg_grid_triangle(g, agg_grid_triangle_number(g)) ;
      tri[0] = (i+0)*nv + j     ;
      tri[1] = (i+1)*nv + j     ;
      tri[2] = (i+1)*nv + j + 1 ;
      agg_grid_triangle_number(g) ++ ;
      tri = agg_grid_triangle(g, agg_grid_triangle_number(g)) ;
      tri[0] = (i+0)*nv + j     ;
      tri[1] = (i+1)*nv + j + 1 ;
      tri[2] = (i+0)*nv + j + 1 ;
      agg_grid_triangle_number(g) ++ ;
    }
  }
  
  return 0 ;
}

gint agg_grid_linear(agg_grid_t *g,
		     gdouble umin, gdouble umax, gint nu, agg_spacing_t su,
		     gdouble vmin, gdouble vmax, gint nv, agg_spacing_t sv)

{
  gint i, j, np, nt, *tri ;
  
  agg_grid_topology(g) = AGG_GRID_LINEAR ;
  agg_grid_point_number(g) = np = 0 ;

  for ( i = 0 ; i < nu ; i ++ ) {
    for ( j = 0 ; j < nv ; j ++ ) {
      agg_grid_point_u(g, np) = agg_spacing_eval(umin, umax, nu, su, i) ;
      agg_grid_point_v(g, np) = agg_spacing_eval(vmin, vmax, nv, sv, j) ;
      np ++ ;
    }
  }

  agg_grid_point_number(g) = np ;
  agg_grid_triangle_number(g) = nt = 0 ;

  for ( i = 0 ; i < nu - 1 ; i ++ ) {
    for ( j = 0 ; j < nv - 1 ; j ++ ) {
      tri = agg_grid_triangle(g, nt) ;
      tri[0] = (i+0)*nv + j + 0 ;
      tri[1] = (i+1)*nv + j + 0 ;
      tri[2] = (i+0)*nv + j + 1 ;
      nt ++ ;
      tri = agg_grid_triangle(g, nt) ;
      tri[0] = (i+1)*nv + j + 0 ;
      tri[1] = (i+1)*nv + j + 1 ;
      tri[2] = (i+0)*nv + j + 1 ;
      nt ++ ;
    }
  }

  agg_grid_triangle_number(g) = nt ;
  g->interp_area = agg_grid_interp_area_linear ;
  
  return 0 ;
}

gint agg_grid_tube(agg_grid_t *g, gint nu, gdouble vmin, gdouble vmax, gint nv)

{
  gint np, nt, i, j, *tri ;
  gdouble t[256], thmin, thmax ;

  g_assert(nv < 256) ;
  
  agg_grid_topology(g) = AGG_GRID_TUBE ;
  agg_grid_point_number(g) = np = 0 ;

  thmin = linear_to_angle(vmin) ;
  thmax = linear_to_angle(vmax) ;

  for ( i = 0 ; i < nv ; i ++ ) {
    /* t[i] = 2.0*M_PI*i/(nv-1) ; */
    t[i] = thmin + (thmax - thmin)*i/(nv-1) ;
    if ( t[i] > M_PI ) t[i] = (1.0 + cos(t[i]))*0.5 ;
    else t[i] = -(1.0 + cos(t[i]))*0.5 ;
  }
  
  for ( i = 0 ; i < nu ; i ++ ) {
    for ( j = 0 ; j < nv ; j ++ ) {
      agg_grid_point_u(g, np) = (gdouble)i/(nu-1) ;
      agg_grid_point_v(g, np) = t[j] ;
      np ++ ;
    }
  }

  agg_grid_point_number(g) = np ;
  agg_grid_triangle_number(g) = nt = 0 ;

  for ( i = 0 ; i < nu - 1 ; i ++ ) {
    for ( j = 0 ; j < nv - 1 ; j ++ ) {
      tri = agg_grid_triangle(g, nt) ;
      tri[0] = (i+0)*nv + j + 0 ;
      tri[1] = (i+1)*nv + j + 0 ;
      tri[2] = (i+0)*nv + j + 1 ;
      nt ++ ;
      tri = agg_grid_triangle(g, nt) ;
      tri[0] = (i+1)*nv + j + 0 ;
      tri[1] = (i+1)*nv + j + 1 ;
      tri[2] = (i+0)*nv + j + 1 ;
      nt ++ ;
    }
  }

  agg_grid_triangle_number(g) = nt ;
  g->interp_area = agg_grid_interp_area_tube ;
  
  return 0 ;
}

gint agg_grid_interp_area_linear(agg_grid_t *g, gdouble *uv,
				 gdouble s, gdouble t,
				 gdouble *u, gdouble *v)

{
  gdouble L[3] ;

  L[0] = 1.0 - s - t ; L[1] = s ; L[2] = t ;

  *u = L[0]*uv[0] + L[1]*uv[2] + L[2]*uv[4] ;
  *v = L[0]*uv[1] + L[1]*uv[3] + L[2]*uv[5] ;

  return 0 ;
}

gint agg_grid_interp_area_tube(agg_grid_t *g, gdouble *uv,
			       gdouble s, gdouble t,
			       gdouble *u, gdouble *v)

{
  gdouble L[3], vt[3] ;

  L[0] = 1.0 - s - t ; L[1] = s ; L[2] = t ;

  vt[0] = linear_to_angle(uv[1]) ;
  vt[1] = linear_to_angle(uv[3]) ;
  vt[2] = linear_to_angle(uv[5]) ;
  
  *v = L[0]*vt[0] + L[1]*vt[1] + L[2]*vt[2] ;

  *u = L[0]*uv[0] + L[1]*uv[2] + L[2]*uv[4] ;
  
  /*v is an angle which needs reconverting to linear coordinate*/
  if ( *v > M_PI ) *v = 0.5*(1.0 + cos(*v)) ;
  else *v = -0.5*(1.0 + cos(*v)) ;
  
  return 0 ;
}

gint agg_grid_area_interpolate(agg_grid_t *g, gdouble *uv,
			       gdouble s, gdouble t,
			       gdouble *u, gdouble *v)

{

  return g->interp_area(g, uv, s, t, u, v) ;
  
  return 0 ;
}

gint agg_grid_line_interpolate(agg_grid_t *g, gdouble *uv,
			       gdouble s,
			       gdouble *u, gdouble *v)

{

  return g->interp_line(g, uv, s, u, v) ;
  
  return 0 ;
}

gint agg_grid_element_interpolate(agg_grid_t *g, gint i,
				  gdouble s, gdouble t,
				  gdouble *u, gdouble *v)

{
  gdouble uv[6] ;
  gint *tri ;

  tri = agg_grid_triangle(g,i) ;
  uv[0] = agg_grid_point_u(g,tri[0]) ;
  uv[1] = agg_grid_point_v(g,tri[0]) ;
  uv[2] = agg_grid_point_u(g,tri[1]) ;
  uv[3] = agg_grid_point_v(g,tri[1]) ;
  uv[4] = agg_grid_point_u(g,tri[2]) ;
  uv[5] = agg_grid_point_v(g,tri[2]) ;

  return g->interp_area(g, uv, s, t, u, v) ;
  
  return 0 ;
}

static gint set_grid_spherical(agg_grid_t *g, gchar *expr[], gdouble *p,
			       gint np)

{
  gint refine ;

  if ( np != 1 ) {
    fprintf(stderr,
	    "spherical grid requires one argument (refinement level)\n") ;
    return 1 ;
  }

  if ( !check_constant("refinement level", expr[0]) ) return 1 ;
  if ( !check_int("refinement level", p[0], &refine) ) return 1 ;

  return agg_grid_spherical(g, refine) ;
}

static gint set_grid_hemispherical(agg_grid_t *g, gchar *expr[], gdouble *p,
				   gint np)

{
  gint refine ;

  if ( np != 1 ) {
    fprintf(stderr,
	    "hemispherical grid requires one argument (refinement level)\n") ;
    return 1 ;
  }

  if ( !check_constant("refinement level", expr[0]) ) return 1 ;
  if ( !check_int("refinement level", p[0], &refine) ) return 1 ;

  return agg_grid_hemispherical(g, refine) ;
}

static gint set_grid_linear(agg_grid_t *g, gchar *expr[], gdouble *p, gint np)

{
  gint ns, nt ;
  agg_spacing_t ss, st ;
  gdouble smin, smax, tmin, tmax ;
  
  /*
   * linear grid can be defined as follows:
   *  ns sdist nt tdist
   *  ns smin smax sdist nt tdist
   *  ns sdist nt tmin tmax tdist
   *  ns smin smax sdist nt tmin tmax tdist
   * 
   */

  smin =   0.0 ; smax = 1.0 ;
  tmin =  -1.0 ; tmax = 1.0 ;
  ns = 0 ; nt = 0 ;
  
  if ( np != 4 && np != 8 ) {
    fprintf(stderr,
	    "linear grid requires four or eight arguments\n"
	    /* "  (# points in s) (distribution in s) " */
	    /* "(# points in t) (distribution in t)\n" */
	    ) ;
    return 1 ;
  }

  if ( np == 4 ) {
    if ( !check_constant("number of points in s", expr[0]) ) return 1 ;
    if ( !check_positive("number of points in s", p[0], NULL) ) return 1 ;
    if ( !check_int("number of points in s", p[0], &ns) ) return 1 ;
    if ( expr[1] == NULL ) {
      fprintf(stderr, "point spacing in s must be a string") ;
      return 1 ;
    }
    if ( agg_parser_constant_parse(expr[1], &ss) != 0 ) {
      fprintf(stderr, "cannot parse s spacing %s\n", expr[1]) ;
      return 1 ;
    }

    if ( !check_constant("number of points in t", expr[2]) ) return 1 ;
    if ( !check_positive("number of points in t", p[2], NULL) ) return 1 ;
    if ( !check_int("number of points in t", p[2], &nt) ) return 1 ;

    if ( expr[3] == NULL ) {
      fprintf(stderr, "point spacing in t must be a string") ;
      return 1 ;
    }
    if ( agg_parser_constant_parse(expr[3], &st) != 0 ) {
      fprintf(stderr, "cannot parse t spacing %s\n", expr[3]) ;
      return 1 ;
    }
  }

  if ( np == 8 ) {
    if ( !check_constant("number of points in s", expr[0]) ) return 1 ;
    if ( !check_positive("number of points in s", p[0], NULL) ) return 1 ;
    if ( !check_int("number of points in s", p[0], &ns) ) return 1 ;

    if ( !check_constant("smin", expr[1]) ) return 1 ;
    smin = p[1] ;
    if ( !check_constant("smax", expr[2]) ) return 1 ;
    smax = p[2] ;

    if ( expr[3] == NULL ) {
      fprintf(stderr, "point spacing in s must be a string") ;
      return 1 ;
    }
    if ( agg_parser_constant_parse(expr[3], &ss) != 0 ) {
      fprintf(stderr, "cannot parse s spacing %s\n", expr[3]) ;
      return 1 ;
    }

    if ( !check_constant("number of points in t", expr[4]) ) return 1 ;
    if ( !check_positive("number of points in t", p[4], NULL) ) return 1 ;
    if ( !check_int("number of points in t", p[4], &nt) ) return 1 ;

    if ( !check_constant("tmin", expr[5]) ) return 1 ;
    tmin = p[5] ;
    if ( !check_constant("tmax", expr[6]) ) return 1 ;
    tmax = p[6] ;
    if ( expr[7] == NULL ) {
      fprintf(stderr, "point spacing in t must be a string") ;
      return 1 ;
    }
    if ( agg_parser_constant_parse(expr[7], &st) != 0 ) {
      fprintf(stderr, "cannot parse t spacing %s\n", expr[7]) ;
      return 1 ;
    }    
  }
  
  return agg_grid_linear(g,
			 smin, smax, ns, ss,
			 tmin, tmax, nt, st) ;
}

static gint set_grid_tube(agg_grid_t *g, gchar *expr[], gdouble *p, gint np)

{
  gint ns, nt ;
  gdouble tmin, tmax ;

  tmin = -1 ; tmax = 1 ;
  if ( np != 2 && np != 4) {
    fprintf(stderr,
	    "tubular grid requires two or four arguments arguments "
	    "(# points in s) (# points in t)\n") ;
    return 1 ;
  }

  if ( !check_constant("number of points in s", expr[0]) ) return 1 ;
  if ( !check_positive("number of points in s", p[0], NULL) ) return 1 ;
  if ( !check_int("number of points in s", p[0], &ns) ) return 1 ;

  if ( !check_constant("number of points in t", expr[1]) ) return 1 ;
  if ( !check_positive("number of points in t", p[1], NULL) ) return 1 ;
  if ( !check_int("number of points in t", p[1], &nt) ) return 1 ;

  if ( np == 4 ) {
    if ( !check_constant("tmin", expr[2]) ) return 1 ;
    tmin = p[2] ;
    if ( !check_constant("tmax", expr[3]) ) return 1 ;
    tmax = p[3] ;
  }

  return agg_grid_tube(g, ns, tmin, tmax, nt) ;
}

static gint set_grid_adaptive(agg_grid_t *g, gchar *expr[], gdouble *p, gint np)

{
  agg_adaptive_grid_data_t data ;
  
  if ( np != 5 ) {
    fprintf(stderr, "adaptive grid requires eight arguments:\n"
	    "  (nsmax, sinterp, ntmax, tinterp, amin, amax, lmin, lmax)\n") ;
    return 1 ;
  }
  
  if ( !check_constant("minimum element area", expr[0]) ) return 1 ;
  if ( !check_non_negative("minimum element area", p[0],
			   &(agg_adaptive_grid_amin(&data))) ) return 1 ;

  if ( !check_constant("maximum element area", expr[1]) ) return 1 ;
  if ( !check_positive("maximum element area", p[1],
		       &(agg_adaptive_grid_amax(&data))) ) return 1 ;

  if ( !check_constant("minimum edge length", expr[2]) ) return 1 ;
  if ( !check_non_negative("minimum edge length", p[2],
			   &(agg_adaptive_grid_lmin(&data))) ) return 1 ;

  if ( !check_constant("maximum edge length", expr[3]) ) return 1 ;
  if ( !check_positive("maximum edge length", p[3],
		       &(agg_adaptive_grid_lmax(&data))) ) return 1 ;

  if ( !check_constant("maximum edge length", expr[4]) ) return 1 ;
  if ( !check_positive("maximum edge length", p[4],
		       &(agg_adaptive_grid_ncos(&data))) ) return 1 ;
  agg_adaptive_grid_ncos(&data) =
    cos(agg_adaptive_grid_ncos(&data)*M_PI/180.0) ;
  
  return agg_grid_adaptive(g, &data) ;

  return 0 ;
}

gint agg_grid_parse(agg_grid_t *g, gchar *name, gchar *expr[],
		    gdouble *p, gint np)

{
  gint (* grid_func)(agg_grid_t *, gchar **, gdouble *, gint) ;
  gint i ;
  gpointer grids[] = {
    "adaptive", set_grid_adaptive,
    "hemispherical", set_grid_hemispherical,
    "linear", set_grid_linear,
    "spherical", set_grid_spherical,
    "tube", set_grid_tube,
    NULL, NULL} ;

  for ( i = 0 ; grids[2*i+0] != NULL ; i ++ ) {
    if ( strcmp(name, (gchar *)(grids[2*i+0])) == 0 ) {
      grid_func = grids[2*i+1] ;
      return grid_func(g, expr, p, np) ;
    }
  }
  
  if ( strcmp(name, "invert") == 0 ) {
    if ( agg_grid_topology(g) != AGG_GRID_NONE ) {
      fprintf(stderr, "%s: set grid inversion before choosing grid topology\n",
	      __FUNCTION__) ;
      return 1 ;
    }
    
    agg_grid_invert(g) = TRUE ;

    return 0 ;
  }
  
  return 1 ;
}
