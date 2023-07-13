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

static const struct {
  gchar *name ;
  agg_patch_mapping_t map ;
} mapping_list[] =
  {
    {"bilinear",  AGG_PATCH_BILINEAR},
    {"spherical", AGG_PATCH_SPHERICAL},
    {"tubular",   AGG_PATCH_TUBULAR},
    {NULL,        -1}
  } ;


/** 
 * @{ 
 *
 * @ingroup patches
 */

agg_patch_t *agg_patch_new(gint nstmax)

{
  agg_patch_t *P ;

  P = (agg_patch_t *)g_malloc0(sizeof(agg_patch_t)) ;

  P->st = (gdouble *)g_malloc0(AGG_PATCH_POINT_SIZE*nstmax*sizeof(gdouble)) ;
  agg_patch_point_number(P) = 0 ;
  agg_patch_point_number_max(P) = nstmax ;

  agg_patch_mapping(P) = AGG_PATCH_BILINEAR ;
  agg_patch_wrap_s(P) = FALSE ;
  agg_patch_wrap_t(P) = FALSE ;
  agg_patch_invert(P) = FALSE ;
  
  return P ;
}

static gint agg_bilinear_map(agg_patch_t *P, gdouble s, gdouble t,
			     gdouble *u, gdouble *v)
{
  /* gdouble w0, w1, w2, w3 ; */

  /* /\*basic bilinear mapping for now */
  /*   https://en.wikipedia.org/wiki/Bilinear_interpolation */
  /* *\/ */
  /* w0 = (1.0-s)*(1.0-t) ; */
  /* w1 = (    s)*(1.0-t) ; */
  /* w2 = (    s)*(    t) ; */
  /* w3 = (1.0-s)*(    t) ; */
  
  /* *u = */
  /*   agg_patch_corner_u(P,0)*w0 + */
  /*   agg_patch_corner_u(P,1)*w1 + */
  /*   agg_patch_corner_u(P,2)*w2 + */
  /*   agg_patch_corner_u(P,3)*w3 ; */
  /* *v = */
  /*   agg_patch_corner_v(P,0)*w0 + */
  /*   agg_patch_corner_v(P,1)*w1 + */
  /*   agg_patch_corner_v(P,2)*w2 + */
  /*   agg_patch_corner_v(P,3)*w3 ; */

  g_assert_not_reached() ; /*untested code*/
  
  *u = s ;
  *v = -1.0 + 2.0*t ;
  
  return 0 ;
}

static gint agg_spherical_map(agg_patch_t *P, gdouble s, gdouble t,
			      gdouble *u, gdouble *v)

{
  *u = 0.5*(1.0-cos(M_PI*s)) ;

  if ( t < 0.5 ) {
    *v = -0.5*(1.0 + cos(2.0*M_PI*t)) ;
  } else {
    *v =  0.5*(1.0 + cos(2.0*M_PI*t)) ;
  }
  
  return 0 ;
}

static gint agg_tubular_map(agg_patch_t *P, gdouble s, gdouble t,
			    gdouble *u, gdouble *v)

{
  *u = s ;

  if ( t < 0.5 ) {
    *v = -0.5*(1.0 + cos(2.0*M_PI*t)) ;
  } else {
    *v =  0.5*(1.0 + cos(2.0*M_PI*t)) ;
  }
  
  return 0 ;
}

/** 
 * Apply mapping in a patch to find parametric variables on surface.
 * 
 * @param P patch whose mapping is to be applied;
 * @param s variable on \a P, \f$0\leq s \leq 1\f$;
 * @param t variable on \a P, \f$0\leq t \leq 1\f$;
 * @param u surface parametric variable;
 * @param v surface parametric variable.
 * 
 * @return 0 on success.
 */

gint agg_patch_map(agg_patch_t *P, gdouble s, gdouble t,
		   gdouble *u, gdouble *v)

{
  if ( s < 0 && agg_patch_wrap_s(P) ) { s = 1.0 + s ; }
  if ( s > 1 && agg_patch_wrap_s(P) ) { s = s - 1.0; }
  if ( t < 0 && agg_patch_wrap_t(P) ) { t = 1.0 + t ; }
  if ( t > 1 && agg_patch_wrap_t(P) ) { t = t - 1.0; }

  if ( s < 0 || s > 1 )
    g_error("%s: s (%lg) out of range", __FUNCTION__, s) ;
  if ( t < 0 || t > 1 )
    g_error("%s: t (%lg) out of range", __FUNCTION__, t) ;

  if ( agg_patch_mapping(P) == AGG_PATCH_BILINEAR ) {
    return agg_bilinear_map(P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_SPHERICAL ) {
    return agg_spherical_map(P, s, t, u, v) ;
  }

  if ( agg_patch_mapping(P) == AGG_PATCH_TUBULAR ) {
    return agg_tubular_map(P, s, t, u, v) ;
  }
  
  g_error("%s: unrecognized mapping %d", __FUNCTION__, agg_patch_mapping(P)) ;
  
  return 0 ;
}

static gint wrapping_parse(gchar *str, gboolean *wrap)

{
  if ( strcmp("wrap", str) == 0 ) { *wrap = TRUE ; return 0 ; }
  if ( strcmp("nowrap", str) == 0 ) { *wrap = FALSE ; return 0 ; }

  return -1 ;
}

/** 
 * Parse string data for an ::agg_patch_t.
 * 
 * @param P on exit, is set according to input variables;
 * @param p array of variables for parameters setting \a P;
 * @param np number of variables in \a p.
 * 
 * @return 0 on success.
 */

gint agg_patch_parse(agg_patch_t *P, agg_variable_t *p, gint np)

{
  gboolean wrap_s, wrap_t, invert ;
  gint i ;

  if ( np < 3 ) {
    g_error("%s: three parameters required for mapping", __FUNCTION__) ;
  }

  invert = FALSE ;
  if ( wrapping_parse(agg_variable_definition(&(p[1])), &wrap_s) != 0 ) {
    g_error("%s: unrecognized mapping \"%s\"",
	    __FUNCTION__, agg_variable_definition(&(p[1]))) ;
  }

  if ( wrapping_parse(agg_variable_definition(&(p[2])), &wrap_t) != 0 ) {
    g_error("%s: unrecognized mapping \"%s\"",
	    __FUNCTION__, agg_variable_definition(&(p[2]))) ;
  }

  if ( np > 3 ) {
    if ( strcmp(agg_variable_definition(&(p[3])), "invert") == 0 ) {
      invert = TRUE ;
    } else {
      g_error("%s: unrecognized inversion parameter \"%s\"",
	      __FUNCTION__, agg_variable_definition(&(p[3]))) ;
    }
  }
  
  for ( i = 0 ; mapping_list[i].name != NULL ; i ++ ) {
    if ( strcmp(mapping_list[i].name, agg_variable_definition(&(p[0]))) == 0 ) {
      agg_patch_mapping(P) = mapping_list[i].map ;
      agg_patch_wrap_s(P)  = wrap_s ;
      agg_patch_wrap_t(P)  = wrap_t ;
      agg_patch_invert(P)  = invert ;
      return 0 ;
    }
  }

  g_error("%s: unrecognized mapping \"%s\"", __FUNCTION__,
	  agg_variable_definition(&(p[0]))) ;
  
  return 0 ;
}


/**
 * @}
 */
