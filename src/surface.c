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

/**
 * @{
 *  @ingroup surfaces
 */

/** 
 * Allocate a new ::agg_surface_t
 * 
 * @param nsmax maximum number of sections in surface.
 * 
 * @return newly allocated surface.
 */

agg_surface_t *agg_surface_new(gint nsmax)

{
  agg_surface_t *S ;
  gint oumax, olmax ;
  gdouble *buf ;
  gint i ;

  oumax = 32 ; olmax = 32 ;
  
  S = (agg_surface_t *)g_malloc0(sizeof(agg_surface_t)) ;

  S->s = (agg_section_t *)g_malloc0(nsmax*sizeof(agg_section_t)) ;
  S->u = (gdouble *)g_malloc0(2*nsmax*sizeof(gdouble)) ;
  S->w = &(S->u[nsmax]) ;
  agg_surface_section_number(S) = 0 ;
  agg_surface_section_number_max(S) = nsmax ;

  agg_surface_transform(S) = agg_transform_new(32) ;
  
  /*allocate memory for the section expansion coefficients*/
  buf = (gdouble *)g_malloc0(nsmax*(oumax+olmax+2)*sizeof(gdouble)) ;

  for ( i = 0 ; i < nsmax ; i ++ ) {
    S->s[i].cu = &(buf[i*(oumax+olmax+2)        ]) ;
    S->s[i].cl = &(buf[i*(oumax+olmax+2)+oumax+1]) ;
    S->s[i].oumax = oumax ; 
    S->s[i].olmax = olmax ; 
  }

  /*default values*/
  agg_surface_umin(S) = 0.0 ;
  agg_surface_umax(S) = 1.0 ;

  return S ;
}

/** 
 * Add a section to a surface
 * 
 * Add a section to a surface at a given value of parameter
 * \f$u\f$. The section data are copied into the surface, so the
 * section variable \a s can be reused afterwards. 
 * 
 * @param S an ::agg_surface_t;
 * @param s an ::agg_section_t;
 * @param u parameter value where section is to be added.
 * 
 * @return 0 on success.
 */

gint agg_surface_section_add(agg_surface_t *S, agg_section_t *s,
			     gdouble u)

{
  gint i ;
  
  if ( agg_surface_section_number(S) >= agg_surface_section_number_max(S) )
    g_error("%s: maximum number of sections (%d) in surface exceeded",
	    __FUNCTION__, agg_surface_section_number_max(S)) ;

  if ( agg_surface_section_number(S) != 0 ) {
    /*check order of new section against existing one*/
    agg_section_t *st = agg_surface_section(S,0) ;
    if ( agg_section_order_upper(st) != agg_section_order_upper(s) ||
	 agg_section_order_lower(st) != agg_section_order_lower(s) )
      g_error("%s: all sections on a surface must have same expansion order",
	      __FUNCTION__) ;
  }

  i = agg_surface_section_number(S) ;

  agg_section_copy(agg_surface_section(S,i), s) ;
  agg_surface_u(S,i) = u ;
  agg_surface_section_number(S) ++ ;
  
  return 0 ;
}

/** 
 * Generate interpolation weights for sections on a surface
 * 
 * Generate interpolation weights for sections on a surface, the
 * preprocessing stage for Berrut and Trefethen's barycentric Lagrangian
 * interpolation, https://dx.doi.org/10.1137/S0036144502417715 
 *
 * @param S an ::agg_surface_t containing a number of sections.
 * 
 * @return 0 on success.
 */

gint agg_surface_weights_make(agg_surface_t *S)

{
  gint j, k ;

  for ( j = 0 ; j < agg_surface_section_number(S) ; j ++ ) {
    agg_surface_weight(S,j) = 1.0 ;
    for ( k = 0 ; k < j ; k ++ ) {
      agg_surface_weight(S,j) *= agg_surface_u(S,j) - agg_surface_u(S,k) ;
    }
    for ( k = j+1 ; k < agg_surface_section_number(S) ; k ++ ) {
      agg_surface_weight(S,j) *= agg_surface_u(S,j) - agg_surface_u(S,k) ;
    }
    agg_surface_weight(S,j) = 1.0/agg_surface_weight(S,j) ;
  }
  
  return 0 ;
}

/** 
 * Interpolate section data from a surface into a local section
 * 
 * Interpolate section data from a surface into a local section, using
 * Berrut and Trefethen's barycentric Lagrangian interpolation,
 * https://dx.doi.org/10.1137/S0036144502417715
 * 
 * @param S an ::agg_surface_t containing a number of sections;
 * @param u parameter at which to interpolate section;
 * @param s on exit, contains interpolated data for section at \a u.
 * 
 * @return 0 on success, 1 if \a S is an empty surface.
 */

gint agg_surface_section_interp(agg_surface_t *S, gdouble u, agg_section_t *s)

{
  agg_section_t *su ;
  gint i, j ;
  gdouble w, wtot ;

  /*empty surface*/
  if ( agg_surface_section_number(S) == 0 ) return 1 ;

  su = agg_surface_section(S, 0) ;
  if ( agg_section_order_upper_max(s) < agg_section_order_upper(su) )
    g_error("%s: upper surface expansion order too high (%d) for section",
	    __FUNCTION__, agg_section_order_upper(su)) ;
  if ( agg_section_order_lower_max(s) < agg_section_order_lower(su) )
    g_error("%s: lower surface expansion order too high (%d) for section",
	    __FUNCTION__, agg_section_order_lower(su)) ;

  /*constant cross-section*/
  if ( agg_surface_section_number(S) == 1 ) {
    agg_section_copy(s, su) ;
    return 0 ;
  }

  agg_section_close(s) = agg_section_close(su) ;

  for ( j = 0 ; j <= agg_section_order_upper(su) ; j ++ ) {
    agg_section_coefficient_upper(s,j) = 0.0 ; 
  }
  for ( j = 0 ; j <= agg_section_order_lower(su) ; j ++ ) {
    agg_section_coefficient_lower(s,j) = 0.0 ; 
  }

  agg_section_eta_left(s) = agg_section_eta_right(s) = 0.0 ;
  agg_section_trailing_edge_upper(s) =
    agg_section_trailing_edge_lower(s) = 0.0 ;

  wtot = 0.0 ;
  for ( i = 0 ; i < agg_surface_section_number(S) ; i ++ ) {
    if ( agg_surface_weight(S,i) == 0.0 ) {
      g_error("%s: surface interpolation weight %d is zero; "
	      "are the weights initialized?", __FUNCTION__, i) ;
    }
    su = agg_surface_section(S, i) ;
    if ( ABS(u - agg_surface_u(S,i)) < 1e-12 ) {
      agg_section_copy(s, su) ;
      return 0 ;
    }
    w = agg_surface_weight(S,i)/(u - agg_surface_u(S,i)) ;
    wtot += w ;

    for ( j = 0 ; j <= agg_section_order_upper(su) ; j ++ ) {
      agg_section_coefficient_upper(s,j) +=
	w*agg_section_coefficient_upper(su,j) ;
    }
    for ( j = 0 ; j <= agg_section_order_lower(su) ; j ++ ) {
      agg_section_coefficient_lower(s,j) +=
	w*agg_section_coefficient_lower(su,j) ;
    }
    
    agg_section_eta_left(s) += w*agg_section_eta_left(su) ;
    agg_section_eta_right(s) += w*agg_section_eta_right(su) ;
    agg_section_trailing_edge_upper(s) +=
      w*agg_section_trailing_edge_upper(su) ;
    agg_section_trailing_edge_lower(s) +=
      w*agg_section_trailing_edge_lower(su) ;    
  }

  for ( j = 0 ; j <= agg_section_order_upper(su) ; j ++ ) {
    agg_section_coefficient_upper(s,j) /= wtot ; 
  }
  for ( j = 0 ; j <= agg_section_order_lower(su) ; j ++ ) {
    agg_section_coefficient_lower(s,j) /= wtot ; 
  }

  agg_section_eta_left(s) /= wtot ;
  agg_section_eta_right(s) /= wtot ;
  agg_section_trailing_edge_upper(s) /= wtot ;
  agg_section_trailing_edge_lower(s) /= wtot ;
  agg_section_order_upper(s) = agg_section_order_upper(su) ;
  agg_section_order_lower(s) = agg_section_order_lower(su) ;

  return 0 ;
}

/** 
 * Allocate a workspace for evaluation of surfaces
 * 
 * @return newly allocated workspace.
 */

agg_surface_workspace_t *agg_surface_workspace_new(void)

{
  agg_surface_workspace_t  *w ;

  w = (agg_surface_workspace_t *)g_malloc(sizeof(agg_surface_workspace_t)) ;

  w->s  = agg_section_new(64, 64) ;
  w->ds = agg_section_new(64, 64) ;
  
  return w ;
}

/** 
 * Evaluate a point on a surface
 * 
 * Evaluate a point on a surface. The surface must have been
 * initialized with ::agg_surface_weights_make, and its transform
 * initialized with ::agg_expression_data_compile and
 * ::agg_transform_expressions_compile, to set up the internal
 * variables and expressions used for point transformation. The
 * section data for \a S are interpolated at \f$u\f$ and the point on
 * the section evaluated at \f$v\f$. The point is then transformed
 * using the surface transformation.
 * 
 * @param S an ::agg_surface_t;
 * @param u surface parameter \f$u\f$;
 * @param v surface parameter \f$v\f$;
 * @param x on exit contains \f$\mathbf{x}(u,v)\f$;
 * @param w workspace allocated using ::agg_surface_workspace_new
 * 
 * @return 0 on success.
 */

gint agg_surface_point_eval(agg_surface_t *S, gdouble u, gdouble v,
			    gdouble *x, agg_surface_workspace_t *w)

{
  agg_transform_t *T ;
  agg_variable_t *var ;
  agg_section_t *s ;
  gdouble y[3] ;
  
  T = agg_surface_transform(S) ;
  s = w->s ;
  
  /*check that the first two variables in T are u and v*/
  var = agg_transform_variable(T, 0) ;
  if ( strcmp(agg_variable_name(var), "u") != 0 ) {
    g_error("%s: first variable in transform must be \"u\"", __FUNCTION__) ;
  } else {
    agg_variable_value(var) = u ;
  }
  /* var = agg_transform_variable(T, 1) ; */
  /* if ( strcmp(agg_variable_name(var), "v") != 0 ) { */
  /*   g_error("%s: second variable in transform must be \"v\"", __FUNCTION__) ; */
  /* } else { */
  /*   agg_variable_value(var) = v ; */
  /* } */

  agg_transform_variables_eval(T) ;  

  agg_surface_section_interp(S, u, s) ;

  y[0] = fabs(v) ;
  y[1] = agg_section_eval(s, v) ;
  y[2] = 0.0 ;

  agg_transform_apply(T, y, x) ;

  /* agg_transform_axes(agg_surface_axes(S), x, x) ; */
  
  return 0 ;
}

/** 
 * Estimate derivatives on a surface
 * 
 * @param S surface on which derivatives are to be evaluated;
 * @param u surface parameter;
 * @param v surface parameter;
 * @param x on exit contains \f$\mathbf{x}(u,v)\f$;
 * @param xu on exit contains \f$\partial\mathbf{x}/\partial u\f$;
 * @param xv on exit contains \f$\partial\mathbf{x}/\partial v\f$;
 * @param w workspace for surface point evaluation.
 * 
 * @return 0 on success.
 */

gint agg_surface_point_diff(agg_surface_t *S, gdouble u, gdouble v,
			    gdouble *x, gdouble *xu, gdouble *xv,
			    agg_surface_workspace_t *w)

/*
 * should do this analytically some day
 */
  
{
  agg_transform_t *T ;
  agg_variable_t *var ;
  gdouble ee ;
  
  T = agg_surface_transform(S) ;
  ee = 1e-6 ;
  
  /*check that the first two variables in T are u and v*/
  var = agg_transform_variable(T, 0) ;
  if ( strcmp(agg_variable_name(var), "u") != 0 ) {
    g_error("%s: first variable in transform must be \"u\"", __FUNCTION__) ;
  } else {
    agg_variable_value(var) = u ;
  }
  var = agg_transform_variable(T, 1) ;
  if ( strcmp(agg_variable_name(var), "v") != 0 ) {
    g_error("%s: second variable in transform must be \"v\"", __FUNCTION__) ;
  } else {
    agg_variable_value(var) = v ;
  }

  agg_surface_point_eval(S, u, v, x, w) ;

  if ( u < agg_surface_umin(S) - ee ) {
    agg_surface_point_eval(S, u+ee, v, xu, w) ;
    xu[0] =  (xu[0] - x[0])/ee ;
    xu[1] =  (xu[1] - x[1])/ee ;
    xu[2] =  (xu[2] - x[2])/ee ;
  } else {
    agg_surface_point_eval(S, u-ee, v, xu, w) ;
    xu[0] = -(xu[0] - x[0])/ee ;
    xu[1] = -(xu[1] - x[1])/ee ;
    xu[2] = -(xu[2] - x[2])/ee ;
  }

  if ( v < 1 - ee ) {
    agg_surface_point_eval(S, u, v+ee, xv, w) ;
    xv[0] =  (xv[0] - x[0])/ee ;
    xv[1] =  (xv[1] - x[1])/ee ;
    xv[2] =  (xv[2] - x[2])/ee ;
  } else {
    agg_surface_point_eval(S, u, v-ee, xv, w) ;
    xv[0] = -(xv[0] - x[0])/ee ;
    xv[1] = -(xv[1] - x[1])/ee ;
    xv[2] = -(xv[2] - x[2])/ee ;
  }
  
  return 0 ;
}

/** 
 * Parse a string specifying a surface grid type
 * 
 * @param str grid string to parse
 * 
 * @return ::agg_grid_t corresponding to \a str, or AGG_GRID_UNDEFINED. 
 */

agg_grid_t agg_grid_parse(char *str)

{
  gint i ;
  static const struct {
    char *name ;
    agg_grid_t grid ;
  } grid_list[] = {
    {"regular",    AGG_GRID_REGULAR},
    {"triangle",   AGG_GRID_TRIANGLE},
    {"sphere-ico", AGG_GRID_SPHERE_ICO},
    {"sphere-uv", AGG_GRID_SPHERE_UV},
    {"hemisphere-ico", AGG_GRID_HEMISPHERE_ICO},
    {"hemisphere-uv", AGG_GRID_HEMISPHERE_UV},
    {NULL,       AGG_GRID_UNDEFINED}
  } ;

  for ( i = 0 ; grid_list[i].name != NULL ; i ++ ) {
    if ( strcmp(grid_list[i].name, str) == 0 )
      return grid_list[i].grid ;
  }

  return AGG_GRID_UNDEFINED ;
}

/**
 *  @}
 */
