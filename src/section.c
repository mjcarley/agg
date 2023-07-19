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

typedef void (*section_parse_func_t)(agg_section_t *s,
				     agg_variable_t *p, gint np) ;

void _agg_aerofoil_parse(agg_section_t *s, agg_variable_t *p, gint np) ;
void _agg_circle_parse(agg_section_t *s, agg_variable_t *p, gint np) ;

static const struct {
  gchar *name ;
  gint np ;
  section_parse_func_t func ;
} parse_data[] = {
  {"aerofoil", 3, _agg_aerofoil_parse},
  {"circle",   0, _agg_circle_parse  },
  {NULL, 0, NULL}
} ;

/**
 * @{
 *
 * @ingroup sections
 *
 */

/** 
 * Allocate a new ::agg_section_t
 * 
 * @param oumax maximum order of upper surface expansion;
 * @param olmax maximum order of lower surface expansion.
 * 
 * @return newly allocated ::agg_section_t
 */

agg_section_t *agg_section_new(gint oumax, gint olmax)

{
  agg_section_t *s ;

  s = (agg_section_t *)g_malloc0(sizeof(agg_section_t)) ;
  s->cu = (gdouble *)g_malloc0((oumax+olmax+2)*sizeof(gdouble)) ;
  s->cl = &(s->cu[oumax+1]) ;
  agg_section_order_upper_max(s) = oumax ;
  agg_section_order_lower_max(s) = olmax ;

  return s ;
}

static gdouble agg_section_ellipse_eval(agg_section_t *s, gdouble x)

{
  gdouble y, *c, C ;
  gint i, order ;
  
  if ( x < 0.0 ) {
    order = agg_section_order_lower(s) ;
    c = &(agg_section_coefficient_lower(s,0)) ;
    x = fabs(x) ;
  } else {
    order = agg_section_order_upper(s) ;
    c = &(agg_section_coefficient_upper(s,0)) ;
  }

  y = 0.0 ;
  C = pow(x, agg_section_eta_left(s))*pow(1.0-x, agg_section_eta_right(s)) ;
  for ( i = 0 ; i <= order ; i ++ ) {
    y += c[i]*agg_bernstein_basis_eval(order, i, x)*C ;
  }
  
  return y ;
}

static gdouble agg_section_aerofoil_eval(agg_section_t *s, gdouble x)

{
  gdouble y, *c, C, yte ;
  gint i, order ;
  
  if ( x < 0.0 ) {
    order = agg_section_order_lower(s) ;
    c = &(agg_section_coefficient_lower(s,0)) ;
    yte = agg_section_trailing_edge_lower(s) ;
    x = fabs(x) ;
  } else {
    order = agg_section_order_upper(s) ;
    c = &(agg_section_coefficient_upper(s,0)) ;
    yte = agg_section_trailing_edge_upper(s) ;
  }

  y = x*yte ;
  C = pow(x, agg_section_eta_left(s))*pow(1.0-x, agg_section_eta_right(s)) ;
  for ( i = 0 ; i <= order ; i ++ ) {
    y += c[i]*agg_bernstein_basis_eval(order, i, x)*C ;
  }

  return y ;
}

/** 
 * Evaluate a section geometry
 * 
 * @param s a ::agg_section_t containing the section data;
 * @param x evaluation point \f$-1\leq x\leq 1\f$. 
 * 
 * @return \f$y=f(x)\f$ on success. An error is reported and execution
 * ceases if \a x is out of range or the section type is undefined.
 */

gdouble agg_section_eval(agg_section_t *s, gdouble x)

{
  if ( agg_section_type(s) == AGG_SECTION_UNDEFINED )
    g_error("%s: undefined section type", __FUNCTION__) ;

  if ( x < -1.0 || x > 1.0 ) 
    g_error("%s: input parameter x (%lg) out of range (-1,1)",
	    __FUNCTION__, x) ;

  if ( agg_section_type(s) == AGG_SECTION_ELLIPSE ) {
    return agg_section_ellipse_eval(s, x) ;
  }

  if ( agg_section_type(s) == AGG_SECTION_AEROFOIL ) {
    return agg_section_aerofoil_eval(s, x) ;
  }
  
  g_assert_not_reached() ;
  
  return 0.0 ;
}

/** 
 * Set a section to a circle of radius 1/2 centred at \f$(0,1/2)\f$,
 * as in Kulfan, B. (2010). Recent extensions and applications of the
 * `CST' universal parametric geometry representation method. The
 * Aeronautical Journal, 114(1153), 157-176
 * https://doi.org/10.1017/S0001924000003614
 * 
 * @param s ::agg_section_t to be initialized;
 * 
 * @return 0 on success. 
 */

gint agg_section_set_circle(agg_section_t *s)

{
  agg_section_type(s) = AGG_SECTION_ELLIPSE ;

  /*basic elliptical section*/
  agg_section_eta_left(s) = 0.5 ;
  agg_section_eta_right(s) = 0.5 ;

  agg_section_order_upper(s) = 0 ;
  agg_section_coefficient_upper(s,0) = 1.0 ;
  agg_section_order_lower(s) = 0 ;
  agg_section_coefficient_lower(s,0) = -1.0 ;
  
  return 0 ;
}

/** 
 * Set a section to a simple aerofoil shape of specified leading edge
 * exponent, maximum thickness, and trailing edge displacement, as in
 * Kulfan, B. (2010). Recent extensions and applications of the `CST'
 * universal parametric geometry representation method. The
 * Aeronautical Journal, 114(1153), 157-176
 * https://doi.org/10.1017/S0001924000003614
 * 
 * @param s ::agg_section_t to be initialized;
 * @param eta leading edge exponent \f$\eta_{1}\f$;
 * @param th aerofoil thickness;
 * @param yte trailing edge displacement (half trailing edge thickness). 
 * 
 * @return 0 on success.
 */

gint agg_section_set_aerofoil(agg_section_t *s, gdouble eta,
			      gdouble th, gdouble yte)

{
  agg_section_type(s) = AGG_SECTION_AEROFOIL ;

  /*basic elliptical section*/
  agg_section_eta_left(s) = eta ;
  agg_section_eta_right(s) = 1.0 ;

  agg_section_order_upper(s) = 0 ;
  agg_section_coefficient_upper(s,0) = th ;
  agg_section_order_lower(s) = 0 ;
  agg_section_coefficient_lower(s,0) = -th ;

  agg_section_trailing_edge_upper(s) =  yte ;
  agg_section_trailing_edge_lower(s) = -yte ;

  if ( yte != 0.0 )
    agg_section_close(s) = TRUE ;
  else
    agg_section_close(s) = FALSE ;
  
  return 0 ;
}

/** 
 * Copy one surface section into another
 * 
 * @param dest destination section;
 * @param src source section.
 * 
 * @return 0 on success, with \a src copied into \a dest.
 */

gint agg_section_copy(agg_section_t *dest, agg_section_t *src)

{
  gint i ;
  
  g_assert(dest != NULL) ;
  g_assert(src != NULL) ;

  if ( agg_section_order_upper_max(dest) < agg_section_order_upper(src) )
    g_error("%s: dest does not have enough space (%d) for upper coefficients",
	    __FUNCTION__, agg_section_order_upper(src)) ;
  
  if ( agg_section_order_lower_max(dest) < agg_section_order_lower(src) )
    g_error("%s: dest does not have enough space (%d) for lower coefficients",
	    __FUNCTION__, agg_section_order_lower(src)) ;

  agg_section_type(dest) = agg_section_type(src) ;
  agg_section_close(dest) = agg_section_close(src) ;
  agg_section_eta_left(dest) = agg_section_eta_left(src) ;
  agg_section_eta_right(dest) = agg_section_eta_right(src) ;

  agg_section_order_upper(dest) = agg_section_order_upper(src) ;
  for ( i = 0 ; i <= agg_section_order_upper(dest) ; i ++ ) {
    agg_section_coefficient_upper(dest,i) = 
      agg_section_coefficient_upper(src,i) ;
  }
  agg_section_order_lower(dest) = agg_section_order_lower(src) ;
  for ( i = 0 ; i <= agg_section_order_lower(dest) ; i ++ ) {
    agg_section_coefficient_lower(dest,i) = 
      agg_section_coefficient_lower(src,i) ;
  }

  agg_section_trailing_edge_upper(dest) =
    agg_section_trailing_edge_upper(src) ;
  agg_section_trailing_edge_lower(dest) =
    agg_section_trailing_edge_lower(src) ;
  
  return 0 ;
}

void _agg_aerofoil_parse(agg_section_t *s, agg_variable_t *p, gint np)

{
  gdouble eta, th, yte ;

  if ( agg_variable_definition(&(p[0])) != NULL )
    g_error("%s: eta must be a numerical constant", __FUNCTION__) ;
  eta = agg_variable_value(&(p[0])) ;
  if ( agg_variable_definition(&(p[1])) != NULL )
    g_error("%s: thickness must be a numerical constant", __FUNCTION__) ;
  th = agg_variable_value(&(p[1])) ;
  if ( agg_variable_definition(&(p[2])) != NULL )
    g_error("%s: trailing edge thickness must be a numerical constant",
	    __FUNCTION__) ;
  yte = agg_variable_value(&(p[2])) ;

  agg_section_set_aerofoil(s, eta, th, yte) ;
  
  return ;
}

void _agg_circle_parse(agg_section_t *s, agg_variable_t *p, gint np)

{
  if ( np > 0 ) 
    g_error("%s: circular section has no parameters", __FUNCTION__) ;

  agg_section_set_circle(s) ;
  
  return ;
}

/** 
 * Parse an ::agg_section_t description given as a name and list of
 * variable parameters
 * 
 * @param s ::agg_section_t to set;
 * @param name name of section;
 * @param p array of ::agg_variable_t containing parameters for section;
 * @param np number of entries in \a p.
 * 
 * @return 0 on success.
 */

gint agg_section_parse(agg_section_t *s, gchar *name,
		       agg_variable_t *p, gint np)

{
  section_parse_func_t func ;
  gint i, nparam ;

  for ( (i = 0), (func = NULL) ; parse_data[i].name != NULL ; i ++ ) {
    if ( strcmp(parse_data[i].name, name) == 0 ) {
      func = parse_data[i].func ;
      nparam = parse_data[i].np ;
    }
  }

  if ( func == NULL ) {
    g_error("%s: unrecognized section name \"%s\"", __FUNCTION__, name) ;
  }
  
  if ( np != nparam ) {
    g_error("%s: section \"%s\" requires %d parameters (not %d)",
	    __FUNCTION__, name, nparam, np) ;
  }

  func(s, p, np) ;
  
  return 0 ;
}

/** 
 * Write a list of available section types to output
 * 
 * @param f output file stream.
 * 
 * @return 0 on success.
 */

gint agg_sections_list(FILE *f)

{
  gint i ;

  for ( i = 0 ; parse_data[i].name != NULL ; i ++ ) {
    fprintf(stderr, "  %s(%d parameters)\n",
	    parse_data[i].name, parse_data[i].np) ;
  }

  return 0 ;
}


static gdouble agg_section_ellipse_diff(agg_section_t *s, gdouble x)

{
  gdouble dy, *c, C, dC, sgn, n1, n2 ;
  gint i, order ;
  
  if ( x < 0.0 ) {
    order = agg_section_order_lower(s) ;
    c = &(agg_section_coefficient_lower(s,0)) ;
    x = fabs(x) ;
    sgn = -1 ;
  } else {
    order = agg_section_order_upper(s) ;
    c = &(agg_section_coefficient_upper(s,0)) ;
    sgn = 1 ;
  }

  dy = 0.0 ;
  n1 = agg_section_eta_left(s) ;
  n2 = agg_section_eta_right(s) ;
  C  = pow(x, n1)*pow(1.0-x, n2) ;
  dC = n1*pow(x, n1-1)*pow(1.0-x, n2) - n2*pow(x, n1)*pow(1.0-x, n2-1) ;

  for ( i = 0 ; i <= order ; i ++ ) {
    dy += sgn*c[i]*(agg_bernstein_derivative_eval(order, i, x)*C +
		    agg_bernstein_basis_eval(order, i, x)*dC) ;
  }
  
  return dy ;
}

static gdouble agg_section_aerofoil_diff(agg_section_t *s, gdouble x)

{
  gdouble dy, *c, C, dC, yte, sgn, n1, n2 ;
  gint i, order ;
  
  if ( x < 0.0 ) {
    order = agg_section_order_lower(s) ;
    c = &(agg_section_coefficient_lower(s,0)) ;
    yte = agg_section_trailing_edge_lower(s) ;
    x = fabs(x) ;
    sgn = -1 ;
  } else {
    order = agg_section_order_upper(s) ;
    c = &(agg_section_coefficient_upper(s,0)) ;
    yte = agg_section_trailing_edge_upper(s) ;
    sgn = 1 ;
  }

  dy = sgn*yte ;
  n1 = agg_section_eta_left(s) ;
  n2 = agg_section_eta_right(s) ;
  C  = pow(x, n1)*pow(1.0-x, n2) ;
  dC = n1*pow(x, n1-1)*pow(1.0-x, n2) - n2*pow(x, n1)*pow(1.0-x, n2-1) ;
  for ( i = 0 ; i <= order ; i ++ ) {
    dy += sgn*c[i]*(agg_bernstein_derivative_eval(order, i, x)*C +
		    agg_bernstein_basis_eval(order, i, x)*dC) ;
  }

  return dy ;
}

/** 
 * Derivative of a section curve
 * 
 * @param s a ::agg_section_t containing the section data;
 * @param x evaluation point \f$-1\leq x\leq 1\f$. 
 * 
 * @return \f$y'=\mathrm{d}f(x)/\mathrm{d}x\f$ on success. An error is
 * reported and execution ceases if \a x is out of range or the
 * section type is undefined.
 */

gdouble agg_section_diff(agg_section_t *s, gdouble x)

{
  if ( agg_section_type(s) == AGG_SECTION_UNDEFINED )
    g_error("%s: undefined section type", __FUNCTION__) ;

  /*code not tested for more complex geometries yet*/
  g_assert(agg_section_order_lower(s) < 1) ;
  g_assert(agg_section_order_upper(s) < 1) ;
  
  if ( x < -1.0 || x > 1.0 ) 
    g_error("%s: input parameter x (%lg) out of range (-1,1)",
	    __FUNCTION__, x) ;

  if ( agg_section_type(s) == AGG_SECTION_ELLIPSE ) {
    return agg_section_ellipse_diff(s, x) ;
  }

  if ( agg_section_type(s) == AGG_SECTION_AEROFOIL ) {
    return agg_section_aerofoil_diff(s, x) ;
  }
  
  g_assert_not_reached() ;
  
  return 0.0 ;
}

/**
 * @}
 */

