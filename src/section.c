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

#include <agg.h>

#include <blaswrap.h>

#include "agg-private.h"

typedef void (*section_parse_func_t)(agg_section_t *s,
				     agg_variable_t *p, gint np) ;

void _agg_aerofoil_parse(agg_section_t *s, agg_variable_t *p, gint np) ;
void _agg_circle_parse(agg_section_t *s, agg_variable_t *p, gint np) ;
void _agg_ellipse_parse(agg_section_t *s, agg_variable_t *p, gint np) ;

static const struct {
  gchar *name ;
  gint np ;
  section_parse_func_t func ;
} parse_data[] = {
  {"aerofoil", 3, _agg_aerofoil_parse},
  {"circle",   0, _agg_circle_parse  },
  {"ellipse",  1, _agg_ellipse_parse},
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
 * Set a section to an ellipse of maximum vertical thickness \a th,
 * centred at \f$(0,1/2)\f$, as in Kulfan, B. (2010). Recent
 * extensions and applications of the `CST' universal parametric
 * geometry representation method. The Aeronautical Journal,
 * 114(1153), 157-176 https://doi.org/10.1017/S0001924000003614
 * 
 * @param s ::agg_section_t to be initialized;
 * @param th section "thickness".
 * 
 * @return 0 on success. 
 */

gint agg_section_set_ellipse(agg_section_t *s, gdouble th)

{
  /*basic elliptical section*/
  agg_section_eta_left(s) = 0.5 ;
  agg_section_eta_right(s) = 0.5 ;

  agg_section_order_upper(s) = 0 ;
  agg_section_coefficient_upper(s,0) = th ;
  agg_section_order_lower(s) = 0 ;
  agg_section_coefficient_lower(s,0) = -th ;

  agg_section_trailing_edge_upper(s) = 0 ;
  agg_section_trailing_edge_lower(s) = 0 ;
  
  return 0 ;
}

/** 
 * Set a section to a circle of radius 1/2 centred at \f$(0,1/2)\f$,
 * as in Kulfan, B. (2010). Recent extensions and applications of the
 * `CST' universal parametric geometry representation method. The
 * Aeronautical Journal, 114(1153), 157-176
 * https://doi.org/10.1017/S0001924000003614
 * 
 * @param s ::agg_section_t to be initialized.
 * 
 * @return 0 on success. 
 */

gint agg_section_set_circle(agg_section_t *s)

{
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
  agg_section_set_circle(s) ;
  
  return ;
}

void _agg_ellipse_parse(agg_section_t *s, agg_variable_t *p, gint np)

{
  gdouble th ;
  
  if ( agg_variable_definition(&(p[0])) != NULL )
    g_error("%s: th must be a numerical constant", __FUNCTION__) ;
  th = agg_variable_value(&(p[0])) ;

  agg_section_set_ellipse(s, th) ;
  
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
  gdouble dy, *c, C, dC, yte, sgn, n1, n2 ;
  gint i, order ;
  
  /*code not tested for more complex geometries yet*/
  g_assert(agg_section_order_lower(s) < 1) ;
  g_assert(agg_section_order_upper(s) < 1) ;
  
  if ( x < -1.0 || x > 1.0 ) 
    g_error("%s: input parameter x (%lg) out of range (-1,1)",
	    __FUNCTION__, x) ;

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
 * Write a section to file as a list of points \f$(|x|, y(x))\f$,
 * \f$-\leq x\leq 1\f$.
 * 
 * @param f output file stream;
 * @param s the section to write;
 * @param T if not NULL, a transform applied to points before writing;
 * @param npts number of points on section.
 * 
 * @return 0 on success.
 */

gint agg_section_write(FILE *f, agg_section_t *s, agg_transform_t *T,
		       gint npts)

{
  gint i ;
  gdouble t, x[3] ;

  for ( i = 0 ; i < npts ; i ++ ) {
    t = -1 + 2.0*(gdouble)i/(npts-1) ;
    x[0] = ABS(t) ;
    x[1] = agg_section_eval(s, t) ;
    if ( T != NULL ) agg_transform_apply(T, x, x) ;
    fprintf(f, "%lg %lg\n", x[0], x[1]) ;
  }
  
  return 0 ;
}

/** 
 * Write a section to file as a list of points \f$(|x|, y(x))\f$,
 * \f$-\leq x\leq 1\f$, using a C printf format string
 * 
 * @param f output file stream;
 * @param s the section to write;
 * @param T if not NULL, a transform applied to points before writing;
 * @param fstr format string for output;
 * @param estr if not NULL, format string used for final point of output;
 * @param npts number of points on section.
 * 
 * @return 0 on success.
 */

gint agg_section_format_write(FILE *f, agg_section_t *s, agg_transform_t *T,
			      gchar *fstr, gchar *estr, gint npts)

{
  gint i ;
  gdouble t, x[3] ;

  for ( i = 0 ; i < npts ; i ++ ) {
    t = -1 + 2.0*(gdouble)i/(npts-1) ;
    x[0] = ABS(t) ;
    x[1] = agg_section_eval(s, t) ;
    if ( T != NULL ) agg_transform_apply(T, x, x) ;
    fprintf(f, fstr, x[0], x[1]) ;
  }

  i = npts-1 ; 
  t = -1 + 2.0*(gdouble)i/(npts-1) ;
  x[0] = ABS(t) ;
  x[1] = agg_section_eval(s, t) ;
  if ( T != NULL ) agg_transform_apply(T, x, x) ;
  
  if ( estr == NULL ) 
    fprintf(f, fstr, x[0], x[1]) ;
  else
    fprintf(f, estr, x[0], x[1]) ;
  
  return 0 ;
}

/** 
 * Fit CST form to section data.
 * 
 * @param s ::agg_section_t to contain fitted section data;
 * @param xu \f$x\f$ coordinates of points on upper surface;
 * @param xustr stride of data in \a xu;
 * @param yu \f$y\f$ coordinates of points on upper surface;
 * @param yustr stride of data in \a yu;
 * @param npu number of points on upper surface;
 * @param xl \f$x\f$ coordinates of points on lower surface;
 * @param xlstr stride of data in \a xu;
 * @param yl \f$y\f$ coordinates of points on lower surface;
 * @param ylstr stride of data in \a yu;
 * @param npl number of points on lower surface;
 * @param n1 exponent at \f$x=0\f$;
 * @param n2 exponent at \f$x=1\f$;
 * @param nu number of coefficients to compute in upper surface expansion;
 * @param nl number of coefficients to compute in lower surface expansion;
 * 
 * @return 0 on success.
 */

gint agg_section_fit(agg_section_t *s,
		     gdouble *xu, gint xustr, gdouble *yu, gint yustr, gint npu,
		     gdouble *xl, gint xlstr, gdouble *yl, gint ylstr, gint npl,
		     gdouble n1, gdouble n2, gint nu, gint nl)

{
  gdouble yteU, yteL, A[2048], work[4096], b[256], C ;
  gint i, j, m, n, i1 = 1, lwork, info, ldb ;

  yteU = yteL = 0 ;
  for ( i = 0 ; i < npu ; i ++ ) {
    if ( xu[i*xustr] == 1.0 ) yteU = yu[i*yustr] ;
  }
  for ( i = 0 ; i < npl ; i ++ ) {
    if ( xl[i*xlstr] == 1.0 ) yteL = yl[i*ylstr] ;
  }

  /*form equations for upper surface*/
  for ( i = 0 ; i < npu ; i ++ ) {
    C = pow(xu[i*xustr], n1)*pow(1.0-xu[i*xustr], n2) ;
    agg_bernstein_basis(nu, xu[i*xustr], work, NULL) ;
    for ( j = 0 ; j <= nu ; j ++ ) {
      A[j*npu + i] = work[j]*C ;
    }
    b[i] = yu[i*yustr] - xu[i*xustr]*yteU ;
  }

  /*solve least squares problem (note this is transposed in Fortran sense)*/
  m = npu ; n = nu + 1 ; lwork = 4096 ; ldb = MAX(m, n) ;
  dgels_("N", &m, &n, &i1, A, &m, b, &ldb, work, &lwork, &info) ;
  for ( i = 0 ; i <= nu ; i ++ ) agg_section_coefficient_upper(s,i) = b[i] ;

  /*form equations for lower surface*/
  for ( i = 0 ; i < npl ; i ++ ) {
    C = pow(xl[i*xlstr], n1)*pow(1.0-xl[i*xlstr], n2) ;
    agg_bernstein_basis(nl, xl[i*xlstr], work, NULL) ;
    for ( j = 0 ; j <= nl ; j ++ ) {
      A[j*npl + i] = work[j]*C ;
    }
    b[i] = yl[i*ylstr] - xl[i*xlstr]*yteL ;
  }

  /*solve least squares problem (note this is transposed in Fortran sense)*/
  m = npl ; n = nl + 1 ; lwork = 4096 ; ldb = MAX(m, n) ;
  dgels_("N", &m, &n, &i1, A, &m, b, &ldb, work, &lwork, &info) ;
  for ( i = 0 ; i <= nl ; i ++ ) agg_section_coefficient_lower(s,i) = b[i] ;  
  
  agg_section_order_upper(s) = nu ;
  agg_section_order_lower(s) = nl ;
  agg_section_eta_left(s) = n1 ;
  agg_section_eta_right(s) = n2 ;
  agg_section_trailing_edge_upper(s) = yteU ;  
  agg_section_trailing_edge_lower(s) = yteL ;  
  
  return 0 ;
}

/**
 * @}
 */

