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

#include "binomials.h"

/** 
 * @ingroup bernstein
 *
 * Evaluate a Bernstein basis monomial, 
 * \f$B_{r}^{(n)}(x)=\binom{n}{r}x^{r}(1-x)^{n-r}\f$
 * 
 * @param n degree of Bernstein basis;
 * @param r mononomial from basis to evaluate;
 * @param x argument of monomial
 * 
 * @return \f$B_{r}^{(n)}(x)\f$
 */

gdouble agg_bernstein_basis_eval(gint n, gint r, gdouble x)

{
  gdouble S ;

  S = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;
  
  return S ;
}

gdouble agg_bernstein_derivative_eval(gint n, gint r, gdouble x)

{
  gdouble dS ;

  if ( n == 0 ) return 0 ;
  
  if ( r == 0 ) {
    dS = -_binomial(n,r)*(n-r)*pow(1-x,n-r-1) ;
    return dS ;
  }

  if ( r == n ) {
    dS = _binomial(n,r)*r*pow(x,r-1) ;
    return dS ;
  }

  dS = _binomial(n,r)*(r*pow(x,r-1)*pow(1-x,n-r) -
		       (n-r)*pow(x, r)*pow(1-x,n-r-1)) ;
  
  return dS ;
}

gint agg_bernstein_basis(gint n, gdouble x, gdouble *S, gdouble *dS)

{
  gint r ;

  if ( dS == NULL ) {
    for ( r = 0 ; r <= n ; r ++ ) 
      S[r] = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;

    return 0 ;
  }

  r = 0 ;
  S[r] = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;
  dS[r] = -_binomial(n,r)*(n-r)*pow(1-x,n-r-1) ;
  r = n ;
  S[r] = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;
  dS[r] = _binomial(n,r)*r*pow(x,r-1) ;
  for ( r = 1 ; r < n ; r ++ ) {
    S[r] = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;
    dS[r] = _binomial(n,r)*(r*pow(x,r-1)*pow(1-x,n-r) -
			    (n-r)*pow(x,r)*pow(1-x,n-r-1)) ;
  }
  
  return 0 ;
}
