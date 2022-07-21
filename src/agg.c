/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2021, 2022 Michael Carley
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


gdouble agg_bernstein_basis_eval(gint n, gint r, gdouble x)

{
  gdouble S ;

  S = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;
  
  return S ;
}

gint agg_bernstein_basis(gint n, gdouble x, gdouble *S, gdouble *dS)

{
  gint r ;

  if ( dS == NULL ) {
    for ( r = 0 ; r <= n ; r ++ ) 
      S[r] = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;

    return 0 ;
  }

  g_assert_not_reached() ; /*derivatives not implemented yet*/
  
  for ( r = 0 ; r <= n ; r ++ ) {
    S[r] = _binomial(n,r)*pow(x,r)*pow(1-x,n-r) ;
  }
  
  return 0 ;
}

agg_workspace_t *agg_workspace_alloc(gint ns)

{
  agg_workspace_t *w ;

  w = (agg_workspace_t *)g_malloc0(sizeof(agg_workspace_t)) ;

  agg_workspace_shape(w) = agg_shape_alloc(ns) ;
  agg_workspace_local_transform(w) = agg_local_transform_alloc() ;
  
  return w ;
}
