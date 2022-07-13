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

agg_cut_t *agg_cut_alloc(void)

{
  agg_cut_t *h ;

  h = (agg_cut_t *)g_malloc0(sizeof(agg_cut_t)) ;
  agg_cut_start(h, 0) = -1.0 ;
  agg_cut_end  (h, 0) =  1.0 ;
  h->nm = 1 ;
  
  return h ;
}

gint agg_cut_add(agg_cut_t *h, gdouble smin, gdouble smax)

{
  gint i, j ;

  g_assert(agg_cut_number(h) > 0) ;
  g_assert(smax >= smin) ;

  if ( smin == agg_cut_start(h,0) ) {
    agg_cut_start(h,0) = smax ;
    return 0 ;
  }

  if ( smax == agg_cut_end(h,agg_cut_number(h)-1) ) {
    agg_cut_end(h,agg_cut_number(h)-1) = smin ;
    return 0 ;
  }
  
  for ( i = 0 ; i < agg_cut_number(h) ; i ++ ) {
    if ( agg_cut_start(h,i) <= smin && agg_cut_end(h,i) >= smax ) break ;
  }

  /*check we have not tried to make an invalid cut*/
  g_assert(i < agg_cut_number(h)) ;
  
  /*shuffle the mappings up one*/
  for ( j = agg_cut_number(h) ; j > i ; j -- ) {
    h->sh[2*j+0] = h->sh[2*(j-1)+0] ;
    h->sh[2*j+1] = h->sh[2*(j-1)+1] ;
  }

  h->sh[2* i   +1] = smin ; 
  h->sh[2*(i+1)+0] = smax ; 

  agg_cut_number(h) ++ ;
  
  return 0 ;
}
