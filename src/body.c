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

agg_body_t *agg_body_alloc(void)

{
  agg_body_t *b ;

  b = (agg_body_t *)g_malloc0(sizeof(agg_body_t)) ;

  agg_body_distribution_number(b) = 0 ;
  
  return b ;
}

gint agg_body_distribution_add(agg_body_t *b, agg_distribution_t *d,
			       gchar *name)

{
  g_assert(agg_body_distribution_number(b) <
	   AGG_BODY_DISTRIBUTION_NUMBER_MAX) ;

  b->names[agg_body_distribution_number(b)] = g_strdup(name) ;
  b->d[agg_body_distribution_number(b)] = d ;
  
  agg_body_distribution_number(b) ++ ;

  return 0 ;
}

gint agg_body_distributions_list(FILE *f, agg_body_t *b)

{
  gint i ;

  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) {
    fprintf(f, "%s\n", b->names[i]) ;
  }
  
  return 0 ;
}

gint agg_body_distribution_locate_u(agg_body_t *b, gdouble u)

{
  gint i ;
  agg_distribution_t *d ;
  
  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) {
    d = agg_body_distribution(b, i) ;
    if ( agg_distribution_parameter_min(d) <= u &&
	 u <= agg_distribution_parameter_max(d) ) return i ;
  }    
  
  return -1 ;
}
