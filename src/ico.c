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

#include "agg-private.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

/*numbering starts at 1 to allow for negative sign to indicate orientation*/
static gint edges_unwrap[] = {
  0, 0,
  0, 1,  /*1*/
  0, 2,
  0, 6,
  0, 7,
  0, 9,   /*5*/
  0, 15,
  1, 2,
  1, 3,
  1, 7,
  1, 8,  /*10*/
  2, 3,
  2, 5,
  2, 9,
  3, 4,
  3, 5,  /*15*/
  3, 8,
  3, 14,
  4, 5,
  4, 10,
  4, 11,  /*20*/
  4, 13,
  5, 9,
  5, 11,
  6, 9,
  6, 10,  /*25*/
  6, 11,
  6, 12,
  7, 8,
  7, 15,
  8, 14,  /*30*/
  9, 11,
  10, 11,
  10, 12,
  10, 13,
  12, 13} ;   /*35*/

static gint tri_unwrap[] = {
-7, -1, 2,
7, 11, -8,
-18, -14, 15,
-30, -16, 17,
29, -6, 4,
24, -5, 3,
-32, -19, 20,
32, -26, 25,
-22, -12, 13,
22, 31, -23,
-28, -9, 10,
35, -34, 33,
12, -15, -11,
-10, 8, 16,
-13, -2, 5,
9, -4, 1,
-31, -24, 26,
-33, -25, 27,
23, -20, 18,
34, -21, 19,
} ;

static gint edges_wrap[] = {
  0, 0,
  0, 1,
  0, 2,
  0, 6,
  0, 7,
  0, 9,
  1, 2,
  1, 3,
  1, 7,
  1, 8,
  2, 3,
  2, 5,
  2, 9,
  3, 4,
  3, 5,
  3, 8,
  4, 5,
  4, 8,
  4, 10,
  4, 11,
  5, 9,
  5, 11,
  6, 7,
  6, 9,
  6, 10,
  6, 11,
  7, 8,
  7, 10,
  8, 10,
  9, 11,
  10, 11} ;

static gint tri_wrap[] = {
  -6, -1,  2,
  6, 10, -7,
  -16, -13, 14,
  17, -15, 13,
  -22, -3,  4,
  23, -5,  3,
  -30, -18, 19,
  30, -25, 24,
  -20, -11, 12,
  20, 29, -21,
  -26, -8,   9,
  26, 28, -27,
  11, -14,  -10,
  -9, 7, 15,
  -12, -2, 5,
  8, -4, 1,
  -29, -23, 25,
   27, -24, 22,
   21, -19, 16,
   -28, -17,  18
} ;

static void add_point(gdouble *x, gdouble *th, gdouble *ph)

{
  gdouble r ;

  r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) ;
  *ph = acos(x[2]/r) ;
  *th = atan2(x[1], x[0]) ;

  if ( *th < 0 ) (*th) += 2.0*M_PI ;

  return ;
}

gint agg_sphere_ico_base(gdouble *th, gint tstr,
			 gdouble *ph, gint pstr,
			 gint *edges, gint edstr,
			 gint *elements, gint elstr,
			 gint *np, gint *ne, gint *nt,
			 gboolean convert,
			 gboolean unwrap)

/*
 * nodes and triangles of regular icosahedron on a unit sphere
 * 
 * 12 nodes, 30 edges, 20 triangular surface elements
 * 
 * https://danielsieger.com/blog/2021/01/03/generating-platonic-solids.html
 * 
 */
  
{
  gdouble x[3], phi, a, b, r ;
  gint i, *ed, *tri ;

  /* g_assert(unwrap) ; /\*standard version not implemented (and may well */
  /* 		       not be)*\/ */
  
  phi = 0.5*(1.0 + sqrt(5.0)) ;
  a = 1.0 ; b = 1.0/phi ;
  r = sqrt(a*a + b*b) ;

  a /= r ; b /= r ;
  
  (*np) = 0 ;
  x[0] =  0 ; x[1] =  b ; x[2] = -a ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] =  b ; x[1] =  a ; x[2] =  0 ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] = -b ; x[1] =  a ; x[2] =  0 ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] =  0 ; x[1] =  b ; x[2] =  a ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] =  0 ; x[1] = -b ; x[2] =  a ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] = -a ; x[1] =  0 ; x[2] =  b ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] =  0 ; x[1] = -b ; x[2] = -a ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] =  a ; x[1] =  0 ; x[2] = -b ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] =  a ; x[1] =  0 ; x[2] =  b ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] = -a ; x[1] =  0 ; x[2] = -b ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] =  b ; x[1] = -a ; x[2] =  0 ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;
  x[0] = -b ; x[1] = -a ; x[2] =  0 ;
  add_point(x, &(th[(*np)*tstr]), &(ph[(*np)*pstr])) ; (*np) ++ ;

  *ne = 31 ; *nt = 20 ;
  ed = edges_wrap ;
  tri = tri_wrap ;
  
  for ( i = 0 ; i < (*np) ; i ++ ) {
    ph[i*pstr] = 1.0 - ph[i*pstr]/M_PI ;
    th[i*tstr] = 1.0 - th[i*tstr]*0.5/M_PI ;
  }

  if ( unwrap ) {
    /*no sense in unwrapping points in (\theta,\phi) space*/
    g_assert(convert) ;
    /*extra points to wrap around patch boundaries in parametric space*/
    ph[(*np)*pstr] = ph[7*pstr] ;
    th[(*np)*tstr] = 0 ;
    (*np) ++ ;
    ph[(*np)*pstr] = ph[8*pstr] ;
    th[(*np)*tstr] = 0 ;
    (*np) ++ ;
    ph[(*np)*pstr] = ph[4*pstr]       ;
    th[(*np)*tstr] = th[4*pstr] + 1.0 ;
    (*np) ++ ;
    ph[(*np)*pstr] = ph[6*pstr]       ;
    th[(*np)*tstr] = th[6*pstr] + 1.0 ;
    (*np) ++ ;

    *ne = 36 ; *nt = 20 ;
    ed = edges_unwrap ;
    tri = tri_unwrap ;
  }

  for ( i = 1 ; i < (*ne) ; i ++ ) {
    edges[i*edstr+0] = ed[2*i+0] ;
    edges[i*edstr+1] = ed[2*i+1] ;
  }

  for ( i = 0 ; i < (*nt) ; i ++ ) {
    elements[i*elstr+0] = tri[i*3+0] ;
    elements[i*elstr+1] = tri[i*3+1] ;
    elements[i*elstr+2] = tri[i*3+2] ;
  }

  return 0 ;
}
