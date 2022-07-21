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

static const struct {
  gchar *shape ;
  gint np ;
  gboolean defaults ;
} shape_parse_list[] =
  {
   {"ellipse",    4, TRUE},
   {"naca",       4, FALSE },
   /* {"root",       11, FALSE}, */
   {NULL, 0, FALSE}
  } ;

#define M411  0
#define M412  1
#define M413  2
#define M414  3
#define M421  4
#define M422  5
#define M423  6
#define M424  7
#define M431  8
#define M432  9
#define M433 10
#define M434 11
#define M441 12
#define M442 13
#define M443 14
#define M444 15

static void invert4x4(gdouble *Ai, gdouble *A)

{
  gdouble det ;

  det = 
    A[M411]*A[M422]*A[M433]*A[M444] + A[M411]*A[M423]*A[M434]*A[M442] + 
    A[M411]*A[M424]*A[M432]*A[M443] +
    A[M412]*A[M421]*A[M434]*A[M443] + A[M412]*A[M423]*A[M431]*A[M444] + 
    A[M412]*A[M424]*A[M433]*A[M441] +
    A[M413]*A[M421]*A[M432]*A[M444] + A[M413]*A[M422]*A[M434]*A[M441] + 
    A[M413]*A[M424]*A[M431]*A[M442] +
    A[M414]*A[M421]*A[M433]*A[M442] + A[M414]*A[M422]*A[M431]*A[M443] + 
    A[M414]*A[M423]*A[M432]*A[M441] -
    A[M411]*A[M422]*A[M434]*A[M443] - A[M411]*A[M423]*A[M432]*A[M444] -
    A[M411]*A[M424]*A[M433]*A[M442] -
    A[M412]*A[M421]*A[M433]*A[M444] - A[M412]*A[M423]*A[M434]*A[M441] -
    A[M412]*A[M424]*A[M431]*A[M443] -
    A[M413]*A[M421]*A[M434]*A[M442] - A[M413]*A[M422]*A[M431]*A[M444] -
    A[M413]*A[M424]*A[M432]*A[M441] -
    A[M414]*A[M421]*A[M432]*A[M443] - A[M414]*A[M422]*A[M433]*A[M441] -
    A[M414]*A[M423]*A[M431]*A[M442] ;

  det = 1.0/det ;

  Ai[M411] = 
    A[M422]*A[M433]*A[M444] + A[M423]*A[M434]*A[M442] + 
    A[M424]*A[M432]*A[M443] -
    A[M422]*A[M434]*A[M443] - A[M423]*A[M432]*A[M444] - 
    A[M424]*A[M433]*A[M442] ;
  Ai[M412] = 
    A[M412]*A[M434]*A[M443] + A[M413]*A[M432]*A[M444] + 
    A[M414]*A[M433]*A[M442] -
    A[M412]*A[M433]*A[M444] - A[M413]*A[M434]*A[M442] - 
    A[M414]*A[M432]*A[M443] ;
  Ai[M413] = 
    A[M412]*A[M423]*A[M444] + A[M413]*A[M424]*A[M442] + 
    A[M414]*A[M422]*A[M443] -
    A[M412]*A[M424]*A[M443] - A[M413]*A[M422]*A[M444] - 
    A[M414]*A[M423]*A[M442] ;
  Ai[M414] = 
    A[M412]*A[M424]*A[M433] + A[M413]*A[M422]*A[M434] + 
    A[M414]*A[M423]*A[M432] -
    A[M412]*A[M423]*A[M434] - A[M413]*A[M424]*A[M432] - 
    A[M414]*A[M422]*A[M433] ;
  Ai[M421] = 
    A[M421]*A[M434]*A[M443] + A[M423]*A[M431]*A[M444] + 
    A[M424]*A[M433]*A[M441] -
    A[M421]*A[M433]*A[M444] - A[M423]*A[M434]*A[M441] - 
    A[M424]*A[M431]*A[M443] ;
  Ai[M422] = 
    A[M411]*A[M433]*A[M444] + A[M413]*A[M434]*A[M441] + 
    A[M414]*A[M431]*A[M443] -
    A[M411]*A[M434]*A[M443] - A[M413]*A[M431]*A[M444] - 
    A[M414]*A[M433]*A[M441] ;
  Ai[M423] = 
    A[M411]*A[M424]*A[M443] + A[M413]*A[M421]*A[M444] +
    A[M414]*A[M423]*A[M441] -
    A[M411]*A[M423]*A[M444] - A[M413]*A[M424]*A[M441] - 
    A[M414]*A[M421]*A[M443] ;
  Ai[M424] = 
    A[M411]*A[M423]*A[M434] + A[M413]*A[M424]*A[M431] + 
    A[M414]*A[M421]*A[M433] -
    A[M411]*A[M424]*A[M433] - A[M413]*A[M421]*A[M434] - 
    A[M414]*A[M423]*A[M431] ;
  Ai[M431] = 
    A[M421]*A[M432]*A[M444] + A[M422]*A[M434]*A[M441] + 
    A[M424]*A[M431]*A[M442] -
    A[M421]*A[M434]*A[M442] - A[M422]*A[M431]*A[M444] - 
    A[M424]*A[M432]*A[M441] ;
  Ai[M432] = 
    A[M411]*A[M434]*A[M442] + A[M412]*A[M431]*A[M444] +
    A[M414]*A[M432]*A[M441] -
    A[M411]*A[M432]*A[M444] - A[M412]*A[M434]*A[M441] - 
    A[M414]*A[M431]*A[M442] ;
  Ai[M433] = 
    A[M411]*A[M422]*A[M444] + A[M412]*A[M424]*A[M441] +
    A[M414]*A[M421]*A[M442] -
    A[M411]*A[M424]*A[M442] - A[M412]*A[M421]*A[M444] - 
    A[M414]*A[M422]*A[M441] ;
  Ai[M434] = 
    A[M411]*A[M424]*A[M432] + A[M412]*A[M421]*A[M434] +
    A[M414]*A[M422]*A[M431] -
    A[M411]*A[M422]*A[M434] - A[M412]*A[M424]*A[M431] - 
    A[M414]*A[M421]*A[M432] ;
  Ai[M441] = 
    A[M421]*A[M433]*A[M442] + A[M422]*A[M431]*A[M443] + 
    A[M423]*A[M432]*A[M441] -
    A[M421]*A[M432]*A[M443] - A[M422]*A[M433]*A[M441] - 
    A[M423]*A[M431]*A[M442] ;
  Ai[M442] = 
    A[M411]*A[M432]*A[M443] + A[M412]*A[M433]*A[M441] + 
    A[M413]*A[M431]*A[M442] -
    A[M411]*A[M433]*A[M442] - A[M412]*A[M431]*A[M443] - 
    A[M413]*A[M432]*A[M441] ;
  Ai[M443] = 
    A[M411]*A[M423]*A[M442] + A[M412]*A[M421]*A[M443] + 
    A[M413]*A[M422]*A[M441] -
    A[M411]*A[M422]*A[M443] - A[M412]*A[M423]*A[M441] - 
    A[M413]*A[M421]*A[M442] ;
  Ai[M444] = 
    A[M411]*A[M422]*A[M433] + A[M412]*A[M423]*A[M431] + 
    A[M413]*A[M421]*A[M432] -
    A[M411]*A[M423]*A[M432] - A[M412]*A[M421]*A[M433] - 
    A[M413]*A[M422]*A[M431] ;

  Ai[M411] *= det ; Ai[M412] *= det ; Ai[M413] *= det ; Ai[M414] *= det ;
  Ai[M421] *= det ; Ai[M422] *= det ; Ai[M423] *= det ; Ai[M424] *= det ;
  Ai[M431] *= det ; Ai[M432] *= det ; Ai[M433] *= det ; Ai[M434] *= det ;
  Ai[M441] *= det ; Ai[M442] *= det ; Ai[M443] *= det ; Ai[M444] *= det ;

  return ;
}

gdouble agg_shape_eval(agg_shape_t *s, gdouble x, gint lim)

{
  gdouble y ;

  switch ( s->type ) {
  default: g_assert_not_reached() ; break ;
  case AGG_SHAPE_AEROFOIL: y = agg_shape_aerofoil_eval(s, x) ; break ;
  case AGG_SHAPE_ELLIPSE:  y = agg_shape_ellipse_eval(s, x, lim) ; break ;
  }
  
  return y ;
}

gdouble agg_shape_aerofoil_eval(agg_shape_t *s, gdouble x)

{
  gdouble y, sgn = 1.0, sc, S[32] ;
  gint i, ns ;

  g_assert(s->nb == 0) ;
  if ( x < 0.0 ) { x = -x ; sgn = -1 ; }
  /*evaluate thickness*/
  sc = pow(x, s->n1[0])*pow(1.0-x, s->n2[0]) ;

  ns = s->i[1] - s->i[0] - 1 ;
  /* ns = MAX(ns, s->i[3] - s->i[2]) ; */

  /*trailing edge thickness*/
  y = sgn*x*s->s[s->i[1]] ;
  agg_bernstein_basis(ns, x, S, NULL) ;

  for ( i = s->i[0] ; i < s->i[1] ; i ++ ) y += sgn*sc*s->s[i]*S[i-s->i[0]] ;

  /*camber*/
  for ( i = s->i[2] ; i < s->i[3] ; i ++ ) y += sc*s->s[i]*S[i-s->i[2]] ;  
  
  return y ;
}

gdouble agg_shape_ellipse_eval(agg_shape_t *s, gdouble x, gint b)


/*
 * evaluate shape s at |x| on break b
 *
 * if b < 0, locate interval for x
 */
  
{
  gdouble y, sgn = 1.0, sc, S[32] ;
  gint i, ns, i0, i1 ;

  if ( b < 0 ) {
    /*locate breakpoint*/
    for ( b = 0 ;
	  b < agg_shape_break_number(s) &&
	    !(agg_shape_break_lower(s,b) <= x &&
	      agg_shape_break_upper(s,b) >= x ) ;
	  b ++ ) ;
  }

  g_assert(b < agg_shape_break_number(s)) ;
  i0 = s->i[b] ; i1 = s->i[b+1] ;
  if ( x < 0.0 ) { x = -x ; }
  
  sc = pow(x, s->n1[b])*pow(1.0-x, s->n2[b]) ;

  ns = i1 - i0 - 1 ;

  agg_bernstein_basis(ns, x, S, NULL) ;

  y = 0.0 ;
  for ( i = i0 ; i < i1 ; i ++ ) y += sgn*sc*s->s[i]*S[i-i0] ;

  return y ;
}

agg_shape_t *agg_shape_alloc(gint nsmax)

{
  agg_shape_t *s ;

  s = (agg_shape_t *)g_malloc0(sizeof(agg_shape_t)) ;

  s->nsmax = 2*nsmax+4 ;  
  s->s = (gdouble *)g_malloc0((s->nsmax)*sizeof(gdouble)) ;
  s->b[0] = s->b[1] = -1 ;
  
  return s ;
}

gint agg_shape_naca(agg_shape_t *s, gdouble t)

{
  s->type = AGG_SHAPE_AEROFOIL ;
  s->closed = TRUE ;

  /*aerofoil: i[0] start of thickness distribution
   *          i[1] trailing edge thickness term
   *          i[2] start of camber
   *          i[3] end of camber
   */

  s->i[0] = 0 ; 
  s->i[1] = 1 ; 
  s->i[2] = 2 ; 
  s->i[3] = 2 ; 

  s->s[0] = t ;
  s->s[1] = 0.002 ;
  
  s->n1[0] = 0.5 ; s->n2[0] = 1.0 ;

  s->nb = 1 ;
  s->b[0] = -1.0 ; s->b[1] = 1.0 ; 

  return 0 ;
}

gint agg_shape_ellipse(agg_shape_t *s)

{
  s->type = AGG_SHAPE_ELLIPSE ;
  s->closed = TRUE ;

  s->n1[0] = 0.5 ; s->n2[0] = 0.5 ;  
  s->n1[1] = 0.5 ; s->n2[1] = 0.5 ;  

  /*ellipse: i[0] start of upper surface
   *         i[1] start of lower surface
   *         i[2] end of lower surface
   */

  /*one breakpoint, at 0*/
  s->nb = 2 ;
  s->b[0] = -1.0 ; s->b[1] =  0.0 ; 
  s->b[2] =  0.0 ; s->b[3] =  1.0 ; 
  
  s->i[0] = 0 ;
  s->i[1] = 1 ;
  s->i[2] = 2 ;

  s->s[0] = -1.0 ;
  s->s[1] =  1.0 ;
  
  return 0 ;
}

gint agg_shape_fit(agg_shape_t *s,
		   gdouble *xu, gdouble *yu, gint nu,
		   gdouble *xl, gdouble *yl, gint nl,
		   gdouble n1, gdouble n2, gdouble d,
		   gint n, gboolean closed, gdouble *work)

{
  gint ns, i, j, i1 = 1, lwork, info, lda, wsize ;
  gdouble *A, *b, S[32] ;

  ns = nu + nl ; lda = 2*n + 2 ;

  g_assert(2*n + 1 < s->nsmax) ;

  /*workspace size*/
  wsize = 2*n+2 + nu + nl + (2*n+2)*(nu+nl) + nu + nl ;
  memset(work, 0, wsize*sizeof(gdouble)) ;
  
  lwork = lda + ns ;
  A = &(work[lwork]) ;
  b = &(A[lda*ns]) ;
  
  for ( i = 0 ; i < nu ; i ++ ) {
    b[i] = yu[i] - xu[i]*d ;
    agg_bernstein_basis(n, xu[i], S, NULL) ;

    for ( j = 0 ; j <= n ; j ++ ) {
      A[i*lda + j] = S[j]*pow(xu[i],n1)*pow(1.0-xu[i],n2) ;
    }
  }
  for ( i = 0 ; i < nl ; i ++ ) {
    b[nu+i] = yl[i] - (-xl[i]*d) ;
    agg_bernstein_basis(n, xl[i], S, NULL) ;

    for ( j = 0 ; j <= n ; j ++ ) {
      A[(nu+i)*lda + n+1 + j] = S[j]*pow(xl[i],n1)*pow(1.0-xl[i],n2) ;
    }
  }

  /*least-squares solve using LAPACK*/
  dgels_("T", &lda, &ns, &i1, A, &lda, b, &ns, work, &lwork, &info) ;

  s->i[0] = 0 ;
  for ( i = 0 ; i <= n ; i ++ ) s->s[i] = 0.5*(b[i] - b[n+1+i]);
  s->i[1] = s->i[0] + n+1 ;
  s->s[s->i[1]] = d ;
  
  s->i[2] = s->i[1]+1 ;  
  for ( i = 0 ; i <= n ; i ++ ) s->s[s->i[2]+i] = 0.5*(b[i] + b[n+1+i]);
  s->i[3] = s->i[2] + n + 1 ;
  
  s->n1[0] = n1 ; s->n2[0] = n2 ; 

  s->closed = closed ;
  
  return 0 ;
}

static void root_fairing(gdouble n1, gdouble n2,
			 gdouble h, gdouble s, gdouble d, gdouble b,
			 gdouble *x, gdouble *y, gint np)

{
  gdouble A[16], Ai[16], q[4], r[4] ;
  gint i ;

  g_assert(s > 0.0 && s < 0.5) ;
  
  /*generate surface centred on 0 (not 1/2)*/
  A[0] = 1.0 ; A[1] =  s ; A[2] =  s*s ; A[3] =  s*s*s ;
  A[4] = 0.0 ; A[5] = 1.0 ; A[6] = 2.0*s ; A[7] = 3.0*s*s ;
  /*shoulder point*/
  A[8] = 1.0 ; A[9] = 0.5+d ; A[10] =  A[9]*(0.5+d) ;
  A[11] = A[10]*(0.5+d) ;
  
  A[12] = 0.0 ; A[13] = 1.0 ; A[14] = 2.0*A[9] ; A[15] = 3.0*A[10] ;

  invert4x4(Ai, A) ;
  /*note s is currently centred on zero, so shift to match Kulfan
    coordinates*/
  r[0] = pow(0.5+s,n1)*pow(0.5-s,n2)*h ;
  r[1] = r[0]/(0.5+s)/(0.5-s)*(n1 - (n1+n2)*(0.5+s)) ;
  r[2] = b ;
  r[3] = 0.0 ;

  q[0] = Ai[0]*r[0] + Ai[1]*r[1] + Ai[2]*r[2] + Ai[3]*r[3] ;
  q[1] = Ai[4]*r[0] + Ai[5]*r[1] + Ai[6]*r[2] + Ai[7]*r[3] ;
  q[2] = Ai[8]*r[0] + Ai[9]*r[1] + Ai[10]*r[2] + Ai[11]*r[3] ;
  q[3] = Ai[12]*r[0] + Ai[13]*r[1] + Ai[14]*r[2] + Ai[15]*r[3] ;

  /*build the curve symmetrically about 0.0*/
  i = 0 ;
  x[i] = 0.0 ; y[i] = pow(0.5+x[i], n1)*pow(0.5-x[i], n2)*h ;
  for ( i = 1 ; i < np ; i ++ ) {
    x[i] = agg_spacing_eval(0.0, 0.5+d, np, AGG_SPACING_LINEAR, i) ;
    if ( x[i] < s ) {
      y[i] = pow(0.5+x[i], n1)*pow(0.5-x[i], n2)*h ;
    } else {
      y[i] = q[0] + x[i]*(q[1] + x[i]*(q[2] + x[i]*q[3])) ;
    }
    x[i+np-1] = -x[i] ; y[i+np-1] = y[i] ;
  }

  /*remap to (0,1)*/
  for ( i = 0 ; i < 2*np-1 ; i ++ ) {
    x[i] = (x[i] + 0.5 + d)/(1+2*d) ;
    y[i] /= 1.0+2*d ;
    /* fprintf(stdout, "%lg %lg\n", x[i], y[i]) ; */
  }
  
  return ;
}

gint agg_shape_root_fit(agg_shape_t *s,
			gdouble n1, gdouble n2,
			gdouble hu, gdouble su, gdouble du, gdouble bu,
			gdouble hl, gdouble sl, gdouble dl, gdouble bl,
			gint n, gdouble *work)

/*
 * construct a "wing root" cross section 
 * 
 * (u,l) (upper, lower)
 *
 * h(u,l): height of main curve
 * s(u,l): breakaway point from curve
 * d(u,l): displacement of "shoulder"
 * b(u,l): height of "shoulder"
 * 
 */
  
{
  gdouble x[129], y[129] ;
  gint np ;

  s->nb = 0 ; s->b[0] = -1 ; s->i[0] = 0 ;
  s->type = AGG_SHAPE_ELLIPSE ;
  s->closed = FALSE ;

  g_assert(du == dl) ; /*transformation not checked yet for du =/= dl*/
  
  np = 64 ;
  root_fairing(n1, n2, hl, sl, dl, bl, x, y, np) ;
  agg_shape_interval_add(s, -1, 0, x, y, 2*np-1, FALSE, 0, 0, n, work) ;

  root_fairing(n1, n2, hu, su, du, bu, x, y, np) ;
  agg_shape_interval_add(s, 0, 1, x, y, 2*np-1, FALSE, 0, 0, n, work) ;
  
  return 0 ;
}

gint agg_shape_interval_add(agg_shape_t *s,
			    gdouble x1, gdouble x2,
			    gdouble *x, gdouble *y, gint np,
			    gboolean upper, 
			    gdouble n1, gdouble n2,
			    gint n, gdouble *work)

/*
 * add a piecewise section to a shape, covering the interval (x1,x2),
 * fitted to the np points in (x,y)
 */
  
{
  gint i, j, i0, i1 = 1, lwork, info, lda, wsize ;
  gdouble *A, *b, S[32] ;

  g_assert(s->nb < 6) ;
  g_assert(x1 <= 1) ;
  g_assert(x2 <= 1) ;
  g_assert(x2 >= x1) ;
  /*find the last breakpoint*/
  i = s->nb ;
  s->n1[i] = n1 ; s->n2[i] = n2 ; 
  g_assert(s->b[2*i+1] <= x1) ;

  /*locate start of coefficients for the new interval*/
  i0 = s->i[i] ;

  /*set up least squares fit*/
  lda = n + 1 ;
  wsize = n+1 + np + (n+1)*(np) + np ;
  memset(work, 0, wsize*sizeof(gdouble)) ;
  lwork = lda + np ;
  A = &(work[lwork]) ;
  b = &(A[lda*np]) ;
  
  for ( i = 0 ; i < np ; i ++ ) {
    b[i] = y[i] ;
    agg_bernstein_basis(n, x[i], S, NULL) ;

    for ( j = 0 ; j <= n ; j ++ ) {
      A[i*lda + j] = S[j]*pow(x[i],n1)*pow(1.0-x[i],n2) ;
    }
  }

  /*least-squares solve using LAPACK*/
  dgels_("T", &lda, &np, &i1, A, &lda, b, &np, work, &lwork, &info) ;

  /*update the bookkeeping information*/
  for ( i = 0 ; i <= n ; i ++ ) {
    s->s[i0+i] = b[i] ;
  }
  
  i = s->nb ;
  s->i[i+1] = i0+n+1 ;
  s->b[2*i+0] = x1 ;
  s->b[2*i+1] = x2 ;
  
  s->nb ++ ;

  return 0 ;
}

gint agg_shape_parse(agg_shape_t *s, gchar *type, gdouble *p, gint np)

/*
 * type: string identifier
 * p:    parameters (all doubles)
 * np:   number of parameters (checked for transform)
 */
  
{
  gint i ;
  gdouble work[4096] ;
  
  for ( i = 0 ; shape_parse_list[i].shape != NULL ; i ++ ) {
    if ( !strcmp(shape_parse_list[i].shape, type) ) break ;
  }
       
  if ( shape_parse_list[i].shape == NULL ) return -1 ;

  /*check number of parameters*/
  if ( np != shape_parse_list[i].np && !(shape_parse_list[i].defaults) ) {
    g_error("%s: shape \"%s\" should have %d parameters, not %d",
	    __FUNCTION__, type, shape_parse_list[i].np, np) ;
  }

  if ( i == 0 ) {
    /*ellipse*/
    agg_shape_ellipse(s) ;

    /*modify settings if defaults are overridden*/
    if ( np > 0 ) s->s[0]  = p[0] ; 
    if ( np > 1 ) s->s[1]  = p[1] ; 
    if ( np > 2 ) s->n1[0] = p[2] ; 
    if ( np > 3 ) s->n2[0] = p[3] ; 
    
    return 0 ;
  }

  if ( i == 1 ) {
    /*naca four*/
    agg_shape_fit_naca_four(s, (gint)p[0], p[1], p[2], p[3], 64, 64, work) ;

    return 0 ;
  }
  
  return 0 ;
}

gint agg_shape_copy(agg_shape_t *dest, agg_shape_t *src)

{
  gint i ;
  
  dest->type = src->type ;
  dest->closed = src->closed ;

  memcpy(dest->n1, src->n1, 8*sizeof(gdouble)) ;
  memcpy(dest->n2, src->n2, 8*sizeof(gdouble)) ;
  /* dest->n1 = src->n1 ; */
  /* dest->n2 = src->n2 ; */

  dest->nb = src->nb ;
  for ( i = 0 ; i < src->nb+2 ; i ++ ) dest->b[i] = src->b[i] ;
  memcpy(dest->i, src->i, 8*sizeof(gint)) ;
  
  g_assert(dest->nsmax >= src->nsmax) ;
  memcpy(dest->s, src->s, src->nsmax*sizeof(gdouble)) ;
  
  return 0 ;
}
