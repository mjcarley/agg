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

agg_chebyshev_t *agg_chebyshev_new(gint Nmax, gint nc)

{
  agg_chebyshev_t *C ;
  
  C = (agg_chebyshev_t *)g_malloc0(sizeof(agg_chebyshev_t)) ;

  C->f = (gdouble *)g_malloc0((Nmax+1)*nc*sizeof(gdouble)) ;

  C->Nmax = (Nmax+1)*nc ;
  
  return C ;
}

gint agg_chebyshev_eval(agg_chebyshev_t *C, gdouble x, gdouble *f)

/*
 * Battles and Trefethen, An extension of Matlab to continuous
 * functions and operators, https://dx.doi.org/10.1137/S1064827503430126
 * 
 */
  
{
  gint i, j, N, inter, nc ;
  gdouble t, ti, wsum, m1pi, xmin, xmax, *F, c, s, dc, ds ;

  for ( j = 0 ; j < agg_chebyshev_component_number(C) ; j ++ ) {
    f[j] = 0.0 ;
  }
  if ( x < C->xinter[0] ) return 0 ;
  if ( x > C->xinter[C->ninter] ) return 0 ;
  
  /*locate interval*/
  for ( inter = 0 ; inter < agg_chebyshev_interval_number(C) ; inter ++ ) {
    if ( agg_chebyshev_interval_xmin(C,inter) <= x &&
	 agg_chebyshev_interval_xmin(C,inter+1) >= x ) break ;
  }

  nc = agg_chebyshev_component_number(C) ;
  xmin = agg_chebyshev_interval_xmin(C,inter  ) ;
  xmax = agg_chebyshev_interval_xmin(C,inter+1) ;
  F = &(C->f[agg_chebyshev_interval_start(C,inter)*nc]) ;
  N = agg_chebyshev_interval_order(C,inter) ;

  if ( x == xmax ) {
    for ( j = 0 ; j < nc ; j ++ ) {
      f[j] = F[0*nc+j] ;
      f[j] = f[j]*sqrt(fabs(x)) + C->x0[j] ;
    }
    return 0 ;
  }
  if ( x == xmin ) {
    for ( j = 0 ; j < agg_chebyshev_component_number(C) ; j ++ ) {
      f[j] = F[N*nc+j] ;
      f[j] = f[j]*sqrt(fabs(x)) + C->x0[j] ;
    }
    return 0 ;
  }
  
  t = -1.0 + 2.0*(x - xmin)/(xmax - xmin) ;

  dc = cos(M_PI/N) ; ds = sin(M_PI/N) ; 
  c = 1.0 ; s = 0.0 ;
  
  m1pi = 1 ; i = 0 ;
  ti = c ;
  for ( j = 0 ; j < nc ; j ++ ) f[j] = 0.5*m1pi*F[i*nc+j]/(t-ti) ;
  wsum = 0.5*m1pi/(t-ti) ; m1pi = -m1pi ;

  for ( i = 1 ; i < N ; i ++ ) {
    ti = c ; c = dc*ti - ds*s ; s = dc*s  + ds*ti ; ti = c ;
    for ( j = 0 ; j < nc ; j ++ ) f[j] += m1pi*F[i*nc+j]/(t-ti) ;
    wsum += m1pi/(t-ti) ; m1pi = -m1pi ;
  }

  i = N ; 
  ti = c ; c = dc*ti - ds*s ; s = dc*s  + ds*ti ; ti = c ;
  for ( j = 0 ; j < nc ; j ++ ) f[j] += 0.5*m1pi*F[i*nc+j]/(t-ti) ;
  wsum += 0.5*m1pi/(t-ti) ;

  for ( j = 0 ; j < nc ; j ++ ) {
    f[j] /= wsum ;
    f[j] = f[j]*sqrt(fabs(x)) + C->x0[j] ;
  }
  
  return 0 ;
}

static gint interval_add_points(agg_chebyshev_t *C,
				agg_surface_t *S, gdouble u,
				gint N, gint inter,
				agg_surface_workspace_t *w)

{
  gint nc, i, j, ns, s[8] ;
  gdouble v, vmax, vmin, *F, *x0, vs[8], vc, x[3], wj, wc ;

  nc = agg_chebyshev_component_number(C) ;

  x0 = C->x0 ;
  agg_surface_point_eval(S, u, 0, x0, w) ;

  F = &(C->f[agg_chebyshev_interval_start(C,inter)*nc]) ;
  vmin = agg_chebyshev_interval_xmin(C,inter  ) ;
  vmax = agg_chebyshev_interval_xmin(C,inter+1) ;
  ns = 0 ;
  for ( i = 0 ; i <= N ; i ++ ) {
    v = cos(M_PI*i/N) ;
    v = 0.5*(vmax - vmin)*(1.0 + v) + vmin ;
    agg_surface_point_eval(S, u, v, &(F[3*i]), w) ;
    F[3*i+0] -= x0[0] ; 
    F[3*i+1] -= x0[1] ; 
    F[3*i+2] -= x0[2] ;
    if ( v != 0 ) {
      F[3*i+0] /= sqrt(fabs(v)) ;
      F[3*i+1] /= sqrt(fabs(v)) ;
      F[3*i+2] /= sqrt(fabs(v)) ;
    } else {
      /*solve for the value of the scaled function at the singularity*/
      s[ns] = i ; vs[ns] = v ; ns ++ ;
      vs[ns] = cos(M_PI*i/N) ;
      F[3*i+0] = F[3*i+1] = F[3*i+2] = 0 ;
    }
  }

  if ( ns == 0 ) return N + 1 ;

  g_assert(ns == 1) ;

  /*reference point on curve for setting of free variable*/
  i = s[0] ;
  vc = cos((N/2+0.5)*M_PI/N) ;
  vc = 0.5*(vmax - vmin)*(1.0 + vc) + vmin ;
  agg_surface_point_eval(S, u, vc, x, w) ;
  x[0] -= x0[0] ; x[1] -= x0[1] ; x[2] -= x0[2] ; 
  x[0] /= sqrt(fabs(vc)) ; x[1] /= sqrt(fabs(vc)) ; x[2] /= sqrt(fabs(vc)) ;
  vc = cos((N/2+0.5)*M_PI/N) ;
  /* vc = cos((i+0.5)*M_PI/N) ; */

  wj = 0.5 ; j = 0 ;
  v = cos(M_PI*j/N) ;
  if ( j != i ) {
    F[3*i+0] += wj/(vc - v)*(x[0] - F[3*j+0]) ;
    F[3*i+1] += wj/(vc - v)*(x[1] - F[3*j+1]) ;
    F[3*i+2] += wj/(vc - v)*(x[2] - F[3*j+2]) ;
  } else {
    wc = wj ;
  }

  wj = -1 ;
  for ( j = 1 ; j < N ; j ++ ) {
    v = cos(M_PI*j/N) ;
    if ( j != i ) {
      F[3*i+0] += wj/(vc - v)*(x[0] - F[3*j+0]) ;
      F[3*i+1] += wj/(vc - v)*(x[1] - F[3*j+1]) ;
      F[3*i+2] += wj/(vc - v)*(x[2] - F[3*j+2]) ;
    } else {
      wc = wj ;
    }
    wj *= -1 ;
  }

  j = N ; wj *= 0.5 ;
  
  v = cos(M_PI*j/N) ;
  if ( j != i ) {
    F[3*i+0] += wj/(vc - v)*(x[0] - F[3*j+0]) ;
    F[3*i+1] += wj/(vc - v)*(x[1] - F[3*j+1]) ;
    F[3*i+2] += wj/(vc - v)*(x[2] - F[3*j+2]) ;
  } else {
    wc = wj ;
  }

  F[3*i+0] *= (vc - vs[0])/wc ;
  F[3*i+1] *= (vc - vs[0])/wc ;
  F[3*i+2] *= (vc - vs[0])/wc ;

  F[3*i+0] += x[0] ; F[3*i+1] += x[1] ; F[3*i+2] += x[2] ; 
  
  return N + 1 ;
}

gint agg_chebyshev_surface_section(agg_chebyshev_t *C,
				   agg_surface_t *S, gdouble u,
				   gint N, agg_surface_workspace_t *w)

{
  gint i, inter, nc ;

  nc = 3 ;
  agg_chebyshev_component_number(C) = nc ;
  C->ninter = 8 ;

  C->inter[0] = 0 ; C->xinter[0] = -1 ;
  for ( i = 1 ; i <= C->ninter ; i ++ ) {
    C->xinter[i] = -1 + 2.0*i/(C->ninter) ;
    C->inter[i]  = C->inter[i-1] + N + 1 ;
  }
  
  g_assert(C->ninter*(N+1)*3 <= C->Nmax) ;

  for ( inter = 0 ; inter < C->ninter ; inter ++ ) {
    C->inter[i+1] =
      C->inter[i] + interval_add_points(C, S, u, N, inter, w) ;
  }

  return 0 ;
}  

static gdouble interval_error(agg_chebyshev_t *C, gint inter,
			      agg_surface_t *S, gdouble u,
			      gint nv,
			      agg_surface_workspace_t *w)

{
  gdouble emax, vmin, vmax, v, x[3], y[3] ;
  gint i ;

  emax = 0.0 ;

  vmin = agg_chebyshev_interval_xmin(C, inter  ) ;
  vmax = agg_chebyshev_interval_xmin(C, inter+1) ;

  for ( i = 0 ; i <= nv ; i ++ ) {
    v = vmin + (vmax - vmin)*i/nv ;
    agg_surface_point_eval(S, u, v, x, w) ;    
    agg_chebyshev_eval(C, v, y) ;
    emax = MAX(emax, fabs(x[0] - y[0])) ;
    emax = MAX(emax, fabs(x[1] - y[1])) ;
    emax = MAX(emax, fabs(x[2] - y[2])) ;
  }
  
  return emax ;
}
  
gint agg_chebyshev_surface_section_adaptive(agg_chebyshev_t *C,
					    agg_surface_t *S, gdouble u,
					    gint N, gdouble tol, gdouble dmin,
					    agg_surface_workspace_t *w)

{
  gint inter, nc ;
  gdouble emax ;
  
  nc = 3 ;
  agg_chebyshev_component_number(C) = nc ;

  C->ninter = 0 ;
  C->inter[0] = 0 ; C->xinter[0] = -1 ;

  do { 
    inter = C->ninter ;
    agg_chebyshev_interval_number(C) ++ ;
    if ( agg_chebyshev_interval_number(C) >=
	 AGG_CHEBYSHEV_INTERVAL_NUMBER_MAX ) {
      g_error("%s: maximum number of Chebyshev intervals (%d) exceeded",
	      __FUNCTION__, AGG_CHEBYSHEV_INTERVAL_NUMBER_MAX) ;
    }
    C->xinter[inter+1] = 1 ;
    C->inter[inter+1] = C->inter[inter] + N + 1 ;
    
    while ( C->xinter[inter+1] - C->xinter[inter] > dmin ) {
      interval_add_points(C, S, u, N, inter, w) ;
      emax = interval_error(C, inter, S, u, 16, w) ;
      if ( emax <= tol ) break ;
      C->xinter[inter+1] = C->xinter[inter] +
	0.5*(C->xinter[inter+1] - C->xinter[inter]) ;
    }
  } while ( C->xinter[inter+1] != 1 ) ;

  /*check the last interval was filled*/
  inter = agg_chebyshev_interval_number(C) - 1 ;
  if ( C->xinter[inter+1] - C->xinter[inter] < dmin ) {
    interval_add_points(C, S, u, N, inter, w) ;
  }
  
  return 0 ;
}  

gdouble agg_chebyshev_interval_shortest(agg_chebyshev_t *C)

{
  gdouble min ;
  gint i ;

  min = G_MAXDOUBLE ;
  for ( i = 0 ; i < agg_chebyshev_interval_number(C) ; i ++ ) {
    min = MIN(min, agg_chebyshev_interval_xmin(C,i+1) -
	      agg_chebyshev_interval_xmin(C,i)) ;	      
  }
  
  return min ;
}
