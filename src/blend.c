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

agg_surface_blend_t *agg_surface_blend_new(void)

{
  agg_surface_blend_t *B ;

  B = (agg_surface_blend_t *)g_malloc0(sizeof(agg_surface_blend_t)) ;
  
  return B ;
}

gint agg_hermite_eval(gdouble s, gdouble *H)

{
  H[0] = s*s*(2*s-3) + 1 ;
  H[1] = 1.0 - H[0] ;
  H[2] = s*(s-1)*(s-1) ;
  H[3] = s*s*(s-1) ;
  
  return 0 ;
}

gint agg_surface_blend_evaluate(agg_surface_blend_t *B,
				gdouble s, gdouble t, gdouble *x,
				agg_surface_workspace_t *w)

/*
 * Blended surface method of Blending parametric surfaces Daniel
 * J. Filip, ACM Transactions on Graphics, Volume 8, Issue 301 July
 * 1989, pp164-173 https://doi.org/10.1145/77055.77057
 */
  
{
  agg_surface_t *S1, *S2 ;
  agg_patch_t *P1, *P2 ;
  agg_patch_clipping_t *c1, *c2 ;
  gdouble xu[3], xv[3], u, v, Z1[3], Z2[3], K[3], C1[3], C2[3] ;
  gdouble V1[3], V2[3], N1[3], N2[3], T1[3], T2[3], g, l, wt ;
  gdouble S1s[3], S1t[3], S2s[3], S2t[3] ;
  gdouble s1, t1, s2, t2, tmp1, tmp2, H[4] ;
  gint i ;
  gboolean use1, use2 ;
  
  g_assert(B->S[0] != NULL) ; g_assert(B->S[1] != NULL) ;
  g_assert(B->P[0] != NULL) ; g_assert(B->P[1] != NULL) ;

  wt = 1.0 ;
  
  S1 = agg_surface_blend_surface(B,0) ;
  S2 = agg_surface_blend_surface(B,1) ;
  P1 = agg_surface_blend_patch(B,0) ;
  P2 = agg_surface_blend_patch(B,1) ;
  i = B->ic[0] ;
  g_assert(i < agg_patch_clipping_number(P1)) ;
  c1 = agg_patch_clipping(P1, i) ;
  i = B->ic[1] ;
  g_assert(i < agg_patch_clipping_number(P2)) ;
  c2 = agg_patch_clipping(P2, i) ;

  /*coordinates on surfaces to be blended*/
  agg_patch_clip_eval(c1, t, &s1, &t1) ;
  agg_patch_clip_eval(c2, t, &s2, &t2) ;

  use1 = use2 = TRUE ;
  
  agg_patch_map(P1, s1, t1, &u, &v) ;
  agg_surface_point_diff(S1, u, v, C1, xu, xv, w) ;
  agg_patch_surface_diff(P1, s1, t1, xu, xv, S1s, S1t) ;
  if ( agg_vector_length(S1s) == 0 || agg_vector_length(S1t) == 0 )
    use1 = FALSE ;
  agg_vector_cross(Z1, S1s, S1t) ;

  agg_patch_map(P2, s2, t2, &u, &v) ;
  agg_surface_point_diff(S2, u, v, C2, xu, xv, w) ;
  agg_patch_surface_diff(P2, s2, t2, xu, xv, S2s, S2t) ;
  if ( agg_vector_length(S2s) == 0 || agg_vector_length(S2t) == 0 )
    use2 = FALSE ;
  agg_vector_cross(Z2, S2s, S2t) ;

  agg_vector_diff(K, C2, C1) ;

  tmp1 = agg_vector_scalar(Z1,Z1) ;
  tmp2 = agg_vector_scalar(K,Z1) ;
  V1[0] = tmp1*K[0] - tmp2*Z1[0] ;
  V1[1] = tmp1*K[1] - tmp2*Z1[1] ;
  V1[2] = tmp1*K[2] - tmp2*Z1[2] ;
  tmp1 = agg_vector_scalar(Z2,Z2) ;
  tmp2 = agg_vector_scalar(K,Z2) ;
  V2[0] = tmp1*K[0] - tmp2*Z2[0] ;
  V2[1] = tmp1*K[1] - tmp2*Z2[1] ;
  V2[2] = tmp1*K[2] - tmp2*Z2[2] ;

  tmp1 = agg_vector_length(V1) ;
  tmp2 = agg_vector_length(V2) ;
  N1[0] = V1[0]/tmp1 ; N1[1] = V1[1]/tmp1 ; N1[2] = V1[2]/tmp1 ; 
  N2[0] = V2[0]/tmp2 ; N2[1] = V2[1]/tmp2 ; N2[2] = V2[2]/tmp2 ; 

  /*check for valid tangent vectors*/
  if ( use1 ) {
    g = agg_vector_length(K) + agg_vector_scalar(N1, K) ;
    l = 2.0*agg_vector_scalar(K,K)/g ;
  } else {
    g = agg_vector_length(K) + agg_vector_scalar(N2, K) ;
    l = 2.0*agg_vector_scalar(K,K)/g ;
  }
  
  T1[0] = wt*l*N1[0] ; T1[1] = wt*l*N1[1] ; T1[2] = wt*l*N1[2] ; 
  T2[0] = wt*l*N2[0] ; T2[1] = wt*l*N2[1] ; T2[2] = wt*l*N2[2] ; 
  
  agg_hermite_eval(s, H) ;

  x[0] = H[0]*C1[0] + H[1]*C2[0] ;
  x[1] = H[0]*C1[1] + H[1]*C2[1] ;
  x[2] = H[0]*C1[2] + H[1]*C2[2] ;

  if ( use1 ) {
    x[0] += H[2]*T1[0] ;
    x[1] += H[2]*T1[1] ;
    x[2] += H[2]*T1[2] ;
  }

  if ( use2 ) {
    x[0] += H[3]*T2[0] ;
    x[1] += H[3]*T2[1] ;
    x[2] += H[3]*T2[2] ;
  }
  
  return 0 ;
}