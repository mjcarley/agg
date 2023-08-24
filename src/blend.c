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

/**
 * @{
 *  @ingroup blends
 */

/** 
 * Allocate a new blended surface
 * 
 * 
 * @return newly allocated ::agg_surface_blend_t
 */

agg_surface_blend_t *agg_surface_blend_new(void)

{
  agg_surface_blend_t *B ;

  B = (agg_surface_blend_t *)g_malloc0(sizeof(agg_surface_blend_t)) ;
  
  return B ;
}

/** 
 * Evaluate Hermite smoothing polynomials
 * 
 * @param s argument of polynomials;
 * @param H on exit contains \f$H_{i}(s)\f$, \f$1\leq i\leq4\f$.
 * 
 * @return 0 on success.
 */

gint agg_hermite_eval(gdouble s, gdouble *H)

{
  H[0] = s*s*(2*s-3) + 1 ;
  H[1] = 1.0 - H[0] ;
  H[2] = s*(s-1)*(s-1) ;
  H[3] = s*s*(s-1) ;
  
  return 0 ;
}

/** 
 * Evaluate a surface blend, using method of Blending parametric
 * surfaces Daniel J. Filip, ACM Transactions on Graphics, Volume 8,
 * Issue 301, July 1989, pp164-173 https://doi.org/10.1145/77055.77057
 * 
 * @param B surface blend, containing details of rail curves and surface
 * to join;
 * @param s parameter \f$s\f$, \f$0\leq s\leq 1\f$
 * @param t parameter \f$t\f$, \f$0\leq t\leq 1\f$
 * @param x on exit, contains surface position at \f$(s,t)\f$;
 * @param w workspace for surface evaluation.
 * 
 * @return 
 */

gint agg_surface_blend_evaluate(agg_surface_blend_t *B,
				gdouble s, gdouble t, gdouble *x,
				agg_surface_workspace_t *w)
  
{
  agg_surface_t *S1, *S2 ;
  agg_patch_t *P1, *P2 ;
  agg_curve_t *c1, *c2 ;
  gdouble Z1[3], Z2[3], K[3], C1[3], C2[3] ;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
  gdouble V1[3], V2[3], N1[3], N2[3], T1[3], T2[3] ;
#pragma GCC diagnostic pop
  gdouble g1, l1, w1, g2, l2, w2 ;
  gdouble S1s[3], S1t[3], S2s[3], S2t[3] ;
  gdouble s1, t1, s2, t2, tmp1, tmp2, H[4] ;
  gint h1, h2 ;
  gboolean use1, use2 ;
  
  w1 = 1.0 ; w2 = 1.0 ;
  
  S1 = agg_surface_blend_surface(B,0) ; g_assert(S1 != NULL) ;
  S2 = agg_surface_blend_surface(B,1) ; g_assert(S2 != NULL) ;
  P1 = agg_surface_blend_patch(B,0) ;   g_assert(P1 != NULL) ;
  P2 = agg_surface_blend_patch(B,1) ;   g_assert(P2 != NULL) ;
  h1 = agg_surface_blend_hole(B,0) ;
  c1 = agg_patch_hole(P1,h1) ;
  h2 = agg_surface_blend_hole(B,1) ;
  c2 = agg_patch_hole(P2,h2) ;  
  /*coordinates on surfaces to be blended*/
  if ( h1 > 1 ) {
    if ( agg_surface_blend_curve_reverse(B,0) )
      agg_curve_eval(c1, 1.0-t, &s1, &t1) ;
    else
      agg_curve_eval(c1,     t, &s1, &t1) ;
  } else {
    g_assert(h1 == 0) ;
    /*cut-off ends are mapped with s==0 (handled internally by the
      patch mapping)*/
    s1 = 0 ; t1 = (agg_surface_blend_curve_reverse(B,0) ? 1.0-t : t) ;
  }
  if ( h2 > 1 ) {
    if ( agg_surface_blend_curve_reverse(B,1) )
      agg_curve_eval(c2, 1.0-t, &s2, &t2) ;
    else
      agg_curve_eval(c2,     t, &s2, &t2) ;
  } else {
    g_assert(h2 == 0) ;
    /*cut-off ends are mapped with s==0 (handled internally by the
      patch mapping)*/
    s2 = 0 ; t2 = (agg_surface_blend_curve_reverse(B,1) ? 1.0-t : t) ;
  }

  use1 = use2 = TRUE ;
  
  agg_patch_point_diff(S1, P1, s1, t1, C1, S1s, S1t, w) ;
  
  if ( agg_vector_length(S1s) == 0 || agg_vector_length(S1t) == 0 )
    use1 = FALSE ;
  agg_vector_cross(Z1, S1s, S1t) ;

  agg_patch_point_diff(S2, P2, s2, t2, C2, S2s, S2t, w) ;
  if ( agg_vector_length(S2s) == 0 || agg_vector_length(S2t) == 0 )
    use2 = FALSE ;
  agg_vector_cross(Z2, S2s, S2t) ;

  agg_vector_diff(K, C2, C1) ;

  tmp1 = agg_vector_scalar(K,Z1)/agg_vector_scalar(Z1,Z1) ;
  V1[0] = K[0] - tmp1*Z1[0] ;
  V1[1] = K[1] - tmp1*Z1[1] ;
  V1[2] = K[2] - tmp1*Z1[2] ;

  tmp1 = agg_vector_scalar(K,Z2)/agg_vector_scalar(Z2,Z2) ;
  V2[0] = K[0] - tmp1*Z2[0] ;
  V2[1] = K[1] - tmp1*Z2[1] ;
  V2[2] = K[2] - tmp1*Z2[2] ;
  
  tmp1 = agg_vector_length(V1) ; tmp2 = agg_vector_length(V2) ;
  N1[0] = V1[0]/tmp1 ; N1[1] = V1[1]/tmp1 ; N1[2] = V1[2]/tmp1 ;
  N2[0] = V2[0]/tmp2 ; N2[1] = V2[1]/tmp2 ; N2[2] = V2[2]/tmp2 ;

  /*
   * if a surface is cut on a constant parameter line, we use a
   * different definition of the tangent vector
   *
   * don't know why these need to be negative, but they do (sign
   * convention for projection in Filip's paper?)
   */
  /* if ( agg_patch_clipping_type(c1) == AGG_CLIP_CONSTANT_S ) { */
  /* if ( h1 == 0 ) { */
  /*   N1[0] = -S1s[0] ; N1[1] = -S1s[1] ; N1[2] = -S1s[2] ; */
  /* } */
  /* /\* if ( agg_patch_clipping_type(c2) == AGG_CLIP_CONSTANT_S ) { *\/ */
  /* if ( h2 == 0 ) { */
  /*   N2[0] = -S2s[0] ; N2[1] = -S2s[1] ; N2[2] = -S2s[2] ; */
  /* } */

  /*check for valid tangent vectors*/
  g1 = agg_vector_length(K) + agg_vector_scalar(N1, K) ;
  l1 = 2.0*agg_vector_scalar(K,K)/g1 ;
  g2 = agg_vector_length(K) + agg_vector_scalar(N2, K) ;
  l2 = 2.0*agg_vector_scalar(K,K)/g2 ;

  T1[0] = w1*l1*N1[0] ; T1[1] = w1*l1*N1[1] ; T1[2] = w1*l1*N1[2] ;
  T2[0] = w2*l2*N2[0] ; T2[1] = w2*l2*N2[1] ; T2[2] = w2*l2*N2[2] ;
  
  agg_hermite_eval(s, H) ;

  x[0] = H[0]*C1[0] + H[1]*C2[0] ;
  x[1] = H[0]*C1[1] + H[1]*C2[1] ;
  x[2] = H[0]*C1[2] + H[1]*C2[2] ;
  
  if ( !use1 || !use2 ) return 0 ;
  
  /* if ( use1 ) { */
  /*   x[0] += H[2]*T1[0] ; */
  /*   x[1] += H[2]*T1[1] ; */
  /*   x[2] += H[2]*T1[2] ; */
  /* } */

  /* if ( use2 ) { */
  /*   x[0] += H[3]*T2[0] ; */
  /*   x[1] += H[3]*T2[1] ; */
  /*   x[2] += H[3]*T2[2] ; */
  /* } */
  
  return 0 ;
}

/**
 *  @}
 */
