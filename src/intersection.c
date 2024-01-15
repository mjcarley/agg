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

#include <blaswrap.h>

#include "agg-private.h"

#include "hefsi.h"

/** 
 * @{ 
 *
 * @ingroup intersections
 */

static gint hefsi_func(gdouble s, gdouble t, gdouble *x, gpointer data)

{
  gpointer *hdata = data ;
  agg_surface_t *S = hdata[0] ;
  agg_patch_t *P = hdata[1] ;
  agg_surface_workspace_t *w = hdata[2] ;
  gdouble u, v ;

  agg_patch_map(P, s, t, &u, &v) ;  
  
  agg_surface_point_eval(S, u, v, x, w) ;

  return 0 ;
}

/** 
 * Allocate a new surface intersection structure. 
 * 
 * @param nstmax maximum number of points on intersection curve.
 * 
 * @return newly allocated ::agg_intersection_t.
 */

agg_intersection_t *agg_intersection_new(gint nstmax)

{
  agg_intersection_t *inter ;

  inter = (agg_intersection_t *)g_malloc(sizeof(agg_intersection_t)) ;

  inter->st = (gdouble *)g_malloc0(AGG_INTERSECTION_DATA_SIZE*nstmax*
				   sizeof(gdouble)) ;

  agg_intersection_point_number(inter) = 0 ;
  agg_intersection_point_number_max(inter) = nstmax ;
  
  return inter ;
}


static void fit_patch_end_cut(hefsi_workspace_t *wh, gint idx,
			      gdouble del, gint N, gdouble *C)

/*
 * fit a Fourier series to the trimming curve on a patch
 */

{
  gdouble A[16384], t, s, rhs[2048], work[1024] ;
  gint i, j, m, n, lwork, info, i1 = 1, ldb ;
  GSList *il ;
  hefsi_segment_t *seg ;

  /*matrix size*/
  n = 2*N + 1 ;
  m = 0 ;
  for ( i = 0 ; i < wh->c->len ; i ++ ) {
    for ( il = hefsi_workspace_curve(wh,i) ; il != NULL ;
	  il = il->next ) {
      j = GPOINTER_TO_INT(il->data) ;
      seg = hefsi_workspace_segment(wh,j) ;
      /*parameter values at this intersection*/
      s = seg->x1[3+2*idx+0] ; t = seg->x1[3+2*idx+1] ;
      rhs[m] = s + del ;
      A[m*n+0] = 1.0 ;
      for ( j = 1 ; j <= N ; j ++ ) {
	A[m*n+2*j-1] = cos(2.0*M_PI*j*t) ;
	A[m*n+2*j-0] = sin(2.0*M_PI*j*t) ;
      }
      m ++ ;
    }
  }

  lwork = 1024 ; info = 0 ; ldb = m ;
  dgels_("T", &n, &m, &i1, A, &n, rhs, &ldb, work, &lwork, &info) ;

  for ( i = 0 ; i <= n ; i ++ ) C[i] = rhs[i] ;
  
  return ;
}

static void fit_patch_hole(hefsi_workspace_t *wh, gint idx,
			   gdouble del, gdouble *C)

/*
 * fit an ellipse to a hole in a patch; for now this just uses the
 * bounding box, will implement 
 *
 * Least-squares orthogonal distances fitting of circle, sphere,
 * ellipse, hyperbola, and parabola, Sung Joon Ahn, Wolfgang Rauh,
 * Hans-Juergen Warnecke
 */

{
  gdouble s0, t0, smin, smax, tmin, tmax, s, t ;
  gint i, j, m ;
  GSList *il ;
  hefsi_segment_t *seg ;

  /*find hole centre*/
  s0 = t0 = 0 ; m = 0 ;
  smin = tmin =  G_MAXDOUBLE ;
  smax = tmax = -G_MAXDOUBLE ;
  for ( i = 0 ; i < wh->c->len ; i ++ ) {
    for ( il = hefsi_workspace_curve(wh,i) ; il != NULL ;
	  il = il->next ) {
      j = GPOINTER_TO_INT(il->data) ;
      seg = hefsi_workspace_segment(wh,j) ;
      /*parameter values at this intersection*/
      s = seg->x1[3+2*idx+0] ;
      t = seg->x1[3+2*idx+1] ;

      smin = MIN(smin, s) ; tmin = MIN(tmin, t) ;
      smax = MAX(smax, s) ; tmax = MAX(tmax, t) ;
      
      s0 += s ; t0 += t ;
      m ++ ;
    }
  }

  s0 /= m ; t0 /= m ;
  /*ellipse centre*/
  C[0] = s0 ; C[1] = t0 ;
  /*semi-major and semi-minor axes*/
  C[2] = MAX(smax-s0,s0-smin) + del ;
  C[3] = MAX(tmax-t0,t0-tmin) + del ;

  return ;
}

static gint fit_patch_trim(agg_patch_t *P, hefsi_workspace_t *wh, gint idx,
			   gdouble del, gint N)

/*
 * fit a Fourier series to the trimming curve on a patch
 */

{
  gint nh ;
  agg_curve_t *c ;
  
  if ( agg_patch_mapping(P) == AGG_PATCH_HEMISPHERICAL ) {
    /*
     * this should really be a check to see if a whole end is cut off
     * 0 <= t <= 1, and which end (this assumes min)
     */
    c = agg_patch_curve_smin(P) ; nh = 0 ;
    g_assert(2*N+1 < AGG_CURVE_DATA_SIZE) ;
    fit_patch_end_cut(wh, idx, del, N, c->data) ;

    agg_curve_type(c) = AGG_CURVE_FOURIER ;
    agg_curve_order(c) = N ;

    return nh ;
  }

  nh = agg_patch_hole_number(P) ;
  c = agg_patch_hole(P,nh) ;
  fit_patch_hole(wh, idx, del, c->data) ;
  agg_curve_type(c) = AGG_CURVE_ELLIPSE ;
  agg_patch_hole_number(P) ++ ;

  /*index of hole on patch, for use in blended surface*/
  return nh ;
}

gboolean agg_surface_patch_trim(agg_surface_t *S1, agg_patch_t *P1, gdouble d1,
				agg_surface_t *S2, agg_patch_t *P2, gdouble d2,
				agg_surface_blend_t *B,
				agg_surface_workspace_t *w)

{
  hefsi_surface_t *h1, *h2 ;
  gpointer data1[] = {S1, P1, w} ;
  gpointer data2[] = {S2, P2, w} ;
  hefsi_workspace_t *wh ;
  hefsi_segment_t *seg1 ;
  gint dmin, dmax, i, j ;
  gdouble scale, tol, n1[3], n2[3] ;
  agg_curve_t *c1, *c2 ;
  GSList *il ;
  
  dmin = 6 ; dmax = 10 ; scale = 18/16.0 ; tol = 1e-6 ;
  
  wh = hefsi_workspace_new() ;

  h1 = hefsi_surface_new(hefsi_func, data1, FALSE, 0, 1, 0, 1) ;
  h2 = hefsi_surface_new(hefsi_func, data2, FALSE, 0, 1, 0, 1) ;

  hefsi_surface_initialize(h1, dmin, dmax) ;
  hefsi_set_bounding_boxes(h1, scale) ;
  hefsi_surface_initialize(h2, dmin, dmax) ;
  hefsi_set_bounding_boxes(h2, scale) ;
  
  hefsi_surface_intersections(h1, h2, tol, wh) ;

  if ( wh->c->len == 0 ) return FALSE ;
  for ( i = 0 ; i < wh->c->len ; i ++ ) {
    for ( il = hefsi_workspace_curve(wh,i) ; il != NULL ;
	  il = il->next ) {
      j = GPOINTER_TO_INT(il->data) ;
      seg1 = hefsi_workspace_segment(wh,j) ;
      for ( j = 0 ; j < 7 ; j ++ ) {
	seg1->x1[j] = 	0.5*(seg1->x1[j] + seg1->x2[j]) ;
      }
      agg_patch_st_correct(P1, &(seg1->x1[3])) ;
      agg_patch_st_correct(P2, &(seg1->x1[5])) ;
    }
  }

  agg_surface_blend_hole(B,0) = fit_patch_trim(P1, wh, 0, d1, 8) ;
  agg_surface_blend_hole(B,1) = fit_patch_trim(P2, wh, 1, d2, 16) ;
  agg_patch_blend(P1,agg_surface_blend_hole(B,0)) = B ;
  agg_patch_blend(P2,agg_surface_blend_hole(B,1)) = B ;
  agg_surface_blend_surface(B,0) = S1 ; 
  agg_surface_blend_surface(B,1) = S2 ; 
  agg_surface_blend_patch(B,0) = P1 ; 
  agg_surface_blend_patch(B,1) = P2 ; 
  agg_surface_blend_curve_reverse(B,0) = FALSE ;
  agg_surface_blend_curve_reverse(B,1) = FALSE ;

  /*check orientation of patch curves and reverse one if necessary*/
  c1 = &(P1->curves[agg_surface_blend_hole(B,0)]) ;
  c2 = &(P2->curves[agg_surface_blend_hole(B,1)]) ;
  agg_curve_plane_normal(c1, S1, P1, n1, w) ;
  agg_curve_plane_normal(c2, S2, P2, n2, w) ;

  /* fprintf(stderr, "(%lg,%lg,%lg).(%lg,%lg,%lg) = %lg\n", */
  /* 	  n1[0], n1[1], n1[2], n2[0], n2[1], n2[2], */
  /* 	  agg_vector_scalar(n1,n2)) ; */

  if ( agg_vector_scalar(n1,n2) < 0 ) {
    /*reverse one rail curve*/
    if ( agg_surface_blend_hole(B,0) == 0 ||
	 agg_surface_blend_hole(B,0) == 1 ) {
      agg_surface_blend_curve_reverse(B,1) = TRUE ;
    } else {
      agg_surface_blend_curve_reverse(B,0) = TRUE ;
    }
  }
  
  return TRUE ;
}

/** 
 * Write an intersection curve to file in the form
 *
 * [\f$x\f$ \f$y\f$ \f$z\f$ \f$s_{1}\f$ \f$t_{1}\f$ \f$s_{2}\f$ \f$t_{2}\f$]
 *
 * where \f$(x,y,z)\f$ are physical coordinates and
 * \f$(s_{i},t_{i})\f$ are parametric coordinates on the respective
 * surface patches (not the surfaces proper)
 * 
 * @param f output file stream;
 * @param inter intersection curve to write;
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */

gint agg_intersection_curve_write(FILE *f, agg_intersection_t *inter,
				  agg_surface_workspace_t *w)

{
  gdouble x[3], u, v ;
  gint i ;

  for ( i = 0 ; i < agg_intersection_point_number(inter) ; i ++ ) {
    agg_patch_map(agg_intersection_patch1(inter),
		  agg_intersection_point_s1(inter,i),
		  agg_intersection_point_t1(inter,i),
		  &u, &v) ;
    agg_surface_point_eval(agg_intersection_surface1(inter), u, v, x, w) ;
    fprintf(f, "%lg %lg %lg %lg %lg %lg %lg\n",
	    x[0], x[1], x[2],
	    agg_intersection_point_s1(inter,i),
	    agg_intersection_point_t1(inter,i),
	    agg_intersection_point_s2(inter,i),
	    agg_intersection_point_t2(inter,i)) ;
  }
  
  return 0 ;
}


/** 
 * Correct out of range values of \f$(s,t)\f$ to lie in the range
 * \f$(0,1)\times(0,1)\f$ where coordinate wrapping is permitted
 * 
 * @param P an ::agg_patch_t;
 * @param st an \f$(s,t)\f$ parameter pair on patch \a P.
 * 
 * @return 0 on success.
 */

gint agg_patch_st_correct(agg_patch_t *P, gdouble *st)

{
  if ( st[0] < 0 && agg_patch_wrap_s(P) ) { st[0] = 1.0 + st[0] ; }
  if ( st[0] > 1 && agg_patch_wrap_s(P) ) { st[0] = st[0] - 1.0; }
  if ( st[1] < 0 && agg_patch_wrap_t(P) ) { st[1] = 1.0 + st[1] ; }
  if ( st[1] > 1 && agg_patch_wrap_t(P) ) { st[1] = st[1] - 1.0; }

  return 0 ;
}

/** 
 * Calculate the bounding box of an intersection curve in the two
 * surface patch mappings.
 * 
 * @param inter intersection curve whose bounding boxes are to be set. 
 * 
 * @return 0 on success.
 */

gint agg_intersection_bbox_set(agg_intersection_t *inter)

{
  gint i ;

  agg_intersection_bbox_s1_min(inter) = 
    agg_intersection_bbox_t1_min(inter) = 
    agg_intersection_bbox_s2_min(inter) = 
    agg_intersection_bbox_t2_min(inter) = G_MAXDOUBLE ;
  agg_intersection_bbox_s1_max(inter) = 
    agg_intersection_bbox_t1_max(inter) = 
    agg_intersection_bbox_s2_max(inter) = 
    agg_intersection_bbox_t2_max(inter) = -G_MAXDOUBLE ;

  for ( i = 0 ; i < agg_intersection_point_number(inter) ; i ++ ) {
    if ( agg_intersection_point_s1(inter,i) <
	 agg_intersection_bbox_s1_min(inter) ) {
      agg_intersection_bbox_s1_min(inter) =
	agg_intersection_point_s1(inter,i) ;
      agg_intersection_bbox_s1_min_index(inter) = i ;
    }
    if ( agg_intersection_point_t1(inter,i) <
	 agg_intersection_bbox_t1_min(inter) ) {
      agg_intersection_bbox_t1_min(inter) =
	agg_intersection_point_t1(inter,i) ;
      agg_intersection_bbox_t1_min_index(inter) = i ;
    }
    if ( agg_intersection_point_s2(inter,i) <
	 agg_intersection_bbox_s2_min(inter) ) {
      agg_intersection_bbox_s2_min(inter) =
	agg_intersection_point_s2(inter,i) ;
      agg_intersection_bbox_s2_min_index(inter) = i ;
    }
    if ( agg_intersection_point_t2(inter,i) <
	 agg_intersection_bbox_t2_min(inter) ) {
      agg_intersection_bbox_t2_min(inter) =
	agg_intersection_point_t2(inter,i) ;
      agg_intersection_bbox_t2_min_index(inter) = i ;
    }

    if ( agg_intersection_point_s1(inter,i) >
	 agg_intersection_bbox_s1_max(inter) ) {
      agg_intersection_bbox_s1_max(inter) =
	agg_intersection_point_s1(inter,i) ;
      agg_intersection_bbox_s1_max_index(inter) = i ;
    }
    if ( agg_intersection_point_t1(inter,i) >
	 agg_intersection_bbox_t1_max(inter) ) {
      agg_intersection_bbox_t1_max(inter) =
	agg_intersection_point_t1(inter,i) ;
      agg_intersection_bbox_t1_max_index(inter) = i ;
    }
    if ( agg_intersection_point_s2(inter,i) >
	 agg_intersection_bbox_s2_max(inter) ) {
      agg_intersection_bbox_s2_max(inter) =
	agg_intersection_point_s2(inter,i) ;
      agg_intersection_bbox_s2_max_index(inter) = i ;
    }
    if ( agg_intersection_point_t2(inter,i) >
	 agg_intersection_bbox_t2_max(inter) ) {
      agg_intersection_bbox_t2_max(inter) =
	agg_intersection_point_t2(inter,i) ;
      agg_intersection_bbox_t2_max_index(inter) = i ;
    }    
  }
  
  return 0 ;
}

static gboolean bracketed(gdouble x0, gdouble x1, gdouble x, gdouble *t)

{
  if ( ABS(x1 - x0) > 0.5 ) return FALSE ;
  
  if ( x0 <= x && x <= x1 ) {
    *t = (x - x0)/(x1 - x0) ;
    return TRUE ;
  }

  if ( x1 <= x && x <= x0 ) {
    *t = (x0 - x)/(x0 - x1) ;
    return TRUE ;
  }
  
  
  return FALSE ;
}

static void interp_st(gdouble *st, gint nst, gint c,
		      gdouble *sti, gdouble t)

{
  gint i, j ;
  gdouble ti ;
  
  for ( i = 0 ; i < nst ; i ++ ) {
    if ( bracketed(st[(i+0)*AGG_INTERSECTION_DATA_SIZE+c],
		   st[(i+1)*AGG_INTERSECTION_DATA_SIZE+c],
		   t, &ti) ) {
      for ( j = 0 ; j < 4 ; j ++ ) {
	sti[j] =
	  st[(i+0)*AGG_INTERSECTION_DATA_SIZE+j] + 	  
	  (st[(i+1)*AGG_INTERSECTION_DATA_SIZE+j] -
	   st[(i+0)*AGG_INTERSECTION_DATA_SIZE+j])*ti ;
      }
      return ;
    }
  }

  i = nst - 1 ;
  if ( bracketed(st[(i+0)*AGG_INTERSECTION_DATA_SIZE+c],
		 st[(  0)*AGG_INTERSECTION_DATA_SIZE+c],
		 t, &ti) ) {
    for ( j = 0 ; j < 4 ; j ++ ) {
      sti[j] =
	st[(i+0)*AGG_INTERSECTION_DATA_SIZE+j] + 	  
	(st[( 0)*AGG_INTERSECTION_DATA_SIZE+j] -
	 st[(i+0)*AGG_INTERSECTION_DATA_SIZE+j])*ti ;
      }
    return ;
  }
  
  g_assert_not_reached() ;
  
  return ;
}

/** 
 * Resample an intersection curve to (approximately) evenly spaced
 * points in the mapping patch.
 * 
 * @param inter an intersection curve;
 * @param nsp number of splines in the resampled curve;
 * @param pps number of points per spline in the resampled curve;
 * @param resample on output contains the resampled intersection data;
 * @param w workspace for surface evaluation.
 * 
 * @return 0 on success.
 */

gint agg_intersection_resample(agg_intersection_t *inter,
			       gint nsp, gint pps,
			       agg_intersection_t *resample,
			       agg_surface_workspace_t *w)

{
  gint np, i, c ;
  gdouble tmin, tmax, t, *st, u, v ;
  hefsi_surface_t *h1, *h2 ;
  gpointer data1[] =
    {agg_intersection_surface1(inter),
     agg_intersection_patch1(inter),
     w} ;
  gpointer data2[] =
    {agg_intersection_surface2(inter),
     agg_intersection_patch2(inter),
     w} ;

  agg_intersection_patch1(resample) = agg_intersection_patch1(inter) ;
  agg_intersection_patch2(resample) = agg_intersection_patch2(inter) ;
  agg_intersection_surface1(resample) = agg_intersection_surface1(inter) ;
  agg_intersection_surface2(resample) = agg_intersection_surface2(inter) ;
  
  h1 = hefsi_surface_new(hefsi_func, data1, FALSE, 0, 1, 0, 1) ;
  h2 = hefsi_surface_new(hefsi_func, data2, FALSE, 0, 1, 0, 1) ;

  /*find the limits on (s,t) on each surface*/
  agg_intersection_bbox_set(inter) ;

  /*number of points in resample*/
  np = nsp*(pps-1) + 1 ;

  agg_intersection_point_number(resample) = 0 ;

  if ( agg_intersection_point_number(inter) == 0 ) return 0 ;
  
  st = resample->st ;

  /*select surface to use for parametric sweep*/
  if ( agg_intersection_bbox_t1_min(inter) == 0 &&
       agg_intersection_bbox_t1_max(inter) == 1 ) {
    tmin = 0 ; tmax = 1 ; c = 1 ;
  } else {
    tmin = agg_intersection_bbox_t2_min(inter) ;
    tmax = agg_intersection_bbox_t2_max(inter) ;
    c = 3 ;
  }

  for ( i = 0 ; i < np ; i ++ ) {
    t = tmin + (tmax - tmin)*(gdouble)i/(np-1) ;
    interp_st(inter->st, agg_intersection_point_number(inter), c,
	      &(st[i*AGG_INTERSECTION_DATA_SIZE]), t) ;

    agg_patch_map(agg_intersection_patch1(resample),
		  agg_intersection_point_s1(resample,i),
		  agg_intersection_point_t1(resample,i),
		  &u, &v) ;
    agg_surface_point_eval(agg_intersection_surface1(resample), u, v,
			   agg_intersection_point(resample,i), w) ;
    if ( t != 0 && t != 1 ) {
      hefsi_refine_point(h1, &(st[i*AGG_INTERSECTION_DATA_SIZE+0]),
			 h2, &(st[i*AGG_INTERSECTION_DATA_SIZE+2]),
			 agg_intersection_point(resample,i), 1e-9, 8) ;
    }
  }

  agg_intersection_point_number(resample) = np ;
  
  return 0 ;
}

#if 0
gint agg_intersection_clip(agg_intersection_t *inter,
			   gint c, agg_patch_clip_t cut,
			   agg_patch_clipping_t *clip)

{
  agg_intersection_bbox_set(inter) ;

  agg_patch_clipping_type(clip) = cut ;
  if ( cut == AGG_CLIP_CONSTANT_S ) {
    agg_patch_clipping_data(clip,0) = inter->bbox[4+2*c+0] ;

    return 0 ;
  }

  if ( cut == AGG_CLIP_CONSTANT_T ) {
    g_assert_not_reached() ;
    agg_patch_clipping_data(clip,0) = inter->bbox[4+2*c+1] ;

    return 0 ;
  }

  if ( cut == AGG_CLIP_ELLIPSE ) {
    gdouble th ;
    /*centre of ellipse is centre of bounding box of intersection
      curve*/
    agg_patch_clipping_data(clip,0) =
      (inter->bbox[2*c+0] + inter->bbox[4+2*c+0])/2 ;
    agg_patch_clipping_data(clip,1) =
      (inter->bbox[2*c+1] + inter->bbox[4+2*c+1])/2 ;

    /*angle to corner of bounding box*/
    th = atan2(inter->bbox[4+2*c+1]-agg_patch_clipping_data(clip,1),
	       inter->bbox[4+2*c+0]-agg_patch_clipping_data(clip,0)) ;

    /*semi-major and semi-minor axes*/
    agg_patch_clipping_data(clip,2) =
      (inter->bbox[4+2*c+0]-agg_patch_clipping_data(clip,0))/cos(th) ;
    agg_patch_clipping_data(clip,3) =
      (inter->bbox[4+2*c+1]-agg_patch_clipping_data(clip,1))/sin(th) ;
    
    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}
#endif

/**
 * @}
 */
