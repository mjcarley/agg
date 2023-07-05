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

#include "hefsi.h"

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


static gboolean duplicate_points(hefsi_segment_t *seg1, hefsi_segment_t *seg2)

{
  if ( seg1->x1[3] != seg2->x1[3] ) return FALSE ;
  if ( seg1->x1[4] != seg2->x1[4] ) return FALSE ;
  
  return TRUE ;
}

static gboolean boundary_crossed(agg_patch_t *P1, agg_patch_t *P2, 
				 hefsi_segment_t *s1, hefsi_segment_t *s2,
				 gdouble *st, gdouble *sn)

{
  gdouble wtol, t, a ;
  gboolean crossed ;
  
  wtol = 2e-1 ;
  crossed = FALSE ;
  if ( agg_patch_wrap_t(P1) ) {
    if ( s1->x1[4] > 1.0-wtol && s2->x1[4] < wtol ) {
      fprintf(stderr, "HELLO1\n") ;
      t = s1->x1[4] - 1.0 ;
      a = (0.0 - t)/(s2->x1[4] - t) ;
      
      st[0] = s1->x1[3] + (s2->x1[3] - s1->x1[3])*a ;
      st[1] = 1.0 ;
      st[2] = 0.5*(s1->x1[5] + s2->x1[5]) ;
      st[3] = 0.5*(s1->x1[6] + s2->x1[6]) ;

      sn[0] = st[0] ; sn[1] = 0.0 ; sn[2] = st[2] ; sn[3] = st[3] ;
      crossed = TRUE ;
    }

    if ( s2->x1[4] > 1.0-wtol && s1->x1[4] < wtol ) {
      fprintf(stderr, "HELLO3\n") ;

      t = s1->x1[4] ;
      a = t/(t - (s2->x1[4]-1)) ;
      
      st[0] = s1->x1[3] + (s2->x1[3] - s2->x1[3])*a ;
      st[1] = 0.0 ;
      st[2] = 0.5*(s1->x1[5] + s2->x1[5]) ;
      st[3] = 0.5*(s1->x1[6] + s2->x1[6]) ;

      sn[0] = st[0] ; sn[1] = 1.0 ; sn[2] = st[2] ; sn[3] = st[3] ;
      crossed = TRUE ;
    }
  }

  if ( agg_patch_wrap_t(P2) ) {
    if ( s1->x1[6] > 1.0-wtol && s2->x1[6] < wtol ) {
      fprintf(stderr, "HELLO2\n") ;
      t = s1->x1[6] - 1.0 ;
      a = (0.0 - t)/(s2->x1[6] - t) ;
      
      st[2] = s1->x1[5] + (s2->x1[5] - s1->x1[5])*a ;
      st[3] = 1.0 ;
      st[0] = 0.5*(s1->x1[3] + s2->x1[3]) ;
      st[1] = 0.5*(s1->x1[4] + s2->x1[4]) ;
      sn[2] = st[2] ; sn[3] = 0.0 ; sn[0] = st[0] ; sn[1] = st[1] ;
      crossed = TRUE ;
    }
    if ( s2->x1[6] > 1.0-wtol && s1->x1[6] < wtol ) {
      fprintf(stderr, "HELLO4\n") ;
      t = s1->x1[6] ;
      a = t/(t - (s2->x1[6]-1)) ;
      
      st[2] = s1->x1[5] + (s2->x1[5] - s2->x1[5])*a ;
      st[3] = 0.0 ;
      st[0] = 0.5*(s1->x1[3] + s2->x1[3]) ;
      st[1] = 0.5*(s1->x1[4] + s2->x1[4]) ;

      sn[2] = st[2] ; sn[3] = 1.0 ; sn[0] = st[0] ; sn[1] = st[1] ;
      crossed = TRUE ;
    }
  }
  
  return crossed ;
}

static void densify_intersection(agg_intersection_t *inter,
				 hefsi_surface_t *h1, hefsi_surface_t *h2,
				 hefsi_segment_t *s1, hefsi_segment_t *s2, 
				 gint np, gdouble tol,
				 agg_surface_workspace_t *w)

{
  gint i, j, ns, nst ;
  gdouble t, u, v, st[16], *sti, x[3] ;
  agg_patch_t *P1, *P2 ;
  agg_surface_t *S1 ;

  P1 = agg_intersection_patch1(inter) ;
  P2 = agg_intersection_patch2(inter) ;
  S1 = agg_intersection_surface1(inter) ;
  
  i = 0 ; ns = 1 ;
  st[i*4+0] = s1->x1[3] ; st[i*4+1] = s1->x1[4] ;
  st[i*4+2] = s1->x1[5] ; st[i*4+3] = s1->x1[6] ; 
  if ( boundary_crossed(P1, P2, s1, s2, &(st[1*4]), &(st[2*4])) ) {
    ns += 2 ;
  }
  st[ns*4+0] = s2->x1[3] ; st[ns*4+1] = s2->x1[4] ;
  st[ns*4+2] = s2->x1[5] ; st[ns*4+3] = s2->x1[6] ; 
  ns ++ ;

  for ( i = 0 ; i < np ; i ++ ) {
    t = (gdouble)i/np ;
    nst = agg_intersection_point_number(inter) ;
    sti = &(inter->st[nst*AGG_INTERSECTION_DATA_SIZE]) ;
    for ( j = 0 ; j < 4 ; j ++ ) {
      sti[j] = st[0*4+j] + t*(st[1*4+j] - st[0*4+j]) ;
    }
    g_assert(sti[1] >= 0) ;
    agg_patch_map(P1, sti[0], sti[1], &u, &v) ;
    agg_surface_point_eval(S1, u, v, x, w) ;
    /* hefsi_refine_point(h1, &(sti[0]), h2, &(sti[2]), x, tol, 8) ; */
    agg_patch_st_correct(P1, &(sti[0])) ;
    agg_patch_st_correct(P2, &(sti[2])) ;
    agg_intersection_point_number(inter) ++ ;
  }

  if ( ns < 3 ) return ;
  /*wrapping over a boundary, need to add the last point from the
    previous interval*/  
  nst = agg_intersection_point_number(inter) ;
  sti = &(inter->st[nst*AGG_INTERSECTION_DATA_SIZE]) ;
  for ( j = 0 ; j < 4 ; j ++ ) {
    sti[j] = st[1*4+j] ;
  }
  agg_patch_map(P1, sti[0], sti[1], &u, &v) ;
  agg_surface_point_eval(S1, u, v, x, w) ;
  agg_intersection_point_number(inter) ++ ;

  for ( i = 0 ; i < np ; i ++ ) {
    t = (gdouble)i/np ;
    nst = agg_intersection_point_number(inter) ;
    sti = &(inter->st[nst*AGG_INTERSECTION_DATA_SIZE]) ;
    for ( j = 0 ; j < 4 ; j ++ ) {
      sti[j] = st[2*4+j] + t*(st[3*4+j] - st[2*4+j]) ;
    }
    g_assert(sti[1] >= 0) ;
    agg_patch_map(P1, sti[0], sti[1], &u, &v) ;
    agg_surface_point_eval(S1, u, v, x, w) ;
    /* hefsi_refine_point(h1, &(sti[0]), h2, &(sti[2]), x, tol, 8) ; */
    agg_intersection_point_number(inter) ++ ;
  }
  
  return ;
}

gint agg_surface_patch_intersection(agg_intersection_t *inter,
				    agg_surface_t *S1, agg_patch_t *P1,
				    agg_surface_t *S2, agg_patch_t *P2,
				    agg_surface_workspace_t *w)

{
  hefsi_surface_t *h1, *h2 ;
  gpointer data1[] = {S1, P1, w} ;
  gpointer data2[] = {S2, P2, w} ;
  hefsi_workspace_t *wh ;
  hefsi_segment_t *seg1, *seg2 ;
  gint dmin, dmax, i, j ;
  gdouble scale, tol, u, v ;
  GSList *il ;

  dmin = 6 ; dmax = 10 ; scale = 18/16.0 ; tol = 1e-6 ;
  
  wh = hefsi_workspace_new() ;

  agg_intersection_surface1(inter) = S1 ; 
  agg_intersection_patch1(inter) = P1 ; 
  agg_intersection_surface2(inter) = S2 ; 
  agg_intersection_patch2(inter) = P2 ; 

  h1 = hefsi_surface_new(hefsi_func, data1, FALSE, 0, 1, 0, 1) ;
  h2 = hefsi_surface_new(hefsi_func, data2, FALSE, 0, 1, 0, 1) ;

  hefsi_surface_initialize(h1, dmin, dmax) ;
  hefsi_set_bounding_boxes(h1, scale) ;
  hefsi_surface_initialize(h2, dmin, dmax) ;
  hefsi_set_bounding_boxes(h2, scale) ;
  
  hefsi_surface_intersections(h1, h2, tol, wh) ;

  if ( wh->c->len == 0 ) return 0 ;
  /*set the segment points to their midpoints and refine*/
  for ( il = hefsi_workspace_curve(wh,0) ; il != NULL ;
	il = il->next ) {
    j = GPOINTER_TO_INT(il->data) ;
    seg1 = hefsi_workspace_segment(wh,j) ;
    for ( j = 0 ; j < 7 ; j ++ ) {
      seg1->x1[j] = seg1->x1[j] ;
	/* 0.5*(seg1->x1[j] + seg1->x2[j]) ; */
    }
    
    /* hefsi_refine_point(h1, &(seg1->x1[3]), h2, &(seg1->x1[5]), seg1->x1, */
    /* 		       tol, 8) ; */
    agg_patch_st_correct(P1, &(seg1->x1[3])) ;
    agg_patch_st_correct(P2, &(seg1->x1[5])) ;
  }

  for ( il = hefsi_workspace_curve(wh,0) ; il->next != NULL ;
	il = il->next ) {  
    j = GPOINTER_TO_INT(il->data) ;
    seg1 = hefsi_workspace_segment(wh,j) ;
    j = GPOINTER_TO_INT(il->next->data) ;
    seg2 = hefsi_workspace_segment(wh,j) ;
    if ( !duplicate_points(seg1, seg2) ) {
      densify_intersection(inter, h1, h2, seg1, seg2, 2, tol, w) ;
    }
  }

  for ( i = 0 ; i < agg_intersection_point_number(inter) ; i ++ ) {
    agg_patch_map(P1,
		  agg_intersection_point_s1(inter,i),
		  agg_intersection_point_t1(inter,i),
		  &u, &v) ;
    agg_surface_point_eval(S1, u, v,
			   agg_intersection_point(inter,i),
			   w) ;
  }
  
  /* for ( i = 0 ; i < agg_intersection_point_number(inter) ; i ++ ) { */
  /*   hefsi_refine_point(h1, &(inter->st[i*AGG_INTERSECTION_DATA_SIZE+0]), */
  /* 		       h2, &(inter->st[i*AGG_INTERSECTION_DATA_SIZE+2]), */
  /* 		       x, tol, 1) ; */
  /* } */
  
  return 0 ;
  
  /* /\* agg_intersection_point_number(inter) = nst ; *\/ */

  /* /\*jump between segments in st plane marks start of curve*\/ */
  /* rmax1 = rmax2 = 0 ; */
  /* istart1 = istart2 = iend1 = iend2 = NULL ; */
  /* i = 0 ; ni = 0 ; */
  /* for ( il = hefsi_workspace_curve(wh,i) ; il->next != NULL ; */
  /* 	il = il->next ) { */
  /*   j = GPOINTER_TO_INT(il->data) ; */
  /*   seg1 = hefsi_workspace_segment(wh,j) ; */
  /*   j = GPOINTER_TO_INT(il->next->data) ; */
  /*   seg2 = hefsi_workspace_segment(wh,j) ; */
  /*   if ( (r = segment_distance2(seg1, seg2, 3)) > rmax1 ) { */
  /*     rmax1 = r ; */
  /*     iend1 = il ; */
  /*     istart1 = il->next ; */
  /*   } */
  /*   if ( (r = segment_distance2(seg1, seg2, 5)) > rmax2 ) { */
  /*     rmax2 = r ; */
  /*     iend2 = il ; */
  /*     istart2 = il->next ; */
  /*   } */
  /*   ni ++ ; */
  /* } */

  /* return 0 ; */
}

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

gint agg_intersection_resample(agg_intersection_t *inter,
			       gint nsp, gint pps,
			       agg_intersection_t *resample,
			       agg_surface_workspace_t *w)

{
  gint np, i, c ;
  gdouble tmin, tmax, t, *st, u, v ;

  agg_intersection_patch1(resample) = agg_intersection_patch1(inter) ;
  agg_intersection_patch2(resample) = agg_intersection_patch2(inter) ;
  agg_intersection_surface1(resample) = agg_intersection_surface1(inter) ;
  agg_intersection_surface2(resample) = agg_intersection_surface2(inter) ;
  
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
  }

  agg_intersection_point_number(resample) = np ;
  
  return 0 ;
}
