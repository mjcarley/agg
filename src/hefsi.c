/* This file is part of HEFSI, code for calculating surface intersections
 *
 * Copyright (C) 2023 Michael Carley
 *
 * HEFSI is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.  HEFSI is distributed in the hope
 * that it will be useful, but WITHOUT ANY WARRANTY; without even the
 * implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HEFSI.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include "hefsi.h"

int tri_tri_intersect_with_isectline(double V0[3],double V1[3],double V2[3],
				     double U0[3],double U1[3],double U2[3],
				     int *coplanar, double isectpt1[3],
				     double isectpt2[3]) ;

static void make_normal(gdouble *n, gdouble *xu, gdouble *xv)

{
  gdouble len ;

  g_assert(!(xu[0] == 0 && xu[1] == 0 && xu[2] == 0)) ;
  g_assert(!(xv[0] == 0 && xv[1] == 0 && xv[2] == 0)) ;
  
  hefsi_vector_cross(n, xu, xv) ;

  len = hefsi_vector_length(n) ;
  n[0] /= len ; n[1] /= len ; n[2] /= len ; 

  return ;
}

static gdouble vector_cos_angle(gdouble *x1, gdouble *x2)
  
{
  gdouble c, l1, l2 ;

  c = hefsi_vector_dot(x1,x2) ;
  l1 = hefsi_vector_length(x1) ; 
  l2 = hefsi_vector_length(x2) ; 

  c /= l1*l2 ;
  
  return c ;
}

static void check_edge(gdouble *x1, gdouble *x2, gdouble *x3,
		       gdouble ncos, gboolean *split)

/*
 * check angle between normals and segments for three nodes on
 * isoparametric line
 */

{
  gdouble nc, dx1[3], dx2[3] ;

  /*angle between normals*/
  nc = vector_cos_angle(&(x1[9]), &(x2[9])) ;
  if ( nc < ncos ) { *split = TRUE ; return ; }
  nc = vector_cos_angle(&(x2[9]), &(x3[9])) ;
  if ( nc < ncos ) { *split = TRUE ; return ; }

  hefsi_vector_diff(dx1, x2, x1) ; 
  hefsi_vector_diff(dx2, x3, x2) ; 

  nc = vector_cos_angle(dx1, dx2) ;
  if ( nc < ncos ) { *split = TRUE ; return ; }
  
  return ;
}

static gboolean check_element_split(gint el,
				    hefsi_surface_t *s,
				    gdouble ncos, gdouble ecos)

{
  gdouble u, v, uc, vc ;
  gint i, j, k, edges[] = {0, 1, 2, 3, 0}, nn ;
  gboolean split ;
  hefsi_element_t enew, *eold ;
  hefsi_node_t node ;
  /*
   * generate five new nodes at midpoints of edges and in element
   * centre, and check for flatness
   */
  uc = vc = 0 ;
  nn = hefsi_surface_node_number(s) ;
  eold = hefsi_surface_element(s,el) ;
  for ( k = 0 ; k < 4 ; k ++ ) {
    i = eold->i[edges[k]] ; j = eold->i[edges[k+1]] ;
    u = (hefsi_surface_u(s,i) + hefsi_surface_u(s,j))/2 ;
    v = (hefsi_surface_v(s,i) + hefsi_surface_v(s,j))/2 ;
    hefsi_surface_evaluate(s, u, v,
			   hefsi_node_x(&node),
			   hefsi_node_xu(&node),
			   hefsi_node_xv(&node),
			   hefsi_node_n (&node)) ;
    uc += hefsi_surface_u(s,i) ; vc += hefsi_surface_v(s,i) ;
    hefsi_node_u(&node) = u ; hefsi_node_v(&node) = v ;
    g_array_append_val(s->n,node) ;
  }

  uc *= 0.25 ; vc *= 0.25 ; k = 4 ;
  hefsi_surface_evaluate(s, uc, vc,
			 hefsi_node_x(&node),
			 hefsi_node_xu(&node),
			 hefsi_node_xv(&node),
			 hefsi_node_n (&node)) ;
  hefsi_node_u(&node) = uc ; hefsi_node_v(&node) = vc ;
  g_array_append_val(s->n,node) ;

  /*
   * end of node list now contains new nodes and surface vectors
   */

  split = FALSE ;
  i = eold->i[0] ; j = nn+0 ; k = eold->i[1] ;
  check_edge(hefsi_surface_x(s,i), hefsi_surface_x(s,j),
	     hefsi_surface_x(s,k), ecos, &split) ;
  i = eold->i[1] ; j = nn+1 ; k = eold->i[2] ;
  check_edge(hefsi_surface_x(s,i), hefsi_surface_x(s,j),
	     hefsi_surface_x(s,k), ecos, &split) ;
  i = eold->i[2] ; j = nn+2 ; k = eold->i[3] ;
  check_edge(hefsi_surface_x(s,i), hefsi_surface_x(s,j),
	     hefsi_surface_x(s,k), ecos, &split) ;
  i = eold->i[3] ; j = nn+3 ; k = eold->i[0] ;
  check_edge(hefsi_surface_x(s,i), hefsi_surface_x(s,j),
	     hefsi_surface_x(s,k), ecos, &split) ;

  i = nn+0 ; j = nn+4 ; k = nn+2 ;
  check_edge(hefsi_surface_x(s,i), hefsi_surface_x(s,j),
	     hefsi_surface_x(s,k), ecos, &split) ;
  i = nn+1 ; j = nn+4 ; k = nn+3 ;
  check_edge(hefsi_surface_x(s,i), hefsi_surface_x(s,j),
	     hefsi_surface_x(s,k), ecos, &split) ;
  
  /*check angle between normals*/
  for ( i = 0 ; i < 5 ; i ++ ) {
    for ( j = i+1 ; j < 5 ; j ++ ) {
      if ( vector_cos_angle(hefsi_surface_n(s,nn+i),
			    hefsi_surface_n(s,nn+i)) < ncos ) {
	split = TRUE ;
      }
    }
  }
  
  if ( !split ) return FALSE ;

  enew.s = NULL ;
  
  eold = hefsi_surface_element(s,el) ;
  enew.i[0] = eold->i[0] ; 
  enew.i[1] = nn + 0 ;
  enew.i[2] = nn + 4 ;
  enew.i[3] = nn + 3 ;
  g_array_append_val(s->e, enew) ;

  eold = hefsi_surface_element(s,el) ;
  enew.i[0] = eold->i[1] ; 
  enew.i[1] = nn + 1 ;
  enew.i[2] = nn + 4 ;
  enew.i[3] = nn + 0 ;
  g_array_append_val(s->e, enew) ;
  
  eold = hefsi_surface_element(s,el) ;
  enew.i[0] = eold->i[2] ; 
  enew.i[1] = nn + 2 ;
  enew.i[2] = nn + 4 ;
  enew.i[3] = nn + 1 ;
  g_array_append_val(s->e, enew) ;

  eold = hefsi_surface_element(s,el) ;
  enew.i[0] = eold->i[3] ; 
  enew.i[1] = nn + 3 ;
  enew.i[2] = nn + 4 ;
  enew.i[3] = nn + 2 ;
  g_array_append_val(s->e, enew) ;
  
  return TRUE ;
}

static void surface_add_node(hefsi_surface_t *s, gdouble u, gdouble v)

{
  hefsi_node_t n ; 

  hefsi_surface_evaluate(s, u, v,
			 hefsi_node_x(&n),
			 hefsi_node_xu(&n),
			 hefsi_node_xv(&n),
			 hefsi_node_n (&n)) ;
  hefsi_node_u(&n) = u ; hefsi_node_v(&n) = v ;

  g_array_append_val(s->n,n) ;
  
  return ;
}  

static gboolean split_surface(GNode *node, hefsi_surface_t *s,
			      gdouble ncos, gdouble ecos,
			      gint dmin, gint dmax)

{
  gint i ;
  GNode *new, *c ;
  gboolean split ;
  
  if ( g_node_depth(node) >= dmax ) return TRUE ;

  split = FALSE ;

  i = GPOINTER_TO_INT(node->data) ;
  if ( g_node_depth(node) <  dmin ) {
    /* set ncos > 1 to force split */
    check_element_split(i, s, 1.1, ecos) ;
    split = TRUE ;
  } else {
    if ( check_element_split(i, s, ncos, ecos) ) split = TRUE ;
  }
  if ( !split ) return TRUE ;

  /*add the split elements as child nodes*/
  for ( i = s->e->len - 4 ; i < s->e->len ; i ++ ) {
    new = g_node_new(GINT_TO_POINTER(i)) ;
    g_node_insert(node, -1, new) ;
  }  

  for ( c = node->children ; c != NULL ; c = c->next ) {
    split_surface(c, s, ncos, ecos, dmin, dmax) ; 
  }
  
  return FALSE ;
}

/** 
 * Initialize elements of a ::hefsi_surface_t by recursion to
 * specified level.
 * 
 * @param s a ::hefsi_surface_t;
 * @param dmin minimum depth in element tree to which to subdivide elements.
 * 
 * @return 0 on success.
 */

gint hefsi_surface_initialize(hefsi_surface_t *s, gint dmin, gint dmax)

{
  gint i ;
  hefsi_element_t *e ;
  gdouble ncos, ecos ;

  if ( s->tree != NULL ) {
    g_node_destroy(s->tree) ;
  }
  s->tree = g_node_new(GINT_TO_POINTER(0)) ;

  for ( i = 0 ; i < s->e->len ; i ++ ) {
    memset(hefsi_surface_element(s,i),0,sizeof(hefsi_element_t)) ;
  }
  for ( i = 0 ; i < s->n->len ; i ++ ) {
    memset(hefsi_surface_node(s,i),0,sizeof(hefsi_node_t)) ;
  }

  g_array_set_size(s->n,0) ;
  g_array_set_size(s->e,0) ;
  surface_add_node(s, hefsi_surface_umin(s), hefsi_surface_vmin(s)) ;
  surface_add_node(s, hefsi_surface_umax(s), hefsi_surface_vmin(s)) ;
  surface_add_node(s, hefsi_surface_umax(s), hefsi_surface_vmax(s)) ;
  surface_add_node(s, hefsi_surface_umin(s), hefsi_surface_vmax(s)) ;

  g_array_set_size(s->e, 1) ;
  e = hefsi_surface_element(s,0) ;
  e->i[0] = 0 ; e->i[1] = 1 ; e->i[2] = 2 ; e->i[3] = 3 ;

  ncos = cos(M_PI/12) ; ecos = cos(M_PI/6.0) ;  
  split_surface(s->tree, s, ncos, ecos, dmin, dmax) ;
  
  return 0 ;
}

static gboolean traverse_write_tri(GNode *node, gpointer *tdata)

{
  hefsi_surface_t *s = tdata[HEFSI_TRAVERSE_DATA_SURFACE] ;
  hefsi_element_t *e ;
  FILE *f = tdata[HEFSI_TRAVERSE_DATA_FILE] ;
  gint i ;

  i = GPOINTER_TO_INT(node->data) ;  
  e = hefsi_surface_element(s,i) ;

  fprintf(f, "%d %d %d\n", e->i[0], e->i[1], e->i[2]) ;
  fprintf(f, "%d %d %d\n", e->i[2], e->i[3], e->i[0]) ;
  
  return FALSE ;
}

/** 
 * Write a triangulated ::hefsi_surface_t to file.
 *
 * Output format is:
 * 
 * [number of nodes] [number of triangles]
 * x y z u v (at each node)
 * [index of node 0] [index of node 1] [index of node 2] (for each triangle)
 * 
 * @param f 
 * @param s 
 * 
 * @return 
 */

gint hefsi_write_surface(FILE *f, hefsi_surface_t *s)

{
  gint i ;
  guint nle ;
  gdouble *x ;
  gpointer tdata[HEFSI_TRAVERSE_DATA_SIZE] ;
  
  /*need the number of leaf elements to write surface*/
  nle = g_node_n_nodes(s->tree, G_TRAVERSE_LEAVES) ;
  fprintf(f, "%d %u\n", hefsi_surface_node_number(s), 2*nle) ;
  for ( i = 0 ; i < hefsi_surface_node_number(s) ; i ++ ) {
    x = hefsi_surface_x(s,i) ;
    fprintf(f, "%e %e %e %e %e\n", x[0], x[1], x[2],
	    hefsi_surface_u(s,i), hefsi_surface_v(s,i)) ;
  }

  tdata[HEFSI_TRAVERSE_DATA_SURFACE] = s ;
  tdata[HEFSI_TRAVERSE_DATA_FILE] = f ;
  g_node_traverse(s->tree, G_IN_ORDER, G_TRAVERSE_LEAVES, -1,
		  (GNodeTraverseFunc)traverse_write_tri, tdata) ;
  
  return 0 ;
}

static void element_set_bbox(hefsi_surface_t *s,
			     gint i,
			     /* hefsi_element_t *e, */
			     
			     /* gdouble *nodes, */
			     gdouble scale)

{
  gint j ;
  gdouble *x, dx[3], xb[3] ;
  hefsi_node_t *node ;
  hefsi_element_t *e ;

  e = hefsi_surface_element(s,i) ;
  e->bbox[0] = e->bbox[2] = e->bbox[4] =  G_MAXDOUBLE ;
  e->bbox[1] = e->bbox[3] = e->bbox[5] = -G_MAXDOUBLE ;

  for ( j = 0 ; j < 4 ; j ++ ) {
    node = hefsi_surface_node(s,e->i[j]) ;
    x = hefsi_node_x(node) ;
    e->bbox[0] = MIN(x[0], e->bbox[0]) ;
    e->bbox[1] = MAX(x[0], e->bbox[1]) ;
    e->bbox[2] = MIN(x[1], e->bbox[2]) ;
    e->bbox[3] = MAX(x[1], e->bbox[3]) ;
    e->bbox[4] = MIN(x[2], e->bbox[4]) ;
    e->bbox[5] = MAX(x[2], e->bbox[5]) ;
  }

  dx[0] = (e->bbox[1] - e->bbox[0])/2.0 ;
  dx[1] = (e->bbox[3] - e->bbox[2])/2.0 ;
  dx[2] = (e->bbox[5] - e->bbox[4])/2.0 ;
  xb[0] = (e->bbox[1] + e->bbox[0])/2.0 ;
  xb[1] = (e->bbox[3] + e->bbox[2])/2.0 ;
  xb[2] = (e->bbox[5] + e->bbox[4])/2.0 ;

  e->bbox[0] = xb[0] - scale*dx[0] ; 
  e->bbox[1] = xb[0] + scale*dx[0] ; 
  e->bbox[2] = xb[1] - scale*dx[1] ; 
  e->bbox[3] = xb[1] + scale*dx[1] ; 
  e->bbox[4] = xb[2] - scale*dx[2] ; 
  e->bbox[5] = xb[2] + scale*dx[2] ; 
  
  return ;
}

static gboolean traverse_bboxes(GNode *node, gpointer *tdata)

{
  hefsi_surface_t *s = tdata[HEFSI_TRAVERSE_DATA_SURFACE] ;
  gdouble scale = *((gdouble *)tdata[HEFSI_TRAVERSE_DATA_SCALE]) ;
  hefsi_element_t *e, *ej ;
  gint i, j ;
  GNode *c ;
  
  if ( G_NODE_IS_LEAF(node) ) {
    /*need to find bbox limits from element*/
    i = GPOINTER_TO_INT(node->data) ;
    element_set_bbox(s, i, scale) ;
    return FALSE ;
  }

  /*limits are acquired from child nodes*/
  i = GPOINTER_TO_INT(node->data) ;
  e = hefsi_surface_element(s,i) ;
  e->bbox[0] = e->bbox[2] = e->bbox[4] =  G_MAXDOUBLE ;
  e->bbox[1] = e->bbox[3] = e->bbox[5] = -G_MAXDOUBLE ;

  for ( c = node->children ; c != NULL ; c = c->next ) {
    j = GPOINTER_TO_INT(c->data) ;
    ej = hefsi_surface_element(s,j) ;
    e->bbox[0] = MIN(e->bbox[0], ej->bbox[0]) ;
    e->bbox[1] = MAX(e->bbox[1], ej->bbox[1]) ;
    e->bbox[2] = MIN(e->bbox[2], ej->bbox[2]) ;
    e->bbox[3] = MAX(e->bbox[3], ej->bbox[3]) ;
    e->bbox[4] = MIN(e->bbox[4], ej->bbox[4]) ;
    e->bbox[5] = MAX(e->bbox[5], ej->bbox[5]) ;
  }
  
  return FALSE ;
}

/** 
 * Set bounding boxes on elements of a ::hefsi_surface_t, scaling to
 * allow for elements which extend beyond their nodes.
 * 
 * @param s a ::hefsi_surface_t with elements allocated;
 * @param scale factor by which to scale element bounding boxes.
 * 
 * @return 0 on success.
 */

gint hefsi_set_bounding_boxes(hefsi_surface_t *s, gdouble scale) 

{
  gpointer tdata[HEFSI_TRAVERSE_DATA_SIZE] ;
  tdata[HEFSI_TRAVERSE_DATA_SURFACE] = s ;
  tdata[HEFSI_TRAVERSE_DATA_SCALE] = &scale ;
  
  g_node_traverse(s->tree, G_POST_ORDER, G_TRAVERSE_ALL, -1,
		  (GNodeTraverseFunc)traverse_bboxes, tdata) ;

  return 0 ;
}

static gboolean traverse_write_bbox(GNode *node, gpointer *tdata)

{
  /* hefsi_element_t *e = tdata[HEFSI_TRAVERSE_DATA_ELEMENTS] ; */
  hefsi_surface_t *s = tdata[HEFSI_TRAVERSE_DATA_SURFACE] ;
  hefsi_element_t *e ;
  FILE *f = tdata[HEFSI_TRAVERSE_DATA_FILE] ;
  gint i ;

  i = GPOINTER_TO_INT(node->data) ;  
  e = hefsi_surface_element(s,i) ;

  fprintf(f, "%lg %lg %lg %lg %lg %lg %u\n",
	  e->bbox[0], e->bbox[1], e->bbox[2], 
	  e->bbox[3], e->bbox[4], e->bbox[5], 
	  g_node_depth(node)) ;
  
  return FALSE ;
}

/** 
 * Write bounding boxes of elements of a ::hefsi_surface_t
 * 
 * @param f output file stream;
 * @param s a ::hefsi_surface_t.
 * 
 * @return 0 on success.
 */

gint hefsi_write_bounding_boxes(FILE *f, hefsi_surface_t *s)

{
  gpointer tdata[HEFSI_TRAVERSE_DATA_SIZE] ;

  tdata[HEFSI_TRAVERSE_DATA_SURFACE] = s ;
  tdata[HEFSI_TRAVERSE_DATA_FILE] = f ;
  g_node_traverse(s->tree, G_IN_ORDER, G_TRAVERSE_ALL, -1,
		  (GNodeTraverseFunc)traverse_write_bbox, tdata) ;
  
  return 0 ;
}

/** 
 * Allocate new surface and set basic parameters
 * 
 * @param dfunc surface position and derivatives evaluator;
 * @param fdata data to pass to \a dfunc;
 * @param umin minimum value of \f$u\f$ on surface;
 * @param umax maximum value of \f$u\f$ on surface;
 * @param vmin minimum value of \f$v\f$ on surface;
 * @param vmax maximum value of \f$v\f$ on surface.
 * 
 * @return newly allocated ::hefsi_surface_t
 */

hefsi_surface_t *hefsi_surface_new(gpointer eval, gpointer data,
				   gboolean derivative,
				   gdouble umin, gdouble umax,
				   gdouble vmin, gdouble vmax)

{
  hefsi_surface_t *s ;

  s = (hefsi_surface_t *)g_malloc0(sizeof(hefsi_surface_t)) ;

  hefsi_surface_umin(s) = umin ;
  hefsi_surface_umax(s) = umax ;
  hefsi_surface_vmin(s) = vmin ;
  hefsi_surface_vmax(s) = vmax ;
  hefsi_surface_eval(s) = eval ;
  hefsi_surface_data(s) = data ;
  hefsi_surface_derivative(s) = derivative ;

  s->e = g_array_new(FALSE, TRUE, sizeof(hefsi_element_t)) ;
  s->n = g_array_new(FALSE, TRUE, sizeof(hefsi_node_t)) ;
  
  return s ;
}

/** 
 * Allocate a new ::hefsi_workspace_t for calculating intersections
 * 
 * @return newly allocated workspace.
 */

hefsi_workspace_t *hefsi_workspace_new(void)

{
  hefsi_workspace_t *w ;

  w = (hefsi_workspace_t *)g_malloc0(sizeof(hefsi_workspace_t)) ;

  w->s = g_array_new(FALSE, TRUE, sizeof(hefsi_segment_t)) ;
  w->c = g_ptr_array_new() ;
  
  return w ;
}


/* static void segment_swap_nodes(hefsi_segment_t *s) */

/* { */
/*   gint i ; */
/*   gdouble tmp ; */

/*   for ( i = 0 ; i < HEFSI_SEGMENT_SIZE ; i ++ ) { */
/*     tmp = s->x1[i] ; s->x1[i] = s->x2[i] ; s->x2[i] = tmp ; */
/*   } */
  
/*   return ; */
/* } */

static gboolean segments_connect(hefsi_segment_t *s1, hefsi_segment_t *s2, 
				 gdouble tol)

{
  if ( hefsi_vector_distance(s1->x1, s2->x1) < tol ) return TRUE ;
  if ( hefsi_vector_distance(s1->x2, s2->x1) < tol ) return TRUE ;
  if ( hefsi_vector_distance(s1->x1, s2->x2) < tol ) return TRUE ;
  if ( hefsi_vector_distance(s1->x2, s2->x2) < tol ) return TRUE ;
  
  return FALSE ;
}

static GSList *insert_segment(hefsi_surface_t *s, GSList *c,
			      hefsi_workspace_t *w, gint i, gdouble tol)

{
  gint j ;
  hefsi_segment_t *segi, *segj ;
  /* GSList *jl, *next ; */

  /*list is empty, so initialize it*/
  if ( c == NULL ) {
    return g_slist_prepend(NULL, GINT_TO_POINTER(i)) ;
  }

  segi = hefsi_workspace_segment(w,i) ;
  /*does segment connect to start of list?*/
  j = GPOINTER_TO_INT(c->data) ;
  segj = hefsi_workspace_segment(w,j) ;
  if ( segments_connect(segj, segi, tol) ) {
    return g_slist_prepend(c, GINT_TO_POINTER(i)) ;    
  }

  /* jl = g_slist_last(c) ; */
  /* j = GPOINTER_TO_INT(jl->data) ; */
  /* if ( segments_connect(&(seg[j]), &(seg[i]), tol) ) { */
  /*   return g_slist_append(c, GINT_TO_POINTER(i)) ; */
  /* } */
  /* for ( jl = c ; jl != NULL ; jl = jl->next ) { */
  /*   j = GPOINTER_TO_INT(jl->data) ; */
  /*   if ( segments_connect(&(seg[j]), &(seg[i]), tol) ) { */
  /*     next = jl->next ; */
  /*     jl->next = g_slist_prepend(NULL, GINT_TO_POINTER(i)) ; */
  /*     jl->next->next = next ; */
  /*     return c ; */
  /*   } */
  /* } */

  /* g_assert_not_reached() ; */
  
  return g_slist_append(c, GINT_TO_POINTER(i)) ;
}

static gboolean area_coordinates(gdouble *x1, gdouble *x2, gdouble *x3,
				 gdouble *x,
				 gdouble *L1, gdouble *L2, gdouble *L3)

{
  gdouble L, r12[3], r23[3], r[3] ;

  hefsi_vector_diff(r12, x2, x1) ;
  hefsi_vector_diff(r23, x3, x2) ;
  hefsi_vector_cross(r, r12, r23) ;
  L = hefsi_vector_length(r) ;

  hefsi_vector_diff(r12, x, x2) ;
  hefsi_vector_diff(r23, x, x3) ;
  hefsi_vector_cross(r, r12, r23) ;
  *L1 = hefsi_vector_length(r) ;

  hefsi_vector_diff(r12, x, x3) ;
  hefsi_vector_diff(r23, x, x1) ;
  hefsi_vector_cross(r, r12, r23) ;
  *L2 = hefsi_vector_length(r) ;
  
  hefsi_vector_diff(r12, x, x1) ;
  hefsi_vector_diff(r23, x, x2) ;
  hefsi_vector_cross(r, r12, r23) ;
  *L3 = hefsi_vector_length(r) ;

  if ( fabs(*L1 + *L2 + *L3 - L) > 1e-1*L ) {
    /*intersection point lies outside triangle*/
    /* g_error("%s: area coordinates out of range", __FUNCTION__) ; */
    return FALSE ;
  }

  *L1 /= L ; *L2 /= L ; *L3 /= L ; 
  
  return TRUE ;
}

static void check_intersection(hefsi_surface_t *s1,
			       hefsi_element_t *e1, gint t1,
			       hefsi_surface_t *s2,
			       hefsi_element_t *e2, gint t2,
			       hefsi_workspace_t *w, gdouble tol,
			       gboolean *inter)

{
  gint idx[] = {0, 1, 2, 0, 2, 3}, *i1, *i2, cp, nseg ;
  gdouble ip1[3], ip2[3], L[3] ;
  gdouble *u11, *u12, *u13, *u21, *u22, *u23 ;
  gdouble *x11, *x12, *x13, *x21, *x22, *x23 ;
  hefsi_segment_t seg ;

  i1 = &(idx[3*t1]) ; i2 = &(idx[3*t2]) ;

  x11 = hefsi_surface_x(s1, e1->i[i1[0]]) ;
  x12 = hefsi_surface_x(s1, e1->i[i1[1]]) ;
  x13 = hefsi_surface_x(s1, e1->i[i1[2]]) ;
  x21 = hefsi_surface_x(s2, e2->i[i2[0]]) ;
  x22 = hefsi_surface_x(s2, e2->i[i2[1]]) ;
  x23 = hefsi_surface_x(s2, e2->i[i2[2]]) ;  

  if ( tri_tri_intersect_with_isectline(x11, x12, x13, x21, x22, x23,
					&cp, ip1, ip2) == 0  ||
       cp == 1 ) return ;

  /*ignore segments shorter than the tolerance/2*/
  if ( hefsi_vector_distance(ip1,ip2) < 0.5*tol ) return ;
  
  u11 = &(hefsi_surface_u(s1, e1->i[i1[0]])) ;
  u12 = &(hefsi_surface_u(s1, e1->i[i1[1]])) ;
  u13 = &(hefsi_surface_u(s1, e1->i[i1[2]])) ;
  u21 = &(hefsi_surface_u(s2, e2->i[i2[0]])) ;
  u22 = &(hefsi_surface_u(s2, e2->i[i2[1]])) ;
  u23 = &(hefsi_surface_u(s2, e2->i[i2[2]])) ;

  seg.x1[0] = ip1[0] ; seg.x1[1] = ip1[1] ; seg.x1[2] = ip1[2] ;
  seg.x2[0] = ip2[0] ; seg.x2[1] = ip2[1] ; seg.x2[2] = ip2[2] ; 
  seg.i1 = t1 ; seg.i2 = t2 ; 

  if ( !area_coordinates(x11, x12, x13, ip1, &(L[0]), &(L[1]), &(L[2])) )
    return ;
  seg.x1[3] = L[0]*u11[0] + L[1]*u12[0] + L[2]*u13[0] ;
  seg.x1[4] = L[0]*u11[1] + L[1]*u12[1] + L[2]*u13[1] ;

  if ( !area_coordinates(x11, x12, x13, ip2, &(L[0]), &(L[1]), &(L[2])) )
    return ;
  seg.x2[3] = L[0]*u11[0] + L[1]*u12[0] + L[2]*u13[0] ;
  seg.x2[4] = L[0]*u11[1] + L[1]*u12[1] + L[2]*u13[1] ;

  if ( !area_coordinates(x21, x22, x23, ip1, &(L[0]), &(L[1]), &(L[2])) )
    return ;
  seg.x1[5] = L[0]*u21[0] + L[1]*u22[0] + L[2]*u23[0] ;
  seg.x1[6] = L[0]*u21[1] + L[1]*u22[1] + L[2]*u23[1] ;

  if ( !area_coordinates(x21, x22, x23, ip2, &(L[0]), &(L[1]), &(L[2])) )
    return ;
  seg.x2[5] = L[0]*u21[0] + L[1]*u22[0] + L[2]*u23[0] ;
  seg.x2[6] = L[0]*u21[1] + L[1]*u22[1] + L[2]*u23[1] ;

  /*add the segment to list in the workspace*/
  nseg = w->s->len ;
  g_array_append_val(w->s, seg) ;
  
  e1->s = insert_segment(s1, e1->s, w, nseg, tol) ;
  e2->s = insert_segment(s2, e2->s, w, nseg, tol) ;

  *inter = TRUE ;

  return ;
}

static gboolean locate_intersection(GNode *node2, gpointer *data)

{
  hefsi_surface_t *s1  = data[0] ; 
  hefsi_surface_t *s2  = data[1] ; 
  GNode *node1         = data[2] ;
  hefsi_workspace_t *w = data[3] ;
  gdouble tol         = *((gdouble *)data[4]) ;
  gint *ninter        = data[5] ;
  gboolean inter ;
  GNode *c ;
  hefsi_element_t *e1, *e2 ;
  gint i1, i2 ;

  i1 = GPOINTER_TO_INT(node1->data) ;
  e1 = hefsi_surface_element(s1,i1) ;
  i2 = GPOINTER_TO_INT(node2->data) ;
  e2 = hefsi_surface_element(s2,i2) ;

  if ( hefsi_element_bboxes_disjoint(e1,e2) ) {
    return TRUE ;
  }

  if ( G_NODE_IS_LEAF(node2) ) {
    /*boxes overlap: check for triangle intersections*/
    inter = FALSE ;
    check_intersection(s1, e1, 0, s2, e2, 0, w, tol, &inter) ;
    check_intersection(s1, e1, 0, s2, e2, 1, w, tol, &inter) ;
    check_intersection(s1, e1, 1, s2, e2, 0, w, tol, &inter) ;
    check_intersection(s1, e1, 1, s2, e2, 1, w, tol, &inter) ;

    if ( inter ) (*ninter) ++ ;
    return inter ;
  }

  for ( c = node2->children ; c != NULL ; c = c->next ) {
    locate_intersection(c, data) ; 
  }
  
  return FALSE ;
}

static gboolean traverse_intersections(GNode *node, gpointer *data)

{
  hefsi_surface_t *s2 = data[1] ; 

  /*node is a leaf node of s1*/
  data[2] = node ;

  locate_intersection(s2->tree, data) ;
  
  return FALSE ;
}

static void intersection_boxes(hefsi_surface_t *s1, hefsi_surface_t *s2,
			       gdouble tol, gint *ninter, hefsi_workspace_t *w)

{
  gpointer data[] = {s1, s2, NULL, w, &tol, ninter} ;

  g_node_traverse(s1->tree, G_IN_ORDER, G_TRAVERSE_LEAVES, -1,
		  (GNodeTraverseFunc)traverse_intersections, data) ;  
  return ;
}

static gboolean connect_curves(hefsi_workspace_t *w,
			       GSList **c1, GSList **c2,
			       gdouble tol)

{
  gint first1, last1, first2, last2 ;
  GSList *j ;
  
  first1 = GPOINTER_TO_INT((*c1)->data) ;
  first2 = GPOINTER_TO_INT((*c2)->data) ;

  j = g_slist_last(*c1) ;
  last1 = GPOINTER_TO_INT(j->data) ;
  j = g_slist_last(*c2) ;
  last2 = GPOINTER_TO_INT(j->data) ;

  if ( segments_connect(hefsi_workspace_segment(w,first1),
			hefsi_workspace_segment(w,first2),
			tol) ) {
    *c2 = g_slist_reverse(*c2) ;
    *c1 = g_slist_concat(*c2, *c1) ;
    return TRUE ;
  }

  if ( segments_connect(hefsi_workspace_segment(w,first1),
			hefsi_workspace_segment(w, last2),
			tol) ) {
    *c1 = g_slist_concat(*c2, *c1) ;
    return TRUE ;
  }

  if ( segments_connect(hefsi_workspace_segment(w, last1),
			hefsi_workspace_segment(w,first2),
			tol) ) {
    *c1 = g_slist_concat(*c1, *c2) ;
    return TRUE ;
  }

  if ( segments_connect(hefsi_workspace_segment(w,last1),
			hefsi_workspace_segment(w,last2),
			tol) ) {
    *c2 = g_slist_reverse(*c2) ;
    *c1 = g_slist_concat(*c1, *c2) ;
    return TRUE ;
  }
  
  return FALSE ;
}

static gboolean traverse_curves(GNode *node, gpointer *data)

{
  hefsi_surface_t *s1 = data[0] ; 
  hefsi_workspace_t *w = data[1] ;
  gdouble tol         = *((gdouble *)data[2]) ;
  gint i, j, c, first, last ;
  GSList *il, *curve ;
  hefsi_element_t *e ;
  gboolean connected ;
  
  i = GPOINTER_TO_INT(node->data) ;
  e = hefsi_surface_element(s1,i) ;
  if ( e->s == NULL ) return FALSE ;
  
  for ( il = e->s ; il != NULL ; il = il->next ) {
    j = GPOINTER_TO_INT(il->data) ;
    connected = FALSE ;
    for ( c = 0 ; c < w->c->len ; c ++ ) {
      curve = hefsi_workspace_curve(w, c) ;
      first = GPOINTER_TO_INT(curve->data) ;
      last  = GPOINTER_TO_INT(g_slist_last(curve)->data) ;
      if ( segments_connect(hefsi_workspace_segment(w,    j),
			    hefsi_workspace_segment(w,first),
			    tol) ) {
	hefsi_workspace_curve(w, c) = g_slist_prepend(curve, il->data) ;
	connected = TRUE ;
      }
      if ( segments_connect(hefsi_workspace_segment(w,   j),
			    hefsi_workspace_segment(w,last),
			    tol) ) {
	hefsi_workspace_curve(w, c) = g_slist_append(curve, il->data) ;
	connected = TRUE ;
      }
    }

    if ( !connected ) {
      /*start a new list of segments*/
      g_ptr_array_add(w->c, g_slist_prepend(NULL, il->data)) ;
    }
  }
  
  return FALSE ;
}

static void sort_curves(hefsi_surface_t *s1, gdouble tol, gint kmax,
		 hefsi_workspace_t *w)

{
  gpointer data[] = {s1, w, &tol} ;
  gint i, j, k, ncold ;
  
  /*traverse the leaves of s1, adding segments to curves*/
  g_node_traverse(s1->tree, G_IN_ORDER, G_TRAVERSE_LEAVES, -1,
		  (GNodeTraverseFunc)traverse_curves, data) ;
  
  /*join curves if possible*/
  ncold = 0 ;
  for ( k = 0 ; (k < kmax) && (w->c->len != ncold) ; k ++ ) {
    /* fprintf(stderr, "nc: %d\n", w->c->len) ; */
    ncold = w->c->len ;
    for ( i = 0 ; i < hefsi_workspace_curve_number(w) ; i ++ ) {
      if ( hefsi_workspace_curve(w, i) != NULL ) {
	for ( j = i + 1 ; j < w->c->len ; j ++ ) {
	  if ( hefsi_workspace_curve(w, j) != NULL ) {
	  if ( connect_curves(w,
			      (GSList **)(&(hefsi_workspace_curve(w,i))),
			      (GSList **)(&(hefsi_workspace_curve(w,j))),
			      tol) ) {
	      g_ptr_array_steal_index_fast(w->c, j) ;
	      j -- ;
	    }
	  }
	}
      }
    }
  }

  return ;
}

static void surface_nearest_point(hefsi_surface_t *s, gdouble *x,
				  gdouble *y, gdouble *yu, gdouble *yv,
				  gdouble *uv, gdouble tol)

{
  gdouble dx[3], A[4], b[2], dv[2], det, rxf ;
  gint i ;

  /* if ( uv[0] < 0 ) { */
  /*   uv[0] = 0 ; */
  /*   hefsi_surface_evaluate(s, uv[0], uv[1], y, yu, yv, NULL) ; */
  /*   return ; */
  /* } */
  /* if ( uv[0] > 1 ) { */
  /*   uv[0] = 1 ; */
  /*   hefsi_surface_evaluate(s, uv[0], uv[1], y, yu, yv, NULL) ; */
  /*   return ; */
  /* } */

  hefsi_surface_evaluate(s, uv[0], uv[1], y, yu, yv, NULL) ;
  /*point lies in surface*/
  if ( hefsi_vector_distance(x, y) < tol ) return ;
  for ( i = 0 ; i < 8 ; i ++ ) {
    hefsi_vector_diff(dx, x, y) ;
  
    /*set up equations*/
    A[0] = hefsi_vector_dot(yu, yu) ; 
    A[1] = hefsi_vector_dot(yu, yv) ; 
    A[2] = hefsi_vector_dot(yv, yu) ; 
    A[3] = hefsi_vector_dot(yv, yv) ; 

    if ( isnan(A[0]) ) {
      g_error("%s: NaN in matrix", __FUNCTION__) ;
    }
    
    b[0] = hefsi_vector_dot(dx, yu) ;
    b[1] = hefsi_vector_dot(dx, yv) ;
    hefsi_invert2x2(A, A, &det) ;
    
    dv[0] = A[0]*b[0] + A[1]*b[1] ; 
    dv[1] = A[2]*b[0] + A[3]*b[1] ; 

    if ( uv[0] + dv[0] > 0 && uv[0] + dv[0] < 1 ) {
      rxf = 1.0 ;
    } else {
      if ( uv[0] + dv[0] < 0 ) {
	rxf = -uv[0]/dv[0]*0.5 ;
      } else {
	rxf = (1.0-uv[0])/dv[0]*0.5 ;
      }
    }
    rxf = MIN(ABS(det), ABS(rxf)) ;
    if ( rxf > 1 ) rxf = 0.5 ;
    uv[0] += rxf*dv[0] ; uv[1] += rxf*dv[1] ; 

    /*this should really be done better, but it does handle the
      calculation okay near a singular point*/
    while ( uv[1] > 1 ) uv[1] -= 1 ;
    while ( uv[1] < 0 ) uv[1] += 1 ;
    
    hefsi_surface_evaluate(s, uv[0], uv[1], y, yu, yv, NULL) ;
  }  

  return ;
}

static void plane_intersection_line(gdouble *y1, gdouble *n1,
				    gdouble *y2, gdouble *n2,
				    gdouble *x0, gdouble *s)

/*
 * parametric equation for line of intersection of two planes
 *
 * x = x0 + s*t is intersection of two planes (x-y1).n1 = 0, (x-y2).n2 = 0
 */

{
  gdouble tol = 1e-3, A[4], b[2], det ;
  gint idx[] = {0, 1, 2, 0, 1}, i ;
  
  x0[0] = x0[1] = x0[2] = 0.0 ;
  
  hefsi_vector_cross(s, n1, n2) ;

  b[0] = hefsi_vector_dot(y1, n1) ;
  b[1] = hefsi_vector_dot(y2, n2) ;

  /*set one of the entries in x0 to zero and solve for the other two*/
  for ( i = 0 ; i < 3 ; i ++ ) {
    if ( fabs(n1[idx[i+0]] - n2[idx[i+0]]) > tol &&
	 fabs(n1[idx[i+1]] - n2[idx[i+1]]) > tol ) {
      A[0] = n1[idx[i+0]] ; A[1] = n1[idx[i+1]] ;
      A[2] = n2[idx[i+0]] ; A[3] = n2[idx[i+1]] ;
      hefsi_invert2x2(A, A, &det) ;
      x0[idx[i+0]] = A[0]*b[0] + A[1]*b[1] ; 
      x0[idx[i+1]] = A[2]*b[0] + A[3]*b[1] ;
      return ;
    }
  }

  /*if we get to here, there is an unhandled special case*/
  /* g_assert_not_reached() ; */

  return ;
}

/* static void check_plane_intersection(gdouble *y1, gdouble *n1, */
/* 				     gdouble *y2, gdouble *n2, */
/* 				     gdouble *x0, gdouble *s) */

/* { */
/*   gdouble t, x[3], dx[3], err ; */

/*   t = 0.3 ; */

/*   x[0] = x0[0] + s[0]*t ; */
/*   x[1] = x0[1] + s[1]*t ; */
/*   x[2] = x0[2] + s[2]*t ; */

/*   hefsi_vector_diff(dx, x, y1) ;   */
/*   err = hefsi_vector_dot(dx, n1) ; */
/*   g_assert(fabs(err) < 1e-9) ; */

/*   hefsi_vector_diff(dx, x, y2) ;   */
/*   err = hefsi_vector_dot(dx, n2) ; */
/*   g_assert(fabs(err) < 1e-9) ; */
  
/*   return ; */
/* } */

static void line_nearest_point(gdouble *x0, gdouble *s, gdouble *y,
			       gdouble *x)

/*
 * on exit x is nearest point to y on line x0 + s*t
 */
  
{
  gdouble t, dx[3] ;

  hefsi_vector_diff(dx, y, x0) ;
  t = hefsi_vector_dot(dx, s)/hefsi_vector_dot(s, s) ;

  x[0] = x0[0] + s[0]*t ;
  x[1] = x0[1] + s[1]*t ;
  x[2] = x0[2] + s[2]*t ;
  
  return ;
}

/** 
 * Refine an estimated point of intersection between two surfaces.
 * 
 * Input is estimated intersection point \a x and estimated
 * \f$(u,v)\f$ coordinates on intersecting surfaces. On exit, \a x is
 * refined estimate of intersection location and \f$(u,v)\f$
 * coordinates are refined estimates of parametric location of \a x.
 * 
 * @param s1 a ::hefsi_surface_t;
 * @param uv1 estimated parametric coordinates of \a x on \a s1;
 * @param s2 a ::hefsi_surface_t;
 * @param uv2 estimated parametric coordinates of \a x on \a s2;
 * @param x estimated position of intersection point;
 * @param tol convergence tolerance for refinement;
 * @param imax maximum number of iterations.
 * 
 * @return 0 on success.
 */

gint hefsi_refine_point(hefsi_surface_t *s1, gdouble *uv1,
			hefsi_surface_t *s2, gdouble *uv2,
			double *x, gdouble tol, gint imax)

{
  gdouble y1[3], y2[3], yu1[3], yv1[3], yu2[3], yv2[3] ;
  gdouble dx1[3], dx2[3],  A[9], Ai[9], b[3], det, d1, d2, d, x0[3], s[3] ;
  gint i ;

  for ( i = 0 ; i < imax ; i ++ ) {
    /*check for x coinciding with y1 or y2*/
    surface_nearest_point(s1, x, y1, yu1, yv1, uv1, tol) ;
    if ( isnan(uv2[0]) )
      g_error("%s: NaN in uv2", __FUNCTION__) ;
    surface_nearest_point(s2, x, y2, yu2, yv2, uv2, tol) ;
    if ( isnan(y2[0]) )
      g_error("%s: NaN in y2", __FUNCTION__) ;

    d1 = hefsi_vector_distance(x , y1) ;
    d2 = hefsi_vector_distance(x , y2) ;
    d  = hefsi_vector_distance(y1, y2) ;

    if ( d1 < tol && d2 < tol && d < tol ) {
      /*converged*/
      x[0] = 0.5*(y1[0]+y2[0]) ;
      x[1] = 0.5*(y1[1]+y2[1]) ;
      x[2] = 0.5*(y1[2]+y2[2]) ;
      g_assert(!isnan(x[0])) ;
      return 0 ;
    }
  
    /*check for colinear x, y1, y2 here?*/
  
    /*generate normals to tangent planes*/
    hefsi_vector_cross(&(A[0]), yu1, yv1) ;
    hefsi_vector_cross(&(A[3]), yu2, yv2) ;
  
    if ( d1 < tol || d2 < tol ) {
      /*
       * estimated point lies on one of the surfaces: project surface
       * points to line of intersection of tangent planes
       */
      plane_intersection_line(y1, &(A[0]), y2, &(A[3]), x0, s) ;
      /* check_plane_intersection(y1, &(A[0]), y2, &(A[3]), x0, s) ; */
      line_nearest_point(x0, s, y1, y1) ; 
      line_nearest_point(x0, s, y2, y2) ;
      x[0] = 0.5*(y1[0] + y2[0]) ;
      x[1] = 0.5*(y1[1] + y2[1]) ;
      x[2] = 0.5*(y1[2] + y2[2]) ;
      if ( isnan(x[0]) ) {
	g_error("%s: NaN in intersection", __FUNCTION__) ;
      }
    } else {
      /*normal to plane containing x, y1, and y2*/
      hefsi_vector_diff(dx1, x, y1) ;
      hefsi_vector_diff(dx2, x, y2) ;
      hefsi_vector_cross(&(A[6]), dx1, dx2) ;

      /*intersection of three planes: form right hand side*/
      b[0] = hefsi_vector_dot(&(A[0]), y1) ;
      b[1] = hefsi_vector_dot(&(A[3]), y2) ;
      b[2] = hefsi_vector_dot(&(A[6]), x) ;
      hefsi_invert3x3(Ai, A, &det) ;
    
      /* x[0] = Ai[0]*b[0] + Ai[1]*b[1] + Ai[2]*b[2] ;  */
      /* x[1] = Ai[3]*b[0] + Ai[4]*b[1] + Ai[5]*b[2] ;  */
      /* x[2] = Ai[6]*b[0] + Ai[7]*b[1] + Ai[8]*b[2] ;  */
      hefsi_matvecmul3x3(x,Ai,b) ;
      g_assert(!isnan(x[0])) ;
    }
  }
  if ( isnan(x[0]) || isnan(x[1]) || isnan(x[2]) ) {
    fprintf(stderr, "%s: NaN\n", __FUNCTION__) ;
  }
  
  return 0 ;
}

static void refine_curves(hefsi_surface_t *s1, hefsi_surface_t *s2, gdouble tol,
			  hefsi_workspace_t *w)

{
  gint c, j, imax ;
  GSList *i, *curve ;
  gdouble *x ;

  imax = 8 ;
  for ( c = 0 ; c < w->c->len ; c ++ ) {
    curve = hefsi_workspace_curve(w, c) ;
    for ( i = curve ; i != NULL ; i = i->next ) {
      j = GPOINTER_TO_INT(i->data) ;
      x = hefsi_workspace_segment(w,j)->x1 ;
      hefsi_refine_point(s1, &(x[3]), s2, &(x[5]), x, tol, imax) ;
      x = hefsi_workspace_segment(w,j)->x2 ;
      hefsi_refine_point(s1, &(x[3]), s2, &(x[5]), x, tol, imax) ;
    }
  }
  
  return ;
}

/** 
 * Calculate the intersection(s) of two surfaces.
 * 
 * @param s1 a ::hefsi_surface_t;
 * @param s2 another ::hefsi_surface_t;
 * @param tol tolerance for distance of estimated intersection points 
 * from surfaces;
 * @param w ::hefsi_workspace_t, on exit contains sorted intersection curves.
 * 
 * @return 0 on success.
 */

gint hefsi_surface_intersections(hefsi_surface_t *s1, hefsi_surface_t *s2, 
				 gdouble tol, hefsi_workspace_t *w)

{
  gint i, kmax, ninter ;
  hefsi_element_t *e ;

  for ( i = 0 ; i < hefsi_workspace_curve_number(w) ; i ++ ) {
    g_slist_free((GSList *)(hefsi_workspace_curve(w, i))) ;
    hefsi_workspace_curve(w, i) = NULL ;
  }
  for ( i = 0 ; i < hefsi_workspace_segment_number(w) ; i ++ ) {
    memset(hefsi_workspace_segment(w,i),0,sizeof(hefsi_segment_t)) ;
  }

  for ( i = 0 ; i < s1->e->len ; i ++ ) {
    e = hefsi_surface_element(s1,i) ;
    g_slist_free(e->s) ;
    e->s = NULL ;
  }

  for ( i = 0 ; i < s2->e->len ; i ++ ) {
    e = hefsi_surface_element(s2,i) ;
    g_slist_free(e->s) ;
    e->s = NULL ;
  }
  
  g_array_set_size(w->s, 0) ;
  g_ptr_array_set_size(w->c, 0) ;

  ninter = 0 ;
  intersection_boxes(s1, s2, tol, &ninter, w) ;

  kmax = 8 ;
  
  sort_curves(s1, tol, kmax, w) ;
  refine_curves(s1, s2, tol, w) ;

  return 0 ;
}

static void eval_c1(hefsi_surface_t *s, gdouble u, gdouble v, gdouble *buf)

{
  hefsi_surface_evaluator_dfunc_t func = hefsi_surface_eval(s) ;  
  func(u, v, &(buf[3*0]),  &(buf[3*1]), &(buf[3*2]), s->data) ;

  return ;
}

static void eval_c0(hefsi_surface_t *s, gdouble u, gdouble v, gdouble *buf)

{
  gdouble buv[3], d, sgn ;

  d = 1e-6 ;
  hefsi_surface_evaluator_func_t func = hefsi_surface_eval(s) ;
  func(u, v, &(buf[3*0]), s->data) ;

  sgn = 0 ;
  if ( u > s->umin + d ) {
    func(u-d, v, buv, s->data) ;
    sgn = 1 ;
  } else {
    func(u+d, v, buv, s->data) ;
    sgn = -1 ;
  }

  buf[3] = sgn*(buf[0] - buv[0])/d ;
  buf[4] = sgn*(buf[1] - buv[1])/d ;
  buf[5] = sgn*(buf[2] - buv[2])/d ;

  sgn = 0 ;
  if ( v > s->vmin + d ) {
    func(u, v-d, buv, s->data) ;
    sgn = 1 ;
  } else {
    func(u, v+d, buv, s->data) ;
    sgn = -1 ;
  }

  buf[6] = sgn*(buf[0] - buv[0])/d ;
  buf[7] = sgn*(buf[1] - buv[1])/d ;
  buf[8] = sgn*(buf[2] - buv[2])/d ;
  
  return ;
}

/** 
 * Evaluate surface point and/or derivatives and/or normal. 
 *
 * If any of \a x, \a xu, \a xv, or \a n are NULL, they are not
 * evaluated and no error is reported. If the evaluator for \a s does
 * not return derivatives, any required derivatives are estimated by
 * numerical differencing.
 * 
 * @param s a ::hefsi_surface_t;
 * @param u parametric coordinate on \a s;
 * @param v parametric coordinate on \a s;
 * @param x on exit contains point \f$\mathbf{x}(u,v)\f$;
 * @param xu on exit contains derivative \f$\partial\mathbf{x}\partial u\f$;
 * @param xv on exit contains derivative \f$\partial\mathbf{x}\partial v\f$;
 * @param n on exit contains surface normal at \f$\mathbf{x}(u,v)\f$.
 * 
 * @return 0 on success.
 */

gint hefsi_surface_evaluate(hefsi_surface_t *s, gdouble u, gdouble v,
			    gdouble *x, gdouble *xu, gdouble *xv,
			    gdouble *n)

{
  gdouble buf[9], del, ztol ;

  del = 1e-3 ; ztol = 1e-6 ;
  if (  hefsi_surface_derivative(s) ) {
    eval_c1(s, u, v, buf) ;
  } else {
    eval_c0(s, u, v, buf) ;
  }

  if ( x  != NULL ) {
    x[0] = buf[3*0+0] ; x[1] = buf[3*0+1] ; x[2] = buf[3*0+2] ;
  }
  if ( xu != NULL ) {
    xu[0] = buf[3*1+0] ; xu[1] = buf[3*1+1] ; xu[2] = buf[3*1+2] ;
  }
  if ( xv != NULL ) {
    xv[0] = buf[3*2+0] ; xv[1] = buf[3*2+1] ; xv[2] = buf[3*2+2] ;
  }

  if ( n == NULL ) return 0 ;

  if ( (ABS(buf[3*1+0]) <= ztol && ABS(buf[3*1+1]) <= ztol &&
	ABS(buf[3*1+2]) <= ztol) ||
       (ABS(buf[3*2+0]) <= ztol && ABS(buf[3*2+1]) <= ztol &&
	ABS(buf[3*2+2]) <= ztol) ) {
    /*singular point: wiggle u and v to estimate the normal*/
    if ( u < s->umax - del ) {
      u += del ;
    } else {
      if ( u > s->umin + del ) {
	u -= del ;
      }
    }
    if ( v < s->vmax - del ) {
      v += del ;
    } else {
      if ( v > s->vmin + del ) {
	v -= del ;
      }
    }
    
    if (  hefsi_surface_derivative(s) ) {
      eval_c1(s, u, v, buf) ;
    } else {
      eval_c0(s, u, v, buf) ;
    }
  }
  
  make_normal(n, &(buf[3*1]), &(buf[3*2])) ;
  
 return 0 ;
}

static gint check_segment(hefsi_surface_t *s1, hefsi_surface_t *s2,
			  hefsi_segment_t *s, gdouble tol, gdouble *emax)

{
  gint n ;
  gdouble x1[3], x2[3], r ;

  n = 0 ;
  hefsi_surface_evaluate(s1, s->x1[3], s->x1[4], x1, NULL, NULL, NULL) ;
  hefsi_surface_evaluate(s2, s->x1[5], s->x1[6], x2, NULL, NULL, NULL) ;

  if ( ( r = hefsi_vector_distance(x1, x2)) > tol ) n ++ ;
  *emax = MAX(r, *emax) ;
  
  hefsi_surface_evaluate(s1, s->x2[3], s->x2[4], x1, NULL, NULL, NULL) ;
  hefsi_surface_evaluate(s2, s->x2[5], s->x2[6], x2, NULL, NULL, NULL) ;

  if ( (r = hefsi_vector_distance(x1, x2)) > tol ) n ++ ;
  *emax = MAX(r, *emax) ;
  
  return n ;
}

/** 
 * Check convergence of calculated surface intersection.
 * 
 * @param s1 a ::hefsi_surface_t;
 * @param s2 another ::hefsi_surface_t;
 * @param w a ::hefsi_workspace_t used to calculate the intersections of 
 * \a s1 and \a s2 (using ::hefsi_surface_intersections);
 * @param tol tolerance for intersection;
 * @param emax maximum error in calculated intersections.
 * 
 * @return number of calculated intersection nodes failing
 * intersection check.
 */

gint hefsi_surface_intersection_check(hefsi_surface_t *s1, hefsi_surface_t *s2,
				      hefsi_workspace_t *w, gdouble tol,
				      gdouble *emax)

{
  gint n, c, i  ;
  GSList *il ;
  hefsi_segment_t *s ;

  n = 0 ; *emax = 0.0 ;
  for ( c = 0 ; c < hefsi_workspace_curve_number(w) ; c ++ ) {
    for ( il = hefsi_workspace_curve(w, c) ; il != NULL ; il = il->next ) {
      i = GPOINTER_TO_INT(il->data) ;
      s = hefsi_workspace_segment(w,i) ;
      n += check_segment(s1, s2, s, tol, emax) ;
    }
  }
  
  /*return number of failed nodes*/
  return n ;
}

