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

#ifndef HEFSI_H_INCLUDED
#define HEFSI_H_INCLUDED

#include <glib.h> 

#define hefsi_vector_distance(_x1,_x2)			  \
  (sqrt(((_x1)[0]-(_x2)[0])*((_x1)[0]-(_x2)[0]) +  \
	((_x1)[1]-(_x2)[1])*((_x1)[1]-(_x2)[1]) +  \
	((_x1)[2]-(_x2)[2])*((_x1)[2]-(_x2)[2])))
/*y = x1 - x2*/
#define hefsi_vector_diff(_y,_x1,_x2)			\
  do {							\
    (_y)[0] = (_x1)[0] - (_x2)[0] ;			\
    (_y)[1] = (_x1)[1] - (_x2)[1] ;			\
    (_y)[2] = (_x1)[2] - (_x2)[2] ;			\
  } while (0) 
#define hefsi_vector_dot(_x1,_x2)					\
  ((_x1)[0]*(_x2)[0]+(_x1)[1]*(_x2)[1]+(_x1)[2]*(_x2)[2])
#define hefsi_vector_cross(_y,_x1,_x2)					\
  do {									\
    (_y)[0] = (_x1)[1]*(_x2)[2] - (_x1)[2]*(_x2)[1] ;			\
    (_y)[1] = (_x1)[2]*(_x2)[0] - (_x1)[0]*(_x2)[2] ;			\
    (_y)[2] = (_x1)[0]*(_x2)[1] - (_x1)[1]*(_x2)[0] ;			\
  } while (0)
#define hefsi_vector_length(_x)					\
  (sqrt((_x[0])*(_x[0]) + (_x[1])*(_x[1]) + (_x[2])*(_x[2])))

#define hefsi_matvecmul3x3(_x,_A,_b)					\
  do {									\
    (_x)[0] = hefsi_vector_dot(&((_A)[0]),(_b))	;			\
    (_x)[1] = hefsi_vector_dot(&((_A)[3]),(_b))	;			\
    (_x)[2] = hefsi_vector_dot(&((_A)[6]),(_b))	;			\
  } while (0)


/*
 * inverse of general 3x3 matrix
 * 
 * cannot be done in-place (Ai must be distinct from A)
 */

#define hefsi_invert3x3(_Ai,_A,_det)					\
  do {									\
    (*(_det)) = ((_A)[0]*((_A)[8]*(_A)[4]-(_A)[7]*(_A)[5]) -		\
		 (_A)[3]*((_A)[8]*(_A)[1]-(_A)[7]*(_A)[2]) +		\
		 (_A)[6]*((_A)[5]*(_A)[1]-(_A)[4]*(_A)[2])) ;		\
    (_Ai)[0] =  ((_A)[8]*(_A)[4] - (_A)[7]*(_A)[5])/(*(_det)) ;		\
    (_Ai)[1] = -((_A)[8]*(_A)[1] - (_A)[7]*(_A)[2])/(*(_det)) ;		\
    (_Ai)[2] =  ((_A)[5]*(_A)[1] - (_A)[4]*(_A)[2])/(*(_det)) ;		\
									\
    (_Ai)[3] = -((_A)[8]*(_A)[3] - (_A)[6]*(_A)[5])/(*(_det)) ;		\
    (_Ai)[4] =  ((_A)[8]*(_A)[0] - (_A)[6]*(_A)[2])/(*(_det)) ;		\
    (_Ai)[5] = -((_A)[5]*(_A)[0] - (_A)[3]*(_A)[2])/(*(_det)) ;		\
									\
    (_Ai)[6] =  ((_A)[7]*(_A)[3] - (_A)[6]*(_A)[4])/(*(_det)) ;		\
    (_Ai)[7] = -((_A)[7]*(_A)[0] - (_A)[6]*(_A)[1])/(*(_det)) ;		\
    (_Ai)[8] =  ((_A)[4]*(_A)[0] - (_A)[3]*(_A)[1])/(*(_det)) ;		\
  } while (0)

/*
 * inverse of general 2x2 matrix
 *
 * inv [a b; c d] = [d -b; -c a]/det
 *
 * det = a*d - b*c
 *
 * can be done in-place (Ai == A)
 */

#define hefsi_invert2x2(_Ai,_A,_det)				\
  do {								\
  gdouble _a, _b, _c, _d ;					\
  _a = (_A)[0] ; _b = (_A)[1] ; _c = (_A)[2] ; _d = (_A)[3] ;	\
  (*(_det)) = _a*_d - _b*_c ;					\
  (_Ai)[0] = _d/(*(_det)) ; (_Ai)[1] = -(_b)/(*(_det)) ;	\
  (_Ai)[2] = -_c/(*(_det)) ; (_Ai)[3] = _a/(*(_det)) ;		\
  } while (0)

#define HEFSI_NODE_SIZE 14

typedef struct _hefsi_node_t hefsi_node_t ;
struct _hefsi_node_t {
  gdouble x[HEFSI_NODE_SIZE] ;
} ;

#define hefsi_node_u(_node)   ((_node)->x[0])
#define hefsi_node_v(_node)   ((_node)->x[1])
#define hefsi_node_x(_node)   (&((_node)->x[2]))
#define hefsi_node_xu(_node)  (&((_node)->x[5]))
#define hefsi_node_xv(_node)  (&((_node)->x[8]))
#define hefsi_node_n(_node)   (&((_node)->x[11]))

#define HEFSI_TRAVERSE_DATA_SIZE     3
#define HEFSI_TRAVERSE_DATA_SURFACE  0
#define HEFSI_TRAVERSE_DATA_FILE     1
#define HEFSI_TRAVERSE_DATA_SCALE    2

#define HEFSI_SEGMENT_SIZE 7
    
typedef gint (*hefsi_surface_evaluator_dfunc_t)(gdouble u, gdouble v,
						gdouble *x,
						gdouble *xu,
						gdouble *xv,
						gpointer data) ;
typedef gint (*hefsi_surface_evaluator_func_t)(gdouble u, gdouble v,
					       gdouble *x,
					       gpointer data) ;

typedef struct _hefsi_element_t hefsi_element_t ;
struct _hefsi_element_t {
  gint i[4] ; /*node indices*/
  gdouble bbox[6] ; /*xmin xmax ymin ymax zmin zmax*/
  GSList *s ; /*intersection segments*/
} ;

typedef struct _hefsi_segment_t hefsi_segment_t ;
struct _hefsi_segment_t {
  gint i1, i2 ; /*triangle index on elements*/
  gdouble
  x1[HEFSI_SEGMENT_SIZE],
    x2[HEFSI_SEGMENT_SIZE] ; /*x y z u1 v1 u2 v2*/ 
} ;

typedef struct _hefsi_workspace_t hefsi_workspace_t ;
struct _hefsi_workspace_t {
  GArray *s ;    /*segments*/
  GPtrArray *c ; /*linked lists of intersection curves*/
} ;

#define hefsi_workspace_segment_number(_w)     ((_w)->s->len) 
#define hefsi_workspace_segment(_w,_i)		\
  (&(g_array_index((_w)->s,hefsi_segment_t,(_i))))
#define hefsi_workspace_curve_number(_w)       ((_w)->c->len)
#define hefsi_workspace_curve(_w,_i) g_ptr_array_index((_w)->c,(_i))

#define hefsi_element_bboxes_disjoint(_e1,_e2)	\
  (((_e1)->bbox[1] < (_e2)->bbox[0]) ||((_e2)->bbox[1] < (_e1)->bbox[0]) || \
   ((_e1)->bbox[3] < (_e2)->bbox[2]) ||((_e2)->bbox[3] < (_e1)->bbox[2]) || \
   ((_e1)->bbox[5] < (_e2)->bbox[4]) ||((_e2)->bbox[5] < (_e1)->bbox[4]))

typedef struct _hefsi_surface_t hefsi_surface_t ;
struct _hefsi_surface_t {
  GNode *tree ;
  GArray *e, *n ;
  gdouble umin, umax, vmin, vmax ;
  gboolean derivative ;
  gpointer eval, data ;
} ;

#define hefsi_surface_umin(_s) ((_s)->umin) 
#define hefsi_surface_umax(_s) ((_s)->umax) 
#define hefsi_surface_vmin(_s) ((_s)->vmin) 
#define hefsi_surface_vmax(_s) ((_s)->vmax) 
#define hefsi_surface_eval(_s) ((_s)->eval) 
#define hefsi_surface_data(_s) ((_s)->data) 
#define hefsi_surface_derivative(_s) ((_s)->derivative) 
#define hefsi_surface_element_number(_s) ((_s)->ne)
#define hefsi_surface_element_number_max(_s) ((_s)->nemax)
#define hefsi_surface_node_number(_s) ((_s)->n->len)

#define hefsi_surface_node(_s,_i)		\
  (&(g_array_index((_s)->n,hefsi_node_t,(_i))))
#define hefsi_surface_element(_s,_i)		\
  (&(g_array_index((_s)->e,hefsi_element_t,(_i))))

#define hefsi_surface_u(_s,_i) (hefsi_node_u(hefsi_surface_node((_s),(_i))))
#define hefsi_surface_v(_s,_i) (hefsi_node_v(hefsi_surface_node((_s),(_i))))
#define hefsi_surface_x(_s,_i) (hefsi_node_x(hefsi_surface_node((_s),(_i))))
#define hefsi_surface_xu(_s,_i) (hefsi_node_xu(hefsi_surface_node((_s),(_i))))
#define hefsi_surface_xv(_s,_i) (hefsi_node_xv(hefsi_surface_node((_s),(_i))))
#define hefsi_surface_n(_s,_i)  (hefsi_node_n(hefsi_surface_node((_s),(_i))))

gint hefsi_surface_initialize(hefsi_surface_t *s, gint dmin, gint dmax) ;
gint hefsi_surface_evaluate(hefsi_surface_t *s, gdouble u, gdouble v,
			    gdouble *x, gdouble *xu, gdouble *xv,
			    gdouble *n) ;
gint hefsi_refine_point(hefsi_surface_t *s1, gdouble *uv1,
			hefsi_surface_t *s2, gdouble *uv2,
			double *x, gdouble tol, gint imax) ;
gint hefsi_set_bounding_boxes(hefsi_surface_t *s, gdouble scale) ;
gint hefsi_write_surface(FILE *f, hefsi_surface_t *s) ;
gint hefsi_write_bounding_boxes(FILE *f, hefsi_surface_t *s) ;

hefsi_surface_t *hefsi_surface_new(gpointer eval, gpointer data,
				   gboolean derivative,
				   gdouble umin, gdouble umax,
				   gdouble vmin, gdouble vmax) ;
hefsi_workspace_t *hefsi_workspace_new(void) ;
gint hefsi_surface_intersections(hefsi_surface_t *s1, hefsi_surface_t *s2, 
				 gdouble tol, hefsi_workspace_t *w) ;
gint hefsi_surface_intersection_check(hefsi_surface_t *s1, hefsi_surface_t *s2,
				      hefsi_workspace_t *w, gdouble tol,
				      gdouble *emax) ;

#endif /*HEFSI_H_INCLUDED*/
