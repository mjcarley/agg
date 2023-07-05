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

#ifndef __AGG_H_INCLUDED__
#define __AGG_H_INCLUDED__

/** @file */

/**
 * @{
 * @ingroup variables
 *
 */

/** @typedef agg_variable_t
 *
 * Data structure holding variable names, values, and symbolic definitions
 */

#ifdef DOXYGEN
typedef agg_variable_t ;
#else /*DOXYGEN*/
typedef struct _agg_variable_t agg_variable_t ;
struct _agg_variable_t {
  gchar *name, *def ; /*variable name and its symbolic definition*/
  gpointer eval ;     /*compiled form of definition*/
  gdouble val ;       /*numerical value*/
} ;
#endif /*DOXYGEN*/

#define agg_variable_name(_v)       ((_v)->name)
#define agg_variable_definition(_v) ((_v)->def)
#define agg_variable_evaluator(_v)  ((_v)->eval)
#define agg_variable_value(_v)      ((_v)->val)

/**
 *  @}
 */

/**
 * @{ 
 * @ingroup sections
 * 
 */

/** @typedef agg_section_type_t
 * 
 * Types of section
 */

typedef enum {
  AGG_SECTION_UNDEFINED = 0,	/**< undefined section type */
  AGG_SECTION_AEROFOIL  = 1,	/**< thickness distribution plus camber */
  AGG_SECTION_ELLIPSE   = 2	/**< generic closed curve */
} agg_section_type_t ;

/** @typedef agg_section_t
 * 
 * Data structure holding information for evaluation of a section
 */

#ifdef DOXYGEN
typedef agg_section_t ;
#else /*DOXYGEN*/
typedef struct _agg_section_t agg_section_t ;

struct _agg_section_t {
  agg_section_type_t type ; 
  gboolean close ;
  gdouble nl, nr, *cu, *cl, ytu, ytl ;
  gint ou, ol, oumax, olmax ;
} ;
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/**
 *  @brief A macro that returns the ::agg_section_type_t of an ::agg_section_t
 */
#define agg_section_type(s)                
/**
 * @brief TRUE if an ::agg_section_t needs to be closed
 */
#define agg_section_close(s)               
/**
 * @brief value of \f$\eta_{1}\f$ controlling curvature around \f$x=0\f$
 */
#define agg_section_eta_left(s)            
/**
 * @brief value of \f$\eta_{2}\f$ controlling curvature around \f$x=1\f$
 */
#define agg_section_eta_right(s)           
/**
 * @brief \f$i\f$th coefficient in Bernstein expansion of upper part of a
    section
 */
#define agg_section_coefficient_upper(s,i)
/**
 * @brief \f$i\f$th coefficient in Bernstein expansion of lower part of a
    section
 */
#define agg_section_coefficient_lower(s,i)
/**
 * @brief order of Bernstein expansion in upper part of section
 */
#define agg_section_order_upper(s)         
/**
 * @brief order of Bernstein expansion in lower part of section
 */
#define agg_section_order_lower(s)         
/**
 * @brief maximum order of Bernstein expansion in upper part of section
 */
#define agg_section_order_upper_max(s)     
/**
 * @brief maximum order of Bernstein expansion in lower part of section
 */
#define agg_section_order_lower_max(s)     
/**
 * @brief semi-thickness of upper part of trailing edge on aerofoil section
 */
#define agg_section_trailing_edge_upper(s) 
/**
 * @brief semi-thickness of lower part of trailing edge on aerofoil section
 */
#define agg_section_trailing_edge_lower(s)
#else /*DOXYGEN*/
#define agg_section_type(_s)                 ((_s)->type)
#define agg_section_close(_s)                ((_s)->close)
#define agg_section_eta_left(_s)             ((_s)->nl)
#define agg_section_eta_right(_s)            ((_s)->nr)
#define agg_section_coefficient_upper(_s,_i) ((_s)->cu[(_i)])
#define agg_section_coefficient_lower(_s,_i) ((_s)->cl[(_i)])
#define agg_section_order_upper(_s)          ((_s)->ou)
#define agg_section_order_lower(_s)          ((_s)->ol)
#define agg_section_order_upper_max(_s)      ((_s)->oumax)
#define agg_section_order_lower_max(_s)      ((_s)->olmax)
#define agg_section_trailing_edge_upper(_s)  ((_s)->ytu)
#define agg_section_trailing_edge_lower(_s)  ((_s)->ytl)
#endif /*DOXYGEN*/

/** @typedef agg_expression_data_t
 *
 * Data structure holding internal information for evaluation of expressions
 */

#ifdef DOXYGEN
typedef agg_expression_data_t ;
#else /*DOXYGEN*/
typedef struct _agg_expression_data_t agg_expression_data_t ;
struct _agg_expression_data_t {
  gint ne, nemax ;
  gpointer data, *expr ; /*wrappers for access to tinyexpr variable data type*/
  gchar **defs ;
} ;
#endif /*DOXYGEN*/

/**
 * @}
 */

#define AGG_OPERATOR_PARAMETER_SIZE 8

/**
 * @{
 * @ingroup transforms
 * 
 */

/** @typedef agg_operation_t
 * 
 * Basic transform operations
 */

typedef enum {
  AGG_TRANSFORM_UNDEFINED = 0,	/**< undefined transform */
  AGG_TRANSFORM_ROTATE    = 1,  /**< anti-clockwise rotation about a point */
  AGG_TRANSFORM_SHRINK    = 2,	/**< scaling by contracting towards a point */
  AGG_TRANSFORM_TRANSLATE = 3,	/**< translation in three dimensions */
  AGG_TRANSFORM_SCALE     = 4	/**< scaling by multiplying a constant factor */
} agg_operation_t ;

typedef gint (*agg_transform_operator_func_t)(agg_operation_t op,
					      agg_variable_t *p, gint np,
					      gdouble *xin, gdouble *xout) ;

/** @typedef agg_transform_operator_t
 *
 * Data structure for basic transform operations
 */

#ifdef DOXYGEN
typedef agg_transform_operator_t ;
#else /*DOXYGEN*/
typedef struct _agg_transform_operator_t agg_transform_operator_t ;
struct _agg_transform_operator_t {
  agg_operation_t op ;
  agg_variable_t p[AGG_OPERATOR_PARAMETER_SIZE] ;
  gint np ;
  agg_transform_operator_func_t func ;
} ;
#endif /*DOXYGEN*/

#define agg_transform_operator_operation(_op)        ((_op)->op)
#define agg_transform_operator_parameters(_op)       ((_op)->p)
#define agg_transform_operator_parameter(_op,_i)     (&((_op)->p[(_i)]))
#define agg_transform_operator_parameter_number(_op) ((_op)->np)
#define agg_transform_operator_func(_op)             ((_op)->func)

#define AGG_TRANSFORM_VARIABLE_NUMBER_MAX 128

/** @typedef agg_transform_t
 *
 * Data structure for transformations converting sections into surfaces
 */

#ifdef DOXYGEN
typedef agg_transform_t ;
#else /*DOXYGEN*/
typedef struct _agg_transform_t agg_transform_t ;
struct _agg_transform_t {
  agg_variable_t v[AGG_TRANSFORM_VARIABLE_NUMBER_MAX] ;
  agg_expression_data_t *e ;
  agg_transform_operator_t **op ;
  gint nop, nopmax, nv ;
} ;
#endif /*DOXYGEN*/

#define agg_transform_operator(_T,_i)           ((_T)->op[(_i)])
#define agg_transform_operator_number(_T)       ((_T)->nop)
#define agg_transform_operator_number_max(_T)   ((_T)->nopmax)
#define agg_transform_variable(_T,_i)           (&((_T)->v[(_i)]))
#define agg_transform_variable_number(_T)       ((_T)->nv)

/**
 *  @}
 */

/**
 * @{
 * @ingroup surfaces
 * 
 */

typedef enum {
  AGG_AXES_PX_PY_PZ = 0,  /**< standard, untransformed orientation */
  AGG_AXES_PY_PZ_PX = 1,  /**< (x,y,z) -> (y, z, x)*/
  AGG_AXES_PZ_PX_PY = 2,  /**< (x,y,z) -> (z, x, y)*/
  AGG_AXES_PZ_PY_PX = 3,  /**< (x,y,z) -> (z, y, x)*/
} agg_axes_t ;

/** @typedef agg_surface_t
 *
 * Data structure to hold surfaces formed by transformation of sections
 */

#ifdef DOXYGEN
typedef agg_surface_t ;
#else /*DOXYGEN*/
typedef struct _agg_surface_t agg_surface_t ;
struct _agg_surface_t {
  agg_section_t *s ;
  agg_axes_t axes ;
  gint ns, nsmax ;
  gdouble umin, umax, *u, *w ;
  agg_transform_t *T ;
} ;
#endif /*DOXYGEN*/

#define agg_surface_section(_s,_i)         (&((_s)->s[(_i)]))
#define agg_surface_section_number(_s)     ((_s)->ns)
#define agg_surface_section_number_max(_s) ((_s)->nsmax)
#define agg_surface_u(_s,_i)               ((_s)->u[(_i)])
#define agg_surface_weight(_s,_i)          ((_s)->w[(_i)])
#define agg_surface_umin(_s)               ((_s)->umin)
#define agg_surface_umax(_s)               ((_s)->umax)
#define agg_surface_transform(_s)          ((_s)->T)
#define agg_surface_axes(_s)               ((_s)->axes)

/**
 *  @}
 */

/**
 * @{
 * @ingroup patches
 * 
 */

#define AGG_PATCH_POINT_SIZE 3

/** @typedef agg_patch_mapping_t
 * 
 * Mappings from \f$(s,t)\f$ on a patch to \f$(u,v)\f$ on a surface
 */

typedef enum {
  AGG_PATCH_BILINEAR  = 0, /**< bilinear mapping between corners */
  AGG_PATCH_SPHERICAL = 1, /**< mapping to surface of sphere */
  AGG_PATCH_TUBULAR   = 2  /**< mapping to surface of open cylinder */
} agg_patch_mapping_t ;

#ifdef DOXYGEN
typedef agg_patch_t ;
#else /*DOXYGEN*/
typedef struct _agg_patch_t agg_patch_t ;
struct _agg_patch_t {
  agg_patch_mapping_t map ;
  gdouble *st, uvc[16] ; /*boundary points*/
  gint nc, ci[8],  /*number of corners, corner indices*/
    nst, nstmax ;   /*number of boundary points, maximum number*/
  gboolean swrap, twrap ; /*dealing with s and t out of range*/
} ;

#define agg_patch_point_number(_p)      ((_p)->nst)
#define agg_patch_point_number_max(_p)  ((_p)->nstmax)
#define agg_patch_point_s(_p,_i)        ((_p)->st[AGG_PATCH_POINT_SIZE*(_i)+0])
#define agg_patch_point_t(_p,_i)        ((_p)->st[AGG_PATCH_POINT_SIZE*(_i)+1])
#define agg_patch_point_len(_p,_i)      ((_p)->st[AGG_PATCH_POINT_SIZE*(_i)+2])
#define agg_patch_corner_u(_p,_i)       ((_p)->uvc[2*(_i)+0])
#define agg_patch_corner_v(_p,_i)       ((_p)->uvc[2*(_i)+1])
#define agg_patch_corner_index(_p,_i)   ((_p)->ci[(_i)])
#define agg_patch_corner_s(_p,_i)		\
  ((_p)->st[AGG_PATCH_POINT_SIZE*(agg_patch_corner_index((_p),(_i)))+0])
#define agg_patch_corner_t(_p,_i)		\
  ((_p)->st[AGG_PATCH_POINT_SIZE*(agg_patch_corner_index((_p),(_i)))+1])
#define agg_patch_corner_number(_p)     ((_p)->nc)
#define agg_patch_mapping(_p)           ((_p)->map)
#define agg_patch_wrap_s(_p)            ((_p)->swrap)
#define agg_patch_wrap_t(_p)            ((_p)->twrap)
#define agg_patch_edge_is_cut(_p,_i)		\
  (((_p)->ci[(_i)+1] - (_p)->ci[(_i)]) > 1)

#endif /*DOXYGEN*/

/**
 * @}
 */


/**
 * @{
 * 
 * @ingroup intersections
 */

#define AGG_INTERSECTION_DATA_SIZE 8

#ifdef DOXYGEN
/** @typedef agg_intersection_t
 *
 * Data structure to hold intersections of surfaces
 */
typedef agg_intersection_t ;
#else  /*DOXYGEN*/
typedef struct _agg_intersection_t agg_intersection_t ;
struct _agg_intersection_t {
  agg_surface_t *S1, *S2 ;
  agg_patch_t *P1, *P2 ;
  gdouble *st, bbox[8] ;
  gint nst, nstmax, ibox[8] ;
} ;
#endif /*DOXYGEN*/

#define agg_intersection_surface1(_i) ((_i)->S1)
#define agg_intersection_surface2(_i) ((_i)->S2)
#define agg_intersection_patch1(_i) ((_i)->P1)
#define agg_intersection_patch2(_i) ((_i)->P2)
#define agg_intersection_point_s1(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+0])
#define agg_intersection_point_t1(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+1])
#define agg_intersection_point_s2(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+2])
#define agg_intersection_point_t2(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+3])
#define agg_intersection_parameter(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+4])
#define agg_intersection_point(_i,_j)			\
  (&((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+5]))
#define agg_intersection_point_number(_i) ((_i)->nst)
#define agg_intersection_point_number_max(_i) ((_i)->nstmax)
#define agg_intersection_bbox_s1_min(_i) ((_i)->bbox[0])
#define agg_intersection_bbox_t1_min(_i) ((_i)->bbox[1])
#define agg_intersection_bbox_s2_min(_i) ((_i)->bbox[2])
#define agg_intersection_bbox_t2_min(_i) ((_i)->bbox[3])
#define agg_intersection_bbox_s1_max(_i) ((_i)->bbox[4])
#define agg_intersection_bbox_t1_max(_i) ((_i)->bbox[5])
#define agg_intersection_bbox_s2_max(_i) ((_i)->bbox[6])
#define agg_intersection_bbox_t2_max(_i) ((_i)->bbox[7])
#define agg_intersection_bbox_s1_min_index(_i) ((_i)->ibox[0])
#define agg_intersection_bbox_t1_min_index(_i) ((_i)->ibox[1])
#define agg_intersection_bbox_s2_min_index(_i) ((_i)->ibox[2])
#define agg_intersection_bbox_t2_min_index(_i) ((_i)->ibox[3])
#define agg_intersection_bbox_s1_max_index(_i) ((_i)->ibox[4])
#define agg_intersection_bbox_t1_max_index(_i) ((_i)->ibox[5])
#define agg_intersection_bbox_s2_max_index(_i) ((_i)->ibox[6])
#define agg_intersection_bbox_t2_max_index(_i) ((_i)->ibox[7])

/**
 * @}
 */


/**
 * @{
 * 
 * @ingroup wireframe
 */

#define AGG_WIREFRAME_POINT_SIZE   5
#define AGG_WIREFRAME_ELEMENT_SIZE 5

#ifdef DOXYGEN
/** @typedef agg_wireframe_t
 *
 * Data structure to hold wireframe (patch) representation of surfaces
 */
typedef agg_wireframe_t ;
#else /*DOXYGEN*/
typedef struct _agg_wireframe_t agg_wireframe_t ;

struct _agg_wireframe_t {
  agg_surface_t *S[32] ;
  agg_patch_t *P[32] ;
  gdouble *p ;
  gint np, npmax, *sp, nsp, nspmax, *isp, *e, ne, nemax, *ptags,
    isect[32], ipps[32], insp[32], ni, ns ;
  agg_intersection_t *inter[32] ; 
} ;
#endif /*DOXYGEN*/

#define agg_wireframe_point(_w,_i)  (&((_w)->p[AGG_WIREFRAME_POINT_SIZE*(_i)]))
#define agg_wireframe_point_s(_w,_i)  ((_w)->p[AGG_WIREFRAME_POINT_SIZE*(_i)+3])
#define agg_wireframe_point_t(_w,_i)  ((_w)->p[AGG_WIREFRAME_POINT_SIZE*(_i)+4])
#define agg_wireframe_point_tag(_w,_i) ((_w)->ptags[(_i)])
#define agg_wireframe_point_number(_w)  ((_w)->np)
#define agg_wireframe_point_number_max(_w)  ((_w)->npmax)
#define agg_wireframe_element(_w,_i)		\
  &((_w)->e[AGG_WIREFRAME_ELEMENT_SIZE*(_i)])
#define agg_wireframe_element_number(_w)  ((_w)->ne)
#define agg_wireframe_element_number_max(_w)  ((_w)->nemax)
#define agg_wireframe_spline_number(_w)  ((_w)->nsp)
#define agg_wireframe_spline_number_max(_w)  ((_w)->nspmax)
#define agg_wireframe_spline(_w,_i)  &((_w)->sp[(_w)->isp[(_i)]])
#define agg_wireframe_spline_length(_w,_i)	\
  ((_w)->isp[(_i)+1]-(_w)->isp[(_i)])
#define agg_wireframe_patch(_w,_i)       ((_w)->P[(_i)])
#define agg_wireframe_surface(_w,_i)     ((_w)->S[(_i)])
#define agg_wireframe_surface_number(_w) ((_w)->ns)
#define agg_wireframe_element_size(_e)			\
  ( ((_e)[4] != 0) ? 5 :				\
    (((_e)[3] != 0) ? 4 : 3) )
#define agg_wireframe_intersection(_w,_i)     ((_w)->inter[(_i)])
#define agg_wireframe_intersection_number(_w) ((_w)->ni)
#define agg_wireframe_intersection_points_start(_w,_i) \
  ((_w)->isect[2*(_i)+0])
#define agg_wireframe_intersection_points_end(_w,_i)	\
  ((_w)->isect[2*(_i)+1])
#define agg_wireframe_intersection_points_per_spline(_w,_i)	\
  ((_w)->ipps[(_i)])
#define agg_wireframe_intersection_spline_number(_w,_i)	\
  ((_w)->insp[(_i)])

/**
 * @}
 */

typedef struct _agg_surface_workspace_t agg_surface_workspace_t ;
struct _agg_surface_workspace_t {
  agg_section_t *s ;
} ;

gint agg_bernstein_basis(gint n, gdouble x, gdouble *S, gdouble *dS) ;
gdouble agg_bernstein_basis_eval(gint n, gint r, gdouble x) ;

agg_section_t *agg_section_new(gint oumax, gint olmax) ;
gdouble agg_section_eval(agg_section_t *s, gdouble x) ;
gint agg_section_copy(agg_section_t *dest, agg_section_t *src) ;
gint agg_section_set_circle(agg_section_t *s) ;
gint agg_section_set_aerofoil(agg_section_t *s, gdouble eta,
			      gdouble th, gdouble yte) ;

agg_transform_operator_t *agg_transform_operator_new(void) ;
agg_transform_t *agg_transform_new(gint nopmax) ;
gint agg_transform_variable_add(agg_transform_t *T,
				gchar *var, gchar *def, gdouble val) ;
gint agg_transform_add_global_variables(agg_transform_t *T,
					agg_variable_t *global, gint nglobal) ;
gint agg_transform_expressions_compile(agg_transform_t *T) ;
gint agg_transform_variables_eval(agg_transform_t *T) ;
gint agg_transform_variables_write(FILE *f, agg_transform_t *T,
				   gboolean write_defs) ;
gint agg_transform_operator_add(agg_transform_t *T, agg_operation_t op,
				gdouble *p, gchar **expr, gint np) ;
gint agg_transform_operators_write(FILE *f, agg_transform_t *T) ;
gint agg_transform_apply(agg_transform_t *T, gdouble *xin, gdouble *xout) ;

gint agg_transform_operator_apply(agg_transform_operator_t *op,
				  gdouble *xin, gdouble *xout) ;
gint agg_transform_operator_rotate(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout) ;
gint agg_transform_operator_translate(agg_operation_t op,
				      agg_variable_t *p, gint np,
				      gdouble *xin, gdouble *xout) ;
gint agg_transform_operator_shrink(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout) ;
gint agg_transform_operator_scale(agg_operation_t op,
				  agg_variable_t *p, gint np,
				  gdouble *xin, gdouble *xout) ;
gint agg_transform_axes(agg_axes_t axes, gdouble *xin, gdouble *xout) ;

agg_variable_t *agg_variable_find(agg_variable_t **vars, gint nv,
				  gchar *var) ;
agg_expression_data_t *agg_expression_data_new(gint nemax) ;
gint agg_expression_data_variable_add(agg_expression_data_t *d,
				      agg_variable_t *v) ;
gpointer agg_expression_compile(gchar *e, agg_expression_data_t *d) ;
gdouble agg_expression_eval(gpointer e) ;
gint agg_expression_data_compile(agg_expression_data_t *d) ;
gint agg_expression_data_eval(agg_expression_data_t *d) ;

agg_surface_t *agg_surface_new(gint nsmax) ;
gint agg_surface_section_add(agg_surface_t *S, agg_section_t *s,
			     gdouble u) ;
gint agg_surface_weights_make(agg_surface_t *S) ;
gint agg_surface_section_interp(agg_surface_t *S, gdouble u, agg_section_t *s) ;
gint agg_surface_point_eval(agg_surface_t *S, gdouble u, gdouble v,
			    gdouble *x, agg_surface_workspace_t *w) ;
agg_surface_workspace_t *agg_surface_workspace_new(void) ;

gdouble agg_naca_four(gdouble t, gdouble p, gdouble m, gdouble x) ;

agg_patch_t *agg_patch_new(gint nstmax) ;
gint agg_patch_map(agg_patch_t *p, gdouble s, gdouble t,
		   gdouble *u, gdouble *v) ;
gint agg_patch_boundary_write(FILE *f, agg_patch_t *P, gdouble ds,
			      agg_surface_workspace_t *w) ;
gint agg_patch_edges_parameterize(agg_patch_t *P) ;
gint agg_patch_sample(agg_patch_t *P, gdouble ss, gdouble ts,
		      gdouble *s, gdouble *t) ;
gint agg_patch_edge_interp(agg_patch_t *P, gint e, gdouble p,
			   gdouble *s, gdouble *t) ;
gint agg_patch_st_correct(agg_patch_t *P, gdouble *st) ;

agg_intersection_t *agg_intersection_new(gint nstmax) ;
gint agg_surface_patch_intersection(agg_intersection_t *inter,
				    agg_surface_t *S1, agg_patch_t *P1,
				    agg_surface_t *S2, agg_patch_t *P2,
				    agg_surface_workspace_t *w) ;
gint agg_intersection_curve_write(FILE *f, agg_intersection_t *inter,
				  agg_surface_workspace_t *w) ;
gint agg_intersection_bbox_set(agg_intersection_t *inter) ;
gint agg_intersection_resample(agg_intersection_t *inter,
			       gint nsp, gint pps,
			       agg_intersection_t *resample,
			       agg_surface_workspace_t *w) ;

agg_wireframe_t *agg_wireframe_new(gint npmax, gint nspmax, gint nemax) ;
gint agg_wireframe_surface_make(agg_wireframe_t *w,
				agg_surface_t *S, agg_patch_t *P,
				agg_intersection_t *inter,
				gint nsec, gint nseg, gint pps,
				agg_surface_workspace_t *ws) ;
gint agg_wireframe_write_gmsh(FILE *f, agg_wireframe_t *w,
			      gchar *len, gint offp, gint offsp, gint offs,
			      gboolean opencascade) ;
gint agg_wireframe_spline_ends(agg_wireframe_t *w, gint s, gint *p0, gint *p1) ;
gint agg_wireframe_spline_interp_points(agg_wireframe_t *w, gint s,
					gint p0, gint p1, gint pps,
					agg_surface_workspace_t *ws) ;
gboolean agg_wireframe_splines_connected(agg_wireframe_t *w,
					 gint i0, gint i1, gint *p) ;
gint agg_wireframe_element_add(agg_wireframe_t *w,
			       gint s0, gint s1, gint s2, gint s3) ;
gboolean agg_wireframe_spline_degenerate(agg_wireframe_t *w, gint s) ;
gint agg_wireframe_intersection_add(agg_wireframe_t *w,
				    agg_intersection_t *inter,
				    gint nsp, gint pps,
				    agg_surface_workspace_t *ws) ;
gint agg_wireframe_surface_add(agg_wireframe_t *w,
			       agg_surface_t *S, agg_patch_t *P,
			       gint nsec, gint nseg, gint pps,
			       agg_surface_workspace_t *ws) ;
gint agg_wireframe_spline_from_endpoints(agg_wireframe_t *w,
					 gint p0, gint p1) ;
gint agg_wireframe_surface_point_add(agg_wireframe_t *w, gint surf,
				     gdouble s, gdouble t,
				     agg_surface_workspace_t *ws) ;

#endif /*__AGG_H_INCLUDED__*/

