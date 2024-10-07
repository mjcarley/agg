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

#ifdef DOXYGEN
/**
 *  @brief
 */
#else /*DOXYGEN*/
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/** @typedef agg_curve_type_t
 *
 * Enumerated data type for various curve definitions
 */
#endif /*DOXYGEN*/

typedef enum {
  AGG_CURVE_POLYNOMIAL = 0, /**< polynomial */
  AGG_CURVE_ELLIPSE    = 1, /**< axis-aligned ellipse with arbitrary centre */
  AGG_CURVE_FOURIER    = 2  /**< Fourier series */
} agg_curve_type_t ;

#define AGG_CURVE_DATA_SIZE 35
  
#ifdef DOXYGEN
/** @typedef agg_curve_t
 *
 * Data structure for evaluation of one-dimensional parametric curves
 */

typedef agg_curve_t agg_curve_t ;

/**
 *  @brief ::agg_curve_type_t of an ::agg_curve_t
 */
#define agg_curve_order(c)
/**
 *  @brief order (number of parameters) of an ::agg_curve_t
 */
#define agg_curve_data(c)
/**
 *  @brief data (coefficients, etc) of an ::agg_curve_t
 */
#define agg_curve_data(c)

#else  /*DOXYGEN*/
typedef struct _agg_curve_t agg_curve_t ;
struct _agg_curve_t {
  agg_curve_type_t type ;
  gint n ;
  gdouble data[AGG_CURVE_DATA_SIZE] ;
} ;

#define agg_curve_type(_c)        ((_c)->type)
#define agg_curve_order(_c)       ((_c)->n)
#define agg_curve_data(_c)        ((_c)->data)
#endif  /*DOXYGEN*/

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
  char *name, *def ; /*variable name and its symbolic definition*/
  gpointer eval ;     /*compiled form of definition*/
  gdouble val ;       /*numerical value*/
} ;
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/**
 *  @brief name of an ::agg_variable_t
 */
#define agg_variable_name(v)       
/**
 *  @brief definition of an ::agg_variable_t, as as string, or NULL
 *  for numerical constants
 */
#define agg_variable_definition(v) 
/**
 *  @brief evaluator for the compiled expression for the variable,
 *  or NULL for numerical constants
 */
#define agg_variable_evaluator(v)  
/**
 *  @brief numerical value of the variable (which may change on
 *  upon reevaluation of the expression)
 */
#define agg_variable_value(v)      
#else  /*DOXYGEN*/
#define agg_variable_name(_v)       ((_v)->name)
#define agg_variable_definition(_v) ((_v)->def)
#define agg_variable_evaluator(_v)  ((_v)->eval)
#define agg_variable_value(_v)      ((_v)->val)
#endif /*DOXYGEN*/

/**
 *  @}
 */

/**
 * @{ 
 * @ingroup sections
 * 
 */

/** @typedef agg_section_t
 * 
 * Data structure holding information for evaluation of a section
 */

#ifdef DOXYGEN
typedef agg_section_t ;
#else /*DOXYGEN*/
typedef struct _agg_section_t agg_section_t ;

struct _agg_section_t {
  gboolean close ;
  gdouble nl, nr, *cu, *cl, ytu, ytl ;
  gint ou, ol, oumax, olmax ;
} ;
#endif /*DOXYGEN*/

#ifdef DOXYGEN
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

#define AGG_CHEBYSHEV_INTERVAL_NUMBER_MAX 512

#ifdef DOXYGEN

typedef agg_chebyshev_t ;

#else /*DOXYGEN*/

typedef struct _agg_chebyshev_t agg_chebyshev_t ;
struct _agg_chebyshev_t {
  gint Nmax, nc, ninter, inter[AGG_CHEBYSHEV_INTERVAL_NUMBER_MAX+1] ;
  gdouble xinter[AGG_CHEBYSHEV_INTERVAL_NUMBER_MAX+1], *f ;
  gdouble x0[3] ;
} ;

#define agg_chebyshev_order(_C)              ((_C)->N)
#define agg_chebyshev_component_number(_C)   ((_C)->nc)
/* #define agg_chebyshev_data(_C,_i,_j)				\ */
/*   ((_C)->f[(_i)*agg_chebyshev_component_number((_C))+(_j)]) */
#define agg_chebyshev_interval_number(_C)    ((_C)->ninter)
#define agg_chebyshev_interval_start(_C,_i)  ((_C)->inter[(_i)])
#define agg_chebyshev_interval_xmin(_C,_i)   ((_C)->xinter[(_i)])
#define agg_chebyshev_interval_order(_C,_i)		\
  (((_C)->inter[(_i)+1]) - ((_C)->inter[(_i)])-1)
#define agg_chebyshev_interval_coordinate(_C,_i,_x)			\
  (2.0*((_x)-((_C)->xinter[(_i)]))/((_C)->xinter[(_i)+1]-((_C)->xmin[(_i)]))\
   - 1.0)
  
#define agg_chebyshev_

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
  char **defs ;
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
  AGG_TRANSFORM_UNDEFINED = 0, /**< undefined transform */
  AGG_TRANSFORM_ROTATE    = 1, /**< anti-clockwise rotation about a point */
  AGG_TRANSFORM_SHRINK    = 2, /**< scaling by contracting towards a point */
  AGG_TRANSFORM_TRANSLATE = 3, /**< translation in three dimensions */
  AGG_TRANSFORM_SCALE     = 4, /**< scaling by multiplying a constant factor */
  AGG_TRANSFORM_SCALE_X   = 5, /**< scaling in x multiplying a constant factor */
  AGG_TRANSFORM_SCALE_Y   = 6  /**< scaling in y multiplying a constant factor */
} agg_operation_t ;

typedef gint (*agg_transform_operator_func_t)(agg_operation_t op,
					      agg_variable_t *p, gint np,
					      gdouble *xin, gdouble *xout,
					      gdouble *dxdu, gdouble *dxdv) ;

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
  agg_variable_t
  p[AGG_OPERATOR_PARAMETER_SIZE],
    pu[AGG_OPERATOR_PARAMETER_SIZE],
    pv[AGG_OPERATOR_PARAMETER_SIZE] ;
  gint np ;
  agg_transform_operator_func_t func ;
  gdouble umin, umax ;
} ;
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/**
 *  @brief operation performed by an ::agg_transform_operator_t
 */
#define agg_transform_operator_operation(op)        
/**
 *  @brief parameters of ::agg_variable_t for an ::agg_transform_operator_t
 */
#define agg_transform_operator_parameters(op)       
/**
 *  @brief \f$i\f$th parameter for ::agg_transform_operator_t
 */
#define agg_transform_operator_parameter(op,i)     
/**
 *  @brief number of parameters for the ::agg_transform_operator_t
 */
#define agg_transform_operator_parameter_number(op) 
/**
 *  @brief function called to execute the transform operation
 */
#define agg_transform_operator_func(op)             
/**
 *  @brief lower limit of parameter range for operation
 */
#define agg_transform_operator_umin(op)
/**
 *  @brief upper limit of parameter range for operation
 */
#define agg_transform_operator_umax(op)
#else /*DOXYGEN*/
#define agg_transform_operator_operation(_op)        ((_op)->op)
#define agg_transform_operator_parameters(_op)       ((_op)->p)
#define agg_transform_operator_parameter(_op,_i)     (&((_op)->p[(_i)]))
#define agg_transform_operator_parameter_u(_op,_i)   (&((_op)->pu[(_i)]))
#define agg_transform_operator_parameter_v(_op,_i)   (&((_op)->pv[(_i)]))
#define agg_transform_operator_parameter_number(_op) ((_op)->np)
#define agg_transform_operator_func(_op)             ((_op)->func)
#define agg_transform_operator_umin(_op)             ((_op)->umin)
#define agg_transform_operator_umax(_op)             ((_op)->umax)
#endif /*DOXYGEN*/

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

#ifdef DOXYGEN
/**
 *  @brief \f$i\f$th operation of an ::agg_transform_t
 */
#define agg_transform_operator(T,i)        
/**
 *  @brief number of operations in an ::agg_transform_t
 */
#define agg_transform_operator_number(T)    
/**
 *  @brief maximum number of operations in an  ::agg_transform_t
 */
#define agg_transform_operator_number_max(T)
/**
 *  @brief \f$i\f$th variable in an  ::agg_transform_t
 */
#define agg_transform_variable(T,i)        
/**
 *  @brief number of variables in an  ::agg_transform_t
 */
#define agg_transform_variable_number(T)    
#else /*DOXYGEN*/
#define agg_transform_operator(_T,_i)           ((_T)->op[(_i)])
#define agg_transform_operator_number(_T)       ((_T)->nop)
#define agg_transform_operator_number_max(_T)   ((_T)->nopmax)
#define agg_transform_variable(_T,_i)           (&((_T)->v[(_i)]))
#define agg_transform_variable_number(_T)       ((_T)->nv)
#endif /*DOXYGEN*/

/**
 *  @}
 */

/**
 * @{
 * @ingroup surfaces
 * 
 */

/** @typedef agg_axes_t
 *
 * Transformation between different coordinate systems by reassignment
 * of axes
 */

typedef enum {
  AGG_AXES_UNDEFINED = 0,
  AGG_AXES_PX_PY_PZ = 1,  /**< standard, untransformed orientation */
  AGG_AXES_PY_PZ_PX = 2,  /**< (x,y,z) -> (y, z, x)*/
  AGG_AXES_PZ_PX_PY = 3,  /**< (x,y,z) -> (z, x, y)*/
  AGG_AXES_PZ_PY_PX = 4,  /**< (x,y,z) -> (z, y, x)*/
  AGG_AXES_PX_PY_MZ = 5,  /**< (x,y,z) -> (x, y,-z) (reflection in x-y)*/
} agg_axes_t ;

/** @typedef agg_grid_t
 *
 * Specification of surface grid generation algorithms
 */

typedef enum {
  AGG_GRID_UNDEFINED  = 0,    /**< undefined (generates an error) */
  AGG_GRID_REGULAR    = 1,    /**< regular quadrilaterals in parametric space */
  AGG_GRID_TRIANGLE   = 2,    /**< mesh generated using Triangle code */
  AGG_GRID_SPHERE_ICO = 3,    /**< spherical mesh generated using icosahedron 
				 subdivision */
  AGG_GRID_SPHERE_UV =  4,    /**< spherical mesh generated using regular grid*/
  AGG_GRID_HEMISPHERE_ICO = 5,  /**< hemispherical mesh generated using 
				  icosahedron subdivision */
  AGG_GRID_HEMISPHERE_UV = 6  /**< hemispherical mesh generated using 
				 regular grid */
} agg_grid_t ;

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
  agg_grid_t grid ;
  gint ns, nsmax, nsec, nsp, sub ;
  gdouble umin, umax, *u, *w, area ;
  agg_transform_t *T ;
} ;
#endif /*DOXYGEN*/

#define agg_surface_section(_s,_i)          (&((_s)->s[(_i)]))
#define agg_surface_section_number(_s)      ((_s)->ns)
#define agg_surface_section_number_max(_s)  ((_s)->nsmax)
#define agg_surface_u(_s,_i)                ((_s)->u[(_i)])
#define agg_surface_weight(_s,_i)           ((_s)->w[(_i)])
#define agg_surface_umin(_s)                ((_s)->umin)
#define agg_surface_umax(_s)                ((_s)->umax)
#define agg_surface_transform(_s)           ((_s)->T)
#define agg_surface_axes(_s)                ((_s)->axes)
#define agg_surface_grid(_s)                ((_s)->grid)
#define agg_surface_grid_section_number(_s) ((_s)->nsec)
#define agg_surface_grid_spline_number(_s)  ((_s)->nsp)
#define agg_surface_grid_element_area(_s)   ((_s)->area)
#define agg_surface_grid_subdivision(_s)    ((_s)->sub)

/**
 *  @}
 */

/**
 * @{
 * @ingroup patches
 * 
 */

/* typedef enum { */
/*   AGG_CLIP_UNDEFINED  = 0, */
/*   AGG_CLIP_CONSTANT_S = 1, */
/*   AGG_CLIP_CONSTANT_T = 2, */
/*   AGG_CLIP_ELLIPSE    = 3 */
/* } agg_patch_clip_t ; */

/* typedef struct _agg_patch_clipping_t agg_patch_clipping_t ; */

/* struct _agg_patch_clipping_t { */
/*   agg_patch_clip_t c ; */
/*   gdouble data[4], ornt ; */
/* } ; */

/* #define agg_patch_clipping_type(_c)        ((_c)->c) */
/* #define agg_patch_clipping_data(_c,_i)     ((_c)->data[(_i)]) */
/* #define agg_patch_clipping_orientation(_c) ((_c)->ornt) */

/* #define AGG_PATCH_POINT_SIZE 3 */

/** @typedef agg_patch_mapping_t
 * 
 * Mappings from \f$(s,t)\f$ on a patch to \f$(u,v)\f$ on a surface
 */

typedef enum {
  AGG_PATCH_BILINEAR      = 0, /**< bilinear mapping between corners */
  AGG_PATCH_SPHERICAL     = 1, /**< mapping to surface of sphere */
  AGG_PATCH_HEMISPHERICAL = 2, /**< mapping to surface of hemisphere */
  AGG_PATCH_TUBULAR       = 3  /**< mapping to surface of open cylinder */
} agg_patch_mapping_t ;

#define AGG_PATCH_HOLE_NUMBER_MAX 8

/** @typedef agg_patch_t
 *
 * Data structure for mapping of surfaces
 */
#ifdef DOXYGEN
typedef agg_patch_t ;
#else /*DOXYGEN*/
typedef struct _agg_patch_t agg_patch_t ;
struct _agg_patch_t {
  agg_patch_mapping_t map ;
  agg_curve_t curves[AGG_PATCH_HOLE_NUMBER_MAX+2] ;
  gpointer B[AGG_PATCH_HOLE_NUMBER_MAX+2] ;
  gint nholes ;
  gboolean swrap, twrap,  /*dealing with s and t out of range*/
    invert ; /*invert triangles to maintain correct normal*/
} ;

#define agg_patch_point_number(_P)      ((_P)->nst)
#define agg_patch_point_number_max(_P)  ((_P)->nstmax)
#define agg_patch_mapping(_P)           ((_P)->map)
#define agg_patch_wrap_s(_P)            ((_P)->swrap)
#define agg_patch_wrap_t(_P)            ((_P)->twrap)
#define agg_patch_invert(_P)            ((_P)->invert)
#define agg_patch_hole_number(_P)       ((_P)->nholes)
#define agg_patch_curve_smin(_P)        (&((_P)->curves[0]))
#define agg_patch_curve_smax(_P)        (&((_P)->curves[1]))
#define agg_patch_hole(_P,_i)           (&((_P)->curves[(_i)]))
#define agg_patch_blend(_P,_i)          ((_P)->B[(_i)])

#endif /*DOXYGEN*/

#ifdef DOXYGEN
/**
 * @typedef agg_surface_blend_t
 * 
 * Data type for evaluation of blended surfaces for bridging surface
 * intersections smoothly
 * 
 */

typedef agg_surface_blend_t ;

/**
 * @brief \f$i\f$th surface (i=0,1) joined by surface blend
 */
#define agg_surface_blend_surface(B,i)       
/**
 * @brief patch used to map \f$i\f$th (i=0,1) surface joined by blend
 */
#define agg_surface_blend_patch(B,i)         
/**
 * @brief if TRUE, invert blend to correct orientation of surface elements
 */
#define agg_surface_blend_invert(B)           
/**
 * @brief rail curve for blend on \f$i\f$th surface (i=0,1) joined by blend
 */
#define agg_surface_blend_hole(B,i)          
/**
 * @brief if TRUE, traverse \f$i\f$th (i=0,1) rail curve in opposite direction
 */
#define agg_surface_blend_curve_reverse(B,i) 

#else /*DOXYGEN*/
typedef struct _agg_surface_blend_t agg_surface_blend_t ;

struct _agg_surface_blend_t {
  agg_surface_t *S[2] ;
  agg_patch_t *P[2] ;
  gboolean invert, reverse[2] ;
  gint hole[2] ; /*index of hole or cut on surfaces to be blended*/
} ;

#define agg_surface_blend_surface(_B,_i)       ((_B)->S[(_i)])
#define agg_surface_blend_patch(_B,_i)         ((_B)->P[(_i)])
#define agg_surface_blend_invert(_B)           ((_B)->invert)
#define agg_surface_blend_hole(_B,_i)          ((_B)->hole[(_i)])
#define agg_surface_blend_curve_reverse(_B,_i) ((_B)->reverse[(_i)])
#endif /*DOXYGEN*/

/**
 * @}
 */


/**
 * @{
 * 
 * @ingroup intersections
 */

#define AGG_INTERSECTION_DATA_SIZE 7

#ifdef DOXYGEN
/** @typedef agg_intersection_t
 *
 * Data structure to hold intersections of surfaces
 */
typedef agg_intersection_t ;
#else  /*DOXYGEN*/
typedef struct _agg_intersection_t agg_intersection_t ;
struct _agg_intersection_t {
  agg_surface_t *S[2] ;
  agg_patch_t *P[2] ;
  gdouble *st, bbox[8] ;
  gint nst, nstmax, ibox[8] ;
} ;
#endif /*DOXYGEN*/

#define agg_intersection_surface1(_i) ((_i)->S[0])
#define agg_intersection_surface2(_i) ((_i)->S[1])
#define agg_intersection_patch1(_i) ((_i)->P[0])
#define agg_intersection_patch2(_i) ((_i)->P[1])
#define agg_intersection_point_s(_i,_j,_c)		\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+2*(_c)+0])
#define agg_intersection_point_t(_i,_j,_c)		\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+2*(_c)+1])
#define agg_intersection_point_s1(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+0])
#define agg_intersection_point_t1(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+1])
#define agg_intersection_point_s2(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+2])
#define agg_intersection_point_t2(_i,_j)	\
  ((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+3])
#define agg_intersection_point(_i,_j)			\
  (&((_i)->st[AGG_INTERSECTION_DATA_SIZE*(_j)+4]))
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
 * @ingroup mesh
 */

#define AGG_MESH_POINT_SIZE   5
#define AGG_MESH_ELEMENT_SIZE 5

#ifdef DOXYGEN
/** @typedef agg_mesh_t
 *
 * Data structure for surface meshes
 */
typedef agg_mesh_t ;
#else /*DOXYGEN*/
typedef struct _agg_mesh_t agg_mesh_t ;

struct _agg_mesh_t {
  agg_surface_t *S[32] ;
  agg_patch_t *P[32] ;
  gdouble *p ;
  gint np, npmax, *sp, nsp, nspmax, *isp, *e, ne, nemax, *ptags,
    isect[32], ipps[32], insp[32], ni, ns, nb ;
  agg_intersection_t *inter[32] ;
  agg_surface_blend_t B[32] ;
} ;
#endif /*DOXYGEN*/

#define agg_mesh_point(_w,_i)  (&((_w)->p[AGG_MESH_POINT_SIZE*(_i)]))
#define agg_mesh_point_s(_w,_i)  ((_w)->p[AGG_MESH_POINT_SIZE*(_i)+3])
#define agg_mesh_point_t(_w,_i)  ((_w)->p[AGG_MESH_POINT_SIZE*(_i)+4])
#define agg_mesh_point_tag(_w,_i) ((_w)->ptags[(_i)])
#define agg_mesh_point_number(_w)  ((_w)->np)
#define agg_mesh_point_number_max(_w)  ((_w)->npmax)
#define agg_mesh_element(_w,_i) &((_w)->e[AGG_MESH_ELEMENT_SIZE*(_i)])
#define agg_mesh_element_number(_w)  ((_w)->ne)
#define agg_mesh_element_number_max(_w)  ((_w)->nemax)
#define agg_mesh_spline_number(_w)  ((_w)->nsp)
#define agg_mesh_spline_number_max(_w)  ((_w)->nspmax)
#define agg_mesh_spline(_w,_i)  &((_w)->sp[(_w)->isp[(_i)]])
#define agg_mesh_spline_length(_w,_i)	\
  ((_w)->isp[(_i)+1]-(_w)->isp[(_i)])
#define agg_mesh_patch(_w,_i)         ((_w)->P[(_i)])
#define agg_mesh_surface(_w,_i)       ((_w)->S[(_i)])
#define agg_mesh_surface_number(_w)   ((_w)->ns)
#define agg_mesh_surface_blend(_w,_i) (&((_w)->B[(_i)]))
#define agg_mesh_surface_blend_number(_w) ((_w)->nb)
#define agg_mesh_element_size(_e)			\
  ( ((_e)[4] != 0) ? 5 :				\
    (((_e)[3] != 0) ? 4 : 3) )
#define agg_mesh_intersection(_w,_i)     ((_w)->inter[(_i)])
#define agg_mesh_intersection_number(_w) ((_w)->ni)
#define agg_mesh_intersection_points_start(_w,_i) \
  ((_w)->isect[2*(_i)+0])
#define agg_mesh_intersection_points_end(_w,_i)	\
  ((_w)->isect[2*(_i)+1])
#define agg_mesh_intersection_points_per_spline(_w,_i)	\
  ((_w)->ipps[(_i)])
#define agg_mesh_intersection_spline_number(_w,_i)	\
  ((_w)->insp[(_i)])

/**
 * @}
 */

/**
 * @{
 * 
 * @ingroup body
 */

#ifdef DOXYGEN
/** @typedef agg_body_t
 *
 * Data structure for bodies made up of collections of surfaces
 */
typedef agg_body_t ;
#else /*DOXYGEN*/
typedef struct _agg_body_t agg_body_t ;
struct _agg_body_t {
  agg_variable_t *g ;
  agg_expression_data_t *e ;
  agg_surface_t **S ;
  agg_patch_t **P ;
  char **names ;
  gint ng, ngmax, ns, nsmax ;
} ;

#define agg_body_global(_b,_i)          (&((_b)->g[(_i)]))
#define agg_body_globals(_b)            ((_b)->g)
#define agg_body_global_number(_b)      ((_b)->ng)
#define agg_body_global_number_max(_b)  ((_b)->ngmax)
#define agg_body_surface(_b,_i)         ((_b)->S[(_i)])
#define agg_body_surface_number(_b)     ((_b)->ns)
#define agg_body_surface_number_max(_b) ((_b)->nsmax)
#define agg_body_surface_last(_b)				\
  agg_body_surface((_b),(agg_body_surface_number((_b))-1))
#define agg_body_surface_name(_b,_i)    ((_b)->names[(_i)])
#define agg_body_patch(_b,_i)           ((_b)->P[(_i)])
#define agg_body_patch_last(_b)				\
  agg_body_patch((_b),(agg_body_surface_number((_b))-1))

#endif /*DOXYGEN*/

/**
 * @}
 */

typedef struct _agg_surface_workspace_t agg_surface_workspace_t ;
struct _agg_surface_workspace_t {
  agg_section_t *s ;
} ;

gint agg_bernstein_basis(gint n, gdouble x, gdouble *S, gdouble *dS) ;
gdouble agg_bernstein_basis_eval(gint n, gint r, gdouble x) ;
gdouble agg_bernstein_derivative_eval(gint n, gint r, gdouble x) ;

agg_section_t *agg_section_new(gint oumax, gint olmax) ;
gdouble agg_section_eval(agg_section_t *s, gdouble x) ;
gdouble agg_section_diff(agg_section_t *s, gdouble x) ;
gint agg_section_copy(agg_section_t *dest, agg_section_t *src) ;
agg_section_t *agg_section_duplicate(agg_section_t *s) ;
gint agg_section_set_circle(agg_section_t *s) ;
gint agg_section_set_ellipse(agg_section_t *s, gdouble th) ;
gint agg_section_set_aerofoil(agg_section_t *s, gdouble eta,
			      gdouble th, gdouble yte) ;
gint agg_section_parse(agg_section_t *s, char *name,
		       agg_variable_t *p, gint np) ;
gint agg_section_write(FILE *f, agg_section_t *s, agg_transform_t *T,
		       gint npts) ;
gint agg_section_format_write(FILE *f, agg_section_t *s, agg_transform_t *T,
			      char *fstr, char *estr, gint npts) ;
gint agg_sections_list(FILE *f) ;
gint agg_section_fit(agg_section_t *s,
		     gdouble *xu, gint xustr, gdouble *yu, gint yustr, gint npu,
		     gdouble *xl, gint xlstr, gdouble *yl, gint ylstr, gint npl,
		     gdouble n1, gdouble n2, gint nu, gint nl) ;

agg_transform_operator_t *agg_transform_operator_new(void) ;
agg_transform_t *agg_transform_new(gint nopmax) ;
gint agg_transform_variable_add(agg_transform_t *T,
				char *var, char *def, gdouble val) ;
gint agg_transform_add_global_variables(agg_transform_t *T,
					agg_variable_t *global, gint nglobal) ;
gint agg_transform_expressions_compile(agg_transform_t *T) ;
gint agg_transform_variables_eval(agg_transform_t *T) ;
gint agg_transform_variables_write(FILE *f, agg_transform_t *T,
				   gboolean write_defs) ;
gint agg_transform_operator_add(agg_transform_t *T, agg_operation_t op,
				gdouble umin, gdouble umax,
				gdouble *p,
				char **expr, char **dedu, char **dedv,
				gint np) ;
gint agg_transform_operators_write(FILE *f, agg_transform_t *T) ;
gint agg_transform_apply(agg_transform_t *T, gdouble *xin, gdouble *xout) ;
gint agg_transforms_list(FILE *f) ;

gint agg_transform_operator_apply(agg_transform_operator_t *op,
				  gdouble *xin, gdouble *xout) ;
gint agg_transform_operator_rotate(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv) ;
gint agg_transform_operator_translate(agg_operation_t op,
				      agg_variable_t *p, gint np,
				      gdouble *xin, gdouble *xout,
				      gdouble *dxdu, gdouble *dxdv) ;
gint agg_transform_operator_shrink(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv) ;
gint agg_transform_operator_scale(agg_operation_t op,
				  agg_variable_t *p, gint np,
				  gdouble *xin, gdouble *xout,
				  gdouble *dxdu, gdouble *dxdv) ;
gint agg_transform_operator_xscale(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv) ;
gint agg_transform_operator_yscale(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv) ;
				   
gint agg_transform_axes(agg_axes_t axes, gdouble *xin, gdouble *xout) ;
gint agg_transform_parse(agg_transform_t *T, agg_variable_t *p, gint np) ;
agg_axes_t agg_axes_parse(char *str) ;

gint agg_variable_write(FILE *f, agg_variable_t *v) ;
agg_expression_data_t *agg_expression_data_new(gint nemax) ;
gint agg_expression_data_variable_add(agg_expression_data_t *d,
				      agg_variable_t *v) ;
gpointer agg_expression_compile(char *e, agg_expression_data_t *d) ;
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
gint agg_surface_point_diff(agg_surface_t *S, gdouble u, gdouble v,
			    gdouble *x, gdouble *xu, gdouble *xv,
			    agg_surface_workspace_t *w) ;
agg_surface_workspace_t *agg_surface_workspace_new(void) ;

gdouble agg_naca_four(gdouble t, gdouble p, gdouble m, gdouble x) ;

agg_patch_t *agg_patch_new(void) ;
gint agg_patch_map(agg_patch_t *p, gdouble s, gdouble t,
		   gdouble *u, gdouble *v) ;
gint agg_patch_st_correct(agg_patch_t *P, gdouble *st) ;
gint agg_patch_parse(agg_patch_t *P, agg_variable_t *p, gint np) ;
gint agg_patch_surface_diff(agg_patch_t *P,
			    gdouble s, gdouble t,
			    gdouble *xu, gdouble *xv,
			    gdouble *xs, gdouble *xt) ;
gint agg_patch_point_diff(agg_surface_t *S, agg_patch_t *P,
			  gdouble s, gdouble t,
			  gdouble *x, gdouble *xs, gdouble *xt,
			  agg_surface_workspace_t *w) ;
gint agg_patch_edge_split(agg_patch_t *P,
			  gdouble s0, gdouble t0, gdouble s1, gdouble t1, 
			  gdouble a, gdouble *s, gdouble *t) ;
gint agg_patch_triangle_interp(agg_patch_t *P,
			       gdouble s0, gdouble t0,
			       gdouble s1, gdouble t1, 
			       gdouble s2, gdouble t2, 
			       gdouble u,  gdouble v,
			       gdouble *s, gdouble *t) ;

agg_intersection_t *agg_intersection_new(gint nstmax) ;
gboolean agg_surface_patch_trim(agg_surface_t *S1, agg_patch_t *P1, gdouble d1,
				agg_surface_t *S2, agg_patch_t *P2, gdouble d2,
				agg_surface_blend_t *B,
				agg_surface_workspace_t *w) ;
gint agg_intersection_curve_write(FILE *f, agg_intersection_t *inter,
				  agg_surface_workspace_t *w) ;
gint agg_intersection_bbox_set(agg_intersection_t *inter) ;
gint agg_intersection_resample(agg_intersection_t *inter,
			       gint nsp, gint pps,
			       agg_intersection_t *resample,
			       agg_surface_workspace_t *w) ;
/* gint agg_intersection_clip(agg_intersection_t *inter, */
/* 			   gint c, agg_patch_clip_t cut, */
/* 			   agg_patch_clipping_t *clip) ; */
/* gint agg_patch_clip_eval(agg_patch_clipping_t *c, gdouble u, */
/* 			 gdouble *s, gdouble *t) ; */
/* gint agg_clipping_orient(agg_patch_clipping_t *c1, agg_patch_t *P1, */
/* 			 agg_surface_t *S1, */
/* 			 agg_patch_clipping_t *c2, agg_patch_t *P2, */
/* 			 agg_surface_t *S2, */
/* 			 agg_surface_workspace_t *w) ; */

agg_mesh_t *agg_mesh_new(gint npmax, gint nspmax, gint nemax) ;
gint agg_mesh_surface_make(agg_mesh_t *w,
			   agg_surface_t *S, agg_patch_t *P,
			   agg_intersection_t *inter,
			   gint nsec, gint nseg, gint pps,
			   agg_surface_workspace_t *ws) ;
gint agg_mesh_write_gmsh(FILE *f, agg_mesh_t *w,
			 char *len, gint offp, gint offsp, gint offs,
			 gboolean opencascade) ;
gint agg_mesh_spline_ends(agg_mesh_t *w, gint s, gint *p0, gint *p1) ;
gint agg_mesh_spline_interp_points(agg_mesh_t *w, gint s,
				   gint p0, gint p1, gint pps,
				   agg_surface_workspace_t *ws) ;
gboolean agg_mesh_splines_connected(agg_mesh_t *w,
				    gint i0, gint i1, gint *p) ;
gint agg_mesh_element_add(agg_mesh_t *w, gint s0, gint s1, gint s2, gint s3) ;
gint agg_mesh_surface_add_triangle(agg_mesh_t *msh, gint isurf,
				   char *args, gint pps,
				   agg_surface_workspace_t *w) ;
gint agg_mesh_surface_blend_add(agg_mesh_t *m, gint iB, gint nsec, gint pps,
				agg_surface_workspace_t *w) ;
gint agg_mesh_surface_add_grid(agg_mesh_t *m, gint iS,
			       gint nsec, gint nsp, gint pps,
			       agg_surface_workspace_t *w) ;
gint agg_mesh_spline_from_endpoints(agg_mesh_t *w, gint p0, gint p1) ;
gint agg_mesh_surface_point_add(agg_mesh_t *w, gint surf,
				gdouble s, gdouble t,
				agg_surface_workspace_t *ws) ;
gint agg_mesh_body(agg_mesh_t *m, agg_body_t *b, gint pps,
		   agg_surface_workspace_t *w) ;
gint agg_mesh_element_nodes(agg_mesh_t *m, gint e,
			    gint *nodes, gint *nnodes, gint *s) ;

agg_body_t *agg_body_new(gint ngmax, gint nsmax) ;
gint agg_body_global_add(agg_body_t *b, char *var, char *def, gdouble val) ;
gint agg_body_globals_write(FILE *f, agg_body_t *b) ;
gint agg_body_globals_compile(agg_body_t *b) ;
gint agg_body_globals_eval(agg_body_t *b) ;
gint agg_body_read(agg_body_t *b, char *file, gboolean echo) ;
gint agg_body_surface_add(agg_body_t *b, agg_surface_t *S, agg_patch_t *P) ;
gint agg_body_surfaces_list(FILE *f, agg_body_t *b) ;

agg_grid_t agg_grid_parse(char *str) ;

agg_surface_blend_t *agg_surface_blend_new(void) ;
gint agg_surface_blend_evaluate(agg_surface_blend_t *B,
				gdouble s, gdouble t, gdouble *x,
				agg_surface_workspace_t *w) ;
gint agg_hermite_eval(gdouble s, gdouble *H) ;

gint agg_library_section_add(char *name, char *description,
			     agg_section_t *s) ;
agg_section_t *agg_library_section_lookup(char *name, char **description) ;
gint agg_library_sections_list(FILE *f, gboolean write_description) ;
gint agg_library_section_write(FILE *f, char *name, char *description,
			       agg_section_t *s) ;
gint agg_library_read(FILE *f) ;
gint agg_library_write(FILE *f) ;

gint agg_curve_eval(agg_curve_t *c, gdouble x, gdouble *s, gdouble *t) ;
gboolean agg_curve_point_orientation(agg_curve_t *c, gdouble del,
				     gdouble s, gdouble t) ;
gint agg_curve_plane_normal(agg_curve_t *c, agg_surface_t *S, agg_patch_t *P,
			    gdouble *n, agg_surface_workspace_t *w) ;

agg_chebyshev_t *agg_chebyshev_new(gint Nmax, gint nc) ;
gint agg_chebyshev_eval(agg_chebyshev_t *C, gdouble x, gdouble *f) ;
gint agg_chebyshev_surface_section(agg_chebyshev_t *C,
				   agg_surface_t *S, gdouble u,
				   gint N, agg_surface_workspace_t *w) ;
gint agg_chebyshev_surface_section_adaptive(agg_chebyshev_t *C,
					    agg_surface_t *S, gdouble u,
					    gint N, gdouble tol, gdouble dmin,
					    agg_surface_workspace_t *w) ;
gdouble agg_chebyshev_interval_shortest(agg_chebyshev_t *C) ;

#endif /*__AGG_H_INCLUDED__*/

