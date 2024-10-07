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

typedef struct _agg_affine_t agg_affine_t ;
struct _agg_affine_t {
  gint order, order_max ;
  gpointer *compiled ;
  char **expr ;
  gdouble *matrix ;
} ;

#define agg_affine_order(_A)     ((_A)->order)
#define agg_affine_order_max(_A) ((_A)->order_max)
#define agg_affine_matrix(_A,_order) (&((_A)->matrix[16*(_order)]))
#define agg_affine_matrix_value(_A,_order,_i,_j)	\
  ((_A)->matrix[16*(_order)+4*(_i)+(_j)])
#define agg_affine_expression(_A,_order,_i,_j)	\
  ((_A)->expr[16*(_order)+4*(_i)+(_j)])
#define agg_affine_compiled(_A,_order,_i,_j)	\
  ((_A)->compiled[16*(_order)+4*(_i)+(_j)])

typedef gint (*agg_affine_func_t)(agg_affine_t *A,
				  gdouble *p, char **expr, gint np) ;

typedef enum {
  AGG_TRIANGULATION_UNDEFINED    = 0,
  AGG_TRIANGULATION_GRID_REGULAR = 1
} agg_triangulation_t ;

typedef enum {
  AGG_SAMPLING_UNDEFINED     = 0,
  AGG_SAMPLING_LINEAR        = 1,
  AGG_SAMPLING_COSINE        = 2,
  AGG_SAMPLING_COSINE_DOUBLE = 3  
} agg_sampling_t ;

typedef struct _agg_triangulation_settings_t agg_triangulation_settings_t ;

struct _agg_triangulation_settings_t {
  agg_triangulation_t tri ;
  agg_sampling_t sample_s, sample_t ;
  gint ns, nt, pps ; 
  gdouble area ;
  gboolean invert ;
} ;

#define agg_triangulation_type(_t)                 ((_t)->tri)
#define agg_triangulation_sampling_s(_t)           ((_t)->sample_s)
#define agg_triangulation_sampling_t(_t)           ((_t)->sample_t)
#define agg_triangulation_section_number(_t)       ((_t)->ns)
#define agg_triangulation_section_point_number(_t) ((_t)->nt)
#define agg_triangulation_points_per_spline(_t)    ((_t)->pps)
#define agg_triangulation_invert(_t)               ((_t)->invert)

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
  /* agg_transform_operator_t **op ; */
  agg_affine_t **affine ;
  gint /* nop, nopmax,  */ nv, n_affine, n_affine_max  ;
  gdouble matrix[16] ;
} ;
#endif /*DOXYGEN*/

#ifdef DOXYGEN
/* /\** */
/*  *  @brief \f$i\f$th operation of an ::agg_transform_t */
/*  *\/ */
/* #define agg_transform_operator(T,i)         */
/* /\** */
/*  *  @brief number of operations in an ::agg_transform_t */
/*  *\/ */
/* #define agg_transform_operator_number(T)     */
/* /\** */
/*  *  @brief maximum number of operations in an  ::agg_transform_t */
/*  *\/ */
/* #define agg_transform_operator_number_max(T) */
/**
 *  @brief \f$i\f$th variable in an  ::agg_transform_t
 */
#define agg_transform_variable(T,i)        
/**
 *  @brief number of variables in an  ::agg_transform_t
 */
#define agg_transform_variable_number(T)    
#else /*DOXYGEN*/
/* #define agg_transform_operator(_T,_i)           ((_T)->op[(_i)]) */
/* #define agg_transform_operator_number(_T)       ((_T)->nop) */
/* #define agg_transform_operator_number_max(_T)   ((_T)->nopmax) */
#define agg_transform_variable(_T,_i)           (&((_T)->v[(_i)]))
#define agg_transform_variable_number(_T)       ((_T)->nv)
#define agg_transform_affine(_T,_i)             ((_T)->affine[(_i)])
#define agg_transform_affine_number(_T)         ((_T)->n_affine)
#define agg_transform_affine_number_max(_T)     ((_T)->n_affine_max)
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
  agg_variable_t vmin, vmax ;
  agg_transform_t *T ;
  gboolean close ;
} ;
#endif /*DOXYGEN*/

#define agg_surface_section(_s,_i)          (&((_s)->s[(_i)]))
#define agg_surface_section_number(_s)      ((_s)->ns)
#define agg_surface_section_number_max(_s)  ((_s)->nsmax)
#define agg_surface_u(_s,_i)                ((_s)->u[(_i)])
#define agg_surface_weight(_s,_i)           ((_s)->w[(_i)])
#define agg_surface_umin(_s)                ((_s)->umin)
#define agg_surface_umax(_s)                ((_s)->umax)
#define agg_surface_vmin(_s)                (&((_s)->vmin))
#define agg_surface_vmax(_s)                (&((_s)->vmax))
#define agg_surface_transform(_s)           ((_s)->T)
#define agg_surface_axes(_s)                ((_s)->axes)
#define agg_surface_grid(_s)                ((_s)->grid)
#define agg_surface_grid_section_number(_s) ((_s)->nsec)
#define agg_surface_grid_spline_number(_s)  ((_s)->nsp)
#define agg_surface_grid_element_area(_s)   ((_s)->area)
#define agg_surface_grid_subdivision(_s)    ((_s)->sub)
#define agg_surface_close(_s)               ((_s)->close)

/**
 *  @}
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
  gdouble *p ;
  gint np, npmax, *sp, nsp, nspmax, *isp, *e, ne, nemax, *ptags,
    isect[32], ipps[32], insp[32], ni, ns, nb ;
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
  agg_triangulation_settings_t *settings, settings_default ;
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
#define agg_body_triangulation_settings(_b,_i) (&((_b)->settings[(_i)]))
#define agg_body_triangulation_settings_default(_b) (&((_b)->settings_default))

#endif /*DOXYGEN*/

/**
 * @}
 */

typedef struct _agg_surface_workspace_t agg_surface_workspace_t ;
struct _agg_surface_workspace_t {
  agg_section_t *s, *ds ;
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
gint agg_transform_evaluate(agg_transform_t *T, gdouble u,
			    gint order, gdouble *matrix) ;

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
gint agg_surface_point_diff_numerical(agg_surface_t *S, gdouble u, gdouble v,
				      gdouble *x, gdouble *xu, gdouble *xv,
				      agg_surface_workspace_t *w) ;
agg_surface_workspace_t *agg_surface_workspace_new(void) ;
gint agg_surface_section_diff(agg_surface_t *S, gdouble u,
			      agg_section_t *s, agg_section_t *ds) ;

gdouble agg_naca_four(gdouble t, gdouble p, gdouble m, gdouble x) ;


agg_mesh_t *agg_mesh_new(gint npmax, gint nspmax, gint nemax) ;
gint agg_mesh_write_gmsh(FILE *f, agg_mesh_t *w,
			 char *len, gint offp, gint offsp, gint offs,
			 gboolean opencascade) ;
gint agg_mesh_spline_ends(agg_mesh_t *w, gint s, gint *p0, gint *p1) ;
gint agg_mesh_spline_interp_points(agg_mesh_t *w, gint s,
				   gint p0, gint p1, gint pps,
				   agg_surface_workspace_t *ws) ;
gboolean agg_mesh_splines_connected(agg_mesh_t *w,
				    gint i0, gint i1, gint *p) ;
gboolean agg_mesh_spline_zero_length(agg_mesh_t *m, gint s) ;
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
gint agg_body_surface_add(agg_body_t *b, agg_surface_t *S) ;
gint agg_body_surfaces_list(FILE *f, agg_body_t *b) ;

agg_grid_t agg_grid_parse(char *str) ;

gint agg_hermite_eval(gdouble s, gdouble *H) ;

gint agg_library_section_add(char *name, char *description,
			     agg_section_t *s) ;
agg_section_t *agg_library_section_lookup(char *name, char **description) ;
gint agg_library_sections_list(FILE *f, gboolean write_description) ;
gint agg_library_section_write(FILE *f, char *name, char *description,
			       agg_section_t *s) ;
gint agg_library_read(FILE *f) ;
gint agg_library_write(FILE *f) ;

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

gint agg_transform_affine_add(agg_transform_t *T, agg_affine_t *A) ;
agg_affine_t *agg_affine_new(gint order_max) ;
gint agg_affine_point_transform(gdouble *y, gdouble *A, gdouble *x) ;

gint agg_affine_expression_add(agg_affine_t *A, gint order,
			       gint i, gint j, gdouble val, char *expr) ;
gint agg_affine_expressions_compile(agg_affine_t *A, agg_expression_data_t *e) ;
gint agg_affine_matrices_evaluate(agg_affine_t *A) ;
gint agg_affine_identity(agg_affine_t *A) ;
gint agg_affine_zero(agg_affine_t *A) ;
gint agg_affine_translation(agg_affine_t *A,
			    gdouble *dx, char **expr, gint ns) ;
gint agg_affine_rotation_x(agg_affine_t *A, gdouble *th, char **expr, gint np) ;
gint agg_affine_rotation_y(agg_affine_t *A, gdouble *th, char **expr, gint np) ;
gint agg_affine_rotation_z(agg_affine_t *A, gdouble *th, char **expr, gint np) ;
gint agg_affine_scale(agg_affine_t *A, gdouble *s, char **expr, gint ns) ;
gint agg_affine_parse(agg_affine_t *A, agg_variable_t *p, gint np) ;
gint agg_affine_axes(agg_affine_t *A, agg_axes_t axes) ;

gint agg_mesh_surface_triangulate(agg_mesh_t *m, gint isurf,
				  agg_triangulation_settings_t *settings,
				  agg_surface_workspace_t *w) ;
gint agg_affine_differentiate(agg_affine_t *A, char *var) ;
gint agg_affine_write(FILE *f, agg_affine_t *A) ;
gint agg_affine_list(FILE *f, char *head, char *tail) ;

#endif /*__AGG_H_INCLUDED__*/

