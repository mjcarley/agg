/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2021, 2022 Michael Carley
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

#define AGG_SYMBOL_GLOBAL       G_TOKEN_LAST + 1
#define AGG_SYMBOL_DISTRIBUTION G_TOKEN_LAST + 2
#define AGG_SYMBOL_TRANSFORM    G_TOKEN_LAST + 3
#define AGG_SYMBOL_SHAPE        G_TOKEN_LAST + 4
#define AGG_SYMBOL_AXES         G_TOKEN_LAST + 5
#define AGG_SYMBOL_GRID         G_TOKEN_LAST + 6
#define AGG_SYMBOL_BODY         G_TOKEN_LAST + 7
#define AGG_SYMBOL_MAX          G_TOKEN_LAST + 7

typedef enum
  {
   AGG_SHAPE_AEROFOIL = 0, /**< thickness distribution plus camber */
   AGG_SHAPE_ELLIPSE = 1   /**< generic closed curve */
  } agg_shape_type_t ;

#define AGG_CROWD_BODY_NUMBER_MAX 64

#define AGG_PARSER_PARAMETER_NUMBER_MAX 128
#define AGG_PARSER_PARAMETER_RESERVED_S   0

typedef struct _agg_parser_t agg_parser_t ;

struct _agg_parser_t {
  gpointer vars ; /*variables for tinyexpr*/
  gpointer expr ; /*compiled expressions for tinyexpr*/
  gdouble values[AGG_PARSER_PARAMETER_NUMBER_MAX] ; /*variable values*/
  gboolean isexpr[AGG_PARSER_PARAMETER_NUMBER_MAX] ; /*expression or const*/
  gboolean global_set ; 
  gint nvars ; 
} ;

#define agg_parser_global_set(_p) ((_p)->global_set)

/** 
 * @ingroup types
 * Enumerated type for spacing in discretizations
 *
 */

typedef enum
  {
   AGG_SPACING_LINEAR  = 0,  /**< \f$t\f$ equal intervals over range */
   AGG_SPACING_COSINE  = 1,  /**< \f$\mathrm{sgn}(t)(1-\cos \pi t)/2\f$ */
   AGG_SPACING_HALFCOS = 2,  /**< \f$\mathrm{sgn}(t)(1-\cos \pi t/2)\f$ */
   AGG_SPACING_HALFSIN = 3   /**< \f$\sin \pi t/2\f$ */
  } agg_spacing_t ;

typedef struct _agg_shape_t agg_shape_t ;

struct _agg_shape_t {
  agg_shape_type_t type ;
  gboolean closed ; 
  gdouble *s,                /** < Bernstein polynomial coefficients */
    b[16],                   /** < break points for piecewise fits  */
    n1[8], n2[8] ;           /** < exponents */
  gint nsmax,                /** < maximum number of entries in s */
    i[8],                    /** < indices into s for different shape types */
    nb ;                     /** < number of breakpoints in piecewise fit */
} ;

#define agg_shape_break_number(_s)   ((_s)->nb)
#define agg_shape_break_lower(_s,_i) ((_s)->b[2*(_i)+0])
#define agg_shape_break_upper(_s,_i) ((_s)->b[2*(_i)+1])
/* #define agg_shape_ */

#define AGG_EXPRESSION_LENGTH_MAX 64

typedef struct _agg_expression_t agg_expression_t ;

struct _agg_expression_t {
  gchar def[AGG_EXPRESSION_LENGTH_MAX] ;
  gpointer expr ;
} ;

#define AGG_TRANSFORM_FUNCTION_NUMBER   16
#define AGG_TRANSFORM_PARAMETER_NUMBER 128

typedef gint (*agg_local_transform_func_t)(gdouble *x, gdouble *p,
					   gdouble *y) ;

typedef struct _agg_local_transform_t agg_local_transform_t ;

struct _agg_local_transform_t {
  agg_local_transform_func_t func[AGG_TRANSFORM_FUNCTION_NUMBER] ;
  gint nt, p1[AGG_TRANSFORM_FUNCTION_NUMBER+1] ;
  gdouble p[AGG_TRANSFORM_PARAMETER_NUMBER] ;
  gboolean isexpr[AGG_TRANSFORM_PARAMETER_NUMBER] ;
  agg_expression_t expr[AGG_TRANSFORM_PARAMETER_NUMBER] ;
} ;

#define agg_local_transform_parameter_index(_t,_i,_j)	\
  ((_t)->p1[(_i)] + (_j))

typedef struct _agg_distribution_t agg_distribution_t ;

struct _agg_distribution_t {
  gint
    ns, /**< number of shapes/sections in distribution */
    nsmax, /**< maximum number of shapes/sections in distribution */
    axes[3] ; /**< permutation of axes (coordinate indices) */
  gdouble smin, smax, *ts, *w ;
  gboolean invert ; /**< when meshing reorient triangles to invert normals */
  agg_shape_t **sh ;
  agg_local_transform_t *t ;
  agg_spacing_t sg, sm, ss ; /*spacing for geometry and mesh
			       generation (geometry, mesh, shape)*/
} ;

#define agg_distribution_shape(_d,_i) ((_d)->sh[(_i)])
#define agg_distribution_station_number(_d) ((_d)->ns)
#define agg_distribution_parameter_min(_d) ((_d)->smin)
#define agg_distribution_parameter_max(_d) ((_d)->smax)
#define agg_distribution_invert(_d) ((_d)->invert)

typedef enum
  {
   AGG_GRID_NONE          = 0,
   AGG_GRID_LINEAR        = 1,
   AGG_GRID_SPHERICAL     = 2,
   AGG_GRID_HEMISPHERICAL = 3,
   AGG_GRID_TUBE          = 4
  } agg_grid_topology_t ;

typedef struct _agg_grid_t agg_grid_t ;
typedef gint (*agg_grid_area_interp_func_t)(agg_grid_t *g, gdouble *uv,
					    gdouble s, gdouble t,
					    gdouble *u, gdouble *v) ;
typedef gint (*agg_grid_line_interp_func_t)(agg_grid_t *g, gdouble *uv,
					    gdouble s,
					    gdouble *u, gdouble *v) ;

struct _agg_grid_t {
  agg_grid_topology_t t ;
  agg_grid_area_interp_func_t interp_area ;
  agg_grid_line_interp_func_t interp_line ;
  gint np, npmax, nt, ntmax, *tri ;
  gdouble *uv ;
  gboolean invert ;
} ;

#define agg_grid_topology(_g)            ((_g)->t)
#define agg_grid_point_number(_g)        ((_g)->np)
#define agg_grid_point_number_max(_g)    ((_g)->npmax)
#define agg_grid_triangle_number(_g)     ((_g)->nt)
#define agg_grid_triangle_number_max(_g) ((_g)->ntmax)
#define agg_grid_triangle(_g,_i)         (&((_g)->tri[3*(_i)]))
#define agg_grid_point_u(_g,_i)          ((_g)->uv[2*(_i)+0])
#define agg_grid_point_v(_g,_i)          ((_g)->uv[2*(_i)+1])
#define agg_grid_invert(_g)              ((_g)->invert)

#define AGG_BODY_DISTRIBUTION_NUMBER_MAX 64

typedef struct _agg_body_t agg_body_t ;

struct _agg_body_t {
  gint nd ;
  gchar *names[AGG_BODY_DISTRIBUTION_NUMBER_MAX] ;
  agg_distribution_t *d[AGG_BODY_DISTRIBUTION_NUMBER_MAX] ;
  agg_parser_t *p ;
  agg_grid_t *g ;
} ;

#define agg_body_distribution_number(_b) ((_b)->nd)
#define agg_body_distribution(_b,_i)     (((_b)->d[(_i)]))
#define agg_body_parser(_b)              ((_b)->p)
#define agg_body_grid(_b)              ((_b)->g)

typedef struct _agg_crowd_t agg_crowd_t ;

struct _agg_crowd_t {
  agg_body_t *b[AGG_CROWD_BODY_NUMBER_MAX] ;
  agg_parser_t *p ;
  gint nb ;
} ;

#define agg_crowd_body_number(_c) ((_c)->nb)
#define agg_crowd_body(_c,_i)     ((_c)->b[(_i)])

typedef struct _agg_function_call_t agg_function_call_t ;

struct _agg_function_call_t {
  gchar func[64] ;
  gint na ;
  gpointer expr[16] ;
  gdouble  x[16] ;
  guint c ;
} ;

#define AGG_MESH_DISTRIBUTION_NUMBER_MAX 32

typedef struct _agg_bbox_t agg_bbox_t ;

struct _agg_bbox_t {
  gdouble xmin, xmax, ymin, ymax, zmin, zmax ;
} ;

typedef struct _agg_mesh_t agg_mesh_t ;

struct _agg_mesh_t {
  agg_grid_t *g ;
  agg_distribution_t *dist[AGG_MESH_DISTRIBUTION_NUMBER_MAX] ;
  gint np, npmax, nt, ntmax, nptags, nttags, ndat, *tri, *ptags, ndist ;
  gdouble *x ;
  agg_bbox_t *bboxes ;
} ;

#define agg_mesh_distribution(_m,_i) ((_m)->dist[(_i)])
#define agg_mesh_distribution_number(_m) ((_m)->ndist)
					  
#define agg_mesh_point_number(_m) ((_m)->np)
#define agg_mesh_triangle_number(_m) ((_m)->nt)
#define agg_mesh_triangle(_m,_i) &((_m)->tri[(3+1+(_m)->nttags)*(_i)])
#define agg_mesh_bbox(_m,_i)     &((_m)->bboxes[(_i)])

#define agg_mesh_point(_m,_i)   (&((_m)->x[(_i)*(3+2+(_m)->ndat)]))
#define agg_mesh_point_x(_m,_i) ((_m)->x[(_i)*(3+2+(_m)->ndat)+0])
#define agg_mesh_point_y(_m,_i) ((_m)->x[(_i)*(3+2+(_m)->ndat)+1])
#define agg_mesh_point_z(_m,_i) ((_m)->x[(_i)*(3+2+(_m)->ndat)+2])
#define agg_mesh_point_u(_m,_i) ((_m)->x[(_i)*(3+2+(_m)->ndat)+3])
#define agg_mesh_point_v(_m,_i) ((_m)->x[(_i)*(3+2+(_m)->ndat)+4])
#define agg_mesh_point_data(_m,_i,_j) ((_m)->x[(_i)*(3+2+(_m)->ndat)+5+(_j)])

#define agg_mesh_point_distribution(_m,_i) ((_m)->ptags[(_i)*((_m)->nptags+1)])
#define agg_mesh_triangle_distribution(_m,_i) \
  ((_m)->tri[(3+1+(_m)->nttags)*(_i)+3+0])
#define agg_mesh_distribution(_m,_i)  ((_m)->dist[(_i)])
#define agg_mesh_point_tag(_m,_i,_j)		\
  ((_m)->ptags[(_i)*((_m)->nptags+1)+(_j)+1])
#define agg_mesh_point_tag_number(_m)  ((_m)->nptags)
#define agg_mesh_point_data_number(_m)  ((_m)->ndat)
#define agg_mesh_triangle_tag_number(_m)  ((_m)->nttags)
#define agg_mesh_triangle_tag(_m,_i,_j)		\
  ((_m)->tri[(3+1+(_m)->nttags)*(_i)+3+1+(_j)])
#define agg_mesh_grid(_m) ((_m)->g)


typedef struct _agg_workspace_t agg_workspace_t ;

struct _agg_workspace_t {
  agg_shape_t *s ;
  agg_local_transform_t *t ;
} ;

#define agg_workspace_shape(_w) ((_w)->s)
#define agg_workspace_local_transform(_w) ((_w)->t)

gdouble agg_bernstein_basis_eval(gint n, gint r, gdouble x) ;
gint agg_bernstein_basis(gint n, gdouble x, gdouble *S, gdouble *dS) ;
gdouble agg_shape_eval(agg_shape_t *s, gdouble x, gint lim) ;
gdouble agg_shape_aerofoil_eval(agg_shape_t *s, gdouble x) ;
gdouble agg_shape_ellipse_eval(agg_shape_t *s, gdouble x, gint lim) ;
agg_shape_t *agg_shape_alloc(gint nsmax) ;
gint agg_shape_naca(agg_shape_t *c, gdouble t) ;
gint agg_shape_ellipse(agg_shape_t *s) ;
gint agg_shape_fit(agg_shape_t *c,
		   gdouble *xu, gdouble *yu, gint nu,
		   gdouble *xl, gdouble *yl, gint nl,
		   gdouble n1, gdouble n2, gdouble d,
		   gint n, gboolean closed, gdouble *work) ;
gint agg_shape_interval_add(agg_shape_t *s,
			    gdouble x1, gdouble x2,
			    gdouble *x, gdouble *y, gint np,
			    gboolean upper, 
			    gdouble n1, gdouble n2,
			    gint n, gdouble *work) ;
gint agg_shape_root_fit(agg_shape_t *s,
			gdouble n1, gdouble n2,
			gdouble hu, gdouble su, gdouble du, gdouble bu,
			gdouble hl, gdouble sl, gdouble dl, gdouble bl,
			gint n, gdouble *work) ;

gint agg_shape_fit_naca_four(agg_shape_t *s, gint n,
			     gdouble th, gdouble p, gdouble m,
			     gboolean sharp,
			     gint nu, gint nl, gdouble *work) ;
gint agg_shape_parse(agg_shape_t *s, gchar *type, gdouble *p, gint np) ;
gint agg_shape_copy(agg_shape_t *dest, agg_shape_t *src) ;

gdouble agg_naca_four(gdouble t, gdouble p, gdouble m, gdouble x) ;

gint agg_functions_eval(agg_function_call_t *func, agg_parser_t *p) ;

agg_distribution_t *agg_distribution_alloc(gint nsmax) ;
gint agg_distribution_add_shape(agg_distribution_t *d,
				gdouble t, agg_shape_t *s) ;
gint agg_distribution_interpolation_weights(agg_distribution_t *d) ;
gint agg_distribution_interpolate_shape(agg_distribution_t *d, gdouble s,
					agg_shape_t *sh) ;
gint agg_distribution_point_eval(agg_distribution_t *d,
				 gdouble u, gdouble v,
				 agg_parser_t *p,
				 agg_shape_t *s,
				 gdouble *x) ;
gint agg_distribution_point_normal_eval(agg_distribution_t *d,
					gdouble u, gdouble v,
					agg_parser_t *p,
					agg_shape_t *sh,
					gdouble *x, gdouble *n, gdouble *J) ;

gint agg_distribution_axes_parse(gchar *s, gint *a) ;

gint agg_distribution_mesh(agg_distribution_t *d,
			   gdouble smin, gdouble smax, gint ns,
			   agg_spacing_t ss,
			   gdouble tmin, gdouble tmax, gint nt,
			   agg_spacing_t st,
			   /* agg_shape_t *sh, agg_local_transform_t *T, */
			   agg_parser_t *p,			   
			   agg_mesh_t *m, agg_workspace_t *w) ;
gint agg_body_mesh_grid(agg_body_t *b, agg_grid_t *g,
			agg_mesh_t *m, agg_workspace_t *w) ;

gint agg_mesh_tri_write(FILE *f, agg_mesh_t *m) ;
gint agg_mesh_points_write(FILE *f, agg_mesh_t *m) ;
agg_mesh_t *agg_mesh_alloc(gint np, gint nt, gint nptags, gint nttags,
			   gint ndat) ;
gint agg_mesh_init(agg_mesh_t *m) ;
gint agg_mesh_element_interp(agg_mesh_t *m, gint i,
			     gdouble s, gdouble t,
			     agg_parser_t *p, agg_shape_t *sh,
			     gdouble *x) ;
gint agg_mesh_element_interp_normal(agg_mesh_t *m, gint i,
				    gdouble s, gdouble t,
				    agg_parser_t *p, agg_shape_t *sh,
				    gdouble *x, gdouble *n, gdouble *J) ;
gint agg_mesh_bounding_boxes(agg_mesh_t *m) ;

gdouble agg_spacing_eval(gdouble tmin, gdouble tmax, gint nt,
			 agg_spacing_t s, gint i) ;

agg_parser_t *agg_parser_alloc(void) ;
gint agg_parser_declaration(agg_parser_t *p, gchar *s) ;
gint agg_parser_declarations_write(FILE *f, agg_parser_t *p) ;
gint agg_parser_expressions_evaluate(agg_parser_t *p) ;
gint agg_parser_variable_add(agg_parser_t *p, gchar *s, gdouble x) ;
gint agg_parser_expression_add(agg_parser_t *p, gchar *v, gchar *w) ;
gint agg_parser_function_add(agg_parser_t *p, gchar *v, gpointer func,
			     gint na) ;

GScanner *agg_scanner_alloc(void) ;
gint agg_parser_body_read(GScanner *scanner, agg_parser_t *p, agg_body_t *b) ;
gint agg_parser_constant_parse(gchar *s, guint *v) ;
gint agg_parser_crowd_read(GScanner *scanner, agg_crowd_t *c) ;

agg_body_t *agg_body_alloc(void) ;
gint agg_body_distribution_locate_u(agg_body_t *b, gdouble u) ;
gint agg_body_point_eval(agg_body_t *b, gdouble u, gdouble v, gdouble *x,
			 agg_workspace_t *w) ;
gint agg_body_distributions_list(FILE *f, agg_body_t *b) ;
gint agg_body_distribution_add(agg_body_t *b, agg_distribution_t *d,
			       gchar *name) ;
gint agg_body_parameter_limits(agg_body_t *b, gdouble *umin, gdouble *umax) ;

gint agg_local_transform_init(agg_local_transform_t *t) ;
agg_local_transform_t *agg_local_transform_alloc(void) ;

gint agg_local_transform_rotate(gdouble *x, gdouble *p, gdouble *y) ;
gint agg_local_transform_shift(gdouble *x, gdouble *p, gdouble *y) ;
gint agg_local_transform_shrink(gdouble *x, gdouble *p, gdouble *y) ;
gint agg_local_transform_scale(gdouble *x, gdouble *p, gdouble *y) ;
gint agg_local_transform_parse(agg_local_transform_t *T, gchar *type,
			       gdouble *p, gint np) ;
gint agg_local_transform_apply(agg_local_transform_t *T, gdouble *y) ;
gint agg_local_transform_eval_parameters(agg_local_transform_t *t) ;
gint agg_local_transform_set_parameters(agg_local_transform_t *T,
					agg_local_transform_t *t) ;
gint agg_local_transform_set_expression(agg_local_transform_t *t,
					gint i, gchar *expr,
					agg_parser_t *p) ;

gdouble agg_function_tipright(gdouble eta, gdouble t0, gdouble t) ;
gdouble agg_function_tipleft(gdouble eta, gdouble t0, gdouble t) ;

agg_grid_t *agg_grid_alloc(gint np, gint nt) ;
gint agg_grid_init(agg_grid_t *g) ;
gint agg_grid_square(agg_grid_t *g,
		     gdouble umin, gdouble umax, gint nu, agg_spacing_t su,
		     gdouble vmin, gdouble vmax, gint nv, agg_spacing_t sv) ;
gint agg_grid_spherical(agg_grid_t *g, gint refine) ;
gint agg_grid_hemispherical(agg_grid_t *g, gint refine) ;
gint agg_grid_linear(agg_grid_t *g,
		     gdouble umin, gdouble umax, gint nu, agg_spacing_t su,
		     gdouble vmin, gdouble vmax, gint nv, agg_spacing_t sv) ;
gint agg_grid_tube(agg_grid_t *g, gint nu,
		   gdouble vmin, gdouble vmax, gint nv) ;
gint agg_grid_cone(agg_grid_t *g, gint nu, gint nv) ;
gint agg_grid_area_interpolate(agg_grid_t *g, gdouble *uv,
			       gdouble s, gdouble t,
			       gdouble *u, gdouble *v) ;
gint agg_grid_line_interpolate(agg_grid_t *g, gdouble *uv,
			       gdouble s,
			       gdouble *u, gdouble *v) ;
gint agg_grid_element_interpolate(agg_grid_t *g, gint i,
				  gdouble s, gdouble t,
				  gdouble *u, gdouble *v) ;
gint agg_grid_parse(agg_grid_t *g, gchar *name, gchar *expr[],
		    gdouble *p, gint np) ;

agg_workspace_t *agg_workspace_alloc(gint ns) ;

#endif /*__AGG_H_INCLUDED__*/

