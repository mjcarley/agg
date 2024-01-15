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

static const struct {
  gchar *name ;
  agg_transform_operator_func_t func ;
  agg_operation_t op ;
  gint np ;
} transform_list[] =
  {
    {"undefined", NULL,                             AGG_TRANSFORM_UNDEFINED, 0},
    {"rotate",    agg_transform_operator_rotate,    AGG_TRANSFORM_ROTATE   , 3},
    {"shrink",    agg_transform_operator_shrink,    AGG_TRANSFORM_SHRINK   , 3},
    {"translate", agg_transform_operator_translate, AGG_TRANSFORM_TRANSLATE, 3},
    {"scale",     agg_transform_operator_scale,     AGG_TRANSFORM_SCALE    , 1},
    {"xscale",    agg_transform_operator_xscale,    AGG_TRANSFORM_SCALE_X  , 1},
    {"yscale",    agg_transform_operator_yscale,    AGG_TRANSFORM_SCALE_Y  , 1},
    {NULL,        NULL, -1, -1}
  } ;

static const struct {
  gchar *name ;
  agg_axes_t axes ;
} axes_list[] =
  {
    {"xyz",    AGG_AXES_PX_PY_PZ},
    {"+x+y+z", AGG_AXES_PX_PY_PZ},

    {"yzx",    AGG_AXES_PY_PZ_PX},
    {"+y+z+x", AGG_AXES_PY_PZ_PX},

    {"zxy",    AGG_AXES_PZ_PX_PY},
    {"+z+x+y", AGG_AXES_PZ_PX_PY},
    
    {"zyx",    AGG_AXES_PZ_PY_PX},
    {"+z+y+x", AGG_AXES_PZ_PY_PX},

    {"xy-z",   AGG_AXES_PX_PY_MZ},
    {"+x+y-z", AGG_AXES_PX_PY_MZ},
    
    {NULL,  -1}
  } ;

/*generating derivatives of operator parameter variables*/
#ifdef HAVE_LIBMATHEVAL
#include <matheval.h>
static void parameter_set_derivative(agg_transform_operator_t *tr,
				     agg_variable_t *v,
				     gchar **expr, gchar **de, gint i,
				     gchar *var)
{
  gpointer eval, diff ;

  if ( de != NULL ) {
    if ( de[i] != NULL )
      g_error("%s: overriding derivative evaluation not implemented",
	      __FUNCTION__) ;
  }
  
  if ( expr[i] == NULL ) {
    /*expression is a constant*/
    agg_variable_definition(v) = NULL ;
    agg_variable_value(v) = 0.0 ;
    return ;
  }

  /*find the derivative*/
  eval = evaluator_create(expr[i]) ;
  diff = evaluator_derivative(eval, var) ;
  agg_variable_definition(v) = g_strdup(evaluator_get_string(diff)) ;
  agg_variable_value(v) = 0.0 ;
  evaluator_destroy(eval) ;
  evaluator_destroy(diff) ;
  
  return ;
}
#else  /*HAVE_LIBMATHEVAL*/
static void parameter_set_derivative(agg_transform_operator_t *tr,
				     agg_variable_t *v,
				     gchar **expr, gchar **de, gint i,
				     gchar *var)

{
  /*if analytical differentiation is not available, set derivatives to
    zero*/
  agg_variable_definition(v) = NULL ;
  agg_variable_value(v) = 0.0 ;

  return ;
}
#endif /*HAVE_LIBMATHEVAL*/

/** 
 * @{ 
 *
 * @ingroup transforms
 */

/** 
 * Allocate a new ::agg_transform_operator_t
 * 
 * This function is mainly called as part of setting up a new
 * transform, via ::agg_transform_operator_add
 * 
 * @return a newly allocated ::agg_transform_operator_t
 */

agg_transform_operator_t *agg_transform_operator_new(void)

{
  agg_transform_operator_t *op ;
  gint i ;
  
  op = (agg_transform_operator_t *)g_malloc0(sizeof(agg_transform_operator_t)) ;
  for ( i = 0 ; i < AGG_OPERATOR_PARAMETER_SIZE ; i ++ ) {
    agg_transform_operator_parameter(op,i)->eval = NULL ;
    agg_transform_operator_parameter_u(op,i)->eval = NULL ;
    agg_transform_operator_parameter_v(op,i)->eval = NULL ;
  }
  agg_transform_operator_umin(op) =  G_MAXDOUBLE ;
  agg_transform_operator_umax(op) = -G_MAXDOUBLE ;

  return op ;
}

/** 
 * Allocate a new ::agg_transform_t for transformation of geometries
 * 
 * @param nopmax maximum number of operations in transform
 * 
 * @return newly allocated ::agg_transform_t
 */

agg_transform_t *agg_transform_new(gint nopmax)

{
  agg_transform_t *T ;

  T = (agg_transform_t *)g_malloc0(sizeof(agg_transform_t)) ;
  T->op = (agg_transform_operator_t **)
    g_malloc0(nopmax*sizeof(agg_transform_operator_t *)) ;
  agg_transform_operator_number(T) = 0 ;
  agg_transform_operator_number_max(T) = nopmax ;
  agg_transform_variable_number(T) = 0 ;

  T->e = agg_expression_data_new(128) ;
  T->e->ne = 0 ;
  
  return T ;
}

/** 
 * Add a variable to a transform
 * 
 * @param T an ::agg_transform_t allocated with ::agg_transform_new;
 * @param var variable name;
 * @param def variable definition (may be NULL);
 * @param val variable value (used to set variable constant if \a def is NULL).
 * 
 * @return 0 on success
 */

gint agg_transform_variable_add(agg_transform_t *T,
				gchar *var, gchar *def, gdouble val)

{
  agg_variable_t *v ;
  
  if ( agg_transform_variable_number(T) >=
       AGG_TRANSFORM_VARIABLE_NUMBER_MAX ) {
    g_error("%s: cannot add variable %s; variable storage space (%d) exceeded",
	    __FUNCTION__, var, AGG_TRANSFORM_VARIABLE_NUMBER_MAX) ;
  }

  if ( var == NULL ) {
    g_error("%s: NULL variable name", __FUNCTION__) ;
  }
  
  v = agg_transform_variable(T,agg_transform_variable_number(T)) ;
  agg_variable_name(v) = g_strdup(var) ;
  if ( def == NULL ) agg_variable_definition(v) = NULL ;
  else agg_variable_definition(v) = g_strdup(def) ;
  agg_variable_value(v) = val ;

  agg_transform_variable_number(T) ++ ;
  
  return 0 ;
}

/** 
 * Add a list of global variables to the variables of a transform
 * 
 * @param T an ::agg_transform_t allocated with ::agg_transform_new;
 * @param global array of global variables;
 * @param nglobal number of global variables.
 * 
 * @return 0 on success.
 */

gint agg_transform_add_global_variables(agg_transform_t *T,
					agg_variable_t *global, gint nglobal)

{
  agg_expression_data_t *e ;
  gint i ;
  
  if ( global == NULL || nglobal == 0 ) return 0 ;

  e = T->e ;
  for ( i = 0 ; i < nglobal ; i ++ ) {
    agg_expression_data_variable_add(e, &(global[i])) ;
  }
  
  return 0 ;
}

/** 
 * Compile expressions in a transform for future evaluation.
 * 
 * @param T an ::agg_transform_t allocated with ::agg_transform_new;
 * 
 * @return 0 on success.
 */

gint agg_transform_expressions_compile(agg_transform_t *T)

{
  agg_expression_data_t *e ;
  agg_variable_t *v ;
  agg_transform_operator_t *tr ;
  gint i, j ;
  
  e = T->e ;
  
  for ( i = 0 ; i < agg_transform_variable_number(T) ; i ++ ) {
    v = agg_transform_variable(T,i) ;
    agg_expression_data_variable_add(e, v) ;
  }

  for ( i = 0 ; i < agg_transform_variable_number(T) ; i ++ ) {
    v = agg_transform_variable(T,i) ;
    if ( agg_variable_definition(v) == NULL ) {
      agg_variable_evaluator(v) = NULL ;
    } else {
      agg_variable_evaluator(v) =
	agg_expression_compile(agg_variable_definition(v), e) ;       
    }
  }

  for ( i = 0 ; i < agg_transform_operator_number(T) ; i ++ ) {
    tr = agg_transform_operator(T, i) ;
    for ( j = 0 ; j < agg_transform_operator_parameter_number(tr) ; j ++ ) {
      v = agg_transform_operator_parameter(tr, j) ;
      if ( agg_variable_definition(v) == NULL ) {
	agg_variable_evaluator(v) = NULL ;
      } else {	
	agg_variable_evaluator(v) =
	  agg_expression_compile(agg_variable_definition(v), e) ;
      } 

      v = agg_transform_operator_parameter_u(tr, j) ;
      if ( agg_variable_definition(v) == NULL ) {
	agg_variable_evaluator(v) = NULL ;
      } else {	
	agg_variable_evaluator(v) =
	  agg_expression_compile(agg_variable_definition(v), e) ;
      }

      v = agg_transform_operator_parameter_v(tr, j) ;
      if ( agg_variable_definition(v) == NULL ) {
	agg_variable_evaluator(v) = NULL ;
      } else {	
	agg_variable_evaluator(v) =
	  agg_expression_compile(agg_variable_definition(v), e) ;
      }      
    }
  }
  
  return 0 ;
}

/** 
 * Evaluate variables in a transform. Expressions should have been
 * compiled using ::agg_transform_expressions_compile
 * 
 * @param T an ::agg_transform_t allocated with ::agg_transform_new;
 * 
 * @return 0 on success.
 */

gint agg_transform_variables_eval(agg_transform_t *T)

{
  agg_variable_t *v ;
  agg_transform_operator_t *tr ;
  gint i, j ;

  agg_expression_data_eval(T->e) ;
  for ( i = 0 ; i < agg_transform_variable_number(T) ; i ++ ) {
    v = agg_transform_variable(T,i) ;
    if ( agg_variable_definition(v) != NULL ) {
      agg_variable_value(v) = agg_expression_eval(agg_variable_evaluator(v)) ;
    }
  }

  for ( i = 0 ; i < agg_transform_operator_number(T) ; i ++ ) {
    tr = agg_transform_operator(T, i) ;
    for ( j = 0 ; j < agg_transform_operator_parameter_number(tr) ; j ++ ) {
      v = agg_transform_operator_parameter(tr, j) ;
      if ( agg_variable_evaluator(v) != NULL ) {
	agg_variable_value(v) =
	  agg_expression_eval(agg_variable_evaluator(v)) ;
      }
      v = agg_transform_operator_parameter_u(tr, j) ;
      if ( agg_variable_evaluator(v) != NULL ) {
	agg_variable_value(v) =
	  agg_expression_eval(agg_variable_evaluator(v)) ;
      }
      v = agg_transform_operator_parameter_v(tr, j) ;
      if ( agg_variable_evaluator(v) != NULL ) {
	agg_variable_value(v) =
	  agg_expression_eval(agg_variable_evaluator(v)) ;
      }
    }
  }

  return 0 ;
}

/** 
 * Write list of variables and values for a transform to output
 * 
 * @param f output file stream;
 * @param T an ::agg_transform_t allocated with ::agg_transform_new;
 * @param write_defs if TRUE, output variable definitions for variables
 * that have them.
 * 
 * @return 0 on success.
 */

gint agg_transform_variables_write(FILE *f, agg_transform_t *T,
				   gboolean write_defs)

{
  gint i ;
  agg_variable_t *v ;

  if ( !write_defs ) {
    for ( i = 0 ; i < agg_transform_variable_number(T) ; i ++ ) {
      v = agg_transform_variable(T,i) ;
      fprintf(f, "%s = %lg\n", agg_variable_name(v), agg_variable_value(v)) ;
    }
    return 0 ;
  }
  
  for ( i = 0 ; i < agg_transform_variable_number(T) ; i ++ ) {
    v = agg_transform_variable(T,i) ;
    if ( agg_variable_definition(v) == NULL ) 
      fprintf(f, "%s = %lg\n", agg_variable_name(v), agg_variable_value(v)) ;
    else
      fprintf(f, "%s = %s = %lg\n",
	      agg_variable_name(v), agg_variable_definition(v),
	      agg_variable_value(v)) ;
  }
  
  return 0 ;
}

/** 
 * Add an operation to a transform
 * 
 * @param T an ::agg_transform_t allocated with ::agg_transform_new;
 * @param op a basic transform operation of type ::agg_operation_t;
 * @param umin lower parameter limit for applying transform operation;
 * @param umax upper parameter limit for applying transform operation;
 * @param p array of parameter values to pass to transform calculation;
 * @param expr expressions to be used in evaluating transform parameters;
 * @param dedu derivatives of expressions with respect to \f$u\f$;
 * @param dedv derivatives of expressions with respect to \f$v\f$;
 * @param np number of parameters to pass to transform.
 * 
 * @return 0 on success
 */

gint agg_transform_operator_add(agg_transform_t *T, agg_operation_t op,
				gdouble umin, gdouble umax,
				gdouble *p,
				gchar **expr, gchar **dedu, gchar **dedv,
				gint np)

{
  gint i ;
  agg_transform_operator_t *tr ;
  agg_variable_t *v ;
  
  for ( i = 0 ;
	(transform_list[i].name != NULL) &&
	  (transform_list[i].op != op) ; i ++ ) ;

  if ( transform_list[i].op != op ) {
    g_error("%s: unrecognized transform %d", __FUNCTION__, op) ;
  }

  if ( transform_list[i].np != np ) {
    g_error("%s: transform \"%s\" requires %d parameters, %d supplied",
	    __FUNCTION__,
	    transform_list[i].name, transform_list[i].np, np) ;
  }

  tr = agg_transform_operator_new() ;

  agg_transform_operator_operation(tr) = op ;
  agg_transform_operator_parameter_number(tr) = np ;
  agg_transform_operator_func(tr) = transform_list[i].func ;
  agg_transform_operator_umin(tr) = umin ; 
  agg_transform_operator_umax(tr) = umax ; 

  for ( i = 0 ; i < np ; i ++ ) {
    v = agg_transform_operator_parameter(tr,i) ;
    agg_variable_name(v) = NULL ;
    if ( expr[i] != NULL ) {
      agg_variable_definition(v) = g_strdup(expr[i]) ;
    } else {
      agg_variable_definition(v) = NULL ;
    }
    agg_variable_value(v) = p[i] ;

    v = agg_transform_operator_parameter_u(tr,i) ;
    parameter_set_derivative(tr, v, expr, dedu, i, "u") ;
    v = agg_transform_operator_parameter_v(tr,i) ;
    parameter_set_derivative(tr, v, expr, dedv, i, "v") ;
  }

  agg_transform_operator(T,agg_transform_operator_number(T)) = tr ;
  agg_transform_operator_number(T) ++ ;
  
  return 0 ;
}

static void write_variable(FILE *f, agg_variable_t *v)

{
  if ( agg_variable_definition(v) != NULL ) {
    fprintf(f, "%s", agg_variable_definition(v)) ;
    return ;
  }

  fprintf(f, "%lg", agg_variable_value(v)) ;
  
  return ;
}

/** 
 * Write a list of operators and parameters in a transform to file
 * 
 * @param f file stream for output;
 * @param T the ::agg_transform_t to write.
 * 
 * @return 0 on success.
 */

gint agg_transform_operators_write(FILE *f, agg_transform_t *T)

{
  gint i, j ;
  agg_transform_operator_t *tr ;
  agg_variable_t *v ;

  for ( i = 0 ; i < agg_transform_operator_number(T) ; i ++ ) {
    tr = agg_transform_operator(T,i) ;
    j = agg_transform_operator_operation(tr) ;
    fprintf(f, "  %s(", transform_list[j].name) ;
    for ( j = 0 ; j < agg_transform_operator_parameter_number(tr)-1 ; j ++ ) {
      v = agg_transform_operator_parameter(tr,j) ;
      write_variable(f, v) ;
      fprintf(f, ", ") ;
    }
    if ( agg_transform_operator_parameter_number(tr) != 0 ) {
      j = agg_transform_operator_parameter_number(tr) - 1 ;
      v = agg_transform_operator_parameter(tr,j) ;
      write_variable(f, v) ;
    }
    fprintf(f, ")\n") ;
  }
  
  return 0 ;
}

gint agg_transform_operator_rotate(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv)

{
  gdouble x0, y0, th, C, S, xt[2] ;
  
  g_assert(np == 3) ;
  g_assert(op == AGG_TRANSFORM_ROTATE) ;

  x0 = p[0].val ; y0 = p[1].val ; th = p[2].val ;

  C = cos(th) ; S = sin(th) ;
  /*copy input point in case this is a transformation being done
    in-place and xin == xout*/
  xt[0] = xin[0] ; xt[1] = xin[1] ;
  
  xout[0] = x0 + (xt[0] - x0)*C - (xt[1] - y0)*S ; 
  xout[1] = y0 + (xt[0] - x0)*S + (xt[1] - y0)*C ; 
  
  return 0 ;
}

gint agg_transform_operator_shrink(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv)

{
  gdouble x0, y0, sc ;
  
  g_assert(np == 3) ;
  g_assert(op == AGG_TRANSFORM_SHRINK) ;

  x0 = p[0].val ; y0 = p[1].val ; sc = p[2].val ;

  xout[0] = x0 + (xin[0] - x0)*sc ;
  xout[1] = y0 + (xin[1] - y0)*sc ;
  xout[2] = xin[2] ;
  
  return 0 ;
}

gint agg_transform_operator_translate(agg_operation_t op,
				      agg_variable_t *p, gint np,
				      gdouble *xin, gdouble *xout,
				      gdouble *dxdu, gdouble *dxdv)

{
  gdouble dx, dy, dz ;
  
  g_assert(np == 3) ;
  g_assert(op == AGG_TRANSFORM_TRANSLATE) ;

  dx = p[0].val ; dy = p[1].val ; dz = p[2].val ;

  xout[0] = xin[0] + dx ; 
  xout[1] = xin[1] + dy ; 
  xout[2] = xin[2] + dz ; 

  return 0 ;
}

gint agg_transform_operator_scale(agg_operation_t op,
				  agg_variable_t *p, gint np,
				  gdouble *xin, gdouble *xout,
				  gdouble *dxdu, gdouble *dxdv)

{
  gdouble sc ;
  
  g_assert(np == 1) ;
  g_assert(op == AGG_TRANSFORM_SCALE) ;

  sc = p[0].val ;

  xout[0] = xin[0]*sc ;
  xout[1] = xin[1]*sc ;
  xout[2] = xin[2]*sc ;

  return 0 ;
}

gint agg_transform_operator_xscale(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv)

{
  gdouble sc ;
  
  g_assert(np == 1) ;
  g_assert(op == AGG_TRANSFORM_SCALE_X) ;

  sc = p[0].val ;

  xout[0] = xin[0]*sc ;
  xout[1] = xin[1] ;
  xout[2] = xin[2] ;

  return 0 ;
}

gint agg_transform_operator_yscale(agg_operation_t op,
				   agg_variable_t *p, gint np,
				   gdouble *xin, gdouble *xout,
				   gdouble *dxdu, gdouble *dxdv)

{
  gdouble sc ;
  
  g_assert(np == 1) ;
  g_assert(op == AGG_TRANSFORM_SCALE_Y) ;

  sc = p[0].val ;

  xout[0] = xin[0] ;
  xout[1] = xin[1]*sc ;
  xout[2] = xin[2] ;

  return 0 ;
}

/** 
 * Apply a transform operation to a point.
 * 
 * @param op transform operation (usually part of an ::agg_transform_t);
 * @param xin input point;
 * @param xout on exit contains transformed version of \a xin.
 * 
 * @return 0 on success.
 */
  
gint agg_transform_operator_apply(agg_transform_operator_t *op,
				  gdouble *xin, gdouble *xout)

{
  agg_transform_operator_func_t func ;

  func = agg_transform_operator_func(op) ;
  g_assert(func != NULL) ;

  func(agg_transform_operator_operation(op),
       agg_transform_operator_parameters(op),
       agg_transform_operator_parameter_number(op),
       xin, xout, NULL, NULL) ;
  
  return 0 ;
}

static gboolean parameter_in_range(gdouble umin, gdouble umax, gdouble u)

{
  g_assert(umin < umax) ;
  if ( umax == 1 && u == 1 ) return TRUE ;

  if ( umin <= u && u < umax ) return TRUE ;

  return FALSE ;
}


/** 
 * Apply an axis transform (swap) to a point
 * 
 * @param axes an ::agg_axes_t for the required axis transformation;
 * @param xin input point (can be equal to \a xout);
 * @param xout on exit contains transformed point data.
 * 
 * @return 0 on success.
 */

gint agg_transform_axes(agg_axes_t axes, gdouble *xin, gdouble *xout)

{
  gdouble xt[3] = {xin[0], xin[1], xin[2]} ;
  
  switch ( axes ) {
  case AGG_AXES_PX_PY_PZ:
    xout[0] = xt[0] ; xout[1] = xt[1] ; xout[2] = xt[2] ;
    break ;
  case AGG_AXES_PY_PZ_PX:
    xout[0] = xt[1] ; xout[2] = xt[2] ; xout[2] = xt[0] ;
    break ;
  case AGG_AXES_PZ_PX_PY:
    xout[0] = xt[2] ; xout[1] = xt[0] ; xout[2] = xt[1] ;
    break ;
  case AGG_AXES_PZ_PY_PX:
    xout[0] = xt[2] ; xout[1] = xt[1] ; xout[2] = xt[0] ;
    break ;
  case AGG_AXES_PX_PY_MZ:
    xout[0] = xt[0] ; xout[1] = xt[1] ; xout[2] = -xt[2] ;
    break ;
  default: g_assert_not_reached() ;
  }
  
  return 0 ;
}

/** 
 * Apply a transform to a point
 * 
 * @param T an ::agg_transform_t containing the sequence of operations;
 * @param xin input point (can be the same as \a xout);
 * @param xout on exit contains transformed version of \a xin.
 * 
 * @return 0 on success.
 */

gint agg_transform_apply(agg_transform_t *T, gdouble *xin, gdouble *xout)

{
  gint i ;
  gdouble xt[3], u ;
  agg_variable_t *var ;
  agg_transform_operator_t *op ;
  
  /*empty transform: pass input to output*/
  if ( agg_transform_operator_number(T) == 0 ) {
    xout[0] = xin[0] ; xout[1] = xin[1] ; xout[2] = xin[2] ; 
  }

  /*extract parameter from T*/
  var = agg_transform_variable(T, 0) ;
  g_assert( strcmp(agg_variable_name(var), "u") == 0 ) ;
  u = agg_variable_value(var) ;
  
  xt[0] = xin[0] ; xt[1] = xin[1] ; xt[2] = xin[2] ;
  for ( i = 0 ; i < agg_transform_operator_number(T) ; i ++ ) {
    op = agg_transform_operator(T,i) ;
    if ( parameter_in_range(agg_transform_operator_umin(op),
			    agg_transform_operator_umax(op), u) ) {
      agg_transform_operator_apply(op, xt, xout) ;
      xt[0] = xout[0] ; xt[1] = xout[1] ; xt[2] = xout[2] ;
    }
  }

  return 0 ;
}

/** 
 * Parse a string specifying an axis change on a surface.
 * 
 * @param str axis string to parse
 * 
 * @return ::agg_axes_t corresponding to \a str, or AGG_AXES_UNDEFINED. 
 */

agg_axes_t agg_axes_parse(gchar *str)

{
  gint i ;

  for ( i = 0 ; axes_list[i].name != NULL ; i ++ ) {
    if ( strcmp(axes_list[i].name, str) == 0 )
      return axes_list[i].axes ;
  }

  return AGG_AXES_UNDEFINED ;
}

/** 
 * Parse a transform operation described by a string and add it to an
 * ::agg_transform_t
 * 
 * @param T ::agg_transform_t to have an operation added;
 * @param p array of ::agg_variable_t holding parameters of transform, 
 * including transform name;
 * @param np number of entries in \a p.
 * 
 * @return 0 on success.
 */

gint agg_transform_parse(agg_transform_t *T, agg_variable_t *p, gint np)

{
  gchar *expr[32], *name ;
  gdouble args[32], umin, umax ;
  gint i, i0 ;

  umin = 0 ; umax = 1 ;
  if ( agg_variable_definition(&(p[0])) != NULL ) {
    /*first parameter is transform name, use default umin, umax*/
    name = agg_variable_definition(&(p[0])) ;
    i0 = 1 ; np -- ;
  } else {
    /*first two parameters must be numerical*/
    if ( agg_variable_definition(&(p[1])) != NULL ) {
      g_error("%s: both parameter limits must be specified", __FUNCTION__) ;
    }
    umin = agg_variable_value(&(p[0])) ;
    umax = agg_variable_value(&(p[1])) ;
    g_assert( agg_variable_definition(&(p[2])) != NULL ) ;
    name = agg_variable_definition(&(p[2])) ;
    i0 = 3 ; np -= 3 ;
  }
  
  for ( i = 0 ; i < np ; i ++ ) {
    args[i] = agg_variable_value(&(p[i0+i])) ;
    expr[i] = agg_variable_definition((&p[i0+i])) ;
  }

  for ( i = 0 ;	(transform_list[i].name != NULL) ; i ++ ) {
    if ( strcmp(transform_list[i].name, name) == 0 ) {
      agg_transform_operator_add(T, transform_list[i].op, umin, umax,
				 args, expr, NULL, NULL, np) ;
      return 0 ;
    }
  }

  g_error("%s: operation \"%s\" not recognized", __FUNCTION__, name) ;  
  
  return 0 ;
}

/** 
 * Write a list of available transforms to file
 * 
 * @param f file stream for output.
 * 
 * @return 0 on success.
 */

gint agg_transforms_list(FILE *f)

{
  gint i ;

  for ( i = 1 ;	transform_list[i].name != NULL ; i ++ ) {
    fprintf(f, "  %s(%d paramete%s)\n",
	    transform_list[i].name,
	    transform_list[i].np,
	    (transform_list[i].np == 1 ? "r" : "rs")) ;
  }
  
  return 0 ;
}

/**
 * @}
 */
