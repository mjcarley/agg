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

#include "agg-private.h"

static const struct {
  char *name ;
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

/** 
 * @{ 
 *
 * @ingroup transforms
 */

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

  T->affine = (agg_affine_t **)
    g_malloc0(nopmax*sizeof(agg_affine_t *)) ;
  agg_transform_affine_number(T) = 0 ;
  agg_transform_affine_number_max(T) = nopmax ;

  T->e = agg_expression_data_new(128) ;
  T->e->ne = 0 ;
  
  return T ;
}

/**
 * Add an affine operation to a transform
 *
 * @param T an ::agg_transform_t allocated with ::agg_transform_new;
 * @param A an ::agg_affine_t containing the required operation.
 *
 * @return 0 on success
 */

gint agg_transform_affine_add(agg_transform_t *T, agg_affine_t *A) 

{
  if ( agg_transform_affine_number(T) >=
       agg_transform_affine_number_max(T) ) {
    g_error("%s: not enough space for %d transforms",
	     __FUNCTION__, agg_transform_affine_number(T) + 1) ;
  }

  agg_transform_affine(T,agg_transform_affine_number(T)) = A ;
  agg_transform_affine_number(T) ++ ;
  
  return 0 ;
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
				char *var, char *def, gdouble val)

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
  agg_affine_t *A ;
  gint i ;
  
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

  for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
    A = agg_transform_affine(T, i) ;
    agg_affine_expressions_compile(A, T->e) ;    
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
  agg_affine_t *A ;
  gint i ;

  agg_expression_data_eval(T->e) ;
  for ( i = 0 ; i < agg_transform_variable_number(T) ; i ++ ) {
    v = agg_transform_variable(T,i) ;
    if ( agg_variable_definition(v) != NULL ) {
      agg_variable_value(v) = agg_expression_eval(agg_variable_evaluator(v)) ;
    }
  }

  for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
    A = agg_transform_affine(T,i) ;
    agg_affine_matrices_evaluate(A) ;
  }
  
  agg_matrix_identity_4x4(T->matrix) ;
  for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
    A = agg_transform_affine(T,i) ;
    agg_mat_mat_mul_4(T->matrix, agg_affine_matrix(A,0), T->matrix) ;
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

/* static void write_variable(FILE *f, agg_variable_t *v) */

/* { */
/*   if ( agg_variable_definition(v) != NULL ) { */
/*     fprintf(f, "%s", agg_variable_definition(v)) ; */
/*     return ; */
/*   } */

/*   fprintf(f, "%lg", agg_variable_value(v)) ; */
  
/*   return ; */
/* } */

/* static gboolean parameter_in_range(gdouble umin, gdouble umax, gdouble u) */

/* { */
/*   g_assert(umin < umax) ; */
/*   if ( umax == 1 && u == 1 ) return TRUE ; */

/*   if ( umin <= u && u < umax ) return TRUE ; */

/*   return FALSE ; */
/* } */

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
  /*empty transform: pass input to output*/
  if ( agg_transform_affine_number(T) == 0 ) {
    xout[0] = xin[0] ; xout[1] = xin[1] ; xout[2] = xin[2] ; 
  }

  agg_affine_point_transform(xout, T->matrix, xin) ;
  
  return 0 ;
}

/** 
 * Parse a string specifying an axis change on a surface.
 * 
 * @param str axis string to parse
 * 
 * @return ::agg_axes_t corresponding to \a str, or AGG_AXES_UNDEFINED. 
 */

agg_axes_t agg_axes_parse(char *str)

{
  gint i ;

  for ( i = 0 ; axes_list[i].name != NULL ; i ++ ) {
    if ( strcmp(axes_list[i].name, str) == 0 )
      return axes_list[i].axes ;
  }

  return AGG_AXES_UNDEFINED ;
}

gint agg_transform_evaluate(agg_transform_t *T, gdouble u,
			    gint order, gdouble *matrix)

{
  gint i, j ;
  agg_variable_t *var ;
  agg_affine_t *A ;

  memset(matrix, 0, 16*sizeof(gdouble)) ;
  
  var = agg_transform_variable(T, 0) ;
  if ( strcmp(agg_variable_name(var), "u") != 0 ) {
    g_error("%s: first variable in transform must be \"u\"", __FUNCTION__) ;
  } else {
    agg_variable_value(var) = u ;
  }
  
  agg_transform_variables_eval(T) ;  

  for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
    A = agg_transform_affine(T,i) ;
    if ( agg_affine_order(A) < order ) {
      g_error("%s: affine transform %d does not have enough derivatives",
	      __FUNCTION__, i) ;
    }
    agg_affine_matrices_evaluate(A) ;
  }

  g_assert(order <= 1) ;
  
  if ( order == 0 ) {
    matrix[0] = matrix[5] = matrix[10] = matrix[15] = 1.0 ;
    for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
      A = agg_transform_affine(T,i) ;
      agg_mat_mat_mul_4(matrix, agg_affine_matrix(A,0), matrix) ;
    }

    return 0 ;
  }

  for ( i = 0 ; i < agg_transform_affine_number(T) ; i ++ ) {
    gdouble tmp[] = {
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1} ;
    for ( j = 0 ; j < i ; j ++ ) {
      A = agg_transform_affine(T,j) ;
      agg_mat_mat_mul_4(tmp, agg_affine_matrix(A,0), tmp) ;
    }
    A = agg_transform_affine(T,i) ;
    agg_mat_mat_mul_4(tmp, agg_affine_matrix(A,1), tmp) ;
    for ( j = i + 1 ; j < agg_transform_affine_number(T) ; j ++ ) {
      A = agg_transform_affine(T,j) ;
      agg_mat_mat_mul_4(tmp, agg_affine_matrix(A,0), tmp) ;
    }
    for ( j = 0 ; j < 16 ; j ++ ) matrix[j] += tmp[j] ;
  }
  
  return 0 ;
}

/**
 * @}
 */
