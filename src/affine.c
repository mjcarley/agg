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

#ifdef HAVE_LIBMATHEVAL
#include <matheval.h>
#endif /*HAVE_LIBMATHEVAL*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

static const struct {
  char *name ;
  agg_affine_func_t func ;
  gint npmin, npmax ;
  char *help ;
} transform_list[] =
  {
    {"undefined",  NULL,                    0,  0, ""},
    {"rotate_x",   agg_affine_rotation_x,   1,  1,
     "rotation about x axis"
    },
    {"rotate_y",   agg_affine_rotation_y,   1,  1,
     "rotation about y axis"
    },
    {"rotate_z",   agg_affine_rotation_z,   1,  1,
     "rotation about z axis"
    },
    {"translate",  agg_affine_translation,  3,  3,
     "translation in three dimensions"
    },
    {"scale",      agg_affine_scale,        1,  3,
     "scale point in three dimensions"
    },
    {NULL,         NULL,                   -1, -1, ""}
  } ;

agg_affine_t *agg_affine_new(gint order_max)

{
  agg_affine_t *A ;

  A = (agg_affine_t *)g_malloc0(sizeof(agg_affine_t)) ;

  A->compiled =            g_malloc0(16*(order_max+1)*sizeof(gpointer)) ;
  A->expr     = (char **)  g_malloc0(16*(order_max+1)*sizeof(char *)) ;
  A->matrix   = (gdouble *)g_malloc0(16*(order_max+1)*sizeof(gdouble)) ;

  A->order = 0 ; 
  A->order_max = order_max ;
  
  return A ;
}

/* static void matvec_mul_4(gdouble *y, gdouble *A, gdouble *x) */

/* { */
/*   gint i, j ; */

/*   for ( i = 0 ; i < 4 ; i ++ ) { */
/*     y[i] = 0.0 ; */
/*     for ( j = 0 ; j < 4 ; j ++ ) { */
/*       y[i] += A[i*4+j]*x[j] ; */
/*     } */
/*   } */
  
/*   return ; */
/* } */

gint agg_affine_point_transform(gdouble *y, gdouble *A, gdouble *x)

{
  /*y := A*x*/
  gdouble xt[4], yt[4] ;

  xt[0] = x[0] ; xt[1] = x[1] ; xt[2] = x[2] ; xt[3] = 1.0 ;
  agg_mat_vec_mul_4(yt, A, xt) ;

  y[0] = yt[0] ; y[1] = yt[1] ; y[2] = yt[2] ;
  
  return 0 ;
}

gint agg_affine_expression_add(agg_affine_t *A, gint order,
			       gint i, gint j, gdouble val, char *expr)

{
  if ( order < 0 || order > agg_affine_order_max(A) ) {
    g_error("%s: order (%d) out of range (0--%d)",
	    __FUNCTION__, order, agg_affine_order_max(A)) ;
  }

  if ( expr == NULL ) {
    agg_affine_expression(A,order,i,j) = NULL ;
    A->matrix[16*order+4*i+j] = val ;
  } else {
    agg_affine_expression(A,order,i,j) = g_strdup(expr) ;
  }
  
  return 0 ;
}

gint agg_affine_expressions_compile(agg_affine_t *A, agg_expression_data_t *e)

{
  gint i, j, k ;

  for ( i = 0 ; i <= agg_affine_order(A) ; i ++ ) {
    for ( j = 0 ; j < 4 ; j ++ ) {
      for ( k = 0 ; k < 4 ; k ++ ) {
	if ( agg_affine_expression(A,i,j,k) != NULL ) {
	  agg_affine_compiled(A,i,j,k) =
	    agg_expression_compile(agg_affine_expression(A,i,j,k), e) ;
	} else {
	  agg_affine_compiled(A,i,j,k) = NULL ;
	}
      }
    }
  }
  
  return 0 ;
}
				    
gint agg_affine_matrices_evaluate(agg_affine_t *A)

{
  gint i, j, k ;
  gdouble *matrix ;
  
  for ( i = 0 ; i <= agg_affine_order(A) ; i ++ ) {
    matrix = agg_affine_matrix(A,i) ;

    for ( j = 0 ; j < 4 ; j ++ ) {
      for ( k = 0 ; k < 4 ; k ++ ) {
	if ( agg_affine_compiled(A,i,j,k) != NULL ) {
	  matrix[4*j+k] = agg_expression_eval(agg_affine_compiled(A,i,j,k)) ;
	}
      }
    }
  }

  return 0 ;
}

gint agg_affine_identity(agg_affine_t *A)

{
  gint i, j ;
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    for ( j = 0 ; j < 4 ; j ++ ) {
      agg_affine_expression_add(A, 0, i, j, 0.0, NULL) ;
    }
    agg_affine_expression_add(A, 0, i, i, 1.0, NULL) ;
  }
  
  return 0 ;
}

gint agg_affine_zero(agg_affine_t *A)

{
  gint i, j ;
  
  for ( i = 0 ; i < 4 ; i ++ ) {
    for ( j = 0 ; j < 4 ; j ++ ) {
      agg_affine_expression_add(A, 0, i, j, 0.0, NULL) ;
    }
  }
  
  return 0 ;
}

gint agg_affine_translation(agg_affine_t *A, gdouble *dx, char **expr, gint ns)
  
{
  agg_affine_identity(A) ;
  
  if ( ns != 3 ) {
    g_error("%s: translation requires three arguments not %d",
	    __FUNCTION__, ns) ;
  }

  agg_affine_expression_add(A, 0, 0, 3, dx[0], expr[0]) ;
  agg_affine_expression_add(A, 0, 1, 3, dx[1], expr[1]) ;
  agg_affine_expression_add(A, 0, 2, 3, dx[2], expr[2]) ;
  
  return 0 ;
}

static void rotation_matrix(agg_affine_t *A, gint *idx, gdouble th, char *expr)

/*
 * idx contains matrix row and column indices for cos, -sin, sin, cos
 * in that order
 */
  
{
  char buf[256] ;

  agg_affine_identity(A) ;
  
  if ( expr == NULL ) {
    agg_affine_expression_add(A, 0, idx[0*2+0], idx[0*2+1],  cos(th), NULL) ;
    agg_affine_expression_add(A, 0, idx[1*2+0], idx[1*2+1], -sin(th), NULL) ;
    agg_affine_expression_add(A, 0, idx[2*2+0], idx[2*2+1],  sin(th), NULL) ;
    agg_affine_expression_add(A, 0, idx[3*2+0], idx[3*2+1],  cos(th), NULL) ;

    return ;
  }

  sprintf(buf, "cos(%s)", expr) ;
  agg_affine_expression_add(A, 0, idx[0*2+0], idx[0*2+1], 0.0, buf) ;
  sprintf(buf, "-sin(%s)", expr) ;
  agg_affine_expression_add(A, 0, idx[1*2+0], idx[1*2+1], 0.0, buf) ;
  sprintf(buf, "sin(%s)", expr) ;
  agg_affine_expression_add(A, 0, idx[2*2+0], idx[2*2+1], 0.0, buf) ;
  sprintf(buf, "cos(%s)", expr) ;
  agg_affine_expression_add(A, 0, idx[3*2+0], idx[3*2+1], 0.0, buf) ;

  return ;
}

gint agg_affine_rotation_x(agg_affine_t *A, gdouble *th, char **expr, gint nth)

{
  gint idx[] = {1, 1, 1, 2, 2, 1, 2, 2} ;

  if ( nth != 1 ) {
    g_error("%s: rotation requires one argument (%d supplied)",
	    __FUNCTION__, nth) ;
  }
  
  rotation_matrix(A, idx, th[0], expr[0]) ;
  
  return 0 ;
}

gint agg_affine_rotation_y(agg_affine_t *A, gdouble *th, char **expr, gint nth)

{
  gint idx[] = {0, 0, 2, 0, 0, 2, 2, 2} ;

  if ( nth != 1 ) {
    g_error("%s: rotation requires one argument (%d supplied)",
	    __FUNCTION__, nth) ;
  }

  rotation_matrix(A, idx, th[0], expr[0]) ;
  
  return 0 ;
}

gint agg_affine_rotation_z(agg_affine_t *A, gdouble *th, char **expr, gint nth)

{
  gint idx[] = {0, 0, 0, 1, 1, 0, 1, 1} ;

  if ( nth != 1 ) {
    g_error("%s: rotation requires one argument (%d supplied)",
	    __FUNCTION__, nth) ;
  }

  rotation_matrix(A, idx, th[0], expr[0]) ;
  
  return 0 ;
}

gint agg_affine_scale(agg_affine_t *A, gdouble *s, char **expr, gint ns) 

{
  agg_affine_identity(A) ;

  if ( ns < 1 || ns > 3 ) {
    g_error("%s: affine scaling operation requires between one and three "
	    "arguments not %d",
	    __FUNCTION__, ns) ;
  }

  if ( ns == 1 ) {
    /*one parameter: scale equally on all axes*/
    agg_affine_expression_add(A, 0, 0, 0, s[0], expr[0]) ;
    agg_affine_expression_add(A, 0, 1, 1, s[0], expr[0]) ;
    agg_affine_expression_add(A, 0, 2, 2, s[0], expr[0]) ;

    return 0 ;
  }

  if ( ns == 2 ) {
    /*two parameters: scale on x and y and leave z unscaled*/
    agg_affine_expression_add(A, 0, 0, 0, s[0], expr[0]) ;
    agg_affine_expression_add(A, 0, 1, 1, s[1], expr[1]) ;

    return 0 ;
  }

  /*general case, all three scaling factors specified*/
  agg_affine_expression_add(A, 0, 0, 0, s[0], expr[0]) ;
  agg_affine_expression_add(A, 0, 1, 1, s[1], expr[1]) ;
  agg_affine_expression_add(A, 0, 2, 2, s[2], expr[2]) ;
  
  return 0 ;
}

gint agg_affine_parse(agg_affine_t *A, agg_variable_t *p, gint np)

{
  agg_affine_func_t func ;
  gint i ;
  char *name, *expr[32] ;
  gdouble val[32] ;
  
  if ( agg_variable_definition(&(p[0])) == NULL ) {
    g_error("%s: transform name NULL", __FUNCTION__) ;
  }

  name = agg_variable_definition(&(p[0])) ;

  func = NULL ;
  for ( i = 0 ; transform_list[i].name != NULL ; i ++ ) {
    if ( strcmp(transform_list[i].name, name) == 0 ) {
      if ( np - 1 < transform_list[i].npmin ||
	   np - 1 > transform_list[i].npmax ) {
	g_error("%s: transform \"%s\" requires between %d and %d arguments, "
		"not %d",
		__FUNCTION__, name,
		transform_list[i].npmin, transform_list[i].npmax,
		np) ;
      }
      func = transform_list[i].func ;
    }
  }
  if ( func == NULL ) {
    g_error("%s: transform \"%s\" not recognised", __FUNCTION__, name) ;
  }

  for ( i = 0 ; i < np - 1 ; i ++ ) {
    expr[i] = agg_variable_definition(&(p[i+1])) ;
    val[i]  = agg_variable_value(&(p[i+1])) ;
  }

  func(A, val, expr, np-1) ;
  
  return 0 ;
}

gint agg_affine_axes(agg_affine_t *A, agg_axes_t axes)

{
  agg_affine_zero(A) ;
  agg_affine_expression_add(A, 0, 3, 3, 1.0, NULL) ;

  switch (axes) {
  default:
  case AGG_AXES_UNDEFINED:
    g_error("%s: axis transformation not defined", __FUNCTION__) ;
    break ;
  case AGG_AXES_PX_PY_PZ:
    agg_affine_expression_add(A, 0, 0, 0, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 1, 1, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 2, 2, +1.0, NULL) ;
    break ; 
  case AGG_AXES_PY_PZ_PX:
    g_assert_not_reached() ;
    agg_affine_expression_add(A, 0, 0, 1, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 2, 0, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 2, 0, +1.0, NULL) ;
    break ; 
  case AGG_AXES_PZ_PX_PY:
    agg_affine_expression_add(A, 0, 0, 2, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 1, 0, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 2, 1, +1.0, NULL) ;
    break ; 
  case AGG_AXES_PZ_PY_PX:
    g_assert_not_reached() ;
    agg_affine_expression_add(A, 0, 0, 2, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 1, 1, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 2, 0, +1.0, NULL) ;
    break ; 
  case AGG_AXES_PX_PY_MZ:
    agg_affine_expression_add(A, 0, 0, 0, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 1, 1, +1.0, NULL) ;
    agg_affine_expression_add(A, 0, 2, 2, -1.0, NULL) ;
    break ; 
  }
  
  return 0 ;
}

static void entry_set_derivative(agg_affine_t *A, gint order,
				 gint i, gint j, char *var)

{
  gpointer eval, diff ;

  agg_affine_matrix_value(A,order+1,i,j) = 0.0 ;
  if ( agg_affine_expression(A,order,i,j) == NULL ) {
    /*constant entry*/
    agg_affine_expression(A,order+1,i,j) = NULL ;
    return ;
  }
  
  /*find the derivative*/
  eval = evaluator_create(agg_affine_expression(A,order,i,j)) ;
  diff = evaluator_derivative(eval, var) ;
  agg_affine_expression(A,order+1,i,j) =
    g_strdup(evaluator_get_string(diff)) ;
  evaluator_destroy(eval) ;
  evaluator_destroy(diff) ;

  return ;
}

gint agg_affine_differentiate(agg_affine_t *A, char *var)

{
  gint order, i, j ;

#ifndef HAVE_LIBMATHEVAL
  g_error("%s: symbolic differentiation not implemented (libmatheval not "
	  "available", __FUNCTION__) ;
#else /*HAVE_LIBMATHEVAL*/
  if ( agg_affine_order(A) + 1 > agg_affine_order_max(A) ) {
    g_error("%s: space allocated for orders up to %d, not enough space for %d",
	    __FUNCTION__, agg_affine_order_max(A), agg_affine_order(A) + 1) ;
  }

  order = agg_affine_order(A) + 1 ;

  for ( i = 0 ; i < 4 ; i ++ ) {
    for ( j = 0 ; j < 4 ; j ++ ) {
      entry_set_derivative(A, order-1, i, j, var) ;
    }
  }
  
    agg_affine_order(A) ++ ;
#endif /*HAVE_LIBMATHEVAL*/
  
  return 0 ;
}

gint agg_affine_write(FILE *f, agg_affine_t *A)

{
  gint order, i, j ;

  for ( order = 0 ; order <= agg_affine_order(A) ; order ++ ) {
    for ( i = 0 ; i < 4 ; i ++ ) {
      for ( j = 0 ; j < 4 ; j ++ ) {
	if ( agg_affine_expression(A,order,i,j) == NULL ) {
	  fprintf(f, " %lg ", agg_affine_matrix_value(A,order,i,j)) ;
	} else {
	  fprintf(f, " (%s) ", agg_affine_expression(A,order,i,j)) ;
	}
      }
      fprintf(stderr, "\n") ;
    }
    fprintf(stderr, "\n") ;
  }
  
  return 0 ;
}

gint agg_affine_list(FILE *f, char *head, char *tail, gboolean help)

{
  gint i ;

  if ( !help ) {
    for ( i = 1 ; transform_list[i].name != NULL ; i ++ ) {
      fprintf(f, "%s%s%s", head, transform_list[i].name, tail) ;
    }

    return 0 ;
  }

  for ( i = 1 ; transform_list[i].name != NULL ; i ++ ) {
    fprintf(f, "%s%s%s", head, transform_list[i].name, tail) ;
    fprintf(f, "%s\n", transform_list[i].help) ;
  }
  
  return 0 ;
}
