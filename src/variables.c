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

#include "tinyexpr.h"

/**
 * @{
 * @ingroup variables
 */

/** 
 * Allocate a new ::agg_expression_data_t to contain lists of
 * variables and definitions for parsing and evaluation.
 * 
 * @param nemax maximum number of expressions in ::agg_expression_data_t.
 * 
 * @return newly allocated ::agg_expression_data_t.
 */

agg_expression_data_t *agg_expression_data_new(gint nemax)

{
  agg_expression_data_t *d ;

  d = (agg_expression_data_t *)g_malloc0(sizeof(agg_expression_data_t)) ;
  d->nemax = nemax ; d->ne = 0 ;
  d->data = g_malloc0(nemax*sizeof(te_variable)) ;
  d->expr = g_malloc0(nemax*sizeof(gpointer)) ;
  d->defs = (char **)g_malloc0(nemax*sizeof(char *)) ;
  
  return d ;
}

/** 
 * Add a variable to an ::agg_expression_data_t.
 * 
 * @param d ::agg_expression_data_t to have variable added;
 * @param v a variable to be added to \a d; this may depend on variables
 * already in \a d, but not on variables yet to be added.
 * 
 * @return 0 on success.
 */

gint agg_expression_data_variable_add(agg_expression_data_t *d,
				      agg_variable_t *v)

{
  te_variable *vars ;
  gint i ;
  
  if ( d->ne >= d->nemax )
    g_error("%s: not enough space for %d variables", __FUNCTION__, d->ne) ;

  vars = d->data ;
  i = d->ne ;

  vars[i].name = g_strdup(agg_variable_name(v)) ;
  vars[i].address = &(agg_variable_value(v)) ;
  vars[i].type = TE_VARIABLE ;
  vars[i].context = NULL ;

  if ( agg_variable_definition(v) == NULL )
    d->defs[i] = NULL ;
  else
    d->defs[i] = g_strdup(agg_variable_definition(v)) ;
  
  d->ne ++ ;
  
  return 0 ;
}

/** 
 * Compile a symbolically defined expression, for future evaluation
 * 
 * @param e the expression to be compiled;
 * @param d structure holding variable data.
 * 
 * @return pointer to compiled expression.
 */

gpointer agg_expression_compile(char *e, agg_expression_data_t *d)

{
  te_expr *expr ;
  te_variable *data ;
  gint error ;

  data = d->data ;
  expr = te_compile(e, data, d->ne, &error) ;
  if ( expr == NULL ) {
    g_error("\n  %s: error at character %d, parsing \"%s\"",
	    __FUNCTION__, error, e) ;
  }
  
  return expr ;
}

/** 
 * Evaluate a compiled expression
 * 
 * @param e a pointer to a tinyexpr expression
 * 
 * @return the value of \a e.
 */

gdouble agg_expression_eval(gpointer e)

{
  te_expr *expr ;
  gdouble f ;
  
  expr = (te_expr *)e ;
  f = te_eval(expr) ;
  
  return f ;
}

/** 
 * Compile expressions in an expression data structure, for future
 * evaluation
 * 
 * @param d an ::agg_expression_data_t
 * 
 * @return 0 on success.
 */

gint agg_expression_data_compile(agg_expression_data_t *d)

{
  te_variable *vars ;
  gint i, error ;

  vars = d->data ;
  for ( i = 0 ; i < d->ne ; i ++ ) {
    if ( d->defs[i] != NULL ) {
      d->expr[i] = te_compile(d->defs[i], vars, d->ne, &error) ;
      if ( d->expr[i] == NULL ) {
	g_error("\n  %s: error at character %d, parsing \"%s\"",
		__FUNCTION__, error, d->defs[i]) ;
      }
    }
  }
  
  return 0 ;
}

/** 
 * Evaluate expressions in an expression data structure, making the
 * values available for other calculations.
 * 
 * @param d an ::agg_expression_data_t whose expressions are to be evaluated. 
 * 
 * @return 0 on success.
 */

gint agg_expression_data_eval(agg_expression_data_t *d)

{
  te_variable *vars ;
  gint i ;

  vars = d->data ;
  for ( i = 0 ; i < d->ne ; i ++ ) {
    if ( d->expr[i] != NULL ) {
      *((gdouble *)vars[i].address) = te_eval(d->expr[i]) ;
    }
  }
  
  return 0 ;
}

/** 
 * Write a variable definition to file. If the variable is defined
 * using a string expression, this is written, otherwise the constant
 * numerical value is output.
 * 
 * @param f output file stream;
 * @param v variable to write.
 * 
 * @return 0 on success. 
 */

gint agg_variable_write(FILE *f, agg_variable_t *v)

{
  if ( agg_variable_name(v) != NULL )
    fprintf(f, "\"%s\" ", agg_variable_name(v)) ;

  if ( agg_variable_definition(v) != NULL )
    fprintf(f, "\"%s\"", agg_variable_definition(v)) ;

  if ( agg_variable_definition(v) == NULL )
    fprintf(f, "%lg", agg_variable_value(v)) ;
  
  return 0 ;
}

/**
 * @}
 */
