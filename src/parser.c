/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2022 Michael Carley
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include <glib.h>

#include <blaswrap.h>

#include <agg.h>

#include "tinyexpr.h"

static const struct {
  gchar *name ;
  guint symbol ;
} agg_token_list[] =
  {
   {"global",       AGG_SYMBOL_GLOBAL      ,},
   {"distribution", AGG_SYMBOL_DISTRIBUTION,},
   {"transform",    AGG_SYMBOL_TRANSFORM   ,},
   {"shape",        AGG_SYMBOL_SHAPE       ,},
   {"axes",         AGG_SYMBOL_AXES        ,},
   {"grid",         AGG_SYMBOL_GRID        ,},
   {"body",         AGG_SYMBOL_BODY        ,},
   {NULL, 0}
  } ;

static const struct {
  gchar *name ;
  guint symbol ;
} _agg_constant_list[] =
  {
   {"linear", AGG_SPACING_LINEAR},
   {"cosine", AGG_SPACING_COSINE},
   {"halfcos", AGG_SPACING_HALFCOS},
   {"halfsin", AGG_SPACING_HALFSIN},
   {NULL, 0}
  } ;

/**
 * @{ 
 *
 * @ingroup parser
 *
 */

static gboolean expression_null(gchar *expr, gchar *format, gint i)

{
  if ( expr == NULL ) return TRUE ;
  fprintf(stderr, format, i) ;
  return FALSE ;
}

static gboolean double_is_int(gdouble x)

{
  return ( x == round(x) ) ;
}

static gboolean string_is_numeric(gchar *s)

{
  gint i ;

  for ( i = 0 ; s[i] != '\0' ; i ++ ) {
    if ( !g_ascii_isdigit(s[i]) ) return FALSE ;
  }
  
  return TRUE ;
}

/** 
 *
 * Allocate an ::agg_parser_t used for evaluation of parametric
 * expressions
 * 
 * @return a newly allocated ::agg_parser_t
 */

agg_parser_t *agg_parser_alloc(void)

{
  agg_parser_t *p ;

  p = (agg_parser_t *)g_malloc0(sizeof(agg_parser_t)) ;

  p->vars  = g_malloc0(AGG_PARSER_PARAMETER_NUMBER_MAX*sizeof(te_variable)) ;
  p->expr  = g_malloc0(AGG_PARSER_PARAMETER_NUMBER_MAX*sizeof(te_expr *)) ;
  p->nvars = 0 ;
  p->global_set = FALSE ;

  agg_parser_variable_add(p, "s", 0.0) ;
  agg_parser_function_add(p, "tipleft", agg_function_tipleft, 3) ;
  agg_parser_function_add(p, "tipright", agg_function_tipright, 3) ;
  
  return p ;
}

/** 
 * Add a variable declaration to an ::agg_parser_t
 *
 * The expression \a s is a variable assignment which can set the
 * variable to a numeric constant or to an expression which is
 * evaluated as required.
 * 
 * @param p an ::agg_parser_t;
 * @param s an expression of the form "var = ..."
 * 
 * @return 0 on success
 */

gint agg_parser_declaration(agg_parser_t *p, gchar *s)

{
  gchar **tokens, *v, *w ;
  te_variable *vars = p->vars ;
  te_expr **expr = p->expr ;
  gint err ;
  
  g_assert(p->nvars < AGG_PARSER_PARAMETER_NUMBER_MAX) ;

  /*check this is a definition*/
  if ( g_strstr_len(s, -1, "=") != g_strrstr(s, "=") ) {
    g_error("%s: invalid declaration \"%s\"", __FUNCTION__, s) ;
  }
  
  tokens = g_strsplit(s, "=", 0) ;

  if ( strlen(tokens[0]) == 0 ) return 0 ;
  if ( strlen(tokens[1]) == 0 ) {
    g_error("%s: invalid declaration \"%s\"", __FUNCTION__, s) ;
  }    

  v = g_strstrip(tokens[0]) ;
  if ( g_strstr_len(v, -1, " ") != NULL )
    g_error("%s: variable name \"%s\" contains whitespace",
	    __FUNCTION__, v) ;
  w = g_strstrip(tokens[1]) ;

  vars[p->nvars].name    = g_strdup(v) ;
  vars[p->nvars].address = &(p->values[p->nvars]) ;
  vars[p->nvars].type    = TE_VARIABLE ;

  if ( string_is_numeric(w) ) {
    /*numeric declarations are constants*/
    p->values[p->nvars] = atof(w) ;
    p->isexpr[p->nvars] = FALSE ;
  } else {
    /*need to insert variable declarations in expression list*/
    p->isexpr[p->nvars] = TRUE ;
    expr[p->nvars] = te_compile(w, vars, p->nvars, &err) ;
    if ( expr[p->nvars] == NULL ) {
      g_error("%s: error parsing \"%s\" at character %d\n",
	      __FUNCTION__, w, err-1) ;
    }
  }
  
  p->nvars ++ ;

  return 0 ;
}

gint agg_parser_variable_add(agg_parser_t *p, gchar *v, gdouble x)

{
  te_variable *vars = p->vars ;
  gint i ;
  gboolean incn ;
  
  g_assert(p->nvars < AGG_PARSER_PARAMETER_NUMBER_MAX) ;

  incn = TRUE ;
  for ( i = 0 ; i < p->nvars ; i ++ ) {
    if ( strcmp(v, vars[i].name) == 0 ) {
      incn = FALSE ;
      break ;
    }
  }

  vars[i].type = TE_VARIABLE ;
  p->values[i] = x ;
  p->isexpr[i] = FALSE ;

  if ( incn ) {
    vars[i].address = &(p->values[p->nvars]) ;
    vars[i].name = g_strdup(v) ;
    p->nvars ++ ;
  }

  return 0 ;
}

gint agg_parser_expression_add(agg_parser_t *p, gchar *v, gchar *w)

{
  te_variable *vars = p->vars ;
  te_expr **expr = p->expr ;
  gboolean incn ;
  gint err, i ;
  
  g_assert(p->nvars < AGG_PARSER_PARAMETER_NUMBER_MAX) ;
  fprintf(stderr, "%s: %s\n", __FUNCTION__, v) ;

  incn = TRUE ;
  for ( i = 0 ; i < p->nvars ; i ++ ) {
    if ( strcmp(v, vars[i].name) == 0 ) {
      incn = FALSE ;
      break ;
    }
  }
  
  vars[i].type = TE_VARIABLE ;

  if ( incn ) {
    vars[i].name    = g_strdup(v) ;
    vars[i].address = &(p->values[p->nvars]) ;
  }

  if ( string_is_numeric(w) ) {
    /*numeric declarations are constants*/
    p->values[i] = atof(w) ;
    p->isexpr[i] = FALSE ;
  } else {
    /*need to insert variable declarations in expression list*/
    p->isexpr[i] = TRUE ;
    expr[i] = te_compile(w, vars, p->nvars, &err) ;
    if ( expr[i] == NULL ) {
      g_error("%s: error parsing \"%s\" at character %d\n",
	      __FUNCTION__, w, err-1) ;
    }
  }
  
  if ( incn ) p->nvars ++ ;

  return 0 ;
}

gint agg_parser_function_add(agg_parser_t *p, gchar *v, gpointer func,
			     gint na)

{
  te_variable *vars = p->vars ;
  gint i ;
  
  g_assert(p->nvars < AGG_PARSER_PARAMETER_NUMBER_MAX) ;
  fprintf(stderr, "%s: %s\n", __FUNCTION__, v) ;
  
  for ( i = 0 ; i < p->nvars ; i ++ ) {
    if ( strcmp(v, vars[i].name) == 0 ) {
      g_error("%s: duplicate function %s", __FUNCTION__, v) ;
    }
  }

  vars[p->nvars].name    = g_strdup(v) ;
  vars[p->nvars].address = func ;
  switch ( na ) {
  default:
    g_error("%s: cannot add function with %d arguments", __FUNCTION__, na) ;
    break ;
  case 0: vars[p->nvars].type = TE_FUNCTION0 ; break ;
  case 1: vars[p->nvars].type = TE_FUNCTION1 ; break ;
  case 2: vars[p->nvars].type = TE_FUNCTION2 ; break ;
  case 3: vars[p->nvars].type = TE_FUNCTION3 ; break ;
  case 4: vars[p->nvars].type = TE_FUNCTION4 ; break ;
  case 5: vars[p->nvars].type = TE_FUNCTION5 ; break ;
  case 6: vars[p->nvars].type = TE_FUNCTION6 ; break ;
  case 7: vars[p->nvars].type = TE_FUNCTION7 ; break ;
  }

  p->nvars ++ ;

  return 0 ;
}

gint agg_parser_expressions_evaluate(agg_parser_t *p)

/*
 * set the values of any variables which are defined by expressions
 */
  
{
  gint i ;
  te_expr **expr = p->expr ;

  for ( i = 0 ; i < p->nvars ; i ++ ) {
    if ( p->isexpr[i] ) p->values[i] = te_eval(expr[i]) ;
  }
  
  return 0 ;
}

gint agg_parser_declarations_write(FILE *f, agg_parser_t *p)

{
  te_variable *vars = p->vars ;
  gint i ;

  for ( i = 0 ; i < p->nvars ; i ++ ) {
    fprintf(f, "%s = %lg\n", vars[i].name, *((gdouble *)vars[i].address)) ;
  }    
  
  return 0 ;
}

GScanner *agg_scanner_alloc(void)

{
  GScanner *s ;
  gint i ;
  
  s = g_scanner_new (NULL);

  s->config->symbol_2_token = TRUE ;
  /* convert non-floats (octal values, hex values...) to G_TOKEN_INT */
  s->config->numbers_2_int = TRUE ;
  /* convert G_TOKEN_INT to G_TOKEN_FLOAT */
  s->config->int_2_float = TRUE ;
  
  g_scanner_set_scope(s, 0) ;

  for ( i = 0 ; agg_token_list[i].name != NULL ; i ++ ) {
    g_scanner_scope_add_symbol(s, 0,
			       agg_token_list[i].name,
			       GINT_TO_POINTER(agg_token_list[i].symbol)) ;
  }

  g_assert(agg_token_list[i-1].symbol == AGG_SYMBOL_MAX) ;

  return s ;
}

static guint next_block(GScanner *scanner, agg_parser_t *p)

{
  guint symbol;

  /* expect a valid symbol */
  g_scanner_get_next_token(scanner) ;
  symbol = scanner->token ;
  if ( symbol >= AGG_SYMBOL_GLOBAL &&
       symbol <= AGG_SYMBOL_MAX ) return symbol ;

  /*content outside a block*/
  return G_TOKEN_ERROR ;
}

static guint scanner_read(GScanner *s)

{
  g_scanner_get_next_token(s) ;
  return s->token ;  
}

static guint scan_arguments(GScanner *scanner, agg_parser_t *p,
			    gchar *expr[], gdouble param[],
			    gint *nparam)

/*
 * read an argument list delimited by parentheses and return as
 * strings (for expressions or keywords) or doubles (ints are encoded
 * as doubles)
 */
  
{
  guint symbol ;
  gboolean negate ;
  
  /*look for opening bracket*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_PAREN ) return G_TOKEN_ERROR ;
  symbol = scanner_read(scanner) ;

  *nparam = 0 ;
  while ( symbol != G_TOKEN_RIGHT_PAREN ) {
    negate = FALSE ;

    if ( scanner->token == '-' ) {
      negate = !negate ;
      symbol = scanner_read(scanner) ;
    }

    if ( symbol == G_TOKEN_FLOAT ) {
      if ( negate )
  	param[(*nparam)] = -scanner->value.v_float ;
      else
  	param[(*nparam)] = scanner->value.v_float ;
      expr[(*nparam)] = NULL ;
      (*nparam) ++ ;
    }

    if ( symbol == G_TOKEN_STRING ) {
      expr[(*nparam)] = g_strdup(scanner->value.v_string) ;
      (*nparam) ++ ;
    }
    
    symbol = scanner_read(scanner) ;
    if ( symbol == G_TOKEN_COMMA ) 
      symbol = scanner_read(scanner) ;
  }

  return G_TOKEN_NONE ;
}

static guint scan_global(GScanner *scanner, agg_parser_t *p)

{
  guint symbol;
  gchar v[64] ;
  gboolean negate ;
  
  fprintf(stderr, "scanning global block\n") ;

  if ( agg_parser_global_set(p) ) return G_TOKEN_ERROR ;
  
  if ( scanner_read(scanner) != G_TOKEN_LEFT_CURLY ) return G_TOKEN_ERROR ;

  /* fprintf(stderr, "%u\n", G_TOKEN_CHAR) ; */
  
  do {
    negate = FALSE ;
    symbol = scanner_read(scanner) ;
    if ( isalpha(symbol) ) {
      sprintf(v, "%c", symbol) ;
    } else {
      switch ( symbol ) {
      case G_TOKEN_IDENTIFIER:
	sprintf(v, "%s", scanner->value.v_identifier) ;
	break ;
      case G_TOKEN_RIGHT_CURLY:
	return G_TOKEN_NONE ;
	break ;
      default:
	return G_TOKEN_ERROR ;
	break ;      
      }
    }

    if ( scanner_read(scanner) != G_TOKEN_EQUAL_SIGN ) return G_TOKEN_ERROR ;

    g_scanner_peek_next_token(scanner) ;
    if (scanner->next_token == '-') {
      g_scanner_get_next_token(scanner) ;
      negate = !negate;
    }

    symbol = scanner_read(scanner) ;
    switch ( symbol ) {
    default: g_assert_not_reached() ; break ;
    case G_TOKEN_ERROR: return G_TOKEN_ERROR ; break ;
    case G_TOKEN_FLOAT:
      {
	gdouble w = scanner->value.v_float ;
	if ( negate ) w = -w ;
	/*add a constant declaration to the parser*/
	agg_parser_variable_add(p, v, w) ;
      }
      break ;
    case G_TOKEN_STRING:
      fprintf(stderr, "string: %s\n", scanner->value.v_string) ;
      /*add an expression declaration to the parser*/
      agg_parser_expression_add(p, v, scanner->value.v_string) ;
      break ;
    }
    g_scanner_peek_next_token(scanner) ;
  } while ( symbol != G_TOKEN_ERROR &&
	    symbol != G_TOKEN_EOF &&
	    scanner->next_token != G_TOKEN_RIGHT_CURLY ) ;

  if ( symbol == G_TOKEN_ERROR || symbol == G_TOKEN_EOF ) return symbol ;
  
  if ( scanner->next_token == G_TOKEN_RIGHT_CURLY ) 
    g_scanner_get_next_token(scanner) ;

  agg_parser_global_set(p) = TRUE ;
  
  return G_TOKEN_NONE ;
}

static guint parse_transform(GScanner *scanner,
			     agg_parser_t *p,
			     agg_local_transform_t *t)

{
  guint symbol ;
  gchar tname[32], *expr[32] ;
  gdouble param[32] = {0} ;
  gint nparam, i, j ;
  
  symbol = scan_arguments(scanner, p, expr, param, &nparam) ;
  if ( symbol != G_TOKEN_NONE ) return G_TOKEN_ERROR ;
  
  if ( expr[0] == NULL ) {
    fprintf(stderr,
	    "first argument to transform must be transform name\n") ;
    return G_TOKEN_ERROR ;
  }
  sprintf(tname, "%s", expr[0]) ;

  /*add transform to t*/
  i = t->nt ;
  if ( agg_local_transform_parse(t, tname, &(param[1]), nparam-1) != 0 ) {
    g_error("%s: unrecognized transform %s", __FUNCTION__, tname) ;
  }
  for ( j = 0 ; j < nparam-1 ; j ++ ) {
    if ( expr[j+1] != NULL ) {
      gint idx ;
      idx = agg_local_transform_parameter_index(t, i, j) ;
      agg_local_transform_set_expression(t, idx, expr[j+1], p) ;
    }
  }
  
  return G_TOKEN_NONE ;
}

static guint parse_shape(GScanner *scanner,
			 agg_parser_t *p,
			 agg_function_call_t *f)

{
  guint symbol ;
  gchar *expr[32] ;
  gdouble param[32] = {0} ;
  gint i, nparam, err ;
  
  f->na = 0 ;
  /*look for opening bracket*/
  symbol = scan_arguments(scanner, p, expr, param, &nparam) ;
  if ( symbol != G_TOKEN_NONE ) return G_TOKEN_ERROR ;
  if ( expr[0] == NULL ) {
    fprintf(stderr,
	    "first argument to shape must be shape name\n") ;
    return G_TOKEN_ERROR ;
  }

  sprintf(f->func, "%s", expr[0]) ;
  fprintf(stderr, "  %s\n", f->func) ;
  
  for ( i = 1 ; i < nparam ; i ++ ) {
    if ( expr[i] == NULL ) {
      f->x[f->na] = param[i] ;
      f->expr[f->na] = NULL ;
      f->na ++ ;      
    } else {
      f->expr[f->na] = te_compile(expr[i], p->vars, p->nvars, &err) ; ;
      f->na ++ ;      
    }
  }

  return G_TOKEN_NONE ;
}

static guint parse_axes(GScanner *scanner,
			agg_parser_t *p,
			agg_distribution_t *d)

{
  /*look for opening bracket*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_PAREN ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ;

  if ( agg_distribution_axes_parse(scanner->value.v_identifier,
				   d->axes) != 0 )
    return G_TOKEN_ERROR ;
  
  /*look for closing bracket*/
  if ( scanner_read(scanner) != G_TOKEN_RIGHT_PAREN ) return G_TOKEN_ERROR ;
  
  return G_TOKEN_NONE ;
}

static guint parse_distribution_string(GScanner *scanner,
				       agg_parser_t *p,
				       agg_distribution_t *d)

{
  /*read string*/
  if ( strcmp("invert", scanner->value.v_identifier) == 0 ) {
    agg_distribution_invert(d) = TRUE ;
    return G_TOKEN_NONE ;
  }

  return G_TOKEN_ERROR ;
}

static guint scan_distribution(GScanner *scanner, agg_parser_t *p,
			       agg_body_t *b)

{
  guint symbol;
  agg_distribution_t *d ;
  gchar v[128], *expr[32] ;
  gdouble smin, smax, s, param[32] = {0} ;
  gint ns, ntr, i, nparam ;
  agg_function_call_t shape ;
  
  fprintf(stderr, "scanning distribution: ") ;

  /*initialize transformation list*/
  ntr = 0 ;
  shape.func[0] = '\0' ; shape.c = AGG_SPACING_COSINE ;
  
  /*need the global data set first*/
  if ( !agg_parser_global_set(p) ) return G_TOKEN_ERROR ;
  
  symbol = scan_arguments(scanner, p, expr, param, &nparam) ;
  if ( symbol != G_TOKEN_NONE ) return G_TOKEN_ERROR ;
  
  /*distribution name*/
  if ( expr[0] == NULL ) return G_TOKEN_ERROR ;
  sprintf(v, "%s", expr[0]) ;
  fprintf(stderr, "%s\n", v) ;

  if ( !expression_null(expr[1], "parameter %d must be constant", 1) ) {
    return G_TOKEN_ERROR ;
  }
  smin = param[1] ;

  if ( !expression_null(expr[2], "parameter %d must be constant", 2) ) {
    return G_TOKEN_ERROR ;
  }
  smax = param[2] ;

  /*number of steps in parameter*/  
  if ( !expression_null(expr[3], "parameter %d must be constant", 3) ) {
    return G_TOKEN_ERROR ;
  }
  ns = (gint)param[3] ;
  
  d = agg_distribution_alloc(ns) ;
  agg_body_distribution_add(b, d, v) ;
 
  /*step spacing in generating shapes*/
  i = agg_parser_constant_parse(expr[4], &(d->sg)) ;
  if ( i != 0 )
    g_error("%s: unrecognized constant \"%s\"",
	    __FUNCTION__, scanner->value.v_identifier) ;

  /*curly brackets for the block information*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_CURLY ) return G_TOKEN_ERROR ;

  symbol = G_TOKEN_NONE ;
  do {
    g_scanner_peek_next_token(scanner) ;
    if ( scanner->next_token == G_TOKEN_RIGHT_CURLY ) break ;

    symbol = scanner_read(scanner) ;
    if ( symbol == AGG_SYMBOL_TRANSFORM ) {
      if ( parse_transform(scanner, p, d->t) != G_TOKEN_NONE )
	return G_TOKEN_ERROR ;
      ntr ++ ;
    }

    if ( symbol == AGG_SYMBOL_SHAPE ) {
      if ( parse_shape(scanner, p, &shape) != G_TOKEN_NONE )
	return G_TOKEN_ERROR ;
    }

    if ( symbol == AGG_SYMBOL_AXES ) {
      if ( parse_axes(scanner, p, d) != G_TOKEN_NONE )
	return G_TOKEN_ERROR ;
    }

    /*check for miscellaneous single keywords*/
    if ( symbol == G_TOKEN_STRING || symbol == G_TOKEN_IDENTIFIER ||
	 symbol == G_TOKEN_SYMBOL ) {
      if ( parse_distribution_string(scanner, p, d) != G_TOKEN_NONE )
	return G_TOKEN_ERROR ;
    }
  } while ( symbol != G_TOKEN_ERROR &&
	    symbol != G_TOKEN_EOF &&
	    scanner->next_token != G_TOKEN_RIGHT_CURLY ) ;

  if ( symbol == G_TOKEN_ERROR || symbol == G_TOKEN_EOF ) return symbol ;
  
  if ( scanner->next_token == G_TOKEN_RIGHT_CURLY ) 
    g_scanner_get_next_token(scanner) ;

  /*generate the distribution*/
  d->smin = smin ; d->smax = smax ;

  for ( i = 0 ; i < ns ; i ++ ) {
    /*data for this station*/
    s = agg_spacing_eval(d->smin, d->smax, ns, d->sg, i) ;
    p->values[AGG_PARSER_PARAMETER_RESERVED_S] = s ;
    agg_parser_expressions_evaluate(p) ;

    agg_shape_t *sh = agg_shape_alloc(32) ;
    agg_functions_eval(&(shape), p) ;
    if ( agg_shape_parse(sh, shape.func, shape.x, shape.na) != 0 ) {
      g_error("%s: unrecognized shape %s", __FUNCTION__, shape.func) ;
    }
    agg_distribution_add_shape(d, s, sh) ;
  }

  return G_TOKEN_NONE ;
}

static guint parse_grid(GScanner *scanner,
			agg_parser_t *p,
			agg_grid_t *g)
{
  guint symbol ;
  gchar gname[32], *expr[32] ;
  gdouble param[32] = {0} ;
  gint nparam ;
  
  sprintf(gname, "%s", scanner->value.v_identifier) ;
  fprintf(stderr, "grid: %s\n", gname) ;

  symbol = scan_arguments(scanner, p, expr, param, &nparam) ;
  if ( symbol != G_TOKEN_NONE ) return G_TOKEN_ERROR ;

  if ( agg_grid_parse(g, gname, expr, param, nparam) != 0 )
    return G_TOKEN_ERROR ;
  
  return G_TOKEN_NONE ;
}

static gint scan_grid(GScanner *scanner,
		      agg_parser_t *p, agg_body_t *b)

{
  guint symbol;
  agg_grid_t *g ;
  gint np, nt, nparam ;
  gchar *expr[32] ;
  gdouble param[32] ;

  fprintf(stderr, "scanning grid\n") ;

  symbol = scan_arguments(scanner, p, expr, param, &nparam) ;
  if ( symbol != G_TOKEN_NONE ) return G_TOKEN_ERROR ;

  if ( nparam != 2 ) {
    fprintf(stderr, "grid takes exactly two arguments\n") ;
    return G_TOKEN_ERROR ;
  }
  if ( expr[0] != NULL || expr[1] != NULL ) {
    fprintf(stderr, "grid arguments must be constants\n") ;

    return G_TOKEN_ERROR ;    
  }

  if ( !double_is_int(param[0]) || !double_is_int(param[1]) ) {
    fprintf(stderr, "grid arguments must be integers\n") ;

    return G_TOKEN_ERROR ;    
  }
  
  np = (gint)param[0] ; nt = (gint)param[1] ; 
  fprintf(stderr, "%d points, %d triangles\n", np, nt) ;

  /*set up the body grid*/
  if ( agg_body_grid(b) == NULL ) {
    agg_body_grid(b) = agg_grid_alloc(np, nt) ;
  }
  g = agg_body_grid(b) ;
  if ( np > agg_grid_point_number_max(g) ) {
    g_error("%s: not enough points in grid (%d available, %d requested)",
	    __FUNCTION__, agg_grid_point_number_max(g), np) ;
  }
  if ( nt > agg_grid_triangle_number_max(g) ) {
    g_error("%s: not enough triangles in grid (%d available, %d requested)",
	    __FUNCTION__, agg_grid_triangle_number_max(g), np) ;
  }
  agg_grid_init(g) ;
  
  /*back to default behaviour*/
  scanner->config->int_2_float = TRUE ;

  /*curly brackets for the block information*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_CURLY ) return G_TOKEN_ERROR ;

  symbol = G_TOKEN_NONE ;
  do {
    g_scanner_peek_next_token(scanner) ;
    if ( scanner->next_token == G_TOKEN_RIGHT_CURLY ) break ;

    symbol = scanner_read(scanner) ;

    if ( symbol == G_TOKEN_STRING ||
	 symbol == G_TOKEN_IDENTIFIER )
      symbol = parse_grid(scanner, p, g) ;
  } while ( symbol != G_TOKEN_ERROR &&
	    symbol != G_TOKEN_EOF &&
	    scanner->next_token != G_TOKEN_RIGHT_CURLY ) ;

  if ( symbol == G_TOKEN_ERROR || symbol == G_TOKEN_EOF ) return symbol ;
    
  if ( scanner->next_token == G_TOKEN_RIGHT_CURLY ) 
    g_scanner_get_next_token(scanner) ;
  
  return G_TOKEN_NONE ;
}
		      
/** 
 * Read a body (a collection of distributions) from file
 * 
 * @param scanner lexical scanner;
 * @param p an ::agg_parser_t for evaluation of parametric variables;
 * @param b on exit contains the body information.
 * 
 * @return 0 on success.
 */

gint agg_parser_body_read(GScanner *scanner, agg_parser_t *p, agg_body_t *b)

{
  guint token ;
  
  do {
    token = next_block(scanner, p);
    switch ( token ) {
    case AGG_SYMBOL_GLOBAL:
      token = scan_global(scanner, p) ;
      break ;
    case AGG_SYMBOL_DISTRIBUTION:
      token = scan_distribution(scanner, p, b) ;
      break ;
    case AGG_SYMBOL_GRID:
      token = scan_grid(scanner, p, b) ;
      break ;
    default:
      token = G_TOKEN_ERROR ;
      break ;
    }
  
    g_scanner_peek_next_token (scanner);
  } while ( token == G_TOKEN_NONE &&
	    scanner->next_token != G_TOKEN_EOF &&
	    scanner->next_token != G_TOKEN_RIGHT_CURLY &&
	    scanner->next_token != G_TOKEN_ERROR ) ;
  
  if ( token == G_TOKEN_ERROR ) {
    fprintf(stderr, "parse error at line %d, character %d\n",
	    g_scanner_cur_line(scanner),
	    g_scanner_cur_position(scanner)) ;
    return 1 ;
  }
  
  agg_body_parser(b) = p ;
  
  return 0 ;
}

gint agg_functions_eval(agg_function_call_t *func, agg_parser_t *p)

{
  gint i ;
  te_expr **expr = (te_expr **)(func->expr) ;

  for ( i = 0 ; i < func->na ; i ++ ) {
    if ( expr[i] != NULL ) func->x[i] = te_eval(expr[i]) ;
  }
  
  return 0 ;
}

gint agg_parser_constant_parse(gchar *s, guint *v)

{
  gint i ;

  for ( i = 0 ; _agg_constant_list[i].name != NULL ; i ++ ) {
    if ( !strcmp(s, _agg_constant_list[i].name) ) {
      *v = _agg_constant_list[i].symbol ;
      return 0 ;
    }
  }
  
  return -1 ;
}

/** 
 * Set the numerical value of parameters defined symbolically in a transform
 * 
 * @param t transform whose parameters are to be evaluated.
 * 
 * @return 0 on success
 */

gint agg_local_transform_eval_parameters(agg_local_transform_t *t)

{
  gint i ;
  te_expr *expr ;

  for ( i = 0 ; i < t->p1[t->nt] ; i ++ ) {
    if ( t->isexpr[i] ) {
      expr = t->expr[i].expr ;
      t->p[i] = te_eval(expr) ;
    }
  }
  
  return 0 ;
}

/** 
 * Set the numerical value of parameters in a transform
 * 
 * This function copies \a t into \a T, evaluating expressions as it
 * goes, allowing \a t to be accessed simultaneously in a thread-safe
 * manner.
 *
 * @param T transform whose parameters are to be set;
 * @param t transform containing expressions for parameters.
 * 
 * @return 0 on success
 */

gint agg_local_transform_set_parameters(agg_local_transform_t *T,
					agg_local_transform_t *t)

{
  gint i ;
  te_expr *expr ;

  for ( i = 0 ; i < t->p1[t->nt] ; i ++ ) {
    if ( t->isexpr[i] ) {
      expr = t->expr[i].expr ;
      T->p[i] = te_eval(expr) ;
    } else {
      T->p[i] = t->p[i] ;
    }
  }

  T->nt = t->nt ;
  memcpy(T->p1, t->p1, (AGG_TRANSFORM_FUNCTION_NUMBER+1)*sizeof(gint)) ;
  memcpy(T->func, t->func,
	 AGG_TRANSFORM_FUNCTION_NUMBER*sizeof(agg_local_transform_func_t)) ;

  return 0 ;
}

/** 
 * Insert a symbolic expression for the evaluation of transform parameters
 * 
 * @param t ::agg_local_transform_t whose parameter is to be set;
 * @param i index of parameter in \a t;
 * @param def symbolic expression to add;
 * @param p ::agg_parser_t containing symbol definitions.
 * 
 * @return 0 on success.
 */

gint agg_local_transform_set_expression(agg_local_transform_t *t,
					gint i, gchar *def,
					agg_parser_t *p)

{
  agg_expression_t *e ;
  te_variable *vars = p->vars ;
  te_expr *expr ;
  gint err ;

  e = &(t->expr[i]) ;

  if ( strlen(def) > AGG_EXPRESSION_LENGTH_MAX ) {
    g_error("%s: expression ""%s"" too long (limit is %d characters)",
	    __FUNCTION__, def, AGG_EXPRESSION_LENGTH_MAX) ;
  }

  strcpy(e->def, def) ;

  if ( string_is_numeric(def) ) {
    /*numeric declarations are constants*/
    t->p[i] = atof(def) ;
    t->isexpr[i] = FALSE ;
  } else {
    /*need to insert variable declarations in expression list*/
    expr = te_compile(def, vars, p->nvars, &err) ;
    if ( expr == NULL ) {
      g_error("%s: error parsing \"%s\" at character %d\n",
  	      __FUNCTION__, def, err-1) ;
    }
    t->isexpr[i] = TRUE ;
    e->expr = expr ;
  }
  
  return 0 ;
}

static guint scan_body(GScanner *scanner, agg_parser_t *p, agg_body_t *b)

{
  guint symbol;
  gchar v[128], *expr[32] ;
  gdouble param[32] = {0} ;
  gint nparam ;
  
  fprintf(stderr, "scanning body: ") ;

  /*need the global data set first*/
  if ( !agg_parser_global_set(p) ) return G_TOKEN_ERROR ;
  
  symbol = scan_arguments(scanner, p, expr, param, &nparam) ;
  if ( symbol != G_TOKEN_NONE ) return G_TOKEN_ERROR ;

  /*distribution name*/
  if ( expr[0] == NULL ) return G_TOKEN_ERROR ;
  sprintf(v, "%s", expr[0]) ;
  fprintf(stderr, "%s\n", v) ;

  /*curly brackets for the block information*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_CURLY ) return G_TOKEN_ERROR ;

  /* do { */
    agg_parser_body_read(scanner, p, b) ;

    g_scanner_peek_next_token(scanner) ;
  /* } while ( symbol != G_TOKEN_ERROR && */
  /* 	    symbol != G_TOKEN_EOF && */
  /* 	    scanner->next_token != G_TOKEN_RIGHT_CURLY ) ; */

  if ( symbol == G_TOKEN_ERROR || symbol == G_TOKEN_EOF ) return symbol ;
  
  if ( scanner->next_token == G_TOKEN_RIGHT_CURLY ) 
    g_scanner_get_next_token(scanner) ;

  return G_TOKEN_NONE ;
}

gint agg_parser_crowd_read(GScanner *scanner, agg_crowd_t *c)

{
  guint token ;
  gint body, i ;
  agg_parser_t *p ;
  agg_body_t *b ;
  
  body = -1 ;

  for ( i = 0 ; i < AGG_CROWD_BODY_NUMBER_MAX ; i ++ ) c->b[i] = NULL ;
  c->nb = 0 ;
  
  p = c->p ;
  do {
    token = next_block(scanner, p);
    switch ( token ) {
    case AGG_SYMBOL_GLOBAL:
      token = scan_global(scanner, p) ;
      break ;
    case AGG_SYMBOL_BODY:
      if ( body != -1 ) {
	fprintf(stderr,
		"%s: body cannot be defined inside another body "
		"line %d\n", __FUNCTION__, g_scanner_cur_line(scanner)) ;
	exit(1) ;
      }
      body = c->nb ;
      c->b[body] = b = agg_body_alloc() ;
      token = scan_body(scanner, p, b) ;
      c->nb ++ ;
      body = -1 ;
      break ;
    case AGG_SYMBOL_DISTRIBUTION:
      if ( body == -1 ) {
	fprintf(stderr,
		"line %d: cannot add distribution without defining body\n",
		g_scanner_cur_line(scanner)) ;
	exit(1) ;
      }
      break ;
    case AGG_SYMBOL_GRID:
      if ( body == -1 ) {
	fprintf(stderr, "line %d: cannot add grid without defining body\n",
		g_scanner_cur_line(scanner)) ;
	exit(1) ;
      }
      /* token = scan_grid(scanner, p, b) ; */
      break ;
    default:
      token = G_TOKEN_ERROR ;
      break ;
    }
  
    g_scanner_peek_next_token (scanner);
  } while ( token == G_TOKEN_NONE &&
	    scanner->next_token != G_TOKEN_EOF &&
	    scanner->next_token != G_TOKEN_ERROR ) ;
  
  if ( token == G_TOKEN_ERROR ) {
    fprintf(stderr, "parse error at line %d, character %d\n",
	    g_scanner_cur_line(scanner),
	    g_scanner_cur_position(scanner)) ;
    return 1 ;
  }

  return 0 ;
}


/* @} */
