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

static gboolean string_is_numeric(gchar *s)

{
  gint i ;

  for ( i = 0 ; s[i] != '\0' ; i ++ ) {
    if ( !g_ascii_isdigit(s[i]) ) return FALSE ;
  }
  
  return TRUE ;
}

agg_parser_t *agg_parser_alloc(void)

{
  agg_parser_t *p ;

  p = (agg_parser_t *)g_malloc0(sizeof(agg_parser_t)) ;

  p->vars  = g_malloc0(AGG_PARSER_PARAMETER_NUMBER_MAX*sizeof(te_variable)) ;
  p->expr  = g_malloc0(AGG_PARSER_PARAMETER_NUMBER_MAX*sizeof(te_expr *)) ;
  p->nvars = 0 ;
  p->global_set = FALSE ;

  agg_parser_variable_add(p, "t", 0.0) ;
  agg_parser_function_add(p, "tipleft", agg_function_tipleft, 3) ;
  agg_parser_function_add(p, "tipright", agg_function_tipright, 3) ;
  
  return p ;
}

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

  vars[i].type    = TE_VARIABLE ;
  p->values[i]    = x ;
  p->isexpr[i]    = FALSE ;

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
  
  vars[i].type    = TE_VARIABLE ;

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

gint agg_parser_read_header(FILE *f, agg_parser_t *p)

{
  gchar line[1024], *s ;
  gint i ;
  
  while ( (i = fscanf(f, "%[^\n]s", line)) != EOF ) {
    s = g_strstrip(line) ;
    fprintf(stderr, "%s\n", s) ;
    if ( s[0] != '#' ) {
      agg_parser_declaration(p, s) ;
    }
    fscanf(f, "%*c") ;
  }

  fprintf(stderr, "header ends\n") ;
  
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
  gboolean negate ;
  gchar tname[32], *expr[32] ;
  gdouble param[32] = {0} ;
  gint nparam, i, j ;
  
  /*look for opening bracket*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_PAREN ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ;
  /*name of the transform*/
  sprintf(tname, "%s", scanner->value.v_identifier) ;

  symbol = scanner_read(scanner) ;

  nparam = 0 ;
  while ( symbol != G_TOKEN_RIGHT_PAREN ) {
    if ( symbol != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
    negate = FALSE ;
    
    /*check for negation here*/
    g_scanner_peek_next_token(scanner) ;
    if (scanner->next_token == '-') {
      g_scanner_get_next_token(scanner) ;
      negate = !negate;
    }

    symbol = scanner_read(scanner) ;
    if ( symbol == G_TOKEN_FLOAT ) {
      if ( negate ) 
	param[nparam] = -scanner->value.v_float ;
      else
	param[nparam] = scanner->value.v_float ;
      expr[nparam] = NULL ;
      nparam ++ ;
    }

    if ( symbol == G_TOKEN_STRING ) {
      expr[nparam] = g_strdup(scanner->value.v_string) ;
      nparam ++ ;
    }
    
    symbol = scanner_read(scanner) ;    
  }

  /*add transform to t*/
  i = t->nt ;
  if ( agg_local_transform_parse(t, tname, param, nparam) != 0 ) {
    g_error("%s: unrecognized transform %s", __FUNCTION__, tname) ;
  }
  for ( j = 0 ; j < nparam ; j ++ ) {
    if ( expr[j] != NULL ) {
      gint idx ;
      idx = agg_local_transform_parameter_index(t, i,j) ;
      agg_local_transform_set_expression(t, idx, expr[j], p) ;
    }
  }
  
  return G_TOKEN_NONE ;
}

static guint parse_shape(GScanner *scanner,
			 agg_parser_t *p,
			 agg_function_call_t *f)

{
  guint symbol ;
  gboolean negate ;
  
  f->na = 0 ;
  /*look for opening bracket*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_PAREN ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ;
  sprintf(f->func, "%s", scanner->value.v_identifier) ;

  fprintf(stderr, "  %s\n", f->func) ;
  
  symbol = scanner_read(scanner) ;

  while ( symbol != G_TOKEN_RIGHT_PAREN ) {
    if ( symbol != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
    negate = FALSE ;
    
    /*check for negation here*/
    g_scanner_peek_next_token(scanner) ;
    if (scanner->next_token == '-') {
      g_scanner_get_next_token(scanner) ;
      negate = !negate;
    }

    symbol = scanner_read(scanner) ;
    if ( symbol == G_TOKEN_FLOAT ) {
      if ( negate ) 
	f->   x[f->na] = -scanner->value.v_float ;
      else
	f->   x[f->na] =  scanner->value.v_float ;
      f->expr[f->na] = NULL ;
      f->na ++ ;      
    }
    
    if ( symbol == G_TOKEN_STRING ||
	 symbol == G_TOKEN_IDENTIFIER ||
	 symbol == G_TOKEN_CHAR ) {
      gint err ;
      f->expr[f->na] = te_compile(scanner->value.v_string,
				  p->vars, p->nvars, &err) ; ;
      f->na ++ ;      
    }
    
    symbol = scanner_read(scanner) ;    
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

  /* fprintf(stderr, "  %s\n", scanner->value.v_identifier) ; */

  if ( agg_distribution_axes_parse(scanner->value.v_identifier,
				   d->axes) != 0 )
    return G_TOKEN_ERROR ;
  
  /*look for closing bracket*/
  if ( scanner_read(scanner) != G_TOKEN_RIGHT_PAREN ) return G_TOKEN_ERROR ;
  
  return G_TOKEN_NONE ;
}

static guint parse_string(GScanner *scanner,
			agg_parser_t *p,
			agg_distribution_t *d)

{
  /*read string*/
  /* if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ; */

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
  gchar v[128] ;
  gdouble smin, smax, s ;
  gint ns, ntr, i ;
  agg_function_call_t shape ;
  
  fprintf(stderr, "scanning distribution: ") ;

  /*initialize transformation list*/
  ntr = 0 ;
  shape.func[0] = '\0' ; shape.c = AGG_SPACING_COSINE ;
  
  /*need the global data set first*/
  if ( !agg_parser_global_set(p) ) return G_TOKEN_ERROR ;
  
  /*open brackets*/
  if ( scanner_read(scanner) != G_TOKEN_LEFT_PAREN ) return G_TOKEN_ERROR ;
  
  /*distribution name*/
  if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ;
  sprintf(v, "%s", scanner->value.v_identifier) ;
  fprintf(stderr, "%s\n", v) ;

  /*parameter limits*/
  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_FLOAT ) return G_TOKEN_ERROR ;
  smin = scanner->value.v_float ;
  
  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_FLOAT ) return G_TOKEN_ERROR ;
  smax = scanner->value.v_float ;
  
  /*number of steps in parameter; need to check for ints at this point*/
  scanner->config->int_2_float = FALSE ;
  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_INT   ) return G_TOKEN_ERROR ;
  ns = scanner->value.v_int ;
  
  d = agg_distribution_alloc(ns) ;
  agg_body_distribution_add(b, d, v) ;

  /*number of sections in mesh*/
  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_INT   ) return G_TOKEN_ERROR ;
  i = scanner->value.v_int ;
  agg_distribution_section_number(d) = i ;

  /*number of nodes per section in mesh*/
  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_INT   ) return G_TOKEN_ERROR ;
  i = scanner->value.v_int ;
  agg_distribution_section_node_number(d) = i ;
  
  /*back to default behaviour*/
  scanner->config->int_2_float = TRUE ;

  /*step spacing in generating shapes*/
  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ;
  i = agg_parser_constant_parse(scanner->value.v_identifier, &(d->sg)) ;
  if ( i != 0 )
    g_error("%s: unrecognized constant \"%s\"",
	    __FUNCTION__, scanner->value.v_identifier) ;

  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ;
  i = agg_parser_constant_parse(scanner->value.v_identifier, &(d->sm)) ;
  if ( i != 0 )
    g_error("%s: unrecognized constant \"%s\"",
	    __FUNCTION__, scanner->value.v_identifier) ;

  if ( scanner_read(scanner) != G_TOKEN_COMMA ) return G_TOKEN_ERROR ;
  if ( scanner_read(scanner) != G_TOKEN_STRING ) return G_TOKEN_ERROR ;
  i = agg_parser_constant_parse(scanner->value.v_identifier, &(d->ss)) ;
  if ( i != 0 )
    g_error("%s: unrecognized constant \"%s\"",
	    __FUNCTION__, scanner->value.v_identifier) ;
  
  /*close brackets*/
  if ( scanner_read(scanner) != G_TOKEN_RIGHT_PAREN ) return G_TOKEN_ERROR ;

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
      if ( parse_string(scanner, p, d) != G_TOKEN_NONE )
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

gint agg_parser_body_read(gint fid, GScanner *scanner,
			  agg_parser_t *p, agg_body_t *b)

/*
 * on output:
 * 
 * p contains global variables used in body definition
 *
 * b contains list of distributions, each with transform and shape
 */
  
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
    default:
      token = G_TOKEN_ERROR ;
      break ;
    }
  
    g_scanner_peek_next_token (scanner);
  } while (token == G_TOKEN_NONE &&
	   scanner->next_token != G_TOKEN_EOF &&
	   scanner->next_token != G_TOKEN_ERROR);
  
  if ( token == G_TOKEN_ERROR )
    fprintf(stderr, "parse error at line %d, character %d\n",
	    g_scanner_cur_line(scanner),
	    g_scanner_cur_position(scanner)) ;

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
 * @param t transform whose parameters are to be evaluted.
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

/* @} */
