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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

#include <glib.h>

#include <agg.h>

#include "agg-private.h"

#define PARSER_DATA_SIZE    8
#define PARSER_DATA_SECTION 0

void _agg_global_read(GScanner *scanner, agg_body_t *b, gboolean echo,
		      gpointer data[]) ;
void _agg_surface_read(GScanner *scanner, agg_body_t *b, gboolean echo,
		      gpointer data[]) ;
void _agg_name_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		      gpointer data[]) ;
void _agg_section_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
			gpointer data[]) ;
void _agg_transform_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
			  gpointer data[]) ;
void _agg_patch_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		      gpointer data[]) ;
void _agg_axes_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		     gpointer data[]) ;
void _agg_grid_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		     gpointer data[]) ;
void _agg_invert_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		       gpointer data[]) ;
void _agg_limits_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		       gpointer data[]) ;

typedef void (*block_read_func_t)(GScanner *scanner, agg_body_t *b,
				  gboolean echo, gpointer data[]) ;

static const struct {
  char *id ;
  block_read_func_t func ;
} block_data[] = {
  {"global",    _agg_global_read     },
  {"surface",   _agg_surface_read    },
  {"name",      _agg_name_parse      },
  {"section",   _agg_section_parse   },
  {"transform", _agg_transform_parse },
  {"patch",     _agg_patch_parse     },
  {"axes",      _agg_axes_parse      },
  {"grid",      _agg_grid_parse      },
  {"invert",    _agg_invert_parse    },
  {"limits",    _agg_limits_parse    },
  {NULL, NULL}
} ;

#define token_is_left_bracket(_t)					\
  (((_t) == G_TOKEN_LEFT_PAREN) || ((_t) == G_TOKEN_LEFT_CURLY) ||	\
   ((_t) == G_TOKEN_LEFT_BRACE))
#define token_is_right_bracket(_t)					\
  (((_t)== G_TOKEN_RIGHT_PAREN) || ((_t)== G_TOKEN_RIGHT_CURLY) ||	\
   ((_t)== G_TOKEN_RIGHT_BRACE))
#define token_is_bracket(_t)						\
  ((token_is_left_bracket((_t))) ||(token_is_right_bracket((_t))))
#define token_is_numeric(_t)						\
  ((_t) == '-' || (_t) == G_TOKEN_INT || (_t) == G_TOKEN_FLOAT )

#define tokens_match(_l,_r)			\
  ((((_l) == G_TOKEN_LEFT_PAREN) && ((_r) == G_TOKEN_RIGHT_PAREN)) || \
   (((_l) == G_TOKEN_LEFT_CURLY) && ((_r) == G_TOKEN_RIGHT_CURLY)) || \
   (((_l) == G_TOKEN_LEFT_BRACE) && ((_r) == G_TOKEN_RIGHT_BRACE)))

#define token_read_and_check(_s,_check,_msg)	     \
  do {						     \
    GTokenType _token ;				     \
    (_token) = g_scanner_get_next_token((_s)) ;	     \
    if ( (_token) != (_check) ) {		     \
      g_error("%s: %s on line %u", __FUNCTION__, (_msg),	\
	      g_scanner_cur_line((_s))) ;			\
    }								\
  } while (0)

#define token_read_string(_s,_str)				\
  do {								\
    GTokenType _token ;						\
    (_token) = g_scanner_get_next_token((_s)) ;			\
    if ( (_token) != G_TOKEN_STRING ) {				\
      g_error("%s: string required on line %u, character %u",	\
	      __FUNCTION__, g_scanner_cur_line((_s)),		\
	      g_scanner_cur_position((_s))) ;			\
    }								\
    (_str) = ((_s)->value.v_string) ;				\
  } while (0)

#define token_copy_string(_s,_str)				\
  do {								\
    GTokenType _token ;						\
    (_token) = g_scanner_get_next_token((_s)) ;			\
    if ( (_token) != G_TOKEN_STRING ) {				\
      g_error("%s: string required on line %u, character %u",	\
	      __FUNCTION__, g_scanner_cur_line((_s)),		\
	      g_scanner_cur_position((_s))) ;			\
    }								\
    strcpy((_str),(_s)->value.v_string) ;			\
  } while (0)

/* (_str) = g_strdup((_s)->value.v_string)		\ */
#if 0
/*useful for debugging parser*/
static char *token_string(GTokenType token)

{
  switch (token) {
  case G_TOKEN_EOF:              return "G_TOKEN_EOF" ;
  case G_TOKEN_LEFT_PAREN:       return "G_TOKEN_LEFT_PAREN" ;
  case G_TOKEN_RIGHT_PAREN:      return "G_TOKEN_RIGHT_PAREN" ;
  case G_TOKEN_LEFT_CURLY:       return "G_TOKEN_LEFT_CURLY" ;
  case G_TOKEN_RIGHT_CURLY:      return "G_TOKEN_RIGHT_CURLY" ;
  case G_TOKEN_LEFT_BRACE:       return "G_TOKEN_LEFT_BRACE" ;
  case G_TOKEN_RIGHT_BRACE:      return "G_TOKEN_RIGHT_BRACE" ;
  case G_TOKEN_EQUAL_SIGN:       return "G_TOKEN_EQUAL_SIGN" ;
  case G_TOKEN_COMMA:            return "G_TOKEN_COMMA" ;
  case G_TOKEN_NONE:             return "G_TOKEN_NONE" ;
  case G_TOKEN_ERROR:            return "G_TOKEN_ERROR" ;
  case G_TOKEN_CHAR:             return "G_TOKEN_CHAR" ;
  case G_TOKEN_BINARY:           return "G_TOKEN_BINARY" ;
  case G_TOKEN_OCTAL:            return "G_TOKEN_OCTAL" ;
  case G_TOKEN_INT:              return "G_TOKEN_INT" ;
  case G_TOKEN_HEX:              return "G_TOKEN_HEX" ;
  case G_TOKEN_FLOAT:            return "G_TOKEN_FLOAT" ;
  case G_TOKEN_STRING:           return "G_TOKEN_STRING" ;
  case G_TOKEN_SYMBOL:           return "G_TOKEN_SYMBOL" ;
  case G_TOKEN_IDENTIFIER:       return "G_TOKEN_IDENTIFIER" ;
  case G_TOKEN_IDENTIFIER_NULL:  return "G_TOKEN_IDENTIFIER_NULL" ;
  case G_TOKEN_COMMENT_SINGLE:   return "G_TOKEN_COMMENT_SINGLE" ;
  case G_TOKEN_COMMENT_MULTI:    return "G_TOKEN_COMMENT_MULTI" ;
  case G_TOKEN_LAST:             return "G_TOKEN_COMMENT_LAST" ;
  }
  
  return NULL ;
}  
#endif

static void echo_variable(FILE *f, char *name, char *def, gdouble val)

{
  if ( def != NULL ) {
    fprintf(f, "%s = %s\n", name, def) ;
  } else {
    fprintf(f, "%s = %lg\n", name, val) ;
  }

  return ;
}

static gdouble scan_numeric(GScanner *scanner)

{
  gdouble v, sgn ;
  GTokenType token ;
  
  sgn = 1 ;
  token = g_scanner_get_next_token(scanner) ;
  if ( token == '-' ) {
    sgn = -1 ; 
    token = g_scanner_get_next_token(scanner) ;
  }
  g_assert(token_is_numeric(token)) ;
  /*this requires that the scanner convert int to float*/
  v = sgn*scanner->value.v_float ;

  return v ;

}

static char *scan_string(GScanner *scanner)

{
  char *s ;
  GTokenType token ;
  
  token = g_scanner_get_next_token(scanner) ;
  g_assert(token == G_TOKEN_STRING ) ;
  s = g_strdup(scanner->value.v_identifier) ;

  return s ;
}

static block_read_func_t block_read_func(char *id)

{
  gint i ;

  for ( i = 0 ; block_data[i].id != NULL ; i ++ ) {
    if ( strcmp(block_data[i].id, id) == 0 ) return block_data[i].func ;
  }
  
  return NULL ;
}

void _agg_name_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		     gpointer data[])
		     
{
  gint s ;
  GTokenType token ;
  
  s = agg_body_surface_number(b) - 1 ;

  if ( echo ) {
    fprintf(stderr, "line %u: setting name of surface %d\n",
	    g_scanner_cur_line(scanner), s) ;	    
  }
  
  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  token = g_scanner_get_next_token(scanner) ;
  if ( token != G_TOKEN_STRING ) {
    g_error("%s: string required for name declaration at line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;	    
  }

  agg_body_surface_name(b,s) = g_strdup(scanner->value.v_string) ;
  
  token_read_and_check(scanner, G_TOKEN_RIGHT_PAREN, "missing right bracket") ;

  if ( echo ) {
    fprintf(stderr, "line %u: name of surface %d set to \"%s\"\n",
	    g_scanner_cur_line(scanner), s,
	    agg_body_surface_name(b, s)) ;	    
  }
  
  return ;
}

static void parameter_list_parse(GScanner *scanner,
				 agg_variable_t *params, gint *nparams)

{
  char *def ;
  gdouble val ;

  (*nparams) = 0 ;
  while ( g_scanner_peek_next_token(scanner) != G_TOKEN_RIGHT_PAREN ) {
    def = NULL ; val = 0 ;
    if ( g_scanner_peek_next_token(scanner) == G_TOKEN_COMMA ) {
      g_scanner_get_next_token(scanner) ;
    }
    if ( token_is_numeric(g_scanner_peek_next_token(scanner)) ) {
      val = scan_numeric(scanner) ;
    } else {
      def = scan_string(scanner) ;
    }
    /* fprintf(stderr, "  (%lg, %s)\n", val, def) ; */
    params[(*nparams)].val = val ;
    params[(*nparams)].def = def ;
    params[(*nparams)].name = NULL ;
    (*nparams) ++ ;
  }

  g_scanner_get_next_token(scanner) ;
  
  return ;
}

void _agg_section_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
			gpointer data[])

{
  agg_surface_t *S ;
  gint srf, nparams, i ;
  gdouble u ;
  agg_section_t *s = data[PARSER_DATA_SECTION] ;
  agg_variable_t params[64] ;
  
  srf = agg_body_surface_number(b) - 1 ;
  S = agg_body_surface_last(b) ;
  if ( echo ) {
    fprintf(stderr, "line %u: setting section of surface %d\n",
	    g_scanner_cur_line(scanner), srf) ;
  }
  
  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  parameter_list_parse(scanner, params, &nparams) ;

  if ( nparams < 1 ) {
    g_error("%s: not enough parameters, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  }
  
  /*first parameter is u: if not specified, defaults to 0*/
  i = 0 ; u = 0.0 ;
  if ( agg_variable_definition(&(params[0])) == NULL ) {
    u = agg_variable_value(&(params[0])) ;
    i = 1 ; nparams -- ;
  }

  if ( agg_variable_definition(&(params[i])) == NULL ) {
    g_error("%s: section definition must be a string, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;    
  }
  
  agg_section_parse(s, agg_variable_definition(&(params[i])),
		    &(params[i+1]), nparams-1) ;
  agg_surface_section_add(S, s, u) ;
  
  if ( echo ) {
    fprintf(stderr, "line %u: section of surface %d at u=%lg set to \"%s\"\n",
	    g_scanner_cur_line(scanner), srf, u,
	    agg_variable_definition(&(params[i]))) ;
    fprintf(stderr, "parameters: ") ;
    for ( i = 0 ; i < nparams-1 ; i ++ ) {
      agg_variable_write(stderr, &(params[i])) ;
      fprintf(stderr, ", ") ;
    }      
    agg_variable_write(stderr, &(params[nparams-1])) ;
    fprintf(stderr, "\n") ;
  }

  return ;
}

void _agg_transform_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
			  gpointer data[])

{
  agg_surface_t *S ;
  gint srf, nparams, i ;
  agg_variable_t params[64] ;
  agg_transform_t *T ;
  agg_affine_t *A ;
  
  srf = agg_body_surface_number(b) - 1 ;
  S = agg_body_surface_last(b) ;
  T = agg_surface_transform(S) ;

  if ( echo ) {
    fprintf(stderr, "line %u: setting transform of surface %d\n",
	    g_scanner_cur_line(scanner), srf) ;
  }
  
  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  parameter_list_parse(scanner, params, &nparams) ;

  if ( nparams < 1 ) {
    g_error("%s: not enough parameters, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  }
  
  /*first parameter is name of transform or first two parameters are
    limits*/
  if ( agg_variable_definition(&(params[0])) == NULL &&
       agg_variable_definition(&(params[1])) != NULL ) {
    g_error("%s: transform definition must be a string, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  }

  A = agg_affine_new(1) ;
  agg_affine_parse(A, params, nparams) ;
  agg_affine_differentiate(A, "u") ;
  agg_transform_affine_add(T, A) ;

  /* agg_transform_parse(T, params, nparams) ; */
  
  if ( echo ) {
    fprintf(stderr, "line %u: transform \"%s\" added to surface %d\n",
	    g_scanner_cur_line(scanner),
	    agg_variable_definition(&(params[0])), srf) ;
    fprintf(stderr, "parameters: ") ;
    for ( i = 0 ; i < nparams-1 ; i ++ ) {
      agg_variable_write(stderr, &(params[i])) ;
      fprintf(stderr, ", ") ;
    }      
    agg_variable_write(stderr, &(params[nparams-1])) ;
    fprintf(stderr, "\n") ;
  }

  return ;
}

void _agg_patch_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		      gpointer data[])

{
  gint nparams ;
  agg_variable_t params[64] ;

  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  parameter_list_parse(scanner, params, &nparams) ;
  agg_patch_parse(agg_body_patch_last(b), params, nparams) ;
  
  return ;
}

void _agg_axes_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		     gpointer data[])

{
  gint nparams ;
  agg_variable_t params[64] ;
  agg_surface_t *S ;
  agg_transform_t *T ;
  agg_affine_t *A ;
  agg_axes_t axes ;
  
  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  S = agg_body_surface_last(b) ;
  T = agg_surface_transform(S) ;
  
  parameter_list_parse(scanner, params, &nparams) ;

  if ( nparams > 1 )
    g_error("%s: axes takes one parameter, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;

  if ( agg_variable_definition(&(params[0])) == NULL )
    g_error("%s: axes must be specified by string, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  
  /* agg_surface_axes(S) = agg_axes_parse(agg_variable_definition(&(params[0]))) ; */
  axes = agg_axes_parse(agg_variable_definition(&(params[0]))) ;

  if ( axes == AGG_AXES_UNDEFINED )
    g_error("%s: unrecognized axes \"%s\" on line %u",
	    __FUNCTION__, agg_variable_definition(&(params[0])),
	    g_scanner_cur_line(scanner)) ;
  
  A = agg_affine_new(1) ;
  agg_affine_axes(A, axes) ;
  agg_affine_differentiate(A, "u") ;
  agg_transform_affine_add(T, A) ;

  return ;
}

void _agg_invert_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		       gpointer data[])

{
  gint nparams ;
  agg_variable_t params[64] ;
  agg_patch_t *P ;
  
  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  P = agg_body_patch_last(b) ;
  
  parameter_list_parse(scanner, params, &nparams) ;

  if ( nparams > 0 )
    g_error("%s: invert takes no parameters, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;

  agg_patch_invert(P) = TRUE ;
  
  return ;
}

void _agg_limits_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		       gpointer data[])

{
  gint nparams ;
  agg_variable_t params[64] ;
  agg_surface_t *S ;
  
  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  S = agg_body_surface_last(b) ;
  
  parameter_list_parse(scanner, params, &nparams) ;

  if ( nparams != 2 ) {
    g_error("%s: limits requires two parameters, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  }
    
  if ( agg_variable_definition(&(params[0])) != NULL )
    g_error("%s: lower u limit must be constant, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  
  if ( agg_variable_definition(&(params[1])) != NULL )
    g_error("%s: upper u limit must be constant, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;

  agg_surface_umin(S) = agg_variable_value(&(params[0])) ;
  agg_surface_umax(S) = agg_variable_value(&(params[1])) ;
  
  return ;
}

void _agg_grid_parse(GScanner *scanner, agg_body_t *b, gboolean echo,
		     gpointer data[])

{
  gint nparams ;
  agg_variable_t params[64] ;
  agg_surface_t *S ;
  agg_grid_t grid ;
  
  token_read_and_check(scanner, G_TOKEN_LEFT_PAREN, "missing left bracket") ;

  S = agg_body_surface_last(b) ;
  
  parameter_list_parse(scanner, params, &nparams) ;

  if ( nparams < 2 ) {
    g_error("%s: grid requires at least two parameters, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  }
  if ( agg_variable_definition(&(params[0])) == NULL )
    g_error("%s: grid type must be specified by string, line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;

  if ( (grid = agg_grid_parse(agg_variable_definition(&(params[0])))) ==
       AGG_GRID_UNDEFINED) {
    g_error("%s: unrecognized grid type \"%s\" , line %u",
	    __FUNCTION__,
	    agg_variable_definition(&(params[0])),
	    g_scanner_cur_line(scanner)) ;    
  }

  agg_surface_grid(S) = grid ;

  if ( grid == AGG_GRID_REGULAR ) {
    if ( nparams != 3 ) {
      g_error("%s: regular grid requires two further parameters, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    }
    
    if ( agg_variable_definition(&(params[1])) != NULL )
      g_error("%s: number of grid sections must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    if ( agg_variable_definition(&(params[2])) != NULL )
      g_error("%s: number of grid section splines must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    agg_surface_grid_section_number(S) = (gint)(params[1].val) ;
    agg_surface_grid_spline_number(S) = (gint)(params[2].val) ;
  
    return ;
  }

  if ( grid == AGG_GRID_TRIANGLE ) {
    if ( nparams != 2 ) {
      g_error("%s: triangle grid requires one further parameter, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    }
    
    if ( agg_variable_definition(&(params[1])) != NULL )
      g_error("%s: element area must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    agg_surface_grid_element_area(S) = params[1].val ;
  
    return ;
  }

  if ( grid == AGG_GRID_SPHERE_ICO ) {
    if ( nparams != 2 ) {
      g_error("%s: icosahedron grid requires one further parameter, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    }
    
    if ( agg_variable_definition(&(params[1])) != NULL )
      g_error("%s: surface subdivision must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    agg_surface_grid_subdivision(S) = (gint)(params[1].val) ;

    if ( agg_surface_grid_subdivision(S) < 0 ) 
      g_error("%s: invalid surface subdivision (%lg), line %u",
	      __FUNCTION__, params[1].val, g_scanner_cur_line(scanner)) ;
    
    return ;
  }

  if ( grid == AGG_GRID_SPHERE_UV ) {
    if ( nparams != 3 ) {
      g_error("%s: spherical UV grid requires two further parameters, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    }
    
    if ( agg_variable_definition(&(params[1])) != NULL )
      g_error("%s: number of grid sections must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    if ( agg_variable_definition(&(params[2])) != NULL )
      g_error("%s: number of grid section splines must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    agg_surface_grid_section_number(S) = (gint)(params[1].val) ;
    agg_surface_grid_spline_number(S) = (gint)(params[2].val) ;
  
    return ;
  }

  if ( grid == AGG_GRID_HEMISPHERE_ICO ) {
    if ( nparams != 2 ) {
      g_error("%s: icosahedron grid requires one further parameter, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    }
    
    if ( agg_variable_definition(&(params[1])) != NULL )
      g_error("%s: surface subdivision must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    agg_surface_grid_subdivision(S) = (gint)(params[1].val) ;

    if ( agg_surface_grid_subdivision(S) < 0 ) 
      g_error("%s: invalid surface subdivision (%lg), line %u",
	      __FUNCTION__, params[1].val, g_scanner_cur_line(scanner)) ;
    
    return ;
  }

  if ( grid == AGG_GRID_HEMISPHERE_UV ) {
    if ( nparams != 3 ) {
      g_error("%s: hemispherical UV grid requires two further "
	      "parameters, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    }
    
    if ( agg_variable_definition(&(params[1])) != NULL )
      g_error("%s: number of grid sections must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    if ( agg_variable_definition(&(params[2])) != NULL )
      g_error("%s: number of grid section splines must be constant, line %u",
	      __FUNCTION__, g_scanner_cur_line(scanner)) ;
    
    agg_surface_grid_section_number(S) = (gint)(params[1].val) ;
    agg_surface_grid_spline_number(S) = (gint)(params[2].val) ;
  
    return ;
  }
  
  g_assert_not_reached() ;
  
  return ;
}

void _agg_global_read(GScanner *scanner, agg_body_t *b, gboolean echo,
		      gpointer data[])

{
  GTokenType token ;
  char *name, *def ;
  gdouble val ;
  
  if ( echo )
    fprintf(stderr, "line %u: globals\n", g_scanner_cur_line(scanner)) ;
  token = g_scanner_get_next_token(scanner) ;

  if ( token != G_TOKEN_LEFT_CURLY ) {
    g_error("%s: expected { on line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  }
  
  while ( ( token = g_scanner_peek_next_token(scanner) )
	  != G_TOKEN_RIGHT_CURLY ) {
    token = g_scanner_get_next_token(scanner) ;
    def = name = NULL ; val = 0.0 ;
    if ( scanner->token == G_TOKEN_IDENTIFIER ) {
      /*variable assignment*/
      name = g_strdup(scanner->value.v_identifier) ;
      def = NULL ;
      token = g_scanner_get_next_token(scanner) ;
      if ( token != G_TOKEN_EQUAL_SIGN ) {
	g_error("%s: misformed assignment on line %u",
		__FUNCTION__, g_scanner_cur_line(scanner)) ;
      }
      if ( token_is_numeric(g_scanner_peek_next_token(scanner)) ) {
	val = scan_numeric(scanner) ;
      } else {
	def = scan_string(scanner) ;
      }	
      if ( echo ) echo_variable(stderr, name, def, val) ;
    }

    agg_body_global_add(b, name, def, val) ;    
  }

  if ( echo )
    fprintf(stderr, "line %u: globals read\n", g_scanner_cur_line(scanner)) ;
  g_scanner_get_next_token(scanner) ;
  
  return ;
}

void _agg_surface_read(GScanner *scanner, agg_body_t *b, gboolean echo,
		       gpointer data[])

{
  GTokenType token ;
  agg_surface_t *S ;
  agg_patch_t *P ;
  block_read_func_t func ;
  agg_transform_t *T ;
  
  if ( echo )
    fprintf(stderr, "line %u: surface\n", g_scanner_cur_line(scanner)) ;
  token = g_scanner_get_next_token(scanner) ;

  S = agg_surface_new(64) ;
  P = agg_patch_new() ;
  agg_body_surface_add(b, S, P) ;

  T = agg_surface_transform(S) ;
  agg_transform_add_global_variables(T,
				     agg_body_globals(b),
				     agg_body_global_number(b)) ;
  
  if ( token != G_TOKEN_LEFT_CURLY ) {
    g_error("%s: expected { on line %u",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  }
  
  while ( ( token = g_scanner_peek_next_token(scanner) )
	  != G_TOKEN_RIGHT_CURLY ) {
    token = g_scanner_get_next_token(scanner) ;
    if ( scanner->token == G_TOKEN_IDENTIFIER ) {
      func = block_read_func(scanner->value.v_identifier) ;
      if ( func == NULL ) {
	g_error("%s: unrecognized identifier \"%s\"",
		__FUNCTION__, scanner->value.v_identifier) ;
      } else {
	func(scanner, b, echo, data) ;
      }
    }
  }

  /*sanity check*/
  if ( agg_surface_section_number(S) == 0 )
    g_error("%s: on line %u, surface contains no section data",
	    __FUNCTION__, g_scanner_cur_line(scanner)) ;
  
  if ( echo )
    fprintf(stderr, "line %u: surface read\n", g_scanner_cur_line(scanner)) ;
  g_scanner_get_next_token(scanner) ;

  /*initialize anything that needs doing*/
  /* agg_surface_umin(S) = 0.0 ;  */
  /* agg_surface_umax(S) = 1.0 ;  */

  agg_expression_data_compile(T->e) ;
  agg_transform_expressions_compile(T) ;

  agg_surface_weights_make(S) ;
  
  return ;
}

/**
 * @{
 *  @ingroup body
 */

/** 
 * Allocate and initialize a new ::agg_body_t
 * 
 * @param ngmax maximum number of global variables in body;
 * @param nsmax maximum number of surfaces in body.
 * 
 * @return newly allocated ::agg_body_t
 */

agg_body_t *agg_body_new(gint ngmax, gint nsmax)

{
  agg_body_t *b ;
  
  b = (agg_body_t *)g_malloc0(sizeof(agg_body_t)) ;

  b->g = (agg_variable_t *)g_malloc0(ngmax*sizeof(agg_variable_t)) ;

  b->S = (agg_surface_t **)g_malloc(nsmax*sizeof(agg_surface_t *)) ;
  b->P = (agg_patch_t **)g_malloc(nsmax*sizeof(agg_patch_t *)) ;
  b->names = (char **)g_malloc(nsmax*sizeof(char *)) ;

  b->e = agg_expression_data_new(ngmax) ;
  b->e->ne = 0 ;
  
  agg_body_global_number(b) = 0 ;
  agg_body_global_number_max(b) = ngmax ;
  agg_body_surface_number(b) = 0 ;
  agg_body_surface_number_max(b) = nsmax ;
  
  return b ;
}

/** 
 * Add a global variable to an ::agg_body_t
 * 
 * @param b body to which variable should be added;
 * @param var name of variable;
 * @param def string expression defining variable (can be NULL);
 * @param val value of variable (set constant if \a expr is NULL).
 * 
 * @return 0 on success.
 */

gint agg_body_global_add(agg_body_t *b, char *var, char *def, gdouble val)

{
  gint i ;
  agg_variable_t *g ;
  
  /*check for duplicate variable*/
  for ( i = 0 ; i < agg_body_global_number(b) ; i ++ ) {
    g = agg_body_global(b,i) ;
    if ( strcmp(agg_variable_name(g), var) == 0 ) {
      g_error("%s: variable %s already in body", __FUNCTION__, var) ;
    }
  }

  if ( agg_body_global_number(b) >= agg_body_global_number_max(b) )
    g_error("%s: not enough room for new variable (%d available)",
	    __FUNCTION__, agg_body_global_number_max(b)) ;
  
  i = agg_body_global_number(b) ;
  g = agg_body_global(b,i) ;
  
  agg_variable_name(g) = g_strdup(var) ;
  if ( def == NULL ) agg_variable_definition(g) = NULL ;
  else agg_variable_definition(g) = g_strdup(def) ;
  agg_variable_value(g) = val ;

  agg_body_global_number(b) ++ ;
  
  return 0 ;
}

/** 
 * Write a body's global variables, their definitions (if given), and
 * values to file
 * 
 * @param f output file;
 * @param b ::agg_body_t whose global variables are to be written
 * 
 * @return 0 on success.
 */

gint agg_body_globals_write(FILE *f, agg_body_t *b)

{
  gint i ;
  agg_variable_t *g ;

  for ( i = 0 ; i < agg_body_global_number(b) ; i ++ ) {
    g = agg_body_global(b,i) ;
    fprintf(f, "%s ", agg_variable_name(g)) ;
    if ( agg_variable_definition(g) != NULL )
      fprintf(f, "= %s ", agg_variable_definition(g)) ;
    fprintf(f, "= %lg\n", agg_variable_value(g)) ;
      
  }
  
  return 0 ;
}

/** 
 * Compile expressions for global variables in an ::agg_body_t. This
 * function should be called before evaluating any variables in the
 * body or its surfaces.
 * 
 * @param b an ::agg_body_t containing definitions of global variables.
 * 
 * @return 0 on success.
 */

gint agg_body_globals_compile(agg_body_t *b)

{
  agg_expression_data_t *e ;
  agg_variable_t *v ;
  gint i ;
  
  e = b->e ;
  
  for ( i = 0 ; i < agg_body_global_number(b) ; i ++ ) {
    v = agg_body_global(b,i) ;
    agg_expression_data_variable_add(e, v) ;
  }

  for ( i = 0 ; i < agg_body_global_number(b) ; i ++ ) {
    v = agg_body_global(b,i) ;
    if ( agg_variable_definition(v) == NULL ) {
      agg_variable_evaluator(v) = NULL ;
    } else {
      agg_variable_evaluator(v) =
	agg_expression_compile(agg_variable_definition(v), e) ;       
    }
  }

  return 0 ;
}

/** 
 * Evaluate global variables defined in an ::agg_body_t. This function
 * should be called before evaluating points on the surfaces which
 * make up the body.
 * 
 * @param b an ::agg_body_t containing global variables.
 * 
 * @return 0 on success.
 */

gint agg_body_globals_eval(agg_body_t *b)

{
  agg_variable_t *v ;
  gint i ;

  agg_expression_data_eval(b->e) ;
  for ( i = 0 ; i < agg_body_global_number(b) ; i ++ ) {
    v = agg_body_global(b,i) ;
    if ( agg_variable_definition(v) != NULL ) {
      agg_variable_value(v) = agg_expression_eval(agg_variable_evaluator(v)) ;
    }
  }

  return 0 ;
}

/** 
 * Read an ::agg_body_t from file
 * 
 * @param b body to which data in file are to be added;
 * @param file name of file in .agg format;
 * @param echo if TRUE, echo progress to stderr.
 * 
 * @return 0 on success.
 */

gint agg_body_read(agg_body_t *b, char *file, gboolean echo)

{
  GScanner *scanner ;
  GTokenType token, delim[64] ;
  gint fd, ndelim ;
  block_read_func_t func ;
  gpointer data[PARSER_DATA_SIZE] ;
  
  scanner = g_scanner_new(NULL) ;
  scanner->config->int_2_float = TRUE ;
  scanner->config->scan_identifier_1char = TRUE ;
  
  data[PARSER_DATA_SECTION] = agg_section_new(64, 64) ;
  
  fd = open(file, O_RDONLY) ;
  if ( fd < 0 ) g_error("%s: cannot open file %s", __FUNCTION__, file) ;

  g_scanner_input_file(scanner, fd) ;

  ndelim = 0 ;
  while ( ( token = g_scanner_get_next_token(scanner) ) != G_TOKEN_EOF ) {
    if ( token_is_left_bracket(token) ) {
      delim[ndelim] = token ; ndelim ++ ;
    }
    if ( token_is_right_bracket(token) ) {
      if ( ndelim > 0 && tokens_match(delim[ndelim-1], token)) {
	ndelim -- ;
      } else {
	g_error("%s: mismatched brackets on line %u, character %u, "
		"of file %s\n",
		__FUNCTION__,
		g_scanner_cur_line(scanner), 
		g_scanner_cur_position(scanner), file) ;
      }
    }
    
    if ( scanner->token == G_TOKEN_IDENTIFIER ) {
      func = block_read_func(scanner->value.v_identifier) ;
      if ( func == NULL ) {
	g_error("%s: unrecognized block \"%s\" at line %u",
		__FUNCTION__, scanner->value.v_identifier,
		g_scanner_cur_line(scanner)) ;
      } else {
	func(scanner, b, echo, data) ;
      }
    }
  }
  
  close(fd) ;
  
  g_scanner_destroy(scanner) ;
  
  return 0 ;
}

/** 
 * Add a surface to a body
 * 
 * @param b ::agg_body_t to which surface is to be added;
 * @param S an ::agg_surface_t;
 * @param P ::agg_patch_t for mapping of surface parameters.
 * 
 * @return 0 on success.
 */

gint agg_body_surface_add(agg_body_t *b, agg_surface_t *S, agg_patch_t *P)

{
  gint ns ;

  ns = agg_body_surface_number(b) ;
  if ( ns >= agg_body_surface_number_max(b) ) {
    g_error("%s: not enough space for surface (%d allocated)",
	    __FUNCTION__, agg_body_surface_number(b)) ;
  }

  agg_body_surface(b,ns) = S ;
  agg_body_patch(b,ns) = P ;

  agg_body_surface_number(b) ++ ;
  
  return 0 ;
}

/** 
 * Write the names of a body's surfaces to file
 * 
 * @param f output file stream;
 * @param b an ::agg_body_t whose surfaces are to be listed.
 * 
 * @return 0 on success.
 */

gint agg_body_surfaces_list(FILE *f, agg_body_t *b)

{
  gint i ;

  for ( i = 0 ; i < agg_body_surface_number(b) ; i ++ ) {
    fprintf(f, "surface %d: \"%s\"\n", i, agg_body_surface_name(b,i)) ;
  }
  
  return 0 ;
}

/**
 *  @}
 */
