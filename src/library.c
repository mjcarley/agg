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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif /*HAVE_CONFIG_H*/

static struct {
  gchar *description, *name ;
  agg_section_t *s ;
} section_list[1024] ;
static gint section_number = 0 ;

agg_section_t *agg_library_section_lookup(gchar *name, gchar **description)

{
  gint i ;

  for ( i = 0 ; i < section_number ; i ++ ) {
    if ( strcmp(section_list[i].name, name) == 0 ) {
      if ( description != NULL ) {
	*description = section_list[i].description ;
      }
      return section_list[i].s ;
    }
  }
  
  return NULL ;
}

gint agg_library_section_add(gchar *name, gchar *description,
			     agg_section_t *s)

{
  agg_section_t *sec ;

  sec = agg_library_section_lookup(name, NULL) ;

  if ( sec != NULL ) {
    g_error("%s: section \"%s\" already in library", __FUNCTION__, name) ;
  }

  sec = agg_section_duplicate(s) ;
  section_list[section_number].name = g_strdup(name) ;
  if ( description != NULL ) 
    section_list[section_number].description = g_strdup(description) ;
  else 
    section_list[section_number].description = NULL ;
  section_list[section_number].s = sec ;

  section_number ++ ;
  
  return 0 ;
}

gint agg_library_sections_list(FILE *f, gboolean write_description)

{
  gint i ;

  if ( !write_description ) { 
    for ( i = 0 ; i < section_number ; i ++ ) {
      fprintf(f, "%s\n", section_list[i].name) ;
    }
    return 0 ;
  }
  
  for ( i = 0 ; i < section_number ; i ++ ) {
    fprintf(f,
	    "%s:\n"
	    "  %s\n\n",
	    section_list[i].name, section_list[i].description) ;
  }

  return 0 ;
}

gint agg_library_section_write(FILE *f, gchar *name, gchar *description,
			       agg_section_t *s)

{

  return 0 ;
}

gint agg_library_read(FILE *f)

{

  return 0 ;
}

gint agg_library_write(FILE *f)

{

  return 0 ;
}
