/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2022, 2023 Michael Carley
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

#include <glib.h>

#include <blaswrap.h>

#include <agg.h>

#include "agg-private.h"

static gboolean points_coincide(gdouble *x, gdouble *y)

{
  return ((x[0] == y[0]) && (x[1] == y[1]) && (x[2] == y[2])) ;
}

/** 
 * Write the points of an ::agg_mesh_t to file
 * 
 * Write points of a mesh to file. Output format is:
 *
 * (number of points) (number of data elements per point) (number of tags)
 *
 * (point index) \f$x\f$ \f$y\f$ \f$z\f$ \f$u\f$ \f$v\f$ (data elements) (tags)
 *
 * with \f$(x,y,z)\$ physical coordinates, \f$(u,v)\f$ parametric
 * surface coordinates.
 * 
 * @param f file pointer for output;
 * @param m mesh containing triangulation.
 * 
 * @return 0 on success.
 */

gint agg_mesh_points_write(FILE *f, agg_mesh_t *m)

{
  gint i, j ;
  gdouble *x ;
  
  /*number of points, number of point data elements, number of point tags*/
  fprintf(f, "%d %d %d\n",
	  agg_mesh_point_number(m),
	  agg_mesh_point_data_number(m),
	  agg_mesh_point_tag_number(m)) ;
  for ( i = 0 ; i < agg_mesh_point_number(m) ; i ++ ) {
    x = &(agg_mesh_point_x(m,i)) ;
    fprintf(f, "%d", i) ;
    for ( j = 0 ; j < 3+2+agg_mesh_point_data_number(m) ; j ++ ) 
      fprintf(f, " %lg", x[j]) ;
    fprintf(f, " %d", agg_mesh_point_distribution(m,i)) ;
    for ( j = 0 ; j < agg_mesh_point_data_number(m) ; j ++ )
      fprintf(f, " %d", agg_mesh_point_tag(m,i,j)) ;
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

/** 
 * Write the triangular elements of an ::agg_mesh_t to file
 *
 * Write triangulation to file. Output format is
 *
 * (number of triangles) (number of tags per element)
 *
 * (index 1) (index 2) (index 3) (data tags)
 * 
 * @param f file pointer for output;
 * @param m mesh containing triangulation.
 * 
 * @return 0 on success.
 */

gint agg_mesh_tri_write(FILE *f, agg_mesh_t *m)

{
  gint i, j, *tr ;

  /*number of triangles, number of triangle tags*/
  fprintf(f, "%d %d\n",
	  agg_mesh_triangle_number(m),
	  agg_mesh_triangle_tag_number(m)) ;
  for ( i = 0 ; i < agg_mesh_triangle_number(m) ; i ++ ) {
    tr = agg_mesh_triangle(m,i) ;
    for ( j = 0 ; j < 3 + 1 + agg_mesh_triangle_tag_number(m) ; j ++ )
      fprintf(f, " %d", tr[j]) ;
    fprintf(f, "\n") ;
    /* fprintf(f, "%d %d %d", tr[0], tr[1], tr[2]) ; */
    
  }
  
  return 0 ;
}

agg_mesh_t *agg_mesh_alloc(gint np, gint nt, gint nptags, gint nttags,
			   gint ndat)

{
  agg_mesh_t *m ;

  m = (agg_mesh_t *)g_malloc0(sizeof(agg_mesh_t)) ;

  m->npmax = np ;
  m->ntmax = nt ;
  m->nptags = nptags ;
  m->nttags = nttags ;
  m->ndat = ndat ;
  
  /*3 vector for position; 2 values for surface coordinates;
    additional data*/
  m->x = (gdouble *)g_malloc0(np*(3+2+ndat)*sizeof(gdouble)) ;
  /*extra tag for the distribution index*/
  m->ptags = (gint *)g_malloc0(np*(nptags+1)*sizeof(gint)) ;
  m->tri = (gint *)g_malloc0(nt*(3+1+nttags)*sizeof(gint)) ;
  m->bboxes = (agg_bbox_t *)g_malloc0(nt*sizeof(agg_bbox_t)) ;
  
  agg_mesh_init(m) ;
  
  return m ;
}

gint agg_mesh_init(agg_mesh_t *m)

{
  m->np = 0 ;
  m->nt = 0 ;
  m->ndist = 0 ;
  
  return 0 ;
}

gint agg_mesh_element_interp(agg_mesh_t *m, gint i,
			     gdouble s, gdouble t,
			     agg_parser_t *p, agg_shape_t *sh,
			     gdouble *x)

{
  gdouble L[3], u, v ;
  gint *tri, id ;
  agg_distribution_t *d ;
  
  /*linear shape functions*/
  L[0] = 1.0 - s - t ; L[1] = s ; L[2] = t ;

  /*select triangular element*/
  tri = agg_mesh_triangle(m, i) ;
  id = agg_mesh_triangle_distribution(m, i) ;
  d = agg_mesh_distribution(m, id) ;
  
  u = L[0]*agg_mesh_point_u(m,tri[0]) +
    L[1]*agg_mesh_point_u(m,tri[1]) +
    L[2]*agg_mesh_point_u(m,tri[2]) ;
  v = L[0]*agg_mesh_point_v(m,tri[0]) +
    L[1]*agg_mesh_point_v(m,tri[1]) +
    L[2]*agg_mesh_point_v(m,tri[2]) ;

  agg_distribution_point_eval(d, u, v, p, sh, x) ;
  
  return 0 ;
}

gint agg_mesh_element_interp_normal(agg_mesh_t *m, gint i,
				    gdouble s, gdouble t,
				    agg_parser_t *p, agg_shape_t *sh,
				    gdouble *x, gdouble *n, gdouble *J)

{
  gdouble L[3], u, v ;
  gint *tri, id ;
  agg_distribution_t *d ;
  
  /*linear shape functions*/
  L[0] = 1.0 - s - t ; L[1] = s ; L[2] = t ;

  /*select triangular element*/
  tri = agg_mesh_triangle(m, i) ;
  id = agg_mesh_triangle_distribution(m, i) ;
  d = agg_mesh_distribution(m, id) ;
  
  u = L[0]*agg_mesh_point_u(m,tri[0]) +
    L[1]*agg_mesh_point_u(m,tri[1]) +
    L[2]*agg_mesh_point_u(m,tri[2]) ;
  v = L[0]*agg_mesh_point_v(m,tri[0]) +
    L[1]*agg_mesh_point_v(m,tri[1]) +
    L[2]*agg_mesh_point_v(m,tri[2]) ;

  agg_distribution_point_normal_eval(d, u, v, p, sh, x, n, J) ;
  
  return 0 ;
}

gint agg_body_mesh_grid(agg_body_t *b, agg_grid_t *g,
			agg_mesh_t *m, agg_workspace_t *w)

{
  gint i, idx, np, ntri, xstr, tstr, *ptags, ptstr ;
  gint *tr, *ttags, ttstr, *tg ;
  gdouble *x, u, v ;
  agg_parser_t *p ;

  g_assert(g != NULL) ;
  g_assert((p = agg_body_parser(b)) != NULL) ;

  /* if ( g->init_func != NULL ) { */
  /*   g->init_func(g, b, m, w) ; */
  /* } */
  
  np = agg_mesh_point_number(m) ;
  ntri = agg_mesh_triangle_number(m) ;
  xstr = 3 + 2 + m->ndat ;
  tstr = 3 + 1 + m->nttags ;
  ptags = &(agg_mesh_point_distribution(m,np)) ;
  ptstr = agg_mesh_point_tag_number(m)+1 ;

  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) 
    agg_mesh_distribution(m, i) = agg_body_distribution(b, i) ;
  agg_mesh_distribution_number(m) = agg_body_distribution_number(b) ;

  x = &(agg_mesh_point_x(m,np)) ;
  for ( i = idx = 0 ; i < agg_grid_point_number(g) ; i ++ ) {
    u = agg_grid_point_u(g,i) ;
    v = agg_grid_point_v(g,i) ;

    agg_body_point_eval(b, u, v, &(x[idx*xstr]), w) ;

    x[idx*xstr+3] = u ;
    x[idx*xstr+4] = v ;
    ptags[idx*ptstr] = 0 ;
    idx ++ ;    
  }
  
  m->np += idx ;
  tr = agg_mesh_triangle(m,ntri) ;
  ttags = &(agg_mesh_triangle_distribution(m,ntri)) ;
  ttstr = agg_mesh_triangle_tag_number(m)+3+1 ;

  for ( i = 0 ; i < agg_grid_triangle_number(g) ; i ++ ) {
    tg = agg_grid_triangle(g, i) ;
    tr[i*(tstr)+0] = tg[0] + np ;
    tr[i*(tstr)+1] = tg[1] + np ;
    tr[i*(tstr)+2] = tg[2] + np ;
    ttags[i*ttstr] = agg_mesh_distribution_number(m) ;
    ntri ++ ;
  }

  agg_mesh_triangle_number(m) = ntri ;
  agg_mesh_grid(m) = g ;

  /* if ( g->refine_func != NULL ) { */
  /*   g->refine_func(g, b, m, w) ; */
  /* }   */
  
  return 0 ;
}

static void bbox_limits(agg_bbox_t *b, gdouble *x1, gdouble *x2, gdouble *x3)

{
  b->xmin = MIN(x1[0], MIN(x2[0], x3[0])) ;
  b->xmax = MAX(x1[0], MAX(x2[0], x3[0])) ;
  b->ymin = MIN(x1[1], MIN(x2[1], x3[1])) ;
  b->ymax = MAX(x1[1], MAX(x2[1], x3[1])) ;
  b->zmin = MIN(x1[2], MIN(x2[2], x3[2])) ;
  b->zmax = MAX(x1[2], MAX(x2[2], x3[2])) ;
  
  return ;
}

gint agg_mesh_bounding_boxes(agg_mesh_t *m)

{
  gint i, *t ;

  for ( i = 0 ; i < agg_mesh_triangle_number(m) ; i ++ ) {
    t = agg_mesh_triangle(m, i) ;
    bbox_limits(&(m->bboxes[i]),
		agg_mesh_point(m,t[0]),
		agg_mesh_point(m,t[1]),
		agg_mesh_point(m,t[2])) ;
  }
  
  return 0 ;
}
