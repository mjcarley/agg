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

#include <glib.h>

#include <blaswrap.h>

#include <agg.h>

#include "agg-private.h"

static gboolean points_coincide(gdouble *x, gdouble *y)

{
  return ((x[0] == y[0]) && (x[1] == y[1]) && (x[2] == y[2])) ;
}

static void add_elements(gint *tr, gint tstr, gint *ttags, gint ttstr, gint id,
			 gdouble *x, gint xstr, gint i0, gint i1,
			 gint idx0, gint d, gint *ntr)

{
  gint j0, j1, j2, j3 ;
  
  /*check for degenerate quadrilateral*/
  j0 = idx0 ; j1 = idx0 + d ; j2 = idx0 + d + 1 ; j3 = idx0 + 1 ;
  
  if ( !points_coincide(&(x[j0*xstr]), &(x[j1*xstr])) &&
       !points_coincide(&(x[j1*xstr]), &(x[j2*xstr])) ) {
    tr[tstr*(*ntr) +  0] = j0 ;
    tr[tstr*(*ntr) + i0] = j1 ;
    tr[tstr*(*ntr) + i1] = j2 ;
    ttags[(*ntr)*ttstr] = id ;
    (*ntr) ++ ;
  }

  if ( !points_coincide(&(x[j2*xstr]), &(x[j3*xstr])) &&
       !points_coincide(&(x[j3*xstr]), &(x[j0*xstr])) ) {
    tr[tstr*(*ntr) +  0] = j0 ;
    tr[tstr*(*ntr) + i0] = j2 ;
    tr[tstr*(*ntr) + i1] = j3 ;
    ttags[(*ntr)*ttstr] = id ;
    (*ntr) ++ ;
  }
  
  return ;
}

static gint distribution_mesh(agg_distribution_t *d, gint id,
			      gdouble smin, gdouble smax, gint ns,
			      agg_spacing_t ss,
			      gdouble tmin, gdouble tmax, gint nt,
			      agg_spacing_t st,
			      agg_shape_t *sh, agg_local_transform_t *T,
			      agg_parser_t *p,
			      gdouble *x, gint xstr, gint *ptags, gint ptstr,
			      gint *tr, gint tstr, gint *ttags, gint ttstr,
			      gint *idx0)
  
{
  gint idx, i, j, k, ntr, *axes, i0, i1 ;
  gdouble s, t, y[3] ;
  gint nb, breaks[16] ;
  
  idx = 0 ; axes = d->axes ;
  for ( i = 0 ; i < ns ; i ++ ) {
    s = agg_spacing_eval(smin, smax, ns, ss, i) ;
    agg_distribution_interpolate_shape(d, s, sh) ;
    
    p->values[AGG_PARSER_PARAMETER_RESERVED_S] = s ;
    agg_parser_expressions_evaluate(p) ;
    agg_local_transform_eval_parameters(d->t) ;

    nb = agg_shape_break_number(sh) ;
    /*initialize list of break points on cross-section/shape*/
    breaks[0] = 0 ;
    for ( j = 0 ; j < nb ; j ++ ) {
      tmin = agg_shape_break_lower(sh, j) ;
      tmax = agg_shape_break_upper(sh, j) ;
      for ( k = 0 ; k < nt ; k ++ ) {
	t = agg_spacing_eval(tmin, tmax, nt, st, k) ;
	y[0] = fabs(t) ;
	y[1] = agg_shape_eval(sh, t, j) ;
	y[2] = 0.0 ;
	agg_local_transform_apply(d->t, y) ;
	x[idx*xstr+0] = y[SIGN(axes[0])*axes[0]-1]*SIGN(axes[0]) ;
	x[idx*xstr+1] = y[SIGN(axes[1])*axes[1]-1]*SIGN(axes[1]) ;
	x[idx*xstr+2] = y[SIGN(axes[2])*axes[2]-1]*SIGN(axes[2]) ;
	x[idx*xstr+3] = s ;
	x[idx*xstr+4] = t ;
	ptags[idx*ptstr] = id ;
	idx ++ ;
      }
      breaks[j+1] = breaks[j]+nt ;
    }
  }

  ntr = 0 ;

  if ( agg_distribution_invert(d) ) {
    i0 = 2 ; i1 = 1 ;
  } else {
    i0 = 1 ; i1 = 2 ;
  }
  for ( i = 0 ; i < ns-1 ; i ++ ) {
    gint off = i*breaks[nb] + *idx0 ;
    for ( j = 0 ; j < nb ; j ++ ) {
      for ( k = breaks[j] ; k < breaks[j+1]-1 ; k ++ ) {
	add_elements(tr, tstr, ttags, ttstr, id, x, xstr,
		     i0, i1, off+k, breaks[nb], &ntr) ;
      }
    }
  }

  *idx0 += idx ;
  
  return ntr ;
}

gint agg_distribution_mesh(agg_distribution_t *d,
			   gdouble smin, gdouble smax, gint ns,
			   agg_spacing_t ss,
			   gdouble tmin, gdouble tmax, gint nt,
			   agg_spacing_t st,
			   agg_shape_t *sh, agg_local_transform_t *T,
			   agg_parser_t *p,
			   agg_mesh_t *m)

{
  gint xstr, tstr, np, ntri ;

  xstr = 3 + 2 + m->ndat ;
  tstr = 3 + 1 + m->nttags ;
  np = agg_mesh_point_number(m) ;
  ntri = agg_mesh_triangle_number(m) ;

  agg_mesh_distribution(m, agg_mesh_distribution_number(m)) = d ;
  ntri = distribution_mesh(d, agg_mesh_distribution_number(m),
			   smin, smax, ns, ss,
			   tmin, tmax, nt, st,
			   sh, T, p,
			   &(agg_mesh_point_x(m,np)), xstr,
			   &(agg_mesh_point_distribution(m,np)),
			   agg_mesh_point_tag_number(m)+1,
			   agg_mesh_triangle(m,ntri), tstr,
			   &(agg_mesh_triangle_distribution(m,ntri)),
			   agg_mesh_triangle_tag_number(m)+3+1,
			   &np) ;
  m->np = np ;
  m->nt += ntri ;
  agg_mesh_distribution_number(m) ++ ;
  
  return 0 ;
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

gint agg_body_mesh_spherical(agg_body_t *b,
			     gdouble smin, gdouble smax, gint ns,
			     agg_spacing_t ss,
			     gdouble tmin, gdouble tmax, gint nt,
			     agg_spacing_t st,
			     agg_shape_t *sh, agg_local_transform_t *T,
			     agg_parser_t *p,
			     agg_mesh_t *m)

{
  agg_distribution_t *d ;
  gint i, j, idx, *axes, np, ntri, xstr, tstr, *ptags, ptstr ;
  gint *tr, *ttags, ttstr ;
  gdouble *x, y[3] ;
  gint faces[] = {0,  11,  5,
		  0 ,  5,  1,
		  0 ,  1,  7,
		  0 ,  7, 10,
		  0 , 10, 11,
		  1 ,  5,  9,
		  5 , 11,  4,
		  11, 10,  2,
		  10,  7,  6,
		  7 ,  1,  8,
		  3 ,  9,  4,
		  3 ,  4,  2,
		  3 ,  2,  6,
		  3 ,  6,  8,
		  3 ,  8,  9,
		  4 ,  9,  5,
		  2 ,  4, 11,
		  6 ,  2, 10,
		  8 ,  6,  7,
		  9 ,  8,  1} ;
  gdouble uv[] = {5.0000000000000000e-01,  2.3713444394043320e-01,
		  5.0000000000000000e-01,  7.6286555605956674e-01,
		  5.0000000000000000e-01, -2.3713444394043320e-01,
		  5.0000000000000000e-01, -7.6286555605956674e-01,
		  9.2532540417601994e-01, -5.0000000000000000e-01,
		  9.2532540417601994e-01,  5.0000000000000000e-01,
		  7.4674595823980006e-02, -5.0000000000000000e-01,
		  7.4674595823980006e-02,  5.0000000000000000e-01,
		  2.3713444394043320e-01,  1.0000000000000000e+00,
		  7.6286555605956674e-01,  1.0000000000000000e+00,
		  2.3713444394043320e-01, -1.1102230246251565e-16,
		  7.6286555605956674e-01,  0.0000000000000000e+00} ;

  np = agg_mesh_point_number(m) ;
  ntri = agg_mesh_triangle_number(m) ;
  xstr = 3 + 2 + m->ndat ;
  tstr = 3 + 1 + m->nttags ;
  ptags = &(agg_mesh_point_distribution(m,np)) ;
  ptstr = agg_mesh_point_tag_number(m)+1 ;

  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) 
    agg_mesh_distribution(m, i) = agg_body_distribution(b, i) ;
			/* agg_mesh_distribution_number(m)) = d ; */
  agg_mesh_distribution_number(m) = agg_body_distribution_number(b) ;

  idx = 0 ;
  x = &(agg_mesh_point_x(m,np)) ;
  for ( i = 0 ; i < 12 ; i ++ ) {
    j = agg_body_distribution_locate_u(b, uv[2*i+0]) ;
    d = agg_body_distribution(b, j) ;
    axes = d->axes ;
    agg_distribution_interpolate_shape(d, uv[2*i+0], sh) ;

    p->values[AGG_PARSER_PARAMETER_RESERVED_S] = uv[2*i+0] ;
    agg_parser_expressions_evaluate(p) ;
    agg_local_transform_eval_parameters(d->t) ;
    
    y[0] = fabs(uv[2*i+1]) ;
    y[1] = agg_shape_eval(sh, uv[2*i+1], -1) ;
    y[2] = 0.0 ;
    agg_local_transform_apply(d->t, y) ;
    x[idx*xstr+0] = y[SIGN(axes[0])*axes[0]-1]*SIGN(axes[0]) ;
    x[idx*xstr+1] = y[SIGN(axes[1])*axes[1]-1]*SIGN(axes[1]) ;
    x[idx*xstr+2] = y[SIGN(axes[2])*axes[2]-1]*SIGN(axes[2]) ;
    x[idx*xstr+3] = uv[2*i+0] ;
    x[idx*xstr+4] = uv[2*i+1] ;
    ptags[idx*ptstr] = j ;
    /* agg_mesh_distribution_number(m) ; */
    idx ++ ;    
  }
  
  m->np += idx ;
  tr = agg_mesh_triangle(m,ntri) ;
  ttags = &(agg_mesh_triangle_distribution(m,ntri)) ;
  ttstr = agg_mesh_triangle_tag_number(m)+3+1 ;

  ntri = 0 ;
  for ( i = 0 ; i < 20 ; i ++ ) {
    tr[i*(tstr)+0] = faces[3*i+0] + np ;
    tr[i*(tstr)+1] = faces[3*i+1] + np ;
    tr[i*(tstr)+2] = faces[3*i+2] + np ;
    ttags[i*ttstr] = agg_mesh_distribution_number(m) ;
    ntri ++ ;
  }
  
  m->nt += ntri ;
  agg_mesh_distribution_number(m) ++ ;
  
  return 0 ;
}

gint agg_body_mesh_grid(agg_body_t *b, agg_grid_t *g,
			agg_shape_t *sh, agg_local_transform_t *T,
			agg_parser_t *p,
			agg_mesh_t *m)

{
  agg_distribution_t *d ;
  gint i, j, idx, *axes, np, ntri, xstr, tstr, *ptags, ptstr ;
  gint *tr, *ttags, ttstr, *tg ;
  gdouble *x, y[3], u, v ;

  np = agg_mesh_point_number(m) ;
  ntri = agg_mesh_triangle_number(m) ;
  xstr = 3 + 2 + m->ndat ;
  tstr = 3 + 1 + m->nttags ;
  ptags = &(agg_mesh_point_distribution(m,np)) ;
  ptstr = agg_mesh_point_tag_number(m)+1 ;

  for ( i = 0 ; i < agg_body_distribution_number(b) ; i ++ ) 
    agg_mesh_distribution(m, i) = agg_body_distribution(b, i) ;
  agg_mesh_distribution_number(m) = agg_body_distribution_number(b) ;

  idx = 0 ;
  x = &(agg_mesh_point_x(m,np)) ;
  for ( i = 0 ; i < agg_grid_point_number(g) ; i ++ ) {
    u = agg_grid_point_u(g,i) ;
    v = agg_grid_point_v(g,i) ;
    j = agg_body_distribution_locate_u(b, u) ;
    d = agg_body_distribution(b, j) ;
    axes = d->axes ;
    agg_distribution_interpolate_shape(d, u, sh) ;

    p->values[AGG_PARSER_PARAMETER_RESERVED_S] = u ;
    agg_parser_expressions_evaluate(p) ;
    agg_local_transform_eval_parameters(d->t) ;
    
    y[0] = fabs(v) ;
    y[1] = agg_shape_eval(sh, v, -1) ;
    y[2] = 0.0 ;
    agg_local_transform_apply(d->t, y) ;
    x[idx*xstr+0] = y[SIGN(axes[0])*axes[0]-1]*SIGN(axes[0]) ;
    x[idx*xstr+1] = y[SIGN(axes[1])*axes[1]-1]*SIGN(axes[1]) ;
    x[idx*xstr+2] = y[SIGN(axes[2])*axes[2]-1]*SIGN(axes[2]) ;
    x[idx*xstr+3] = u ;
    x[idx*xstr+4] = v ;
    ptags[idx*ptstr] = j ;
    idx ++ ;    
  }
  
  m->np += idx ;
  tr = agg_mesh_triangle(m,ntri) ;
  ttags = &(agg_mesh_triangle_distribution(m,ntri)) ;
  ttstr = agg_mesh_triangle_tag_number(m)+3+1 ;

  ntri = 0 ;
  for ( i = 0 ; i < agg_grid_triangle_number(g) ; i ++ ) {
    tg = agg_grid_triangle(g, i) ;
    tr[i*(tstr)+0] = tg[0] + np ;
    tr[i*(tstr)+1] = tg[1] + np ;
    tr[i*(tstr)+2] = tg[2] + np ;
    ttags[i*ttstr] = agg_mesh_distribution_number(m) ;
    ntri ++ ;
  }
  
  m->nt += ntri ;
  
  return 0 ;
}
