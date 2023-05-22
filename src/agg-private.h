/* This file is part of AGG, a library for Aerodynamic Geometry Generation
 *
 * Copyright (C) 2021 Michael Carley
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

#ifndef __AGG_PRIVATE_H_INCLUDED__
#define __AGG_PRIVATE_H_INCLUDED__

#define AGG_TRANSFORM_PARAMETER_SCALE             0 
#define AGG_TRANSFORM_PARAMETER_XSCALE            1 
#define AGG_TRANSFORM_PARAMETER_YSCALE            2 
#define AGG_TRANSFORM_PARAMETER_SHIFT             3
#define AGG_TRANSFORM_PARAMETER_PLANE_ROTATE      6
#define AGG_TRANSFORM_PARAMETER_SHRINK            9
#define AGG_TRANSFORM_PARAMETER_      13
/* #define AGG_TRANSFORM_PARAMETER_ */
/* #define AGG_TRANSFORM_PARAMETER_ */

#define plural_character(_i)  ((_i) > 1 ? "s" : "")

#define SIGN(_x) (((_x) < 0) ? -1 : 1)

#define agg_transform_scale_factor(_t)		\
  ((_t)->p[AGG_TRANSFORM_PARAMETER_SCALE])
#define agg_transform_shift_vector(_t)		\
  (&((_t)->p[AGG_TRANSFORM_PARAMETER_SHIFT]))
#define agg_transform_plane_rotation_data(_t)	\
  (&((_t)->p[AGG_TRANSFORM_PARAMETER_PLANE_ROTATE]))
#define agg_transform_shrink_data(_t)	\
  (&((_t)->p[AGG_TRANSFORM_PARAMETER_SHRINK]))

#define agg_vector_cross(_C,_A,_B)				\
  ((_C)[0] = (_A)[1]*(_B)[2] - (_A)[2]*(_B)[1],			\
   (_C)[1] = (_A)[2]*(_B)[0] - (_A)[0]*(_B)[2],			\
   (_C)[2] = (_A)[0]*(_B)[1] - (_A)[1]*(_B)[0])
#define agg_vector_diff(_A,_B,_C)		\
  do {						\
    (_A)[0] = (_B)[0] - (_C)[0] ;		\
    (_A)[1] = (_B)[1] - (_C)[1] ;		\
    (_A)[2] = (_B)[2] - (_C)[2] ;		\
  } while (0)
#define agg_vector_length(_A)				\
  (sqrt(((_A)[0])*((_A)[0])+				\
	((_A)[1])*((_A)[1]) +				\
	((_A)[2])*((_A)[2])))

#define agg_point_copy(_fb,_i,_fe,_j)			\
  do {							\
    (_fb)[2*(_j)+0] = (_fe)[2*(_i)+0] ;			\
    (_fb)[2*(_j)+1] = (_fe)[2*(_i)+1] ;			\
  } while ( 0 ) 

#define agg_point_interp3(_fb,_i,_fe,_L0,_L1,_L2)			\
  do  {									\
  (_fb)[2*(_i)+0] = (_L0)*(_fe)[0] + (_L1)*(_fe)[2] + (_L2)*(_fe)[4] ;	\
  (_fb)[2*(_i)+1] = (_L0)*(_fe)[1] + (_L1)*(_fe)[3] + (_L2)*(_fe)[5] ;	\
} while (0)

#define agg_triangle_divide_loop30(_fe,_fl)				\
  {									\
    agg_point_copy((_fl), 0, (_fe), 0) ;				\
    agg_point_interp3((_fl), 1, (_fe), 0.5, 0.5, 0.0) ;			\
    agg_point_interp3((_fl), 2, (_fe), 0.5, 0.0, 0.5) ;			\
} while (0)

#define agg_triangle_divide_loop31(_fe,_fl)				\
  {									\
    agg_point_copy((_fl), 1, (_fe), 1) ;				\
    agg_point_interp3((_fl), 0, (_fe), 0.5, 0.5, 0.0) ;			\
    agg_point_interp3((_fl), 2, (_fe), 0.0, 0.5, 0.5) ;			\
} while (0)

#define agg_triangle_divide_loop32(_fe,_fl)				\
  {									\
    agg_point_copy((_fl), 2, (_fe), 2) ;				\
    agg_point_interp3((_fl), 0, (_fe), 0.5, 0.0, 0.5) ;			\
    agg_point_interp3((_fl), 1, (_fe), 0.0, 0.5, 0.5) ;			\
} while (0)

#define agg_triangle_divide_loop33(_fe,_fl)				\
  {									\
    agg_point_interp3((_fl), 0, (_fe), 0.5, 0.5, 0.0) ;			\
    agg_point_interp3((_fl), 1, (_fe), 0.0, 0.5, 0.5) ;			\
    agg_point_interp3((_fl), 2, (_fe), 0.5, 0.0, 0.5) ;			\
} while (0)

#define agg_triangle_divide_loop(_i,_fe,_fl)		 \
  {							 \
   switch ((_i)) {					 \
   case 0: agg_triangle_divide_loop30(_fe,_fl) ; break ; \
   case 1: agg_triangle_divide_loop31(_fe,_fl) ; break ; \
   case 2: agg_triangle_divide_loop32(_fe,_fl) ; break ; \
   case 3: agg_triangle_divide_loop33(_fe,_fl) ; break ; \
   }							 \
  } while(0)

gint agg_grid_interp_area_linear(agg_grid_t *g, gdouble *uv,
				 gdouble s, gdouble t,
				 gdouble *u, gdouble *v) ;
gint agg_grid_interp_area_spherical(agg_grid_t *g, gdouble *uv,
				    gdouble s, gdouble t,
				    gdouble *u, gdouble *v) ;
gint agg_grid_interp_line_spherical(agg_grid_t *g, gdouble *uv,
				    gdouble s,
				    gdouble *u, gdouble *v) ;
gint agg_grid_interp_area_hemispherical(agg_grid_t *g, gdouble *uv,
					gdouble s, gdouble t,
					gdouble *u, gdouble *v) ;
gint agg_grid_interp_area_tube(agg_grid_t *g, gdouble *uv,
			       gdouble s, gdouble t,
			       gdouble *u, gdouble *v) ;

int coplanar_tri_tri(double N[3],double V0[3],double V1[3],double V2[3],
                     double U0[3],double U1[3],double U2[3]) ;
int tri_tri_intersect(double V0[3],double V1[3],double V2[3],
                      double U0[3],double U1[3],double U2[3]) ;
int tri_tri_intersect_with_isectline(double V0[3],double V1[3],double V2[3],
				     double U0[3],double U1[3],double U2[3],
				     int *coplanar, double isectpt1[3],
				     double isectpt2[3]) ;
int NoDivTriTriIsect(double V0[3],double V1[3],double V2[3],
                     double U0[3],double U1[3],double U2[3]) ;

#endif /*__AGG_PRIVATE_H_INCLUDED__*/
