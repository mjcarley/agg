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

#ifndef __AGG_PRIVATE_H_INCLUDED__
#define __AGG_PRIVATE_H_INCLUDED__


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

#define agg_vector_scalar(_A,_B)				\
  (((_A)[0])*((_B)[0])+((_A)[1])*((_B)[1])+((_A)[2])*((_B)[2]))

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
#define agg_vector_distance2(_A,_B)		\
  ( ((_A)[0]-(_B)[0])*((_A)[0]-(_B)[0]) +	\
    ((_A)[1]-(_B)[1])*((_A)[1]-(_B)[1]) +	\
    ((_A)[2]-(_B)[2])*((_A)[2]-(_B)[2]) )

#define agg_vector_distance(_A,_B)	\
  (sqrt((agg_vector_distance2(_A,_B))))

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

/*
 * inverse of general 2x2 matrix
 *
 * inv [a b; c d] = [d -b; -c a]/det
 *
 * det = a*d - b*c
 *
 * can be done in-place (Ai == A)
 */

#define agg_invert2x2(_Ai,_A,_det)				\
  do {								\
  gdouble _a, _b, _c, _d ;					\
  _a = (_A)[0] ; _b = (_A)[1] ; _c = (_A)[2] ; _d = (_A)[3] ;	\
  (*(_det)) = _a*_d - _b*_c ;					\
  (_Ai)[0] = _d/(*(_det)) ; (_Ai)[1] = -(_b)/(*(_det)) ;	\
  (_Ai)[2] = -_c/(*(_det)) ; (_Ai)[3] = _a/(*(_det)) ;		\
  } while (0)

#define agg_matrix_identity_4x4(_A)		\
  do {						\
    (_A)[ 0] = 1.0 ; (_A)[ 1] = 0.0 ; (_A)[ 2] = 0.0 ; (_A)[ 3] = 0.0 ;	\
    (_A)[ 4] = 0.0 ; (_A)[ 5] = 1.0 ; (_A)[ 6] = 0.0 ; (_A)[ 7] = 0.0 ;	\
    (_A)[ 8] = 0.0 ; (_A)[ 9] = 0.0 ; (_A)[10] = 1.0 ; (_A)[11] = 0.0 ;	\
    (_A)[12] = 0.0 ; (_A)[13] = 0.0 ; (_A)[14] = 0.0 ; (_A)[15] = 1.0 ;	\
  } while (0)

#define agg_mat_vec_mul_4(_y,_A,_x)		\
  do {						\
    gint _i, _j ;				\
    for ( _i = 0 ; _i < 4 ; _i ++ ) {		\
      (_y)[_i] = 0.0 ;				\
      for ( _j = 0 ; _j < 4 ; _j ++ ) {		\
	(_y)[_i] += (_A)[_i*4+_j]*(_x)[_j] ;	\
      }						\
    }						\
  } while (0)


#define agg_mat_mat_mul_4(_C,_A,_B)				\
  do {								\
    gint _i, _j, _k ;						\
    gdouble _s[16] ;						\
    for ( _i = 0 ; _i < 4 ; _i ++ ) {				\
      for ( _j = 0 ; _j < 4 ; _j ++ ) {				\
	_s[4*_i+_j] = 0.0 ;					\
	for ( _k = 0 ; _k < 4 ; _k ++ ) {			\
	  _s[4*_i+_j] += (_A)[4*_i+_k]*(_B)[4*_k+_j] ;		\
	}							\
      }								\
    }								\
    for ( _i = 0 ; _i < 16 ; _i ++ ) (_C)[_i] = (_s)[_i] ;	\
  } while (0)

#endif /*__AGG_PRIVATE_H_INCLUDED__*/
