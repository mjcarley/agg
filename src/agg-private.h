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

#endif /*__AGG_PRIVATE_H_INCLUDED__*/
