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

/**
 *
 * @{ 
 *
 * @ingroup functions
 *
 */

/** 
 * Smoothing function intended to generate rounded wingtips
 * 
 * Evaluate the function \f$f(t)=1\f$, \f$t<t_{0}\f$,
 * \f$f(t)=p(t-t_{0})(1-t)^{\eta}/(1-t_{0})^{\eta}\f$, \f$t\geq
 * t_{0}\f$, where \f$p(t-t_{0})\f$ is a polynomial selected to make
 * \f$f(t)\f$ \f$C^{2}\f$ continuous.
 * 
 * @param eta \f$\eta\f$ tip radius exponent;
 * @param t0 \f$t_{0}\f$ limit of constant region;
 * @param t \f$t\f$ input variable.
 * 
 * @return \f$f(t)\f$ on success.
 */

gdouble agg_function_tipright(gdouble eta, gdouble t0, gdouble t)

{
  gdouble a, b, c, f, dt ;
  
  if ( t < t0 ) return 1.0 ;

  a = 1.0 ;
  b = eta/(1-t0) ;
  c = 0.5*eta*(eta+1)/(1-t0)/(1-t0) ;

  dt = t - t0 ;
  f = pow((1-t)/(1-t0), eta)*(a + dt*(b+c*dt)) ;
  
  return f ;
}

/** 
 * Smoothing function intended to generate rounded wingtips
 * 
 * Evaluate the function \f$f(t)=1\f$, \f$t>t_{0}\f$,
 * \f$f(t)=p(t_{0}-t)t^{\eta}/t_{0}^{\eta}\f$, \f$t\leq
 * t_{0}\f$, where \f$p(t_{0}-t)\f$ is a polynomial selected to make
 * \f$f(t)\f$ \f$C^{2}\f$ continuous.
 * 
 * @param eta \f$\eta\f$ tip radius exponent;
 * @param t0 \f$t_{0}\f$ limit of constant region;
 * @param t \f$t\f$ input variable.
 * 
 * @return \f$f(t)\f$ on success.
 */

gdouble agg_function_tipleft(gdouble eta, gdouble t0, gdouble t)

{
  gdouble a, b, c, f, dt ;
  
  if ( t > t0 ) return 1.0 ;

  b = eta/t0 ;
  a = 1.0 ;
  c = -0.5*eta*(eta+1)/t0/t0 ;

  dt = t0 - t ;
  f = pow(t/t0, eta)*(a + dt*(b+c*dt)) ;
  
  return f ;
}

/* 
 * @}
 */
