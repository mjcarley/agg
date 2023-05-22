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

#include <glib.h>

#include <blaswrap.h>

/*
 * implementation of Broyden, C. G., A class of methods for solving
 * nonlinear simultaneous equations, Mathematics of Computation,
 * 19:577--593, 196, http://www.jstor.org/stable/2003941
 */

static void broyden_inverse_jacobian(gint (*func)(gdouble *, gint n,
						  gdouble *, gpointer),
				     gint n, gdouble *x, gdouble ee,
				     gpointer data, gdouble *H,
				     gdouble *work)

/*
 * find -J^{-1}, with J estimated Jacobian matrix
 */
  
{
  gdouble *xp, *xm, *fp, *fm ;
  gint one = 1, i, j, ipiv[64], info, lwork ;

  /*initial estimate of Jacobian*/
  xp = work ; xm = &(xp[n]) ;
  fp = &(xm[n]) ; fm = &(fp[n]) ;
  for ( j = 0 ; j < n ; j ++ ) {
    blaswrap_dcopy(n, x, one, xp, one) ;
    blaswrap_dcopy(n, x, one, xm, one) ;
    xp[j] += ee*0.5 ; xm[j] -= ee*0.5 ;
    func(xp, n, fp, data) ;
    func(xm, n, fm, data) ;
    for ( i = 0 ; i < n ; i ++ ) {
      H[i*n+j] = -(fp[i] - fm[i])/ee ;
    }
  }

  lwork = n*n ;
  /*invert approximate Jacobian*/
  dgetrf_(&n, &n, H, &n, ipiv, &info) ;
  dgetri_(&n, H, &n, ipiv, work, &lwork, &info) ;  

  return ;
}

gint broyden_solve(gint (*func)(gdouble *, gint n, gdouble *, gpointer),
		   gint n, gdouble *x, gdouble tol, gint imax,
		   gpointer data, gdouble *work)

{
  gdouble *H, *iwork, *f, *fp1, ee, *p, al, bt, t, *tmp, *t1, *t2, sc, err ;
  gint i, one = 1 ;

  ee = MAX(tol, 1e-6) ;

  H = &(work[0]) ;
  iwork = &(H[n*n]) ;

  f = iwork ;
  fp1 = &(f[n]) ;
  p = &(fp1[n]) ;
  t1 = &(p[n]) ;
  t2 = &(t1[n]) ;

  broyden_inverse_jacobian(func, n, x, ee, data, H, iwork) ;
  func(x, n, f, data) ;
  /*get to iterating*/
  for ( (i = 0), (err = G_MAXDOUBLE) ; (i < imax) && (err > tol) ; i ++ ) {
    al = 1.0 ; bt = 0 ;
    blaswrap_dgemv(FALSE, n, n, al, H, n, f, one, bt, p, one) ;

    /*this could be set using Broyden's Appendix 2*/
    t = 1e-3 ;
    
    blaswrap_daxpy(n, t, p, one, x, one) ;
    if ( func(x, n, fp1, data) != 0 ) return -1 ;

    /*y := f - fp1 (stored in f)*/
    al = -1 ;
    blaswrap_daxpy(n, al, fp1, one, f, one) ;

    /*update H*/
    /*t1 = Hy + tp (noting that y = -f)*/
    al = -1.0 ;
    blaswrap_dgemv(FALSE, n, n, al, H, n, f, one, bt, t1, one) ;
    sc = 1.0/blaswrap_ddot(n, p, one, t1, one) ;
    blaswrap_daxpy(n, t, p, one, t1, one) ;

    /*t2 = p^T*H = (H^T*p)^T*/
    al = 1.0 ;
    blaswrap_dgemv(TRUE, n, n, al, H, n, p, one, bt, t2, one) ;
    /*H = H - sc*t1*t2*/
    al = -sc ; bt = 1.0 ;
    blaswrap_dgemm(FALSE, FALSE, n, n, one, al, t1, one, t2, n, bt, H, n) ;
		   
    /*swap f and fp1 so that f contains f(x)*/
    tmp = f ; f = fp1 ; fp1 = tmp ;

    err = blaswrap_dnrm2(n, f, one) ;
  }

  if ( err <= tol ) return i ;

  return -i ;
}
