function [x, t, u, v, y]=hemisphere()

## http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
##
## modified to return upper half (z >= 0) sphere
  
  t = [0    2   5;
       0    5   1;
       0    1   7;
       0    7   3;
       0    3   2;
       1    5   9;
       5    2   4;
       3    7   6;
       7    1   8;
       4    9   5;
       8    6   7;
       9    8   1] ;

  rt = 0.5*(1.0 + sqrt(5.0)) ;
  R = sqrt(1+rt.^2) ;

  x( 1, 2) = -1.0 ; x( 1, 3) =  rt  ; x( 1, 1) =  0.0 ; 
  x( 2, 2) =  1.0 ; x( 2, 3) =  rt  ; x( 2, 1) =  0.0 ; 
  x( 3, 2) = -rt  ; x(3, 3) =  0.0 ;  x(3, 1) =   1.0 ;
  x( 4, 2) = -rt  ; x(4, 3) =  0.0 ;  x(4, 1) =  -1.0 ;

  x( 5, 2) =  0.0 ; x( 5, 3) =  0.0 ; x( 5, 1) =   R ;
  x( 6, 2) =  0.0 ; x( 6, 3) =  1.0 ; x( 6, 1) =   rt ; 

  x( 7, 2) =  0.0 ; x( 7, 3) =  0.0 ; x( 7, 1) =  -R ;
  x( 8, 2) =  0.0 ; x( 8, 3) =  1.0 ; x( 8, 1) =  -rt ; 
		        		  
  x( 9, 2) =  rt  ; x( 9, 3) =  0.0 ; x( 9, 1) =  -1.0 ; 
  x(10, 2) =  rt  ; x(10, 3) =  0.0 ; x(10, 1) =   1.0 ; 

  x /= 2*R ;

  u = x(:,3)*2 ;

  r = sqrt(1-u.^2) ;
  v = x(:,1)./r + 0.5 ;
  ii = find(x(:,2)<0) ;
  v(ii) = -v(ii) ;

  v(8) = 0 ;
  y = [(abs(v)-0.5).*r+0.0 sign(v).*sqrt(abs(v)).*sqrt(1-abs(v)).*r u/2] ;
  

  return
  
  x(:,1) += 0.5 ;
  
  u = x(:,3)*2 ;

  r = sqrt(1-u.^2) ;

  v = (x(:,1)-0.5)./r + 0.5 ;
  ii = find(x(:,2)<0) ;
  v(ii) = -v(ii) ;

  y = [(abs(v)-0.5).*r+0.5 sign(v).*sqrt(abs(v)).*sqrt(1-abs(v)).*r u/2] ;
