function [x0,y0,C]=fitellipse(x, y)


  x0 = mean(x) ; y0 = mean(y) ;

  r = sqrt((x-x0).^2 + (y-y0).^2) ;
  t = atan2(y-y0, x-x0) ;
  t = t(:) ; r = r(:) ;
  
  A = [cos(t).^2 -2*sin(t).*cos(t) sin(t).^2 1./r.^2] ;

  C = A\ones(size(t)) ;

  
