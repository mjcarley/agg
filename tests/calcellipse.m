function r=calcellipse(C, t)

  f = 1 - C(1)*cos(t).^2 + C(2)*sin(t).*cos(t) - C(3)*sin(t).^2 ;
  r = sqrt(C(4)./f) ;
  
