function x=uvsphere(u, v)


  u = u(:) ; v = v(:) ;
  x = zeros(length(u), 3) ;

  x(:,3) = u ;
  x(:,1) = sqrt(u).*sqrt(1-u).*(2*abs(v)-1) ;
  ii = find(v >= 0) ;
  x(ii,2) = sqrt(u(ii)).*sqrt(1-u(ii)).*sqrt(abs(v(ii))).*sqrt(1-abs(v(ii)))*2 ;
  ii = find(v < 0) ;
  x(ii,2) = -sqrt(u(ii)).*sqrt(1-u(ii)).*sqrt(abs(v(ii))).*sqrt(1-abs(v(ii)))*2 ;

  return
  
  x(:,3) = u ;
  ##x(:,1) = 2*(abs(v)-1/2).*sqrt(u).*sqrt(1-u) ;
  x(:,1) = 2*(abs(v)-0).*sqrt(u).*sqrt(1-u) ;
  [sign(v) sqrt(abs(v)) sqrt(1-abs(v)) sqrt(u) sqrt(1-u)]
  x(:,2) = 2*sign(v).*sqrt(abs(v)).*sqrt(1-abs(v)).*sqrt(u).*sqrt(1-u) ;
