function f=rround(eta, t0, t)

  f = ones(size(t)) ;
  ii = find(t > t0) ;

  b = eta/(1-t0) ;
  a = 1.0 ;
  c = 0.5*eta*(eta+1)/(1-t0)^2 ;

  dt = t(ii) - t0 ;
  f(ii) = ((1-t(ii))/(1-t0)).^eta.*(a + dt.*(b+c*dt)) ;
