function f=lround(eta, t0, t)

  f = ones(size(t)) ;
  ii = find(t < t0) ;

  b = eta/t0 ;
  a = 1.0 ;
  c = 0.5*eta*(eta+1)/t0^2 ;

  dt = t0 - t(ii) ;
  f(ii) = (t(ii)/t0).^eta.*(a + dt.*(b+c*dt)) ;
