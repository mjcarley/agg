t0 = 0.9 ;
n = 0.5 ;

t = linspace(0.8, 1, 65) ;
##t = sort([t0 t]) ;

a = 1 ;
b = n/(1-t0) ;
c = 0.5*n*(n+1)/(1-t0)^2 ;

p = [c b a] ;

f = ones(size(t)) ;
ii = find(t >= t0) ;

f(ii) = ((1-t(ii))/(1-t0)).^n.*polyval(p, t(ii)-t0) ;


t1 = 0.1 ;
t = linspace(0, 0.2, 65) ;
a = 1 ;
b = n/t1 ;
c = -0.5*n*(n+1)/t1^2 ;

p1 = [c b a] ;

t = 0.02243806429580494 ;
t = 0.044864830350514931 ;

f = ones(size(t)) ;
ii = find(t <= t1) ;

f(ii) = (t(ii)/t1).^n.*polyval(p1, t1-t(ii)) ;

