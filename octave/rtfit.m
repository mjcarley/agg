n1 = n2 = 0.5 ;

h = -0.9 ;

s = 0.7 ;
b = -0.1 ;

x = linspace(0, 1, 65) ;

y = x.^n1.*(1-x).^n2*h ;
dy = y*h./x./(1-x).*(n1 - (n1+n2)*x) ;

#s = 1.0 - s ;

ys = s.^n1.*(1-s).^n2*h ;
dys = ys./s./(1-s).*(n1 - (n1+n2)*s) ;

A = [1 s s^2 s^3 ;
     0 1 2*s 3*s^2 ;
#     0 1 1 1 ;
#     0 1 2 3] ;
          1 1 1 1 ;
     0 1 2 3] ;

r = [ys dys b 0]' ;

q = A\r ;

q = flipud(q) ;

xs = linspace(s-0.05, 1, 65) ;
##xs = linspace(0, s+0.05, 65) ;
ys = polyval(q, xs) ;
