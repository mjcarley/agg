## tests of different spacing terms
n = 33 ;

i = (0:(n-1))/(n-1) ;

tmin = -1 ; tmax = 1 ;

## cosine
t = 0.5*sign(i).*(1-cos(pi*i)) ;
t = tmin + (tmax - tmin)*t ;


