function f=fteval(C, th)

  f = C(1)*ones(size(th)) ;
  N = (length(C)-1)/2 ;

  for i=1:N
    f += C(2*i)*cos(i*th) + C(2*i+1)*sin(i*th) ;
  endfor
  
