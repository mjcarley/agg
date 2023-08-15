function C=ffit(th, f, N)

  if 1
  [th,i] = sort(th) ; f = f(i) ;
  th = [th; th(1)+2*pi] ; f = [f; f(1)] ;

  C = trapz(th,f.*cos(0*th))/2/pi ;

  for i=1:N
    C = [C; trapz(th,f.*cos(i*th))/2/pi ; trapz(th,f.*sin(i*th))/2/pi] ;
  endfor
  
  return ;
  endif
  A = zeros(length(f), 2*N+1) ;
  A(:,1) = 1 ;

  for i=1:length(f) ;
    A(i,2:2:end) = cos((1:N)*th(i)) ;
    A(i,3:2:end) = sin((1:N)*th(i)) ;
  endfor

  C = A\f ;
  
  
