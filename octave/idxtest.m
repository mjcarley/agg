np = 15 ;

p = 0:np-1 ;

idx = [] ;

for i=1:length(p)
  for j=1:length(p)
    if ( i ~= j )
      ii = max([p(i) p(j)])*np + min([p(i) p(j)]) ;
      idx = [idx; ii] ;
    endif
  endfor
endfor


