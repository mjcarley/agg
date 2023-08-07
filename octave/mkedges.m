## edges of an icosahedron
## faces (unwrapped)
tri = [3    2    1
       2    3    4
       6    5    4
       15    9    4
       8   16    1
       7   10    1
       12   11    5
       11   12    7
       10    6    3
       6   10   12
       9    8    2
       13   14   11
       3    6    4
       9    2    4
       10    3    1
       2    8    1
       12   10    7
       13   11    7
       6   12    5
       11   14    5] ;
## wrapped
tri = [3, 2, 1,
2, 3, 4,
6, 5, 4,
5, 9, 4,
8, 7, 1,
7, 10, 1,
12, 11, 5,
11, 12, 7,
10, 6, 3,
6, 10, 12,
9, 8, 2,
8, 9, 11,
3, 6, 4,
9, 2, 4,
10, 3, 1,
2, 8, 1,
12, 10, 7,
8, 11, 7,
6, 12, 5,
11, 9, 5] ;

ed = [tri(:,1:2); tri(:,2:3); tri(:,[3 1])] ;
ed = [min(ed'); max(ed')]' ;

ed = unique(ed, "rows") ;

tt = [] ;
for i=1:size(tri,1)
  tmp = [] ;
  if ( tri(i,1) > tri(i,2) )
    s = -1 ; ee = tri(i,[2 1]) ;
  else
    s =  1 ; ee = tri(i,1:2) ;
  endif
  ii = find(ed(:,1) == ee(1) & ed(:,2) == ee(2)) ;
  if ( isempty(ii) ) error("edge missing") ; endif
  tmp = [tmp s*ii] ;
  if ( tri(i,2) > tri(i,3) )
    s = -1 ; ee = tri(i,[3 2]) ;
  else
    s =  1 ; ee = tri(i,2:3) ;
  endif
  ii = find(ed(:,1) == ee(1) & ed(:,2) == ee(2)) ;
  if ( isempty(ii) ) error("edge missing") ; endif
  tmp = [tmp s*ii] ;
  if ( tri(i,3) > tri(i,1) )
    s = -1 ; ee = tri(i,[1 3]) ;
  else
    s =  1 ; ee = tri(i,[3 1]) ;
  endif
  ii = find(ed(:,1) == ee(1) & ed(:,2) == ee(2)) ;
  if ( isempty(ii) ) error("edge missing") ; endif
  tmp = [tmp s*ii] ;

  tt = [tt; tmp] ;
  
endfor

