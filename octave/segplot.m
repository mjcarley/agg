function segplot(st, seg)

  ## plot a triangulation and number segments
  s = reshape(st(seg,1),size(seg,1),2) ;
  t = reshape(st(seg,2),size(seg,1),2) ;

  hold off
  plot(s', t', "k") ;
  hold on

  for i=1:size(seg,1)
    id = seg(i,:) ;
    ss = sum(st(id,1))/2 ;
    tt = sum(st(id,2))/2 ;
    text(ss+0.025, tt+0.025, num2str(i-1)) ;
  endfor
