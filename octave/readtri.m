function [t,x]=readtri(file)

  fid = fopen(file, "r") ;

  dat = fscanf(fid, "%d", 3) ;

  np = dat(1) ; nt = dat(2) ; ps = dat(3) ;

  x = fscanf(fid, "%f", np*ps) ;
  x = reshape(x, ps, np)' ;

  t = fscanf(fid, "%d", nt*3) ;
  t = reshape(t, 3, nt)' ;
  
  fclose(fid) ;
