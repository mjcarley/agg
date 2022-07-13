function [x, t, xt, n]=aggmesh(file)

  ## [X, T] = AGGMESH(FILE): read an AGG mesh from file FILE
  ##

  fid = fopen(file, "r") ;

  dat = fscanf(fid, "%d", 3) ;
  np    = dat(1) ;
  ndat  = dat(2) ;
  ntags = dat(3) ;
  len = 1 + 3 + 2 + ndat + 1 + ntags ;
  dat = fscanf(fid, "%f", len*np) ;

  dat = reshape(dat, len, np)' ;

  x = dat(:,2:end) ;

  dat = fscanf(fid, "%d", 2) ;
  ntri = dat(1) ;
  ntag = dat(2) ;
  len = 3 + 1 + ntag ;
  dat = fscanf(fid, "%d", len*ntri) ;

  t = reshape(dat, len, ntri)' ;
  
  fclose(fid) ;

  ## if required, find the triangle centroids and normals
  xx = reshape(x(t(:,1:3)'+1,1), 3, ntri)' ;
  yy = reshape(x(t(:,1:3)'+1,2), 3, ntri)' ;
  zz = reshape(x(t(:,1:3)'+1,3), 3, ntri)' ;
  
  xt = [mean(xx, 2) mean(yy, 2) mean(zz, 2)] ;

  n = cross([xx(:,2) - xx(:,1) yy(:,2) - yy(:,1) zz(:,2) - zz(:,1)],
	    [xx(:,3) - xx(:,1) yy(:,3) - yy(:,1) zz(:,3) - zz(:,1)]) ;

  len = sqrt(sum(n.^2, 2)) ;
  n = n./[len len len] ;
