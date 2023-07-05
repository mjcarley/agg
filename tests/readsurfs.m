nu = 65 ; nv = 129 ;

dat = load("surface1.dat") ;
x1 = reshape(dat(:,1), nv, nu) ;
y1 = reshape(dat(:,2), nv, nu) ;
z1 = reshape(dat(:,3), nv, nu) ;

dat = load("surface2.dat") ;
x2 = reshape(dat(:,1), nv, nu) ;
y2 = reshape(dat(:,2), nv, nu) ;
z2 = reshape(dat(:,3), nv, nu) ;
