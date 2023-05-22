u0 = mean(uv) ;

r = sqrt((uv(:,1)-u0(1)).^2 + (uv(:,2)-u0(2)).^2) ;
th = atan2(uv(:,2)-u0(2), uv(:,1)-u0(1)) ;

[th,i] = sort(th) ;

uv = uv(i,:) ;


