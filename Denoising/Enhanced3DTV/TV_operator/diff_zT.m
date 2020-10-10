function tv_zT = diff_zT(a, sizeD)
n    = prod(sizeD);
tenZ = reshape(a(1:n), sizeD);
dfz     = diff(tenZ, 1, 3);
dfzT   = zeros(sizeD);
dfzT(:,:,1)     =  tenZ(:,:,end) - tenZ(:,:,1);
dfzT(:,:,2:end) = -dfz;
tv_zT=dfzT(:);