function tv_xT = diff_xT(a, sizeD)
n    = prod(sizeD);
tenX = reshape(a(1: n), sizeD);
dfx     = diff(tenX, 1, 1);
dfxT   = zeros(sizeD);
dfxT(1,:,:) = tenX(end, :, :) - tenX(1, :, :); %
dfxT(2:end,:,:) = -dfx;
tv_xT=dfxT(:);
