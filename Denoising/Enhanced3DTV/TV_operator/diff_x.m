function tv_x = diff_x(x,sizeD)

tenX   = reshape(x, sizeD);
dfx1     = diff(tenX, 1, 1);
dfx      = zeros(sizeD);
dfx(1:end-1,:,:) = dfx1;
dfx(end,:,:)     =  tenX(1,:,:) - tenX(end,:,:);
tv_x=dfx(:);

