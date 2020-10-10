function tv_z = diff_z(x,sizeD)

tenX   = reshape(x, sizeD);
dfz1     = diff(tenX, 1, 3);
dfz      = zeros(sizeD);
dfz(:,:,1:end-1) = dfz1;
dfz(:,:,end)     = tenX(:,:,1) - tenX(:,:,end);
tv_z=dfz(:);
