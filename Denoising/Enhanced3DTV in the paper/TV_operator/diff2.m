function diff_x = diff2(x,sizeD)

tenX   = reshape(x, sizeD);

dfx1     = diff(tenX, 1, 1);
dfy1     = diff(tenX, 1, 2);

dfx      = zeros(sizeD);
dfy      = zeros(sizeD);
dfx(1:end-1,:,:) = dfx1;
dfx(end,:,:)     =  tenX(1,:,:) - tenX(end,:,:);
dfy(:,1:end-1,:) = dfy1;
dfy(:,end,:)     = tenX(:,1,:) - tenX(:,end,:);

diff_x = [dfx(:); dfy(:)];

end

