function tv_y = diff_y(x,sizeD)

tenX   = reshape(x, sizeD);
dfy1     = diff(tenX, 1, 2);
dfy      = zeros(sizeD);
dfy(:,1:end-1,:) = dfy1;
dfy(:,end,:)     = tenX(:,1,:) - tenX(:,end,:);
tv_y=dfy(:);
