function tv_yT = diff_yT(a, sizeD)
n    = prod(sizeD);
tenY = reshape(a(1:n), sizeD);
dfy     = diff(tenY, 1, 2);
dfyT   = zeros(sizeD);
dfyT(:,1,:)     =  tenY(:,end,:) - tenY(:,1,:);
dfyT(:,2:end,:) = -dfy;
tv_yT=dfyT(:);