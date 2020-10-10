function diffT_a = diffT2(a, sizeD)

N    = prod(sizeD); 
tenX = reshape(a(1: N), sizeD);
tenY = reshape(a((N+1):2*N), sizeD);

dfx     = diff(tenX, 1, 1);
dfy     = diff(tenY, 1, 2);

dfxT   = zeros(sizeD);
dfyT   = zeros(sizeD);
dfxT(1,:,:) = tenX(end, :, :) - tenX(1, :, :); %
dfxT(2:end,:,:) = -dfx;
dfyT(:,1,:)     =  tenY(:,end,:) - tenY(:,1,:);
dfyT(:,2:end,:) = -dfy;


diffT_a = dfxT + dfyT;
diffT_a = diffT_a(:);
end