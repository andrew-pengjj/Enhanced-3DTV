function diffT_a = diff_dT2(a, sizeD)

n    = prod(sizeD);
Dxx = reshape(a(1: n), sizeD);
Dyx = reshape(a((n+1):2*n), sizeD);
Dyy = reshape(a((2*n+1):3*n), sizeD);

dfxx     = diff(Dxx, 1, 1);
dfyx     = diff(Dyx, 1, 1);
dfyy     = diff(Dyy, 1, 2);

dfxxT   = zeros(sizeD);
dfyxT   = zeros(sizeD);
dfyyT   = zeros(sizeD);
dfxxT(1,:,:) = Dxx(end, :, :) - Dxx(1, :, :); %
dfxxT(2:end,:,:) = -dfxx;
dfyxT(1,:,:) = Dyx(end, :, :) - Dyx(1, :, :); %
dfyxT(2:end,:,:) = -dfyx;
dfyyT(:,1,:)     =  Dyy(:,end,:) - Dyy(:,1,:);
dfyyT(:,2:end,:) = -dfyy;


diffT_a1 = dfxxT + dfyxT;
diffT_a2 = dfyxT + dfyyT;
diffT_a = [diffT_a1(:);diffT_a2(:)];
end