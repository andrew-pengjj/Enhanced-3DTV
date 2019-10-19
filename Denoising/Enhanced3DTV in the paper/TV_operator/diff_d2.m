function diff_x = diff_d2(tv,sizeD)
%% tv's size is 2*prod(sizeD),containing [dx(:);dy(:)];
%  dx(:)is total variation along the x_axis;
%  dy(:)is total variation along the y_axis;
%  sizeD is size of original image;
n    = prod(sizeD);
DX   = reshape(tv(1:n), sizeD);
DY   = reshape(tv(n+1:2*n),sizeD);
dfxx1     = diff(DX, 1, 1);% construct the dfxx;
dfxx      = zeros(sizeD);
dfxx(1:end-1,:,:) = dfxx1;
dfxx(end,:,:)     =  DX(1,:,:) - DX(end,:,:);

dfyx1     = diff(DY, 1, 1);% construct the dfyx
dfyx      = zeros(sizeD);
dfyx(1:end-1,:,:) = dfyx1;
dfyx(end,:,:)     =  DY(1,:,:) - DY(end,:,:);

dfyy1     = diff(DY, 1, 2);% construct the dfyy
dfyy      = zeros(sizeD);
dfyy(:,1:end-1,:) = dfyy1;
dfyy(end,:,:)     =  DY(:,1,:) - DY(:,end,:);

diff_x = [dfxx(:); dfyx(:);dfyy(:)];

end

