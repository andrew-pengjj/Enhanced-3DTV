function [tenPches] = extractPches(X,rGrid,cGrid,blkSz)
%************************************************************************%
%INPUTs:
% - X : 3-order Tensor
% - rGrid,cGrid : gird for spatial dimension
% - blkSz : block size of each block patch
%
%OUPUTs:
% - tenPches : 4-order Tensor collecting all the block patches.
%
%See also, aggrePches
%E-mail: caowenf2006@163.com 
%*************************************************************************%

sizeD = size(X);

nfrm = sizeD(3);

nblk = numel(rGrid)*numel(cGrid);

tenPches = zeros(blkSz,blkSz,nfrm,nblk);

cnt =0;
for p = 1 : numel(rGrid)
    
    for q = 1 : numel(cGrid)
        
        cnt = cnt +1;
        
        ridx = rGrid(p); cidx = cGrid(q);
        
        tenPches(:,:,:,cnt) = X(ridx:ridx+blkSz-1,cidx:cidx+blkSz-1,:);
    end
end

end

        
        


