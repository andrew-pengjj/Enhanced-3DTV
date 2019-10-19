function [aveX,weig,noAveX] = aggrePches(tenPches,rGrid,cGrid,blkSz,sizeD)
%************************************************************************
%This script can be used to arregate the block patches into a tensor.
%-----------------------------------------------------------------------
%INPUTs:
% - TenPches : the tensor patches
% - rGird : grid for row
% - cGrid : grid for col
% - blkSz : the size of block
% - sizeD : the size of tensor to be constr.
%
%OUTPUTs:
% - aveX : averaged Tensor
% - weig : weights
% - noAveX : not averaged Tensor
%
%see also, extractPches
%***********************************************************************

weig = zeros(sizeD);

noAveX = zeros(sizeD);

cnt =0;
for p = 1 : numel(rGrid)
    
     for q = 1 : numel(cGrid)
        
        cnt = cnt +1;
        
        ridx = rGrid(p); cidx = cGrid(q);
        
        noAveX(ridx : ridx + blkSz-1, cidx : cidx + blkSz -1,:) =...
            noAveX(ridx : ridx + blkSz-1, cidx : cidx + blkSz -1,:) + tenPches(:,:,:,cnt);
        
        weig(ridx : ridx + blkSz -1, cidx : cidx + blkSz -1,:) = ...
                 weig(ridx : ridx + blkSz -1, cidx : cidx + blkSz -1,:) + 1;
    end
    
end

aveX = noAveX./weig;

end


        

