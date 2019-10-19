function Y = extracting3DPatch(X,f)
% X :  h x w x d
% Y:   f^2 x L x d, where L is the numbe of all 3-D patches.
% f:   size of 3D patch

d         = size(X,3);
N         =   size(X,1)-f+1;
M         =   size(X,2)-f+1;
L         =   N*M;
Y         =   zeros(f*f,L,d, 'single'); %% f*f x L x d

k    =  0;
for i  = 1:f
    for j  = 1:f
        
        k    =  k+1;
        blk  =  X(i:end-f+i,j:end-f+j,:);
        blk  = reshape(blk,size(blk,1)*size(blk,2),size(blk,3));
        Y(k,:,:) =  blk;
        
    end
end


return;

end