%% «ÛΩ‚Tucker decomposition 
function [Core, U] = tensorSVD(D,rank)
% performing Higth order SVD
dim        = length(size(D));
sizeD      = size(D);
U          = cell(dim,1);
Core       = D;
if nargin<2
    rank = sizeD;
end
for i = 1:dim
    UnfoD = Unfold( double(D), sizeD, i);
    if size(UnfoD,1)<size(UnfoD,2)
        DxD   = UnfoD*UnfoD';
        [U{i} , ~, ~] = svd(DxD);
        U{i}          = U{i}(:,1:rank(i));
        Core          = my_ttm(Core,{U{1:i-1},U{i}',U{i+1:dim}},i,sizeD,rank,dim);
        sizeD(i)      = rank(i);
    else
        [U , ~, ~]    = svd(UnfoD,'econ');
        U{i}          = U(:,1:rank(i));
        Core          = my_ttm(Core,{U{1:i-1},U{i}',U{i+1:dim}},i,sizeD,rank,dim);
        sizeD(i)      = rank(i);
    end
end

Core = double(Core);
end

