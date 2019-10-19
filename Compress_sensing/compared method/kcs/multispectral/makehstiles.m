function xmat = makehstiles(x)

% This function creates a tiling of a 3-D cube of data for display purposes. 
% Cuts along the first dimension are concatenated into a matrix, which is 
% displayed as an image.
%
% Written by: Marco F. Duarte, Rice University
% Created: December 2006

N = size(x,2);
S = size(x,1);
Sf = cumprod(factor(S));
Sf = Sf(find(Sf <= sqrt(S)));
Sf = [S/Sf(end) Sf(end)];

xmat = zeros(Sf*N);
for i=1:Sf(1),
    for j=1:Sf(2),
        xmat(((i-1)*N)+(1:N),((j-1)*N)+(1:N)) = squeeze(x((i-1)*Sf(2)+j,:,:));
    end
end