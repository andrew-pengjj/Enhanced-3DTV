function img = JointBlocks( blocks, params )

%==========================================================================
% Reconstruct the multi-spectral imagery (MSI) through blocks.
%
% Syntax:
%   img = JointBlocks( blocks, params );
%
% Input arguments:
%   blocks......a stack of 3D (non-decoupled) or 2D (decoupled) array that
%               indicats the blocks.(required)
% 	params......an option structure whose fields are listed as follows:
%       overlap_sz: a 1 x 2 integer vector that indicates the size of
%               overlaps. If overlap_sz is a scalar, overlap size of all
%               dimensions will be the same. (required)
%       block_num: a 1 x 2 integer vector that indicates how many blocks will
%               be extracted. (required)
%       decouple: a logical variable. If it is true, then each spatial
%               slice of blocks is vectorized.(default false)
%       block_sz: a 1 x 2 integer vector that indicates the SPATIAL size of
%               a block. A scalar value of this field is acceptable.
%               (required only for 'decouple' mode)
%       
% REMARK: prod(block_num) must EQUAL TO size(blocks, end). This argument
%   can be obtained from the outputs of function ExtractBlocks.
%
% Output argument:
% 	img.........the reconstructed MSI tensor. Its first and second mode are
%               y axis and x axis respectively.
%
% See also ExtractBlocks
%
% by Yi Peng
%==========================================================================

if isfield(params, 'decouple')
    if islogical(params.decouple) || params.decouple == 0 || params.decouple == 1
        decouple = (params.decouple == true);
    else
        error('Illegal value of logical parameter ''decouple''!');
    end
else
    decouple = false;
end

if decouple
    [L, S, N] = size(blocks);
    if ~isfield(params, 'block_sz')
        error('Spatial size of blocks isrequired for ''decouple'' mode!');
    elseif isscalar(params.block_sz)
        block_sz = repmat(params.block_sz, 1, 2);
    else
        block_sz = [params.block_sz(1), params.block_sz(2)];
    end
    if prod(block_sz) ~= L
        error('Size of blocks doesn''t match the parameter!');
    else
        blocks = reshape(blocks, [block_sz, S, N]);
    end
else
    [H, W, S, N] = size(blocks);
    block_sz = [H, W];
end

if ~isfield(params, 'block_num')
    error('Number of blocks must be specified!\n');
elseif isscalar(params.block_num)
    block_num = repmat(params.block_num, 1, 2);
else
    block_num = [params.block_num(1), params.block_num(2)];
end
if ~(prod(block_num) == N)
    error('Block number does not match!\n');
end

if ~isfield(params, 'overlap_sz')
    error('Size of overlap must be specified!\n');
elseif isscalar(params.block_num)
    overlap_sz = repmat(params.overlap_sz, 1, 2);
else
    overlap_sz = [params.overlap_sz(1), params.overlap_sz(2)];
end

img_sz = [block_num.*(block_sz - overlap_sz) + overlap_sz, S];
mult = zeros(img_sz);
img = zeros(img_sz);
for i = 1:block_num(1)
    for j = 1:block_num(2)
        ii = 1 + (i - 1)*(block_sz(1) - overlap_sz(1));
        jj = 1 + (j - 1)*(block_sz(2) - overlap_sz(2));
        idx = (j-1)*block_num(1) + i;
        mult(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) = ...
            mult(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) + 1;
        img(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) = ...
            img(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :) + blocks(:, :, :, idx);
    end
end
img = img./mult;