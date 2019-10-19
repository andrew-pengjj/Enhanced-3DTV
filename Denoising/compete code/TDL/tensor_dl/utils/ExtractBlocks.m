function  [ blocks, reg_params ] = ExtractBlocks( img, params )

%==========================================================================
% Divide a multi-spectral imagery (MSI) into blocks, each of which is a 3rd
%   order tensor with all bands and a specific spatial domain.
%
% Syntax:
%   [ blocks, reg_params ] = ExtractBlocks( img, params );
%
% Input arguments:
%   img.........the MSI tensor.(required)
%   params......an option structure whose fields are as follows:
%       block_sz: a 1-by-2 integer vector that indicates the size of blocks. 
%               If block_sz is a scalar, size of all dimensions will be the
%               same. (required)
%       overlap_sz: a 1-by-2 integer vector that indicates the size of
%               overlaps. If overlap_sz is a scalar, overlap size of all
%               dimensions will be the same.(required)
%       decouple: a logical variable. If it is true, then each spatial
%               slice of blocks is vectorized.(default false)
%
% Output arguments:
% 	blocks......the extracted block tensor. Its first and second mode
%               are y axis and x axis respectively. Its last mode indicates
%               the various blocks.
% 	reg_params..the regularized parameters including overlap_sz and
%               block_num, which may be used in function JointBlocks.
%
% See also JointBlocks
%
% by Yi Peng
%==========================================================================

sz = size(img);
if ndims(img) ~= 3
    if ndims(img) == 2
        sz(3) = 1;
    else
        error(['Wrong imagery data format!\n' ...
            'Make sure it is a  height x width x nbands tensor.']);
    end
end

if ~exist('params', 'var')
    params = [];
end

if ~isfield(params, 'block_sz')
    error('Block size must be specified!\n');
elseif isscalar(params.block_sz)
    block_sz = repmat(params.block_sz, 1, 2);
else
    block_sz = [params.block_sz(1), params.block_sz(2)];
end
reg_params.block_sz = block_sz;

if ~isfield(params, 'overlap_sz')
    error('Overlap size must be specified!\n');
elseif isscalar(params.overlap_sz)
    overlap_sz = repmat(params.overlap_sz, 1, 2);
else
    overlap_sz = [params.overlap_sz(1), params.overlap_sz(2)];
end
reg_params.overlap_sz = overlap_sz;
block_num = floor((sz(1:2) - overlap_sz)./(block_sz - overlap_sz));
reg_params.block_num = block_num;

if isfield(params, 'decouple')
    if islogical(params.decouple) || params.decouple == 0 || params.decouple == 1
        decouple = (params.decouple == true);
        reg_params.decouple = decouple;
    else
        error('Illegal value of logical parameter ''decouple''!');
    end
else
    decouple = false;
    reg_params.decouple = 0;
end

blocks = zeros([block_sz, sz(3), prod(block_num)]);
for i = 1:block_num(1)
    for j = 1:block_num(2)
        ii = 1 + (i - 1)*(block_sz(1) - overlap_sz(1));
        jj = 1 + (j - 1)*(block_sz(2) - overlap_sz(2));
        idx = (j-1)*block_num(1) + i;
        blocks(:, :, :, idx) = ...
            img(ii:ii+block_sz(1)-1, jj:jj+block_sz(2)-1, :);
    end
end
if decouple
    blocks = reshape(blocks, [block_sz(1)*block_sz(2), sz(3), prod(block_num)]);
end

