% 
% function reconstructed_image = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols)
% 
%   This function performs SPL reconstruction of y using a DCT
%   sparsity basis. Phi gives the projection matrix. The reconstructed
%   image, of size num_rows x num_cols, is returned as
%   reconstructed_image.
%
%   See:
%     S. Mun and J. E. Fowler, "Block Compressed Sensing of Images
%     Using Directional Transforms," to appear in the IEEE
%     International Conference on Image Processing, 2009
%
%   Originally written by SungKwang Mun, Mississippi State University
%   Modified: 10/11/10, norm is used instead of RMS 
%

%
% BCS-SPL: Block Compressed Sensing - Smooth Projected Landweber
% Copyright (C) 2009  James E. Fowler
% http://www.ece.mstate.edu/~fowler
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
%


function reconstructed_image = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols)

[M N] = size(Phi);
block_size = sqrt(N);

Psi = DCT2D_Matrix(block_size);

% Lambda is changed from 6 to 3.

lambda = 3;
TOL = 0.001;
D_prev = 0;
D_prev2 = 0;

num_factor = 0;
max_iterations = 400;

x = Phi' * y;

for i = 1:max_iterations
  [x D] = SPLIteration(y, x, Phi, Psi, ...
      block_size, num_rows, num_cols, lambda);
  if num_factor == 4
      break;
  end
  if (D_prev ~= 0) && (D_prev2 ~=0) && ( abs(D - D_prev)  < TOL || abs(D - D_prev2)  < TOL || abs(D + D_prev2 - 2*D_prev ) < TOL)
    lambda = lambda*0.6;
    num_factor = num_factor + 1;
  end
  D_prev2 = D_prev;
  D_prev = D;
end


reconstructed_image = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x D] = SPLIteration(y, x, Phi, Psi, ...
    block_size, num_rows, num_cols, lambda)

% num_pixels_inblock = block_size*block_size;
% num_blocks = num_cols*num_rows/num_pixels_inblock;

x = col2im(x, [block_size block_size], ...
    [num_rows num_cols], 'distinct'); 
x_hat = wiener2(x, [3, 3]);
x_hat = im2col(x_hat, [block_size block_size], 'distinct');

x_hat = x_hat + Phi' * (y - Phi * x_hat);

x_check = Psi' * x_hat;
%x_check2 = zeros(num_pixels_inblock,num_blocks);

threshold = lambda * sqrt(2 * log(num_rows * num_cols)) * ...
    (median(abs(x_check(:))) / 0.6745);
x_check (abs(x_check(:))<threshold) = 0;
%index = abs(x_check(:))>threshold;
%x_check2(index) = x_check(index);  % hard threshold
% x_check2(index) = sign(x_check(index)).*(abs(x_check(index))-threshold); % soft threshold
x_bar = Psi * x_check;

x = x_bar + Phi' * (y - Phi * x_bar);

% figure(1);imagesc(x2);colormap(gray);axis image;
D = norm(x_hat-x);


%%%%%%%%% Auxilary Functions

function a = col2im(b,block,mat,kind)

if strcmp(kind, 'distinct')
 
  % Find size of padded A.
  mpad = rem(mat(1),block(1)); if mpad>0, mpad = block(1)-mpad; end
  npad = rem(mat(2),block(2)); if npad>0, npad = block(2)-npad; end
  mpad = mat(1)+mpad; npad = mat(2)+npad;
  
  mblocks = mpad/block(1);
  nblocks = npad/block(2);
  aa = mkconstarray(class(b), 0, [mpad npad]);
  x = mkconstarray(class(b), 0, block);
  rows = 1:block(1); cols = 1:block(2);
  for i=0:mblocks-1,
    for j=0:nblocks-1,
      x(:) = b(:,i+j*mblocks+1); 
      aa(i*block(1)+rows,j*block(2)+cols) = x;
    end
  end
  a = aa(1:mat(1),1:mat(2));
else
  error('Images:col2im:unknownBlockType',[deblank(kind),' is a unknown block type']);

end

function b=im2col(a, block, kind)

if strcmp(kind, 'distinct')
    % Pad A if size(A) is not divisible by block.
    [m,n] = size(a);
    mpad = rem(m,block(1)); if mpad>0, mpad = block(1)-mpad; end
    npad = rem(n,block(2)); if npad>0, npad = block(2)-npad; end
    aa = mkconstarray(class(a), 0, [m+mpad n+npad]);
    aa(1:m,1:n) = a;
    
    [m,n] = size(aa);
    mblocks = m/block(1);
    nblocks = n/block(2);
    
    b = mkconstarray(class(a), 0, [prod(block) mblocks*nblocks]);
    x = mkconstarray(class(a), 0, [prod(block) 1]);
    rows = 1:block(1); cols = 1:block(2);
    for i=0:mblocks-1,
        for j=0:nblocks-1,
            x(:) = aa(i*block(1)+rows,j*block(2)+cols);
            b(:,i+j*mblocks+1) = x;
        end
    end
else
    eid = sprintf('Images:%s:internalErrorUnknownBlockType', mfilename);
    msg = sprintf('%s is an unknown block type', kind);
    error(eid, msg);
end

function out = mkconstarray(class, value, size)
out = repmat(feval(class, value), size);


