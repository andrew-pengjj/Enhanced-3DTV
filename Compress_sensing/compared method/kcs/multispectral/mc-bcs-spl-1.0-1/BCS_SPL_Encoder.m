% 
% function y = BCS_SPL_Encoder(current_image, Phi)
% 
%   This function performs BCS projections of each block of
%   current_image. The number of columns of the projection matrix,
%   Phi, determines the size of the blocks into which current_image
%   is partitioned. The projections are returned as the columns of y.
%
%   See:
%     S. Mun and J. E. Fowler, "Block Compressed Sensing of Images
%     Using Directional Transforms," to appear in the IEEE
%     International Conference on Image Processing, 2009
%
%   Originally written by SungKwang Mun, Mississippi State University
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


function y = BCS_SPL_Encoder(current_image, Phi)

[M N] = size(Phi);

block_size = sqrt(N);

[num_rows num_cols] = size(current_image);

x = im2col(current_image, [block_size block_size], 'distinct');

y = Phi * x;
