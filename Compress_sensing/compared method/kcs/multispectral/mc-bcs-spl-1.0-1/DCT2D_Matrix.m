%
% function Psi = DCT2D_Matrix(N)
%
%   This function returns the N^2 x N^2 orthonormal transform matrix
%   associated with the N^2-point DCT.
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


function Psi = DCT2D_Matrix(N)

Psi = zeros(N * N, N * N);

for row = 1:N
  for col = 1:N
    X = zeros(N, N);
    X(row, col) = 1;
    x = idct2(X);
    Psi(:, (row - 1) * N + col) = x(:);
  end
end

