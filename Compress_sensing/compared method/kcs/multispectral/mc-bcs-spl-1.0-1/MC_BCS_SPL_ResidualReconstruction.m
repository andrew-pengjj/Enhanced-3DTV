% 
% function reconstructed_current_image = ...
%     MC_BCS_SPL_ResidualReconstruction(motion_compensated_image, y, Phi, ...
%     num_rows, num_cols)
% 
%   This function performs SPL reconstruction of residual of y and 
%   projected motion compensated frame using a DCT sparsity basis. Once 
%   residual is reconstructed it is added back to motion compensated frame
%   that is reconstuced image. Phi gives the projection matrix. 
%   The reconstructed image, of size num_rows x num_cols, is 
%   returned as reconstructed_image. 
%
%   See MC-BCS-SPL section on:
%     S. Mun and J. E. Fowler, "Residual Reconstruction For Block-
%     Based Compressed Sensing of Video," to appear in the Data
%     Compression Conference, 2011
%
%   Originally written by SungKwang Mun, Mississippi State University
%

%
% MC-BCS-SPL: Motion Compensated Block Compressed Sensing 
%               - Smooth Projected Landweber
% Copyright (C) 2011  James E. Fowler
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

function reconstructed_current_image = ...
    MC_BCS_SPL_ResidualReconstruction(motion_compensated_image, y, Phi, ...
    num_rows, num_cols)

% projecting the compensated image onto random basis
y_mc = BCS_SPL_Encoder(motion_compensated_image, Phi);

% taking the difference to pursue the sparcity.
y_r = y - y_mc;

% reconstrucing the sparse difference error.
x_r = BCS_SPL_DCT_Decoder(y_r, Phi, num_rows, num_cols);

reconstructed_current_image = motion_compensated_image + x_r;

% change the range of the image [0,1] to [0,255] and clipping
reconstructed_current_image(reconstructed_current_image > 255) = 255;
reconstructed_current_image(reconstructed_current_image < 0 ) = 0;
