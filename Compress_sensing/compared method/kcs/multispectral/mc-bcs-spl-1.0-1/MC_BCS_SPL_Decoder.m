% 
% function [reconstructed_current_frame recon_indep] = ...
%     MC_BCS_SPL_Decoder(reference_frame, y, Phi, num_rows, num_cols)
% 
%   This function performs MC-BCS-SPL method to reconstruct current frame
%   based on reference frame. Multihypothesis initialization of independent 
%   reconstruction and residual reconstruction without motino is utilized 
%   for the initial reconstruction of current frame to obtain motion 
%   vectors. After that, residual reconstruction is performed which enables
%   better reconstruction than direct reconstruction because it is sparser.
%   This function returns reconstructed frame by MC-BCS-SPL
%   and independent reconstruction by BCS-SPL. 
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

function [reconstructed_current_frame recon_indep] = ...
    MC_BCS_SPL_Decoder(reference_frame, y, Phi, num_rows, num_cols)

[M N] = size(Phi);
block_size = sqrt(N);

% intial reconstruction
recon_indep = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols);
recon_without_motion = ...
    MC_BCS_SPL_ResidualReconstruction(reference_frame, y, Phi, ...
    num_rows, num_cols);
reconstructed_current_frame = (recon_indep+recon_without_motion)/2;

% iterative motion compensation
accuracy = 'qp';
num_iter = 5;
for i = 1:num_iter    
    motion_compensated_frame = MC_BCS_SPL_MotionCompensation( ...
        reconstructed_current_frame, reference_frame, block_size,accuracy);  
    
    reconstructed_current_frame = MC_BCS_SPL_ResidualReconstruction( ...
        motion_compensated_frame, y, Phi, num_rows, num_cols);
end
