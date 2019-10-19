% 
%function [reconstructed_center_frame recon_indep] = ...
%    MC_BCS_SPL_Decoder_Center(y, previous_frame, next_frame, Phi)
% 
%   This function performs MC-BCS-SPL method to reconstruct center frame
%   based on two reference frame in bidirecional way. Mulityhypothesis 
%   initialization includes independent reconstruction and two residual
%   reconstruction without motion. Two residual reconstruction of motion 
%   compensated frames composited by adjacent frames are averaged as final
%   reconstruction of center frame. 
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

function [reconstructed_center_frame recon_indep] = ...
    MC_BCS_SPL_Decoder_Center(y, previous_frame, next_frame, Phi)

[num_rows num_cols] = size(previous_frame);

[M N]= size(Phi);
block_size = sqrt(N);

% Hypothesis Initialization
recon_indep = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols);
recon_residual_without_motion1 = ...
    MC_BCS_SPL_ResidualReconstruction(previous_frame, y, Phi, num_rows, num_cols);
recon_residual_without_motion2 = ...
    MC_BCS_SPL_ResidualReconstruction(next_frame, y, Phi, num_rows, num_cols);
reconstructed_current_frame = ...
    (recon_indep + recon_residual_without_motion1 + recon_residual_without_motion2)/3;

% Bidirectional iterative motion compensation
num_iter = 5; 
accuracy = 'qp';
reference_frame = previous_frame;
recon_current_frame_forward = ...
    MC_BCS_SPL_IterativeReconstruction(reconstructed_current_frame, ...
    reference_frame, y, Phi, num_iter, block_size, accuracy); 

reference_frame = next_frame;
recon_current_frame_backward = ...
    MC_BCS_SPL_IterativeReconstruction(reconstructed_current_frame, ...
    reference_frame, y, Phi, num_iter, block_size, accuracy); 
% average result
reconstructed_center_frame = ...
    (recon_current_frame_forward + recon_current_frame_backward)/2;

% Auxiliary function
function reconstructed_frame = ...
    MC_BCS_SPL_IterativeReconstruction(current_frame, reference_frame, ...
    y, Phi, num_iter, block_size, accuracy) 

[num_rows num_cols] = size(current_frame);
for i = 1:num_iter
    motion_compensated_image = MC_BCS_SPL_MotionCompensation( ...
        current_frame, reference_frame, block_size, accuracy);
    reconstructed_frame = ...
        MC_BCS_SPL_ResidualReconstruction(motion_compensated_image,...
        y, Phi, num_rows, num_cols);
end