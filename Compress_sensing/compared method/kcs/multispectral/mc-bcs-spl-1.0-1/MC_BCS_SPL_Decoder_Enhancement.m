% 
%function [reconstructed_center_frame recon_indep] = ...
%    MC_BCS_SPL_Decoder_Center(y, previous_frame, next_frame, Phi)
% 
%   This function peforms MC-BCS-SPL reconstruction on center frame 
%   based on the previous and next frame. Emperically, it provide high 
%   quality of reconstruction, which can help find more accurated motion 
%   vectors on previous and next frame. based on the reconstructed center 
%   frame. 
%   This function returns updated reconstructed frame by MC-BCS-SPL
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

function reconstructed_frames = ... 
    MC_BCS_SPL_Decoder_Enhancement(y, reconstructed_frames, Phi)

[M N]= size(Phi);
block_size = sqrt(N);

previous = 1; current = 2; next = 3;

num_iter = 5; 
accuracy = 'qp';
% refining current frame again
reference_frame = reconstructed_frames{previous};
recon_current_frame_forward = ...
    MC_BCS_SPL_IterativeReconstruction(reconstructed_frames{current}, ...
    reference_frame, y{current}, Phi, num_iter, block_size, accuracy); 

reference_frame = reconstructed_frames{next};
recon_current_frame_backward = ...
    MC_BCS_SPL_IterativeReconstruction(reconstructed_frames{current}, ...
    reference_frame, y{current}, Phi, num_iter, block_size, accuracy); 

% average result
reconstructed_frames{current} = ...
    (recon_current_frame_forward + recon_current_frame_backward)/2;

% previous frame enhancement
reference_frame = reconstructed_frames{current};
reconstructed_frames{previous} = ...
    MC_BCS_SPL_IterativeReconstruction(reconstructed_frames{previous}, ...
    reference_frame, y{previous}, Phi, num_iter, block_size, accuracy); 

% next frame enhancement
reconstructed_frames{next} = ...
    MC_BCS_SPL_IterativeReconstruction(reconstructed_frames{next}, ...
    reference_frame, y{next}, Phi, num_iter, block_size, accuracy); 

% subfunction
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
