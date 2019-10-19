% 
% function [reconstructed_frames reconstructed_frames_indep] = ...
%     MC_BCS_SPL_Decoder_GOP(y, Phi_full, GOP_size, num_rows, num_cols)
% 
%   This function performs Forward/Backward reconstruction from the given
%   two key frames to center frame. Center frame is reconstructed by
%   bidirectional reconstruction. Further enhancement is possible with
%   the reconstructed center frame whose quality is emperically higher than
%   other reconstructed frames.
%   See:
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

function [reconstructed_frames reconstructed_frames_indep] = ...
    MC_BCS_SPL_Decoder_GOP(y, Phi_full, GOP_size, num_rows, num_cols)

ind_center_frame = GOP_size/2 + 1;  

disp('< Forward Processing >');
[recon_forward_frames recon_forward_frames_indep] = ...
    MC_BCS_SPL_Decoder_ForwardBackward(y(1:ind_center_frame-1), Phi_full, ...
    num_rows, num_cols);

disp('< Backward Processing >' );
[recon_backward_frames recon_backward_frames_indep]= ...
    MC_BCS_SPL_Decoder_ForwardBackward( y(GOP_size+1:-1:ind_center_frame+1), ...
    Phi_full, num_rows, num_cols);

[M N] = size(y{2});
Phi = Phi_full(1:M,:);

disp('< Bidirectional Processing for Center frame >');
[recon_center_frame recon_center_frame_indep] = ...
    MC_BCS_SPL_Decoder_Center(y{ind_center_frame}, ...
    recon_forward_frames{end}, recon_backward_frames{end}, Phi);

reconstructed_frames = ...
    [recon_forward_frames recon_center_frame recon_backward_frames(end:-1:2)];
reconstructed_frames_indep = ...
    [recon_forward_frames_indep recon_center_frame_indep recon_backward_frames_indep(end:-1:2)];

% Enhancing other frames in reverse order only when GOP size is 8
if GOP_size == 8
    disp('< Enhancement on 3th, 4th, and 5th frames >');
    reconstructed_frames(4:6) = ...
        MC_BCS_SPL_Decoder_Enhancement(y(4:6), reconstructed_frames(4:6), Phi);
    
    disp('< Enhancement on 1st, 2nd, and 3rd frames >');
    reconstructed_frames(2:4) = ...
        MC_BCS_SPL_Decoder_Enhancement(y(2:4), reconstructed_frames(2:4), Phi);
    
    disp('< Enhancement on 5th, 6th, and 7th frames >');
    reconstructed_frames(6:8) = ...
        MC_BCS_SPL_Decoder_Enhancement(y(6:8), reconstructed_frames(6:8), Phi);
end