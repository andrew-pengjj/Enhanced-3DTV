% 
% function [reconstructed_frames reconstructed_frames_indep] = ...
%     MC_BCS_SPL_Decoder_ForwardBackward(y, Phi, num_rows, num_cols)
% 
%   This function performs consecutive reconstruction which starts from
%   the key frame of GOP which normally sampled at high subrate. 
%   In this way, reconstructed current frame will be reference frame for
%   the next frame. This function returns reconstructed frames by MC-BCS-SPL
%   and independent reconstruction by BCS-SPL. 
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
    MC_BCS_SPL_Decoder_ForwardBackward(y, Phi, num_rows, num_cols)

num_frames = length(y);
[M N]= size(y{2});
reconstructed_frames = cell(1,num_frames);
reconstructed_frames_indep = cell(1,num_frames);

% first reference iamge - full reconstruction
reconstructed_reference_frame = BCS_SPL_DCT_Decoder(y{1}, Phi, num_rows, num_cols);
reconstructed_frames{1} = reconstructed_reference_frame;
reconstructed_frames_indep{1} = reconstructed_reference_frame;

Phi_sub = Phi(1:M,:);
for i = 2:num_frames
    [reconstructed_frames{i} reconstructed_frames_indep{i}] = ...
        MC_BCS_SPL_Decoder(reconstructed_frames{i-1}, y{i}, Phi_sub, num_rows, num_cols);
end
