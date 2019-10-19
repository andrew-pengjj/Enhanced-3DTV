% 
% function [y Phi] = MC_BCS_SPL_Encoder(frames, subrate, key_subrate)
% 
%   This function performs BCS projections of each block of
%   the group of frames. Based on assumption that the key frame subrate
%   is larger than the non-key frame subrate, the program returns Phi 
%   whose rows are determined by key frame. The number of columns of 
%   the projection matrix, Phi, is determined by the size of the blocks
%   into which current_image is partitioned. Actual Projection is 
%   performed by BCS_SPL_Encoder frame by frame. The projections are 
%   returned as the columns of y.
%
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

function [y Phi] = MC_BCS_SPL_Encoder(frames, subrate, key_subrate)


[num_rows num_cols] = size(frames{1});
num_frames = length(frames);

block_size = 16;

block_num_rows = floor(num_rows / block_size);
block_num_cols = floor(num_cols / block_size);

if ((num_rows ~= block_num_rows * block_size) || ...
      (num_cols ~= block_num_cols * block_size))
  disp('Frame size not integer multiple of block size');
  return;
end

% calculating M which is subrate of block
N = block_size*block_size;
M = round(key_subrate*N);

% construct random basis measurement matrix M*N which is undersampling overator
% first frame is full-subsampled
Phi = orth(randn(N,M))'; 

y = cell(1,num_frames);

y{1} = BCS_SPL_Encoder(frames{1}, Phi);

% other frame
M = round(subrate*N);
Phi_sub = Phi(1:M,:);

for i = 2 : num_frames
    y{i} = BCS_SPL_Encoder(frames{i}, Phi_sub);
end

y{end} = BCS_SPL_Encoder(frames{end}, Phi);
