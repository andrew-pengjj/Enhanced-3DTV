% 
% function psnr = run_simple_two_frames()
% 
%   This function runs the experiments for MC-BCS-SPL on current frame
%   with given independent reconstruction of reference frame
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
function PSNR_mc = run_simple_two_frames()
addpath(genpath('Sequences'));

sequence_name = 'susie';

key_subrate = 0.7;
subrate = 0.1;

reference_frame = double(imread([ sequence_name '.000.pgm']));
current_frame = double(imread([ sequence_name '.001.pgm']));

[num_rows num_cols] = size(current_frame);

block_size = 16;
N = block_size*block_size;
M = round(key_subrate*N);
Phi = orth(randn(N,M))'; 

x = im2col(reference_frame, [block_size block_size], 'distinct');
y = Phi * x;
reconstructed_reference_frame = BCS_SPL_DCT_Decoder(y, Phi, num_rows, num_cols);

M = round(subrate*N);
Phi = Phi(1:M,:); 
x = im2col(current_frame, [block_size block_size], 'distinct');
y = Phi * x;

disp('< Current frame Reconstruction by MC-BCS-SPL >');
disp(['Reference frame Subrate(M/N): ' num2str(key_subrate)  ... 
    ' , Current frame subrate: ' num2str(subrate) ]);
tic
[reconstructed_current_frame reconstructed_current_frame_indep] =  ...
    MC_BCS_SPL_Decoder(reconstructed_reference_frame, y, Phi, num_rows, num_cols);
toc
PSNR_indep = PSNR(current_frame, reconstructed_current_frame_indep);
PSNR_mc = PSNR(current_frame, reconstructed_current_frame);
disp('< PSNR result Comparison (in dB) >');
disp(['Independent: ' num2str(PSNR_indep) ', MC-BCS-SPL: ' num2str(PSNR_mc) ]);
