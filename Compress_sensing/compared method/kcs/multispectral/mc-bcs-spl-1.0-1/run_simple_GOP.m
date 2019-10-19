% 
% function psnr = run_simple_GOP()
% 
%   This function runs the experiments for MC-BCS-SPL in Table 1 of
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

function [PSNR_mc PSNR_indep] = run_simple_GOP()

addpath(genpath('Sequences'));

% testing parameters
key_subrate = 0.7;
subrate = 0.1;

sequence_name = 'susie';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_num_frames = 24;
GOP_size = 8 ; % the number of frames in a GOP

if floor(GOP_size/2) ~= ceil(GOP_size/2)
    error('Frame number should be even.')
end
if total_num_frames < GOP_size
    error('Total number of frames should be larger than GOP size')
end
if mod(total_num_frames,GOP_size) ~= 0
    error('Total number of frames should be multiple of GOP size')
end

PSNR_mc = zeros(1,total_num_frames);
PSNR_indep = zeros(1,total_num_frames);
for k = 1:GOP_size:total_num_frames;
    frames = readframes(sequence_name, GOP_size, k);    
    [num_rows num_cols] = size(frames{1});
    % encoding
    [y Phi_full] = MC_BCS_SPL_Encoder(frames, subrate, key_subrate);
    
    % decoding
    [reconstructed_frames reconstructed_frames_indep] ...
        = MC_BCS_SPL_Decoder_GOP(y,Phi_full,GOP_size, num_rows, num_cols);
 
    disp('< PSNR result Comparison (in dB) >');
    for i = 1 : GOP_size;
        PSNR_mc(k+i-1) = PSNR(frames{i}, reconstructed_frames{i});
        PSNR_indep(k+i-1) = PSNR(frames{i}, reconstructed_frames_indep{i});
        disp(['Frame # '  num2str(k+i-1) '  Independent: ' num2str(PSNR_indep(k+i-1)) ...
            ', MC-BCS-SPL: ' num2str(PSNR_mc(k+i-1)) ]);
    end

    % save image result of center frame of GOP
    if (k < GOP_size)
        file_name = [sequence_name '.00' num2str(round(GOP_size/2)) '_' ...
            num2str(key_subrate) '_' num2str(subrate) '.pgm'];
        imwrite(uint8(reconstructed_frames{5}),file_name);
    end
end

function frames = readframes(sequence_name, GOP_size, start_sequence)
frames = cell(1,GOP_size);
index = 1;
for i = start_sequence: start_sequence+GOP_size
    if i <= 10
        file_name = [ sequence_name '.00' num2str(i-1) '.pgm'];
    elseif (i>10) && (i<=100)
        file_name = [ sequence_name '.0' num2str(i-1) '.pgm'];
    else
        file_name = [ sequence_name '.' num2str(i-1) '.pgm'];
    end
    frames{index} = double(imread(file_name));
    index = index +1;
end
