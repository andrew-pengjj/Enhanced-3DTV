% 
% function motion_compensated_image = ...
%     MC_BCS_SPL_MotionCompensation(current_image, reference_image, ...
%     block_size, accuracy)
% 
%   This function performs the motion estimation and compensation between
%     current image and reference image. For Linux/Unix environment, mex
%     file is provided. For Windows system, compiled binary file is 
%     provided. 
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

function motion_compensated_image = ...
    MC_BCS_SPL_MotionCompensation(current_image, reference_image, ...
    block_size, accuracy)

window_size = 20;

switch accuracy
  case 'fp'
    motion_accuracy = '-fp';
  case 'hp'
    motion_accuracy = '-hp';
  case 'qp'
    motion_accuracy = '-qp';
  otherwise
    disp('Invalid motion-vector accuracy');
    return;
end

filename = tempname;
mc_filename = [filename '.mc.icp'];
reference_filename = [filename '.ref.icp'];
current_filename = [filename '.cur.icp'];
mvfilename = [filename '.mv.txt'];
icpwrite(current_filename, current_image);
icpwrite(reference_filename, reference_image);

% str = computer('arch');
% if strcmp(str,'win32')
%     unix(['memc_w32.exe -b ' num2str(block_size) ...
%         ' ' motion_accuracy ' -f1 mpeg4.flt -w ' num2str(window_size) ... 
%         ' -mc ' mc_filename ' ' reference_filename ' ' current_filename ...
%         ' ' mvfilename]);
% elseif (strcmp(str,'glnx86') || strcmp(str,'glnxa64'))
    memc('memc','-b',num2str(block_size),motion_accuracy,'-f1', ...
        'mpeg4.flt','-w',num2str(window_size), ...
        '-mc',mc_filename,reference_filename,current_filename,mvfilename);
% else
%     disp(['MEMC does not support this machine: ' str]);
%     disp('Please email to sm655@msstate.edu');
% % end

motion_compensated_image = icpread(mc_filename);
delete(current_filename);
delete(reference_filename);
delete(mc_filename);
delete(mvfilename);
