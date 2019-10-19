function icpwrite(filename, A)
% ICPWRITE(FILENAME, A)
%
% Writes the image A to FILENAME as a QccPack ICP-format
% file. A should be an array of size
%   num_rows x num_cols
%
% Note that the values in A are cast as doubles before being written
% to filename as 32-bit floating-point numbers, as this is the
% requisite format for ICP files.
%

if (ndims(A) ~= 2)
  error('Array must be two-dimensional');
end

fid = fopen(filename , 'w', 'ieee-be');

[num_rows num_cols] = size(A);

fprintf(fid, 'ICP0.50\n');
fprintf(fid, '%d %d\n', [num_cols num_rows ]);
fprintf(fid, '% 16.9e % 16.9e\n', double([min(min(A)) max(max(A))]));

fwrite(fid, reshape(A', (num_rows * num_cols), 1), 'float32');

fclose(fid);
