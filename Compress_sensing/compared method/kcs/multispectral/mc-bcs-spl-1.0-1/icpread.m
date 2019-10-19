function A = icpread(filename)
% A = ICPREAD(FILENAME)
% 
% Reads the QccPack ICP-format file specified by FILENAME into A.
%
% A is returned as a double array of size
%   num_rows x num_cols
%

fid = fopen(filename, 'r', 'ieee-be');
header_value = fgetl(fid);
magic_num = sscanf(header_value, '%[A-Z]');

if (magic_num ~= 'ICP')
  error([filename ' is not an ICP file'])
end

num_cols = fscanf(fid, '%d', 1);
num_rows = fscanf(fid, '%d', 1);
min_val = fscanf(fid, '%f', 1);
max_val = sscanf(fgets(fid), '%f', 1);

data = fread(fid, num_rows * num_cols, 'float32');

A = reshape(data, num_cols, num_rows)';

fclose(fid);
