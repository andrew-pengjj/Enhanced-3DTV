function S = tt_matrix2cellstr(M)
%TT_MATRIX2CELLSTR Convert a matrix to a cell array of strings.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: tt_matrix2cellstr.m,v 1.10 2010/03/19 23:46:31 tgkolda Exp $

fmt = get(0,'FormatSpacing');
format compact
S = evalc('disp(M)');
set(0,'FormatSpacing',fmt)
S = textscan(S,'%s','delimiter','\n','whitespace','');
S = S{1};
end
