function a = double(s)
%DOUBLE Converts a sparse tensor to a dense multidimensional array.
%
%  See also SPTENSOR, SPTENSOR/FULL.
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
% $Id: double.m,v 1.13 2010/03/19 23:46:30 tgkolda Exp $

a = zeros([size(s) 1 1]);
if nnz(s) > 0
    a(tt_sub2ind(size(s),s.subs)) = s.vals;
end
