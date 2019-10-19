function T = permute(T,order)
%PERMUTE Permute tensor dimensions.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The result has the
%   same values of A, but the order of the subscripts needed to access
%   any particular element are rearranged as specified by ORDER.
%
%   See also TENSOR, TENSOR/SIZE, TENSOR/NDIMS, PERMUTE.
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
% $Id: permute.m,v 1.11 2010/03/19 23:46:31 tgkolda Exp $

if ndims(T) ~= numel(order)
  error('Invalid permutation order');
end

% Check for special case of permuting an order-1 object (which has
% no effect but confuses MATLAB's permute command which doesn't
% think that there is such a thing as a 1D-array.
if isequal(order,1)
    return;
end

% Note that permute does error checking on order, so we don't worry
% about it. 
T.data = permute(T.data,order);
T.size = T.size(order);

return;
