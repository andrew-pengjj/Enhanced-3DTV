function Z = mrdivide(X,Y)
%MRDIVIDE Slash right division for tensors.
%
%   MRDIVIDE(A,B) is called for the syntax 'A / B' when A is a tensor and B
%   is a scalar. 
%
%   Example
%   X = tenrand([4 3 2],5);
%   X / 3
%
%   See also TENSOR, TENSOR/RDIVIDE.
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
% $Id: mrdivide.m,v 1.4 2010/03/19 23:46:31 tgkolda Exp $

if isscalar(Y)
    Z = tenfun(@rdivide,X,Y);
    return;
end

error('MRDIVIDE only supports the scalar case for tensors');

