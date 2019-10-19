function Z = times(X,Y)
%TIMES Array multiplication for tensors.
%
%   TIMES(A,B) is called for the syntax 'A .* B' when A or B is a 
%   tensor. A and B must have the same size, unless one is a scalar. A
%   scalar can be multiplied by a tensor of any size.
%
%   See also TENSOR.
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
% $Id: times.m,v 1.10 2010/03/19 23:46:31 tgkolda Exp $

Z = tenfun(@times,X,Y);
