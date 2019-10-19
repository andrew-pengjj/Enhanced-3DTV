function C = ldivide(A,B)
%LDIVIDE Array right division for sparse tensors.
%
%   LDIVIDE(A,B) is called for the syntax 'A .\ B' when A or B is a sparse
%   tensor. A and B must have the same size, unless one is a scalar. 
%
%   See also SPTENSOR.
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
% $Id: ldivide.m,v 1.4 2010/03/19 23:46:30 tgkolda Exp $

C = rdivide(B,A);
