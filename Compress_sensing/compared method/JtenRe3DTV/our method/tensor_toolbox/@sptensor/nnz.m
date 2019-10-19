function a = nnz(t)
%NNZ Number of nonzeros in sparse tensor.
%
%   NNZ(T) is the number of nonzero elements in T.
%
%   See also SPTENSOR, SPTENSOR/FIND.
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
% $Id: nnz.m,v 1.9 2010/03/19 23:46:30 tgkolda Exp $

if isempty(t.subs)
    a = 0;
else
    a = size(t.subs,1);
end
