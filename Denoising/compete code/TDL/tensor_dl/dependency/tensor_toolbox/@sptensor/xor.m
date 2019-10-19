function C = xor(A,B)
%XOR Logical XOR for sptensors.
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
% $Id: xor.m,v 1.7 2010/03/19 23:46:30 tgkolda Exp $

%% Observations for sparse matrix case.
% The result of xor(a,5) is dense!
% The result of xor(a,0) is dense!
% The result of xor(a,full(a)) is dense!
% The result of xor(a,b) is sparse.

%% Case 1: One argument is a scalar
if isscalar(B) || isa(B,'tensor')
    C = xor(full(A),B);
    return;
end
if isscalar(A)
    C = xor(A,full(B));
    return;
end


%% Case 2: Both x and y are tensors of some sort
if ~isequal(size(A),size(B))
    error('Must be tensors of the same size');
end

if isa(A,'sptensor') && isa(B,'sptensor')
    C = sptensor([A.subs; B.subs], 1, size(A), @(x) length(x) == 1);
    return;
end

if isa(B,'tensor')
    Bsubs = find(B ~= 0);
    C = sptensor([A.subs; Bsubs], 1, size(A), @(x) length(x) == 1);
    return;    
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
