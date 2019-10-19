function t = reshape(t,siz)
%RESHAPE Change tensor size.
%   RESHAPE(X,SIZ) returns the tensor whose elements
%   have been reshaped to the appropriate size.
%
%   See also SQUEEZE.
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
% $Id: reshape.m,v 1.3 2010/03/19 23:46:31 tgkolda Exp $

if prod(t.size) ~= prod(siz)
    error('Number of elements cannot change');
end

t.data = reshape(t.data,siz);
t.size = siz;
