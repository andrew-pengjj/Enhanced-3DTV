function z = isequal(x,y)
%ISEQUAL for tensors.
%
%   ISEQUAL(A,B) compares the tensors A and B for equality.
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
% $Id: isequal.m,v 1.4 2010/03/19 23:46:30 tgkolda Exp $

%%
if ~isequal(x.size,y.size) 
    z = false;
elseif isa(x,'tensor') && isa(y,'tensor') 
    z = isequal(x.data,y.data);
elseif isa(y,'sptensor')
    z = isequal(x,full(y));
else
    z = false;
end
