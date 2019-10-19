function U = tocell(X,N)
%TOCELL Convert X to a cell array.
%
%   TOCELL(X) converts X to a cell array, evenly distributing the
%   weight in lambda.
%
%   TOCELL(X,N) absorbs the weights into the Nth factor matrix.
%
%   See also KTENSOR, NORMALIZE.
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
% $Id: tocell.m,v 1.3 2010/03/19 23:46:30 tgkolda Exp $

if exist('N','var')
    X = normalize(X,N);
    U = X.u;
    return;
end

if isequal(X.lambda,ones(size(X.lambda)))
    U = X.u;
    return;
end


lsplit = nthroot(X.lambda,ndims(X));
R = length(X.lambda);
U = X.u;
D = diag(lsplit);
for n = 1:ndims(X)
    U{n} = U{n} * D;
end
        

