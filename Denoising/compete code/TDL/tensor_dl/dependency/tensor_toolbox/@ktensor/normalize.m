function X = normalize(X,N)
%NORMALIZE Normalizes the columns of the factor matrices.
%
%   NORMALIZE(X) normalizes the columns of each factor matrix, absorbing
%   the excess weight into lambda. Also ensures that lambda is positive.  
%
%   NORMALIZE(X,N) absorbs the weights into the Nth factor matrix instead
%   of lambda. (All the lambda values are 1.)
%
%   NORMALIZE(X,0) equally divides the weights across the factor matrices.
%   (All the lambda values are 1.)
%
%   See also KTENSOR, NCOMPONENTS.
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
% $Id: normalize.m,v 1.4 2010/03/19 23:46:30 tgkolda Exp $

%% Ensure that matrices are normalized
for r = 1:length(X.lambda)
    for n = 1:ndims(X)
        tmp = norm(X.u{n}(:,r));
        if (tmp > 0)            
            X.u{n}(:,r) = X.u{n}(:,r) / tmp;
        end
        X.lambda(r) = X.lambda(r) * tmp;        
    end
end

%% Check that all the lambda values are positive
idx = find(X.lambda < 0);
X.u{1}(:,idx) = -1 * X.u{1}(:,idx);
X.lambda(idx) = -1 * X.lambda(idx);

%% Absorb the weight into one factor, if requested
if exist('N','var');
    if (N == 0)
        D = diag(nthroot(X.lambda,ndims(X)));
        X.u = cellfun(@(x) x*D, X.u, 'UniformOutput', false);
        X.lambda = ones(size(X.lambda));
    else
        X.u{N} = X.u{N} * diag(X.lambda);
        X.lambda = ones(size(X.lambda));
    end
end



