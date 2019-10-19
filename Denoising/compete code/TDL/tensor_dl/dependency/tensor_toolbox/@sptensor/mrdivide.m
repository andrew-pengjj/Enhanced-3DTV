function C = mrdivide(A,B)
%MRDIVIDE Slash right division for sparse tensors.
%
%   MRDIVIDE(A,B) is called for the syntax 'A / B' when A is a sparse
%   tensor and B is a scalar.
%
%   Example
%   X = sptenrand([4 3 2],5);
%   X / 3
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
% $Id: mrdivide.m,v 1.4 2010/03/19 23:46:30 tgkolda Exp $

if isscalar(B)
    newsubs = A.subs;
    newvals = A.vals / B;
    if B == 0
        nansubs = setdiff(allsubs(A),newsubs,'rows');
        newsubs = [newsubs; nansubs];
        newvals = [newvals; repmat(NaN,size(nansubs,1),1)];
    end
    C = sptensor(newsubs,newvals,A.size);
    return;
end

error('MRDIVIDE only supports the scalar case for sparse tensors');
