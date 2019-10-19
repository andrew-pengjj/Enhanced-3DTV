function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for tenmat.  
%
%   Examples 
%   X = tenmat(rand(3,4,2),1); 
%   X(1:2,1:2) = ones(2,2); <-- Calls SUBSASGN 
%
%   See also TENMAT, TENMAT/SUBSREF.
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
% $Id: subsasgn.m,v 1.7 2010/03/19 23:46:30 tgkolda Exp $

switch s.type    
    case '()'
        [m n] = size(t.data);
	t.data(s.subs{:}) = b;
        if ~isequal([m n],size(t.data))
            error('Ambiguous change in size')
        end
    otherwise
        error('Invalid assignment for tenmat.')
end


