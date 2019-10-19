function ok = tt_valscheck(vals)
%TT_VALSCHECK Checks for valid values.
%
%  TT_VALSCHECK(S) throws an error if S is not a valid values
%  array, which means that S is a column array.
%
%  X = TT_VALSCHECK(S) returns true if S is a valid and false otherwise.
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
% $Id: tt_valscheck.m,v 1.5 2010/03/19 23:46:31 tgkolda Exp $

if isempty(vals)
    ok = true;
elseif ndims(vals) == 2 && size(vals,2) == 1
    ok = true;
else
    ok = false;
end

if ~ok && nargout == 0
    error('Values must be a column array');
end
