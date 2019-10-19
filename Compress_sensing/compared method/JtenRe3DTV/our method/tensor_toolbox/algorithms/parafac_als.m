function [P,Uinit] = parafac_als(X,R,opts)
%PARAFAC_ALS Deprecated. Use CP_ALS instead.
%   
%   See also CP_ALS.
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
% $Id: parafac_als.m,v 1.19 2010/03/19 23:46:32 tgkolda Exp $

if (nargout == 2) & (nargin == 3)
    [P,Uinit] = cp_als(X,R,opts);  
elseif (nargout == 2) & (nargin == 2)
    [P,Uinit] = cp_als(X,R);  
elseif (nargout == 2) & (nargin == 1)
    [P,Uinit] = cp_als(X);  
elseif (nargout == 1) & (nargin == 3)
    P = cp_als(X,R,opts);  
elseif (nargout == 1) & (nargin == 2)
    P = cp_als(X,R);  
elseif (nargout == 1) & (nargin == 1)
    P = cp_als(X);  
end
    

