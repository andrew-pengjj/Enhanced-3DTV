function [f,g] = cp_wfun(Zdata,W,x,normZsqr,memflag)
%CP_WFUN Computes function and gradient for weighted CP.
%
%   [F,G] = CP_WFUN(Z,W,x,normZsqr) calculates the function and gradient
%   for the function 0.5 * || W .* (Z - ktensor(A)) ||^2 where W is an
%   indicator for missing data (0 = missing, 1 = present), Z is the data
%   tensor that is being fit (assumed that missing entries have already
%   been set to zero), A is a cell array of factor matrices that is created
%   from the vector x, and normZsqr in the norm of Z squared.
%
%   [F,G] = CP_WFUN(Zvals,W,x,normZsqr) is a special version that takes
%   just the nonzeros in Z as calculated by the helper function
%   CP_WFG_SPARSE_SETUP.
%
%   [F,G] = CP_WFUN(....,false) uses a more memory efficient version for
%   the sparse code.
%
%   See also CP_WFG, CP_WFG_SPARSE, CP_WFG_SPARSE_SETUP
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by Brett Bader, Tamara Kolda,
% Evrim Acar, and Daniel Dunlavy.
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: cp_wfun.m,v 1.2 2010/03/19 23:46:32 tgkolda Exp $

%% Convert x to factor matrices (i.e., a cell array).
% Normally we would pass in the data tensor, but we may have a data
% tensor or a data array if we are doing the sparse
% calculation. Therefore, we exploit the fact that W is the same
% size as Z and pass it into the function.
A = cp_vec_to_fac(x,W);

%% Compute the function and gradient
if isa(Zdata,'tensor') || isa(Zdata,'sptensor')
    if ~exist('normZsqr','var')
        normZsqr = norm(Zdata)^2;
    end
    [f,G] = cp_wfg(Zdata,W,A,normZsqr);
else
    if ~exist('normZsqr','var')
        normZsqr = sum(Zdata.^2);
    end
    if ~exist('memflag','var')
        memflag = true;
    end
    [f,G] = cp_wfg_sparse(Zdata,W,A,normZsqr,memflag);
end

%% Convert gradient to a vector
g = fac_to_vec(G);


