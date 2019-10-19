function [f,g] = cp_fun(x,Z,Znormsqr)
%CP_FUN Calculate function and gradient for CP fit function.
%
%  [F,G] = CP_FUN(X,Z) where X is a vector containing the entries of the
%  components of the model and Z is the tensor to be fit.
%
%  See also CP_VEC_TO_FAC, FAC_TO_VEC, CP_FG, CP_OPT
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
% $Id: cp_fun.m,v 1.3 2010/03/19 23:46:32 tgkolda Exp $

%% Convert x to a cell array of matrices
A = cp_vec_to_fac(x,Z);

%% Call cp_fit and cp_gradient using cp_fg
[f,G]=cp_fg(Z,A,Znormsqr);

%% Convert a cell array to a vector
g = fac_to_vec(G);


