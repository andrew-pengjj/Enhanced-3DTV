function [f,G] = cp_fg(Z,A,Znormsqr)
%CP_FG Computes function and gradient of the CP function.
%
%   [F,G] = CP_FG(Z,A) calculates F = (1/2) ||Z - ktensor(A)||^2 where
%   Z is an N-way tensor and A is a ktensor or a cell array with N
%   factor matrices. It also calculates the gradient of the CP fit
%   function where Z is an N-way tensor and A is a ktensor or a
%   cell array with N factor matrices. The result is also a cell
%   array with N factor matrices corresponding to the gradients; in
%   other words, G{n}(:,r) is the partial derivative of the fit
%   function with respect to A{n}(:,r). 
%
%   [F,G] = CP_FG(Z,W,A,NORMZSQR) also passes in the pre-computed
%   norm of Z, which makes the computations faster. 
%
%   See also CP_OPT, CP_FUN.
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
% $Id: cp_fg.m,v 1.3 2010/03/19 23:46:32 tgkolda Exp $


%% Set-up
if ~isa(Z,'tensor') 
    error('Z must be a tensor');
end
N = ndims(Z);

if ~iscell(A) && ~isa(A,'ktensor');
    error('A must be a cell array or ktensor');
end

if isa(A,'ktensor')
    A = tocell(A);
end
R = size(A{1},2);

%% Upsilon and Gamma
Upsilon = cell(N,1);
for n = 1:N
    Upsilon{n} = A{n}'*A{n};
end

Gamma = cell(N,1);
for n = 1:N
    Gamma{n} = ones(R,R);
    for m = [1:n-1,n+1:N]
        Gamma{n} = Gamma{n} .* Upsilon{m};
    end
end


%% Calculation

%F1
if exist('Znormsqr','var')
    f_1 = Znormsqr;
else
    f_1 = norm(Z)^2;
end

%% Calculate gradient and F2
G = cell(N,1);
U = mttkrp(Z,A,1);
V = A{1} .* U;
f_2 = sum(V(:));
G{1} = -U + A{1}*Gamma{1};
for n = 2:N
    U = mttkrp(Z,A,n);
    G{n} = -U + A{n}*Gamma{n};
end

%F3
W = Gamma{1} .* Upsilon{1};
f_3 = sum(W(:));

%SUM
f = 0.5 * f_1 - f_2 + 0.5 * f_3;


