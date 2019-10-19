function [P, P0, output] = cp_wopt(Z,W,R,varargin)
%CP_WOPT Fits a weighted CP model to a tensor via optimization.
%
%   K = CP_WOPT(X,W,R) fits an R-component weighted CANDECOMP/PARAFAC
%   (CP) model to the tensor X, where W is an indicator for missing
%   data (0 = missing, 1 = present). The result K is a ktensor. It is
%   assumed that missing entries of X have been sent to zero (but not
%   that all zeros correspond to missing entries.) The function being
%   optimized is F(K) = 1/2 || W .* (X - K) ||^2.
% 
%   K = CP_WOPT(X,W,R,'param', value,...) specifies additional
%   parameters for the method. Specifically...
%
%   'alg' - Specfies optimization algorithm (default: 'ncg')
%      'ncg'   Nonlinear Conjugate Gradient Method
%      'lbfgs' Limited-Memory BFGS Method
%      'tn'    Truncated Newton
%
%   'init' - Initialization for factor matrices. (default:
%   'random'). This can be a cell array with the initial matrices or
%   one of the following strings:
%      'random' Randomly generated via randn function
%      'nvecs'  Selected as leading left singular vectors of X(n)
%
%   'alg_options' - Parameter settings for selected optimization
%   algorithm. For example, type OPTIONS = NCG('defaults') to get
%   the NCG algorithm options which can then me modified as passed
%   through this function to NCG.
%   
%   'fun' - Specifies the type of implementation (default: 'auto')
%       'auto'           Dense implementation
%       'sparse'         Sparse implementation
%       'sparse_lowmem'  More memory efficient version of sparse implementation
%
%   [K, U0] = CP_WOPT(...) also returns the initial guess.
%
%   [K, U0, OUT] = CP_WOPT(...) also returns a structure with the
%   optimization exit flag, the final relative fit, and the full
%   output from the optimization method.
%
%   See also CP_WFUN, CP_WFG, CP_WFG_SPARSE, CP_WFG_SPARSE_SETUP,
%   CP_VEC_TO_FAC, CP_OPT.
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
% $Id: cp_wopt.m,v 1.5 2010/03/22 18:07:17 tgkolda Exp $

%% Check for POBLANO
if ~exist('poblano_params','file')
    error(['CP_WOPT requires Poblano Toolbox for Matlab. This can be ' ...
           'downloaded at http://software.sandia.gov/trac/poblano.']);
end

%% Set parameters
params = inputParser;
params.addParamValue('alg','ncg', @(x) ismember(x,{'ncg','tn','lbfgs'}));
params.addParamValue('init','random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addParamValue('fun','auto', @(x) ismember(x,{'auto','default','sparse','sparse_lowmem'}));
params.addParamValue('alg_options', '', @isstruct);
params.parse(varargin{:});

%% Set up optimization algorithm
switch (params.Results.alg)
    case 'ncg'
        opthandle = @ncg;
    case 'tn'
        opthandle = @tn;
    case 'lbfgs'
        opthandle = @lbfgs;
end

%% Set up optimization algorithm options
if isempty(params.Results.alg_options)
    options = feval(opthandle, 'defaults');
else
    options = params.Results.alg_options;
end

%% Set up function handle
normZsqr = norm(Z)^2;
funtype = params.Results.fun;

if (isequal(funtype,'auto') && isa(Z,'tensor')) || isequal(funtype,'default')
    funhandle = @(x) cp_wfun(Z,W,x,normZsqr);
else
    if ~isa(Z,'sptensor') || ~isa(W,'sptensor')
        warning('Converting dense tensor to sparse');
        Z = sptensor(Z);
        W = sptensor(W);
    end
    Zvals = cp_wfg_sparse_setup(Z,W);
    fflag = ~isequal(funtype,'sparse_lowmem');
    funhandle = @(x) cp_wfun(Zvals,W,x,normZsqr,fflag);
end
    
%% Initial guess
sz = size(Z);
N = length(sz);

if iscell(params.Results.init)
    P0 = params.Results.init;
elseif strcmpi(params.Results.init,'random')
    P0 = cell(N,1);
    for n=1:N
        P0{n} = randn(sz(n),R);
        for j=1:R
            P0{n}(:,j) = P0{n}(:,j) / norm(P0{n}(:,j));
        end
    end
elseif strcmpi(params.Results.init,'nvecs')
    P0 = cell(N,1);
    for n=1:N
        P0{n} = nvecs(Z,n,R);
    end
else
    error('Initialization type not supported')
end

%% Fit CP using CP_WOPT by ignoring missing entries
out = feval(opthandle, funhandle, fac_to_vec(P0), options);

P  = ktensor(cp_vec_to_fac(out.X,Z));
if nargout > 1
    output.ExitFlag  = out.ExitFlag;
    output.FuncEvals = out.FuncEvals;
    output.f = out.F;
    output.G = cp_vec_to_fac(out.G,W);
end


%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);
% Fix the signs
P = fixsigns(P);
