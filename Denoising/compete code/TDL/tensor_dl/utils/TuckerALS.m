function [T, Uinit] = TuckerALS(X, R, varargin)

%==========================================================================
% TuckerALS is an analogy of standard function tucker_als in tensor
%   toolbox. This function is based on the original code in tucker_als.m,
%   except for some detailed implementation adapted to the tensor
%   dictionary learning.
%
% Syntax:
%   [ T, Uinit ] = TuckerALS(X, R, varargin);
%
% See also tucker_als
%
% by Yi Peng
%==========================================================================

% Extract number of dimensions and norm of X.
N = ndims(X);
normX = norm(X);

%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParamValue('tol',1e-4,@isscalar);
params.addParamValue('maxiters',50,@(x) isscalar(x) & x > 0);
params.addParamValue('dimorder',1:N,@(x) isequal(sort(x),1:N));
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','my_nvecs','eigs'})));
params.addParamValue('printitn',1,@isscalar);
params.parse(varargin{:});

%% Copy from params object
fitchangetol = params.Results.tol;
maxiters = params.Results.maxiters;
dimorder = params.Results.dimorder;
init = params.Results.init;
printitn = params.Results.printitn;

if numel(R) == 1
    R = R * ones(N,1);
end

%% Error checking 
% Error checking on maxiters
if maxiters < 0
    error('OPTS.maxiters must be positive');
end

% Error checking on dimorder
if ~isequal(1:N,sort(dimorder))
    error('OPTS.dimorder must include all elements from 1 to ndims(X)');
end

%% Set up and error checking on initial guess for U.
if iscell(init)
    Uinit = init;
    if numel(Uinit) ~= N
        error('OPTS.init does not have %d cells',N);
    end
    for n = dimorder(2:end);
        if ~isequal(size(Uinit{n}),[size(X,n) R(n)])
            error('OPTS.init{%d} is the wrong size',n);
        end
    end
else
    % Observe that we don't need to calculate an initial guess for the
    % first index in dimorder because that will be solved for in the first
    % inner iteration.
    if strcmp(init,'random')
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            Uinit{n} = rand(size(X,n),R(n));
        end
    elseif strcmp(init,'my_nvecs') || strcmp(init,'eigs') 
        % Compute an orthonormal basis for the dominant
        % Rn-dimensional left singular subspace of
        % X_(n) (1 <= n <= N).
        Uinit = cell(N,1);
        for n = dimorder(2:end)
            fprintf('  Computing %d leading e-vectors for factor %d.\n', ...
                    R(n),n);
            [Uinit{n}, R(n)] = my_nvecs(X,n,R(n));
        end
    else
        error('The selected initialization method is not supported');
    end
end

%% Set up for iterations - initializing U and the fit.
U = Uinit;
fit = 0;

if printitn > 0
    fprintf('\nTucker Alternating Least-Squares:\n');
end

%% Main Loop: Iterate until convergence
for iter = 1:maxiters

    fitold = fit;

    % Iterate over all N modes of the tensor
    for n = dimorder(1:end)
	Utilde = ttm(X, U, -n, 't');

	% Maximize norm(Utilde x_n W') wrt W and
	% keeping orthonormality of W
        [U{n}, R(n)] = my_nvecs(Utilde,n,R(n));
    end

    % Assemble the current approximation
    core = ttm(Utilde, U, n, 't');

    % Compute fit
    normresidual = sqrt( normX^2 - norm(core)^2 );
    fit = 1 - (normresidual / normX); %fraction explained by model
    fitchange = abs(fitold - fit);

    if mod(iter,printitn)==0
        fprintf(' Iter %2d: fit = %e fitdelta = %7.1e\n', iter, fit, fitchange);
    end

    % Check for convergence
    if (iter > 1) && (fitchange < fitchangetol)
        break;
    end

end

T = ttensor(core, U);

end
