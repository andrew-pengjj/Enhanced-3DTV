function sol = prox_TV(b, lambda, param)
% PROX_TV - Total variation proximal operator
%
% sol = prox_TV(y, lambda, param) solves:
%
%   min_{x} ||y - x||_2^2 + lambda * ||x||_{TV}
%
% The input argument param contains the following fields:
%
%   - max_iter: max. nb. of iterations (default: 200)
%
%   - rel_obj: minimum relative change of the objective value (default:
%   1e-4)
%       The algorithm stops if
%           | f(x(t))-f(x(t-1)) | / f(x(t)) < rel_obj,
%       where x(t) is the estimate of the solution at iteration t.
%
%   - verbose: 0 no log, 1 a summary log at convergence, 2 print main
%   steps (default: 1)
%
%
%
% Reference:
% [1] A. Beck and  M. Teboulle, "Fast gradient-based algorithms for
% constrained Total Variation Image Denoising and Deblurring Problems", 
% IEEE Transactions on Image Processing, VOL. 18, NO. 11, 2419-2434, 
% November 2009.
%
%
%
% Author: Gilles Puy
% E-mail: gilles.puy@epfl.ch
% Date: October 15, 2010
%

% Optional input arguments

if nargin<3, param=struct; end

if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end

% Initializations
[r, s] = gradient_op(b*0);
pold = r; qold = s;
told = 1; prev_obj = 0;

% Main iterations
if param.verbose > 1
    fprintf('  Proximal TV operator:\n');
end
for iter = 1:param.max_iter
    
    % Current solution
    sol = b - lambda*div_op(r, s);
    
    % Objective function value
    obj = .5*norm(b(:)-sol(:), 2)^2 + lambda * TV_norm(sol);
    rel_obj = abs(obj-prev_obj)/obj;
    prev_obj = obj;
    
    % Stopping criterion
    if param.verbose>1
        fprintf('   Iter %i, obj = %e, rel_obj = %e\n', ...
            iter, obj, rel_obj);
    end
    if rel_obj < param.rel_obj
        crit_TV = 'TOL_EPS'; break;
    end
    
    % Udpate divergence vectors and project
    [dx, dy] = gradient_op(sol);
    r = r - 1/(8*lambda) * dx; s = s - 1/(8*lambda) * dy;
    weights = max(1, sqrt(abs(r).^2+abs(s).^2));
    p = r./weights; q = s./weights;
    
    % FISTA update
    t = (1+sqrt(4*told^2))/2;
    r = p + (told-1)/t * (p - pold); pold = p;
    s = q + (told-1)/t * (q - qold); qold = q;
    told = t;
    
end

% Log after the minimization
if ~exist('crit_TV', 'var'), crit_TV = 'MAX_IT'; end
if param.verbose >= 1
    fprintf(['  Prox_TV: obj = %e, rel_obj = %e,' ...
        ' %s, iter = %i\n'], obj, rel_obj, crit_TV, iter);
end

end