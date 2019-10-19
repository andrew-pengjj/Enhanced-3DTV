function [x,x_debias,objective,times,debias_start,mses]= ...
         IST(y,A,tau,varargin)
  
% This function solves the convex problem
% arg min_theta = 0.5*|| y - A x ||_2^2 + tau ||x||_1
% using the MM/IST algorithm, described in the paper
% M. Figueiredo, and R. Nowak, "A Bound Optimization Approach 
% to Wavelet-Based Image Deconvolution", IEEE International 
% Conference on Image Processing - ICIP'2005, Genoa, Italy, September  2005.
%
% Copyright (2007): Mario Figueiredo and Robert Nowak
%
%  ===== Required inputs =============
%
%  y: 1D vector or 2D array (image) of observations
%     
%  A: if y and x are both 1D vectors, A can be a 
%     k*n (where k is the size of y and n the size of x)
%     matrix or a handle to a function that computes
%     products of the form A*v, for some vector v.
%     In any other case (if y and/or x are 2D arrays), 
%     A has to be passed as a handle to a function which computes 
%     products of the form A*x; another handle to a function 
%     AT which computes products of the form A'*x is also required 
%     in this case. The size of x is determined as the size
%     of the result of applying AT.
%
%  tau: usually, a non-negative real parameter of the objective 
%       function (see above). It can also be an array, the same 
%       size as x, with non-negative entries; in this  case,
%       the objective function weights differently each element 
%       of x, that is, it becomes
%       0.5*|| y - A x ||_2^2 + tau^T * abs(x)
%
%  ===== Optional inputs =============
%
%  
%  'AT'    = function handle for the function that implements
%            the multiplication by the conjugate of A, when A
%            is a function handle. If A is an array, AT is ignored.
%
%  'StopCriterion' = type of stopping criterion to use
%                    0 = algorithm stops when the relative 
%                        change in the number of non-zero 
%                        components of the estimate falls 
%                        below 'ToleranceA'
%                    1 = stop when the relative 
%                       change in the objective function 
%                       falls below 'ToleranceA'
%                    4 = stop when the objective function 
%                        becomes equal or less than toleranceA.
%                    Default = 3.
%
%  'ToleranceA' = stopping threshold; Default = 0.01
% 
%  'Debias'     = debiasing option: 1 = yes, 0 = no.
%                 Default = 0.
%
%  'ToleranceD' = stopping threshold for the debiasing phase:
%                 Default = 0.0001.
%                 If no debiasing takes place, this parameter,
%                 if present, is ignored.
%
%  'MaxiterA' = maximum number of iterations allowed in the
%               main phase of the algorithm.
%               Default = 1000
%
%  'MiniterA' = minimum number of iterations performed in the
%               main phase of the algorithm.
%               Default = 5
%
%  'MaxiterD' = maximum number of iterations allowed in the
%               debising phase of the algorithm.
%               Default = 200
%
%  'MiniterD' = minimum number of iterations to perform in the
%               debiasing phase of the algorithm.
%               Default = 5
%
%  'Initialization' must be one of {0,1,2,array}
%               0 -> Initialization at zero. 
%               1 -> Random initialization.
%               2 -> initialization with A'*y.
%           array -> initialization provided by the user.
%               Default = 0;
%             
%  'True_x' = if the true underlying x is passed in 
%                this argument, MSE plots are generated.
%
% ===================================================  
% ============ Outputs ==============================
%   x = solution of the main algorithm
%
%   x_debias = solution after the debiasing phase;
%                  if no debiasing phase took place, this
%                  variable is empty, x_debias = [].
%
%   objective = sequence of values of the objective function
%
%   times = CPU time after each iteration
%
%   debias_start = iteration number at which the debiasing 
%                  phase started. If no debiasing took place,
%                  this variable is returned as zero.
%
%   mses = sequence of MSE values, with respect to True_x,
%          if it was given; if it was not given, mses is empty,
%          mses = [].
% ========================================================

% test for number of required parametres
if (nargin-length(varargin)) ~= 3
     error('Wrong number of required parameters');
end

% Set the defaults for the optional parameters
stopCriterion = 3;
tolA = 0.01;
tolD = 0.0001;
debias = 0;
maxiter = 1000;
maxiter_debias = 500;
miniter = 5;
miniter_debias = 5;
init = 0;
compute_mse = 0;
AT = 0;

% Set the defaults for outputs that may not be computed
debias_start = 0;
x_debias = [];
mses = [];

% Read the optional parameters
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch varargin{i}
     case 'StopCriterion'
       stopCriterion = varargin{i+1};
     case 'ToleranceA'       
       tolA = varargin{i+1};
     case 'ToleranceD'
       tolD = varargin{i+1};
     case 'Debias'
       debias = varargin{i+1};
     case 'MaxiterA'
       maxiter = varargin{i+1};
     case 'MaxiterD'
       maxiter_debias = varargin{i+1};
     case 'MiniterA'
       miniter = varargin{i+1};
     case 'MiniterD'
       miniter_debias = varargin{i+1};
     case 'Initialization'
       if prod(size(varargin{i+1})) > 1   % we have an initial x
          init = 33333;    % some flag to be used below
          x = varargin{i+1};
       else 
          init = varargin{i+1};
       end
     case 'True_x'
       compute_mse = 1;
       true = varargin{i+1};
     case 'AT'
       AT = varargin{i+1};
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized option: ''' varargin{i} '''']);
    end;
  end;
end
%%%%%%%%%%%%%%

% if A is a function handle, we have to check presence of AT,
if isa(A, 'function_handle') & ~isa(AT,'function_handle')
   error(['The function handle for transpose of A is missing']);
end 

% if A is a matrix, we find out dimensions of y and x,
% and create function handles for multiplication by A and A',
% so that the code below doesn't have to distinguish between
% the handle/not-handle cases
if ~isa(A, 'function_handle')
   AT = @(x) A'*x;
   A = @(x) A*x;
end
% from this point down, A and AT are always function handles.

% Precompute A'*y since it'll be used a lot
Aty = AT(y);

% Initialization
switch init
    case 0   % initialize at zero, using AT to find the size of x
       x = AT(zeros(size(y)));
    case 1   % initialize randomly, using AT to find the size of x
       x = randn(size(AT(zeros(size(y)))));
    case 2   % initialize x0 = A'*y
       x = Aty; 
    case 33333
       % initial x was given as a function argument; just check size
       if size(A(x)) ~= size(y)
          error(['Size of initial x is not compatible with A']); 
       end
    otherwise
       error(['Unknown ''Initialization'' option']);
end

% now check if tau is an array; if it is, it has to 
% have the same size as x
if prod(size(tau)) > 1
   try,
      dummy = x.*tau;
   catch,
      error(['Parameter tau has wrong dimensions; it should be scalar or size(x)']),
   end
end
      
% if the true x was given, check its size
if compute_mse & (size(true) ~= size(x))  
   error(['Initial x has incompatible size']); 
end

% define the indicator vector or matrix of nonzeros in x
nz_x = (x ~= 0.0);
num_nz_x = sum(nz_x(:));

% Compute and store initial value of the objective function
resid =  y - A(x);
f = 0.5*(resid(:)'*resid(:)) + sum(tau(:).*abs(x(:)));

% start the clock
t0 = cputime;

times(1) = cputime - t0;
objective(1) = f;

if compute_mse
   mses(1) = sum(sum((x-true).^2));
end

cont_outer = 1;
iter = 1;

fprintf(1,'   initial obj=%10.6e, nonzeros=%7d\n', f,num_nz_x);

while cont_outer

  x_previous = x;
  
  x = soft(x_previous + AT(y-A(x_previous)),tau);
  
  % update the nonzero indicator vector
  prev_f = f;
  resid =  y - A(x);
  f = 0.5*(resid(:)'*resid(:)) + sum(tau(:).*abs(x(:)));
  nz_x_prev = nz_x;
  nz_x = (x~=0.0);
  num_nz_x = sum(nz_x(:));
  num_changes_active = sum(nz_x(:)~=nz_x_prev(:));
  
  % compute the objective-based stopping criterion, based on the current and
  % previous iterates.
  criterionObjective = abs(f-prev_f)/(prev_f);
  
  % compute the active-set-based criterion, based on the current and
  % previous iterates.
  if num_nz_x >= 1
    criterionActiveSet = num_changes_active / num_nz_x;
  else
    criterionActiveSet = 1.0;
  end

    
  % print out the various stopping criteria
  fprintf(1,'Iter=%4d, obj=%10.6e, nz=%7d, chg=%6d, cObj=%7.3e\n', iter, f, ...
      num_nz_x, num_changes_active, criterionObjective);

  iter = iter + 1;
  prev_f = f;
  objective(iter) = f;
  times(iter) = cputime-t0;
  
  if compute_mse
    err = true - x;
    mses(iter) = (err(:)'*err(:));
  end
  
  % take no less than miniter and no more than maxiter iterations
  if iter<=miniter
    cont_outer = 1;
  else
    switch stopCriterion
      case 0,
	cont_outer = ((iter <= maxiter) & ...
	    (criterionActiveSet > tolA)|(lambda<1.0));
      case 1,
	cont_outer = ((iter <= maxiter) & ...
	    (criterionObjective > tolA));
      case 3,
	cont_outer = ((iter <= maxiter) & ...
	    (criterionObjective > tolA));
      case 4,
	cont_outer = ((iter <= maxiter) & ...
	    (f > tolA));          
      otherwise,
	cont_outer = ((iter <= maxiter) & ...
	    (criterionLCP > tolA));
    end
  end

end % end of the main loop of the BB-QP algorithm

% Printout results
fprintf(1,'\nFinished the main algorithm!\nResults:\n')
fprintf(1,'||A x - y ||_2 = %10.3e\n',resid(:)'*resid(:))
fprintf(1,'||x||_1 = %10.3e\n',sum(abs(x(:))))
fprintf(1,'Objective function = %10.3e\n',f);
fprintf(1,'Number of non-zero components = %d\n',num_nz_x);
fprintf(1,'CPU time so far = %10.3e\n', times(iter));
fprintf(1,'\n');


% If the 'Debias' option is set to 1, we try to
% remove the bias from the l1 penalty, by applying CG to the 
% least-squares problem obtained by omitting the l1 term 
% and fixing the zero coefficients at zero.
if debias
  fprintf(1,'\n')
  fprintf(1,'Starting the debiasing phase...\n\n')
  
  x_debias = x;
  zeroind = (x_debias~=0); 
  cont_debias_cg = 1;
  debias_start = iter;
  
  % calculate initial residual
  resid = A(x_debias);
  resid = resid-y;
  resid_prev = eps*ones(size(resid));
  
  rvec = AT(resid);
  
  % mask out the zeros
  rvec = rvec .* zeroind;
  rTr_cg = rvec(:)'*rvec(:);
  
  % set convergence threshold for the residual || RW x_debias - y ||_2
  tol_debias = tolD * (rvec(:)'*rvec(:));
  
  % initialize pvec
  pvec = -rvec;
  
  % main loop
  while cont_debias_cg
    
    % calculate A*p = Wt * Rt * R * W * pvec
    RWpvec = A(pvec);      
    Apvec = AT(RWpvec);
    
    % mask out the zero terms
    Apvec = Apvec .* zeroind;
    
    % calculate alpha for CG
    alpha_cg = rTr_cg / (pvec(:)'* Apvec(:));
    
    % take the step
    x_debias = x_debias + alpha_cg * pvec;
    resid = resid + alpha_cg * RWpvec;
    rvec  = rvec  + alpha_cg * Apvec;
    
    rTr_cg_plus = rvec(:)'*rvec(:);
    beta_cg = rTr_cg_plus / rTr_cg;
    pvec = -rvec + beta_cg * pvec;
    
    rTr_cg = rTr_cg_plus;
    
    iter = iter+1;
    
    objective(iter) = 0.5*(resid(:)'*resid(:)) + ...
	sum(tau(:).*abs(x_debias(:)));
    times(iter) = cputime - t0;
    
    if compute_mse
      err = true - x_debias;
      mses(iter) = (err(:)'*err(:));
    end
    
    % in the debiasing CG phase, always use convergence criterion
    % based on the residual (this is standard for CG)
    fprintf(1,' Iter = %5d, debias resid = %13.8e, convergence = %8.3e\n', ...
	iter, resid(:)'*resid(:), rTr_cg / tol_debias);
    cont_debias_cg = ...
	(iter-debias_start <= miniter_debias )| ...
	((rTr_cg > tol_debias) & ...
	(iter-debias_start <= maxiter_debias));
        
  end
  fprintf(1,'\nFinished the debiasing phase!\nResults:\n')
  fprintf(1,'||A x - y ||_2 = %10.3e\n',resid(:)'*resid(:))
  fprintf(1,'||x||_1 = %10.3e\n',sum(abs(x(:))))
  fprintf(1,'Objective function = %10.3e\n',f);
  nz = (x_debias~=0.0);
  fprintf(1,'Number of non-zero components = %d\n',sum(nz(:)));
  fprintf(1,'CPU time so far = %10.3e\n', times(iter));
  fprintf(1,'\n');
end

if compute_mse
  mses = mses/length(true(:));
end

% define the soft threshold function, which is used above.
function y = soft(x,tau)
y = sign(x).*max(abs(x)-tau,0);



