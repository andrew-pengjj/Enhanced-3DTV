% tveq_logbarrier.m
%
% Solve equality constrained TV minimization
% min TV(x)  s.t.  Ax=b.
%
% Recast as the SOCP
% min sum(t) s.t.  ||D_{ij}x||_2 <= t,  i,j=1,...,n
%                  Ax=b
% and use a log barrier algorithm.
%
% Usage:  xp = tveq_logbarrier(x0, A, At, b, lbtol, mu, slqtol, slqmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a K 
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
%
% b - Kx1 vector of observations.
%
% lbtol - The log barrier algorithm terminates when the duality gap <= lbtol.
%         Also, the number of log barrier iterations is completely
%         determined by lbtol.
%         Default = 1e-3.
%
% mu - Factor by which to increase the barrier constant at each iteration.
%      Default = 10.
%
% slqtol - Tolerance for SYMMLQ; ignored if A is a matrix.
%     Default = 1e-8.
%
% slqmaxiter - Maximum number of iterations for SYMMLQ; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = tveq_logbarrier(x0, A, At, b, lbtol, mu, slqtol, slqmaxiter)  

largescale = isa(A,'function_handle'); 

if (nargin < 5), lbtol = 1e-3; end
if (nargin < 6), mu = 10; end
if (nargin < 7), slqtol = 1e-8; end
if (nargin < 8), slqmaxiter = 200; end

newtontol = lbtol;
newtonmaxiter = 50;

N = length(x0);
n = round(sqrt(N));

% create (sparse) differencing matrices for TV
Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);

% starting point --- make sure that it is feasible
if (largescale)
  if (norm(A(x0)-b)/norm(b) > slqtol)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    AAt = @(z) A(At(z));
    [w,cgres] = cgsolve(AAt, b, slqtol, slqmaxiter, 0);
    if (cgres > 1/2)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = At(w);
  end
else
  if (norm(A*x0-b)/norm(b) > slqtol)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    opts.POSDEF = true; opts.SYM = true;
    [w, hcond] = linsolve(A*A', b, opts);
    if (hcond < 1e-14)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = A'*w;
  end  
end
x = x0;
Dhx = Dh*x;  Dvx = Dv*x;
t = (0.95)*sqrt(Dhx.^2 + Dvx.^2) + (0.1)*max(sqrt(Dhx.^2 + Dvx.^2));

% choose initial value of tau so that the duality gap after the first
% step will be about the origial TV
tau = N/sum(sqrt(Dhx.^2+Dvx.^2));

lbiter = ceil((log(N)-log(lbtol)-log(tau))/log(mu));
disp(sprintf('Number of log barrier iterations = %d\n', lbiter));
totaliter = 0;
for ii = 1:lbiter
  
  [xp, tp, ntiter] = tveq_newton(x, t, A, At, b, tau, newtontol, newtonmaxiter, slqtol, slqmaxiter);
  totaliter = totaliter + ntiter;
  
  tvxp = sum(sqrt((Dh*xp).^2 + (Dv*xp).^2));
  disp(sprintf('\nLog barrier iter = %d, TV = %.3f, functional = %8.3f, tau = %8.3e, total newton iter = %d\n', ...
    ii, tvxp, sum(tp), tau, totaliter));
  
  x = xp;
  t = tp;
  
  tau = mu*tau;
  
end
                   