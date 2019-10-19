% tvdantzig_logbarrier.m
%
% Solve the total variation Dantzig program
%
% min_x TV(x)  subject to  ||A'(Ax-b)||_\infty <= epsilon
%
% Recast as the SOCP
% min sum(t) s.t.  ||D_{ij}x||_2 <= t,  i,j=1,...,n
%                  <a_{ij},Ax - b> <= epsilon  i,j=1,...,n
% and use a log barrier algorithm.
%
% Usage:  xp = tvdantzig_logbarrier(x0, A, At, b, epsilon, lbtol, mu, cgtol, cgmaxiter)
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
% epsilon - scalar, constraint relaxation parameter
%
% lbtol - The log barrier algorithm terminates when the duality gap <= lbtol.
%         Also, the number of log barrier iterations is completely
%         determined by lbtol.
%         Default = 1e-3.
%
% mu - Factor by which to increase the barrier constant at each iteration.
%      Default = 10.
%
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%     Default = 1e-8.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = tvdantzig_logbarrier(x0, A, At, b, epsilon, lbtol, mu, cgtol, cgmaxiter)  

largescale = isa(A,'function_handle');

if (nargin < 6), lbtol = 1e-3; end
if (nargin < 7), mu = 10; end
if (nargin < 8), cgtol = 1e-8; end
if (nargin < 9), cgmaxiter = 200; end

newtontol = lbtol;
newtonmaxiter = 50;

N = length(x0);
n = round(sqrt(N));

% create (sparse) differencing matrices for TV
Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);

if (largescale)
  if (norm(A(x0)-b) > epsilon)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    AAt = @(z) A(At(z));
    [w, cgres] = cgsolve(AAt, b, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = At(w);
  end
else
  if (norm(A*x0-b) > epsilon)
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
t = 1.05*sqrt(Dhx.^2 + Dvx.^2) + .01*max(sqrt(Dhx.^2 + Dvx.^2));

% choose initial value of tau so that the duality gap after the first
% step will be about the origial TV
tau = 3*N/sum(sqrt(Dhx.^2+Dvx.^2));

lbiter = ceil((log(3*N)-log(lbtol)-log(tau))/log(mu));
disp(sprintf('Number of log barrier iterations = %d\n', lbiter));
totaliter = 0;
for ii = 1:lbiter
  
  [xp, tp, ntiter] = tvdantzig_newton(x, t, A, At, b, epsilon, tau, newtontol, newtonmaxiter, cgtol, cgmaxiter);
  totaliter = totaliter + ntiter;
  tvxp = sum(sqrt((Dh*xp).^2 + (Dv*xp).^2));
  
  disp(sprintf('\nLog barrier iter = %d, TV = %.3f, functional = %8.3f, tau = %8.3e, total newton iter = %d\n', ...
    ii, tvxp, sum(tp), tau, totaliter));
  
  x = xp;
  t = tp;
  
  tau = mu*tau;
  
end
                   