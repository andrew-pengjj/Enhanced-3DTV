% tvqc_newton.m
%
% Newton algorithm for log-barrier subproblems for TV minimization
% with quadratic constraints.
%
% Usage: 
% [xp,tp,niter] = tvqc_newton(x0, t0, A, At, b, epsilon, tau, 
%                             newtontol, newtonmaxiter, cgtol, cgmaxiter)
%
% x0,t0 - starting points
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
% tau - Log barrier parameter.
%
% newtontol - Terminate when the Newton decrement is <= newtontol.
%
% newtonmaxiter - Maximum number of iterations.
%
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function [xp, tp, niter] = tvqc_newton(x0, t0, A, At, b, epsilon, tau, newtontol, newtonmaxiter, cgtol, cgmaxiter) 

largescale = isa(A,'function_handle'); 

alpha = 0.01;
beta = 0.5;  

N = length(x0);
n = round(sqrt(N));

% create (sparse) differencing matrices for TV
Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);

if (~largescale),  AtA = A'*A;  end;

% initial point
x = x0;
t = t0;
if (largescale), r = A(x) - b;  else  r = A*x - b; end  
Dhx = Dh*x;  Dvx = Dv*x;
ft = 1/2*(Dhx.^2 + Dvx.^2 - t.^2);
fe = 1/2*(r'*r - epsilon^2);
f = sum(t) - (1/tau)*(sum(log(-ft)) + log(-fe));

niter = 0;
done = 0;
while (~done)
  
  if (largescale),  Atr = At(r);  else  Atr = A'*r;  end
  ntgx = Dh'*((1./ft).*Dhx) + Dv'*((1./ft).*Dvx) + 1/fe*Atr;
  ntgt = -tau - t./ft;
  gradf = -(1/tau)*[ntgx; ntgt];
  
  sig22 = 1./ft + (t.^2)./(ft.^2);
  sig12 = -t./ft.^2;
  sigb = 1./ft.^2 - (sig12.^2)./sig22;
  
  w1p = ntgx - Dh'*(Dhx.*(sig12./sig22).*ntgt) - Dv'*(Dvx.*(sig12./sig22).*ntgt);
  if (largescale)
    h11pfun = @(z) H11p(z, A, At, Dh, Dv, Dhx, Dvx, sigb, ft, fe, Atr);
    [dx, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;  tp = t;
      return
    end
    Adx = A(dx);
  else
    H11p =  Dh'*sparse(diag(-1./ft + sigb.*Dhx.^2))*Dh + ...
      Dv'*sparse(diag(-1./ft + sigb.*Dvx.^2))*Dv + ...
      Dh'*sparse(diag(sigb.*Dhx.*Dvx))*Dv + ...
      Dv'*sparse(diag(sigb.*Dhx.*Dvx))*Dh - ...
      (1/fe)*AtA + (1/fe^2)*Atr*Atr';
    opts.POSDEF = true; opts.SYM = true;
    [dx,hcond] = linsolve(H11p, w1p, opts);
    if (hcond < 1e-14)
      disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;  tp = t;
      return
    end
    Adx = A*dx;
  end
  Dhdx = Dh*dx;  Dvdx = Dv*dx;
  dt = (1./sig22).*(ntgt - sig12.*(Dhx.*Dhdx + Dvx.*Dvdx));
  
  % minimum step size that stays in the interior
  aqt = Dhdx.^2 + Dvdx.^2 - dt.^2;   
  bqt = 2*(Dhdx.*Dhx + Dvdx.*Dvx - t.*dt);  
  cqt = Dhx.^2 + Dvx.^2 - t.^2;
  tsols = [(-bqt+sqrt(bqt.^2-4*aqt.*cqt))./(2*aqt); ...
    (-bqt-sqrt(bqt.^2-4*aqt.*cqt))./(2*aqt) ];
  indt = find([(bqt.^2 > 4*aqt.*cqt); (bqt.^2 > 4*aqt.*cqt)] & (tsols > 0));
  aqe = Adx'*Adx;   bqe = 2*r'*Adx;   cqe = r'*r - epsilon^2;
  smax = min(1,min([...
    tsols(indt); ...
    (-bqe+sqrt(bqe^2-4*aqe*cqe))/(2*aqe)
    ]));
  s = (0.99)*smax;
  
  % backtracking line search
  suffdec = 0;
  backiter = 0;
  while (~suffdec)
    xp = x + s*dx;  tp = t + s*dt;
    rp = r + s*Adx;  Dhxp = Dhx + s*Dhdx;  Dvxp = Dvx + s*Dvdx;
    ftp = 1/2*(Dhxp.^2 + Dvxp.^2 - tp.^2);
    fep = 1/2*(rp'*rp - epsilon^2);
    fp = sum(tp) - (1/tau)*(sum(log(-ftp)) + log(-fep));
    flin = f + alpha*s*(gradf'*[dx; dt]);
    suffdec = (fp <= flin);
    s = beta*s;
    backiter = backiter + 1;
    if (backiter > 32)
      disp('Stuck on backtracking line search, returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;  tp = t;
      return
    end
  end
  
  % set up for next iteration
  x = xp; t = tp;
  r = rp;  Dvx = Dvxp;  Dhx = Dhxp; 
  ft = ftp; fe = fep; f = fp;
  
  lambda2 = -(gradf'*[dx; dt]);
  stepsize = s*norm([dx; dt]);
  niter = niter + 1;
  done = (lambda2/2 < newtontol) | (niter >= newtonmaxiter);
  
  disp(sprintf('Newton iter = %d, Functional = %8.3f, Newton decrement = %8.3f, Stepsize = %8.3e', ...
    niter, f, lambda2/2, stepsize));
  if (largescale)
    disp(sprintf('                  CG Res = %8.3e, CG Iter = %d', cgres, cgiter));
  else
    disp(sprintf('                  H11p condition number = %8.3e', hcond));
  end
 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H11p auxiliary function
function y = H11p(v, A, At, Dh, Dv, Dhx, Dvx, sigb, ft, fe, atr)

Dhv = Dh*v;
Dvv = Dv*v;

y = Dh'*((-1./ft + sigb.*Dhx.^2).*Dhv + sigb.*Dhx.*Dvx.*Dvv) + ...
  Dv'*((-1./ft + sigb.*Dvx.^2).*Dvv + sigb.*Dhx.*Dvx.*Dhv) - ...
  1/fe*At(A(v)) + 1/fe^2*(atr'*v)*atr;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
