% tveq_newton.m
%
% Newton algorithm for log-barrier subproblems for TV minimization
% with equality constraints.
%
% Usage: 
% [xp,tp,niter] = tveq_newton(x0, t0, A, At, b, tau, 
%                             newtontol, newtonmaxiter, slqtol, slqmaxiter)
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
% tau - Log barrier parameter.
%
% newtontol - Terminate when the Newton decrement is <= newtontol.
%
% newtonmaxiter - Maximum number of iterations.
%
% slqtol - Tolerance for SYMMLQ; ignored if A is a matrix.
%
% slqmaxiter - Maximum number of iterations for SYMMLQ; ignored
%     if A is a matrix.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function [xp, tp, niter] = tveq_newton(x0, t0, A, At, b, tau, newtontol, newtonmaxiter, slqtol, slqmaxiter) 

largescale = isa(A,'function_handle'); 

alpha = 0.01;
beta = 0.5;  

N = length(x0);
n = round(sqrt(N));
K = length(b);

% create (sparse) differencing matrices for TV
Dv = spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
Dh = spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);

% auxillary matrices for preconditioning
Mdv = spdiags([reshape([ones(n-1,n); zeros(1,n)],N,1) ...
  reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
Mdh = spdiags([reshape([ones(n,n-1) zeros(n,1)],N,1) ...
  reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);
Mmd = reshape([ones(n-1,n-1) zeros(n-1,1); zeros(1,n)],N,1);


% initial point
x = x0;
t = t0;
Dhx = Dh*x;  Dvx = Dv*x;
ft = 1/2*(Dhx.^2 + Dvx.^2 - t.^2);
f = sum(t) - (1/tau)*(sum(log(-ft)));

niter = 0;
done = 0;
while (~done)
  
  ntgx = Dh'*((1./ft).*Dhx) + Dv'*((1./ft).*Dvx);
  ntgt = -tau - t./ft;
  gradf = -(1/tau)*[ntgx; ntgt];
   
  sig22 = 1./ft + (t.^2)./(ft.^2);
  sig12 = -t./ft.^2;
  sigb = 1./ft.^2 - (sig12.^2)./sig22;
  
  w1p = ntgx - Dh'*(Dhx.*(sig12./sig22).*ntgt) - Dv'*(Dvx.*(sig12./sig22).*ntgt);
  wp = [w1p; zeros(K,1)];
  if (largescale)
    % diagonal of H11p
    dg11p = Mdh'*(-1./ft + sigb.*Dhx.^2) + Mdv'*(-1./ft + sigb.*Dvx.^2) + 2*Mmd.*sigb.*Dhx.*Dvx;
    afac = max(dg11p);
    hpfun = @(z) Hpeval(z, A, At, Dh, Dv, Dhx, Dvx, sigb, ft, afac);
    [dxv,slqflag,slqres,slqiter] = symmlq(hpfun, wp, slqtol, slqmaxiter);
    if (slqres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
  else
    H11p =  Dh'*sparse(diag(-1./ft + sigb.*Dhx.^2))*Dh + ...
      Dv'*sparse(diag(-1./ft + sigb.*Dvx.^2))*Dv + ...
      Dh'*sparse(diag(sigb.*Dhx.*Dvx))*Dv + ...
      Dv'*sparse(diag(sigb.*Dhx.*Dvx))*Dh;
    afac = max(diag(H11p));
    Hp = full([H11p afac*A'; afac*A zeros(K)]);
    %keyboard
    opts.SYM = true;
    [dxv, hcond] = linsolve(Hp, wp, opts);
    if (hcond < 1e-14)
      disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;  tp = t;
      return
    end
  end
  dx = dxv(1:N);  
  Dhdx = Dh*dx;  Dvdx = Dv*dx;
  dt = (1./sig22).*(ntgt - sig12.*(Dhx.*Dhdx + Dvx.*Dvdx));
  
  % minimum step size that stays in the interior
  aqt = Dhdx.^2 + Dvdx.^2 - dt.^2;   
  bqt = 2*(Dhdx.*Dhx + Dvdx.*Dvx - t.*dt);  
  cqt = Dhx.^2 + Dvx.^2 - t.^2;
  tsols = [(-bqt+sqrt(bqt.^2-4*aqt.*cqt))./(2*aqt); ...
    (-bqt-sqrt(bqt.^2-4*aqt.*cqt))./(2*aqt) ];
  indt = find([(bqt.^2 > 4*aqt.*cqt); (bqt.^2 > 4*aqt.*cqt)] & (tsols > 0));
  smax = min(1, min(tsols(indt))); 
  s = (0.99)*smax;
  
  % line search
  suffdec = 0;
  backiter = 0;
  while (~suffdec)
    xp = x + s*dx;  tp = t + s*dt;
    Dhxp = Dhx + s*Dhdx;  Dvxp = Dvx + s*Dvdx;
    ftp = 1/2*(Dhxp.^2 + Dvxp.^2 - tp.^2);
    fp = sum(tp) - (1/tau)*(sum(log(-ftp)));
    flin = f + alpha*s*(gradf'*[dx; dt]);
    suffdec = (fp <= flin);
    s = beta*s;
    backiter = backiter + 1;
    if (backiter > 32)
      disp('Stuck backtracking, returning last iterate. (See Section 4 of notes for more information.)');
      xp = x;  tp = t;
      return
    end
  end
  
  % set up for next iteration
  x = xp; t = tp;
  Dvx = Dvxp;  Dhx = Dhxp; 
  ft = ftp; f = fp;
  
  lambda2 = -(gradf'*[dx; dt]);
  stepsize = s*norm([dx; dt]);
  niter = niter + 1;
  done = (lambda2/2 < newtontol) | (niter >= newtonmaxiter);

  disp(sprintf('Newton iter = %d, Functional = %8.3f, Newton decrement = %8.3f, Stepsize = %8.3e', ...
    niter, f, lambda2/2, stepsize));
  if (largescale)
    disp(sprintf('                  SYMMLQ Res = %8.3e, SYMMLQ Iter = %d', slqres, slqiter));
  else
    disp(sprintf('                  H11p condition number = %8.3e', hcond));
  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implicit application of Hessian
function y = Hpeval(z, A, At, Dh, Dv, Dhx, Dvx, sigb, ft, afac)

N = length(ft);
K = length(z)-N;
w = z(1:N);
v = z(N+1:N+K);

Dhw = Dh*w;
Dvw = Dv*w;

y1 = Dh'*((-1./ft + sigb.*Dhx.^2).*Dhw + sigb.*Dhx.*Dvx.*Dvw) + ...
  Dv'*((-1./ft + sigb.*Dvx.^2).*Dvw + sigb.*Dhx.*Dvx.*Dhw) + afac*At(v);
y2 = afac*A(w);

y = [y1; y2];
