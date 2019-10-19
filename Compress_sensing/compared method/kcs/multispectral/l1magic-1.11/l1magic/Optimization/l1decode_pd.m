% l1decode_pd.m
%
% Decoding via linear programming.
% Solve
% min_x  ||b-Ax||_1 .
%
% Recast as the linear program
% min_{x,u} sum(u)  s.t.  -Ax - u + y <= 0
%                          Ax - u - y <= 0
% and solve using primal-dual interior point method.
%
% Usage: xp = l1decode_pd(x0, A, At, y, pdtol, pdmaxiter, cgtol, cgmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a M 
%     vector, or a MxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes an M vector and returns an N vector.
%      If A is a matrix, At is ignored.
%
% y - Mx1 observed code (M > N).
%
% pdtol - Tolerance for primal-dual algorithm (algorithm terminates if
%     the duality gap is less than pdtol).  
%     Default = 1e-3.
%
% pdmaxiter - Maximum number of primal-dual iterations.  
%     Default = 50.
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

function xp = l1decode_pd(x0, A, At, y, pdtol, pdmaxiter, cgtol, cgmaxiter)

largescale = isa(A,'function_handle');   

if (nargin < 5), pdtol = 1e-3; end
if (nargin < 6), pdmaxiter = 50; end
if (nargin < 7), cgtol = 1e-8; end
if (nargin < 8), cgmaxiter = 200; end

N = length(x0);
M = length(y);

alpha = 0.01;
beta = 0.5;
mu = 10;

gradf0 = [zeros(N,1); ones(M,1)];

x = x0;
if (largescale),  Ax = A(x);  else  Ax = A*x;  end
u = (0.95)*abs(y-Ax) + (0.10)*max(abs(y-Ax));

fu1 = Ax - y - u;
fu2 = -Ax + y - u;

lamu1 = -1./fu1;
lamu2 = -1./fu2;

if (largescale), Atv = At(lamu1-lamu2);  else Atv = A'*(lamu1-lamu2);  end

sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*M/sdg;

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [Atv; -lamu1-lamu2];
resnorm = norm([rdual; rcent]);

pditer = 0;
done = (sdg < pdtol)| (pditer >= pdmaxiter);
while (~done)
  
  pditer = pditer + 1;
  
  w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
  
  sig1 = -lamu1./fu1 - lamu2./fu2;
  sig2 = lamu1./fu1 - lamu2./fu2;
  sigx = sig1 - sig2.^2./sig1;
  
  if (largescale)
    w1 = -1/tau*(At(-1./fu1 + 1./fu2));
    w1p = w1 - At((sig2./sig1).*w2);
    h11pfun = @(z) At(sigx.*A(z));
    [dx, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0); 
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    Adx = A(dx);
  else
    w1 = -1/tau*(A'*(-1./fu1 + 1./fu2));
    w1p = w1 - A'*((sig2./sig1).*w2);
    H11p = A'*(sparse(diag(sigx))*A);
    opts.POSDEF = true; opts.SYM = true;
    [dx, hcond] = linsolve(H11p, w1p,opts);
    if (hcond < 1e-14)
      disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    Adx = A*dx;
  end
  
  du = (w2 - sig2.*Adx)./sig1;
  
  dlamu1 = -(lamu1./fu1).*(Adx-du) - lamu1 - (1/tau)*1./fu1;
  dlamu2 = (lamu2./fu2).*(Adx + du) -lamu2 - (1/tau)*1./fu2;
  if (largescale), Atdv = At(dlamu1-dlamu2);  else  Atdv = A'*(dlamu1-dlamu2); end
  
  % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
  indl = find(dlamu1 < 0);  indu = find(dlamu2 < 0);
  s = min([1; -lamu1(indl)./dlamu1(indl); -lamu2(indu)./dlamu2(indu)]);
  indl = find((Adx-du) > 0);  indu = find((-Adx-du) > 0);
  s = (0.99)*min([s; -fu1(indl)./(Adx(indl)-du(indl)); -fu2(indu)./(-Adx(indu)-du(indu))]);
  
  % backtrack 
  suffdec = 0;
  backiter = 0;
  while(~suffdec)
    xp = x + s*dx;  up = u + s*du;
    Axp = Ax + s*Adx;  Atvp = Atv + s*Atdv;
    lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
    fu1p = Axp - y - up;  fu2p = -Axp + y - up;
    rdp = gradf0 + [Atvp; -lamu1p-lamu2p];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    suffdec = (norm([rdp; rcp]) <= (1-alpha*s)*resnorm);
    s = beta*s;
    backiter = backiter + 1;
    if (backiter > 32)
      disp('Stuck backtracking, returning last iterate.  (See Section 4 of notes for more information.)')
      xp = x;
      return
    end
  end
  
  % next iteration
  x = xp;  u = up;
  Ax = Axp;  Atv = Atvp;
  lamu1 = lamu1p;  lamu2 = lamu2p;
  fu1 = fu1p;  fu2 = fu2p;
  
  % surrogate duality gap
  sdg = -(fu1'*lamu1 + fu2'*lamu2);
  tau = mu*2*M/sdg;
  rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
  rdual = rdp;
  resnorm = norm([rdual; rcent]);
  
  done = (sdg < pdtol) | (pditer >= pdmaxiter);
  
  disp(sprintf('Iteration = %d, tau = %8.3e, Primal = %8.3e, PDGap = %8.3e, Dual res = %8.3e',...
    pditer, tau, sum(u), sdg, norm(rdual)));
  if (largescale)
    disp(sprintf('                CG Res = %8.3e, CG Iter = %d', cgres, cgiter));
  else
    disp(sprintf('                  H11p condition number = %8.3e', hcond));
  end  
  
end

