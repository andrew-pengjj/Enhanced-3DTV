function [A_hat E_hat iter] = inexact_alm_rpca_wsnm(D, C, lambda, delta, p, tol, maxIter)

% Jun 2015
% This matlab code implements the inexact augmented Lagrange multiplier 
% method for Robust PCA, with the following objective function:
%
% Initialize A,E,N,\Gamma,\beta
% while ~converged 
%   minimize (inexactly, update A, E and N only once)
%     L(A,E,N,Y,\beta) = C*|A|_{w,p}^{p} + \lambda * |E|_1 - <Y, A+E+N-D> + \beta/2 * |A+E+N-D|_F^2;
%     Y = Y - \beta * (A + E + N - D);
%     \beta = \rho * \beta;
% end
% D - m x n matrix of observations/data (required input)
%
% lambda - weight on sparse error term in the cost function
%
% delta - stand for the Gaussian noise level 
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
% 
%
% Yuan Xie, Jun 2015. Modified mingming chen's RPCA code; 

% addpath PROPACK;
[m n] = size(D);

if nargin < 2
    lambda = 1 / sqrt(m);
end

if nargin < 3
    tol = 1e-7;
elseif tol == -1
    tol = 1e-7;
end

if nargin < 4
    maxIter = 1000;
elseif maxIter == -1
    maxIter = 1000;
end

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros(m, n);
E_hat = zeros(m, n);
N_hat = zeros(m, n);
beta = 1.25/norm_two; % this one can be tuned
beta_bar = beta * 1e7;
rho = 1.2;          % this one can be tuned
d_norm = norm(D, 'fro');

iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 5;

while ~converged       
    iter = iter + 1;
    
    % step 1: Solving N
    temp_Z = (1/beta)*Y + D - A_hat - E_hat;
    norm_temp_Z = norm(temp_Z, 'fro');
    N_hat  = (min(norm_temp_Z, delta)/norm_temp_Z).*temp_Z;
    
    % step 2: Solving E
    temp_T = D - A_hat - N_hat + (1/beta)*Y;
    E_hat = max(temp_T - lambda/beta, 0);
    E_hat = E_hat+min(temp_T + lambda/beta, 0);

    % step 3: Solving A
    if choosvd(n, sv) == 1
        [U S V] = lansvd(D - N_hat - E_hat + (1/beta)*Y, sv, 'L');
    else
        [U S V] = svd(D - N_hat - E_hat + (1/beta)*Y, 'econ');
    end
    
    % step 3a: WSNM via GST
    %%%%%%%%%%%%%
    diagS = diag(S);
    [tempDiagS,svp]=IterativeWSNM(diagS,C*m*sqrt(n)/beta, p);
%     [tempDiagS,svp]=IterativeWSNM(diagS,C*sqrt(m*n)/beta, p);
    A_hat = U(:,1:svp)*diag(tempDiagS(1:svp))*V(:,1:svp)';  
    %%%%%%%%%%%%
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end   
    total_svd = total_svd + 1;
    
    % step 4: Updating Y
    Z = A_hat + E_hat + N_hat - D;  
    Y = Y - beta*Z;
    
    % step 5: Increasing \beta
    beta = min(beta*rho, beta_bar);
        
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end    
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

