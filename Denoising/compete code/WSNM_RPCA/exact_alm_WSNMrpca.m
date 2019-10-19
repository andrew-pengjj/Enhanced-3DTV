function [A_hat E_hat iter] = exact_alm_WSNMrpca(D, C, p, tol, maxIter)
%加权核范数的用来精确RPCA的方法，没有迭代的reweighted情况
%用精确ALM方法求解下面的问题
% Objective function:
% min |A|^{p}_w,p + |E|_1, s.t. A + E = D
% w_i = C*sqrt(m*n)/(SigmaX_i + eps);
% while ~converged 
%   minimize
%     L(A,E,Y,u) = C*|A|_* + |E|_1 + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end


addpath PROPACK;

[m n] = size(D);
p = 0.7;


if nargin < 4
    tol = 1e-6;
elseif tol == -1
    tol = 1e-6;
end

if nargin < 5
    maxIter = 100;
elseif maxIter == -1
    maxIter = 100;
end

% initialize
Y = sign(D);
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf);
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
E_hat = zeros( m, n);
dnorm = norm(D, 'fro');
tolProj = 1e-6 * dnorm;
total_svd = 0;
mu = 0.5/norm_two*sqrt(m) ;% this one can be tuned 0.5
%mu = .5/sqrt(m) ;% this one can be tuned 0.5
rho = 5 ;         % this one can be tuned

iter = 0;
converged = false;
stopCriterion = 1;
sv = 5;
svp = sv;
while ~converged       
    iter = iter + 1;
    
    % solve the primal problem by alternative projection
    primal_converged = false;
    primal_iter = 0;
    sv = sv + round(n * 0.1);
    while primal_converged == false
        
        temp_T = D - A_hat + (1/mu)*Y;
        temp_E = max( temp_T - 1/mu,0) + min( temp_T + 1/mu,0); %这两步是更新E
        
        if choosvd(n, sv) == 1%这是一个节约时间的技巧吧。。
            [U S V] = lansvd(D - temp_E + (1/mu)*Y, sv, 'L');
        else
            [U S V] = svd(D - temp_E + (1/mu)*Y, 'econ');
        end
        
        %%%%%%%%%%%%%%%%%%
        diagS = diag(S);
        [tempDiagS,svp]=IterativeWSNM(diagS,C*sqrt(m*n)/mu,p);
        %[tempDiagS,svp]=IterativeWSNM(diagS,C*sqrt(m)/mu,p);
%         sv
%         svp
        temp_A=U(:,1:svp)*diag(tempDiagS(1:svp))*V(:,1:svp)';          %这里更新A
        %%%%%%%%%%%%%%%%%%
        
        if svp < sv
            sv = min(svp + 1, n); %1
        else
            sv = min(svp + round(0.05*n), n);         
        end


        
        if norm(A_hat - temp_A, 'fro') < tolProj && norm(E_hat - temp_E, 'fro') < tolProj
            primal_converged = true;
        end
        A_hat = temp_A;
        E_hat = temp_E;
        primal_iter = primal_iter + 1;
        total_svd = total_svd + 1;
               
    end
        
    Z = D - A_hat - E_hat;        
    Y = Y + mu*Z;%更新Y
    mu = rho * mu;
    
    %% stop Criterion    
    stopCriterion = norm(Z, 'fro') / dnorm;
    if stopCriterion < tol
        converged = true;
    end    
    
    disp(['Iteration' num2str(iter) ' #svd ' num2str(total_svd) ' r(A) ' num2str(svp)...
        ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
        ' stopCriterion ' num2str(stopCriterion)]);
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
    end
end

if nargin == 5
    fclose(fid);
end


