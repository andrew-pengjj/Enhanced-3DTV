function [L S err] = sparcs(b, r, K, B, sizeD,svdMethod, tol, maxIter, VERBOSE,tol_early)
%%%% 
%Solves the problem
% min | b - A(L+S) |_2    s.t   rank(L) = r, |S|_0 = K
%
% At = Adjoint of A
% svdMethod = {'svds', 'svdlibc','propack'} SVD solvers
% tol, maxIter are used for termination
% verb = verbose. 1 - on, 0 - off
% tol_early.  Terminate early when residual change is less than tol_early
%%%%%%%
A =@(x) B*x(:);
At =@(x) reshape(B'*x(:),sizeD(1),sizeD(2));

if (nargin < 5)
    disp('Not enough arguments');
    L = []; S = []; err=  [];
    return;
end

if (nargin < 6)     svdMethod = 'svds'; end
if (nargin < 7)     tol = 5e-4;         end
if (nargin < 8)     maxIter = 7*(r+1);  end
if (nargin < 9)     VERBOSE = 1;           end
if (nargin < 10)    tol_early = -Inf;  end

%Declare greedy low rank vars
temp = At(b);
rows = size(temp,1); cols = size(temp,2);
N = rows*cols;
Xhat = zeros(rows,cols);
PsihatU = []; PsihatV = [];
err = 1; 

%Declare sparse matrix vars
aa= zeros(N,1); 
s_cosamp = zeros(N,1);
selectAtom = 1;


for iterCnt = 1 : maxIter

    
    %%Do greedy matrix steps
    %Compute new basis for residual, store in Psiprime
    rt = At(b-A(Xhat)-A(s_cosamp));
    switch lower(svdMethod)
        case 'svdlibc'
            [U diagS V] = svdlibc(rt, selectAtom*r);
        case 'propack'
            [U,S,V] = lansvd(rt,selectAtom*r,'L');
            diagS = diag(S);
        case 'svds'
            [U, S, V] = svds(rt, selectAtom*r);
            diagS = diag(S);
        case 'spsvd'
            [U, S, V] = spsvd(rt, selectAtom*r, 1e-2);
            diagS = diag(S);
        otherwise
            [U,S,V] = svd(rt,0);
            diagS = diag(S);
    end
    
    Uprime = U(:, 1:min(selectAtom*r,size(U,2)));
    Vprime = V(:, 1:min(selectAtom*r,size(V,2)));
    
    
    %Merge supports
    PsitildeU = [PsihatU Uprime];
    PsitildeV = [PsihatV Vprime];
    
    %Do least squares estimation to get Xhat
    AP = @(z) APsiUV(z,A,PsitildeU, PsitildeV);
    APt = @(z) APsitUV(z,At,PsitildeU, PsitildeV);
    ALS = @(z) APt(AP(z));
    [alpha, res, iter] = cgsolve(ALS, APt(b-A(s_cosamp)), 1e-6, 100, 0);
    Xtilde = PsitildeU*diag(alpha)*PsitildeV'; 
    
    %Update Psihat
    switch lower(svdMethod)
        case 'svdlibc'
            [U diagS V] = svdlibc(Xtilde, r);
        case 'propack'
            [U,S,V] = lansvd(Xtilde,r,'L');
            diagS = diag(S);
        case 'svds'
            [U, S, V] = svds(Xtilde, r, 'L');
            diagS = diag(S);
        otherwise
            [U,S,V] = svd(Xtilde,0);
            diagS = diag(S);
    end
    PsihatU = U; PsihatV = V;
    Xhat = U*S*V';
    
    %%Do sparse matrix steps
    rr = b - A(s_cosamp) - A(Xhat);
    proxy = At(rr);proxy = proxy(:);
    %---Estimate support
    [tmp,ww]= sort(abs(proxy),'descend');
    tt= union(find(ne(s_cosamp,0)),ww(1:(selectAtom*K)));
    
    % Preparation for cg_solve
    PP_tt = @(z) A_I(A,z,tt,N);
    PP_transpose_tt = @(z) A_I_transpose(At,z,tt);
    qq = PP_transpose_tt(b-A(Xhat));
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    %Pseudo-inverse
    [w, res, iter2] = cgsolve(PPtranspose_PP_tt, qq,1e-6, 100, 0);
    bb= 0*s_cosamp; bb(tt)= w;
    %---Prune
    kk = iterCnt;
    [tmp,ww2]= sort(abs(bb),'descend');
    s_cosamp=0*bb;
    s_cosamp(ww2(1:K)) = bb(ww2(1:K));
    aa = s_cosamp;
    
    errTemp = norm(A(Xhat)+A(s_cosamp)-b)/norm(b);
    if (VERBOSE)
        fprintf('Iter: %d Err: %f\n',iterCnt,errTemp);
    end
    
    if ((err-errTemp)/err <= tol_early) %Early termination condition
        err = errTemp;
        if (VERBOSE)
            disp(sprintf('Terminating.'));
        end
        break 
    end
    
   err = errTemp;
   
   if err<tol && iter >3
        break;
   end
   
end

%Get final data and return
L = Xhat;
xcosamp= aa;
S = reshape(xcosamp,size(Xhat,1),size(Xhat,2));
if (VERBOSE)
	disp(sprintf('Final error: %3.6f', errTemp));
end
