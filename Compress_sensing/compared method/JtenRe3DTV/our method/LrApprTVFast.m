function [x,psnrpath] = LrApprTVFast(A, y, sizeD, rk, lam, opts)

%%  check and initilization
if ~exist('opts','var'); opts=[]; end
if isfield(opts,'maxIter'); maxIter = opts.maxIter; else maxIter = 250; end
if isfield(opts,'tol');  tol = opts.tol; else tol = 1e-6; end;
if isfield(opts,'beta'); beta = opts.beta; else beta = ones(3,1)*1e-5/mean(abs(y(:))); beta(2)=beta(2)*10;end
if isfield(opts,'trX');   trX = opts.trX; end;
if isfield(opts,'isCont')
    isCont = opts.isCont; contpar = opts.contpar; contfactor = opts.contfactor;
 else
    isCont = true; contpar = 0.95; contfactor = 1.15;
end
disp(opts);


% init. variables 
N    = prod(sizeD);
h    = sizeD(1);
w    = sizeD(2);
s    = sizeD(3);

x   = zeros(N,1);
p    = zeros(3*N, 1);
v1    = zeros(3*N, 1);
v2    = zeros(size(y));
e = zeros(N,1);
v3   = e;

dimRank    = numel(rk);
if dimRank == 1
    Aty        = A'*y;
    Aty_mat    = reshape(Aty,h*w,s);
    [U,Sig,V]  = MySVD(Aty_mat);
    U          = U(:,1:rk);
    V          = (V(:,1:rk))';
    L          = U*Sig(1:rk,1:rk)*V;
    L=L(:);
else
    Aty        = A'*y;
    Atb_ten    = reshape(Aty,sizeD);
    Ten        = tucker_als(tensor(Atb_ten),rk,'tol',1e-6,'printitn',0);
    L          = double(Ten);
    L=L(:);
end


%% main loop
stopPath     = zeros(maxIter, 1);
gamma    = 1.15;
display  = 1;                                                                                                                                                                                                                                                                                                                                                                                                               
psnrpath      = zeros(maxIter,1);
for iter = 1 : maxIter
    
    fprintf('\n*****************************iter: %d ******************************\n', iter'); 
    
    % x subproblem
    x= myPCG(A,y,L,x,p,v1,v2,v3,e,beta,sizeD);
%     numer=reshape((v3+beta(3)*(e+L)-diffT3(v1,sizeD)+beta(1)*diffT3(p,sizeD)-A'*v2+beta(2)*(A'*y)),sizeD);
%     x  = real(ifftn( fftn(numer) ./ (beta(1)*denom + beta(2)+beta(3)) ) );
%     x  = x(:);
    
    %e subproblem
    e = (beta(3)*(x-L)-v3)/(1+beta(3));
    
    %- [G;U1,U2,U3] - subproblem
    L = x-e-v3/beta(3);
    if numel(rk) == 1 %% low rank matarix     
        L = reshape(L, h*w, s);
        for i_u =  1 : 2
         U = L * V';
         V = pinv(U) * L;
        end
        L = U*V;
        L = L(:);
    else %% low tensor      
        L = reshape(L, sizeD);
        L  = double( tucker_als(tensor(L), rk, 'tol',1e-6,'printitn',0) );
        L  = L(:);
    end
    %- p subproblem
    
    diff_x = diff3(x, sizeD); 
    p      = softthre( diff_x + v1/beta(1), lam/beta(1) );
     
    %- updating multipliers
	v1 = v1 - gamma*beta(1)*(p - diff_x);
	v2 = v2 - gamma*beta(2)*(y - A*x);
    v3 = v3 - gamma*beta(3)*(x - L -e);
    
    if exist('trX','var')
        xten=reshape(x,sizeD);
        psnrpath(iter) = PSNR(trX,xten);
        if display
            fprintf('psnr :%4.4f\n', psnrpath(iter));
        end
    end
    %- terminating the algorithm and saving some ralted information 
    stopCond = norm(A*x - y)/norm(y);
    stopPath(iter) = stopCond;
    fprintf('stopFormulaVal: %4.4e \n', stopCond);
    if (iter> 50) &&  (stopCond < tol ) 
          disp(' !!!stopped by termination rule!!! ');  break;
    end
    
    %- continuation
      if  isCont
            nr1 = norm(p - diff_x, 'fro');
            nr2 = norm(y - A*x, 'fro');
            nr3 = norm(x-e-L,'fro');
            if iter >1 && nr1 > contpar * nr1_pre
                beta(1) = contfactor*beta(1);
            end
            if iter>1 && nr2 > contpar * nr2_pre
                beta(2) = contfactor*beta(2);
            end
            if iter>1 && nr3 > contpar * nr3_pre
            beta(3) = contfactor*beta(3);
            end
            nr1_pre =nr1;    nr2_pre = nr2;  nr3_pre = nr3; 
      end    
    
end

end





