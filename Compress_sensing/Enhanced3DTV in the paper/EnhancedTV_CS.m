%% the model algorithm
%      1/2*||e||_F^2 + \lambda*sum_i||U_i||_1
%         s.t. y=Ax, x = x_1+e, \Delta_i(x1) =U_i * V_i^T 

function [x,x1,e,psnrpath] = EnhancedTV_CS(A, y, sizeD, rank, weight, opts)

%%  check and initilization
if ~exist('opts','var')
    opts=[]; 
end

if isfield(opts,'maxIter')
    maxIter = opts.maxIter; 
else
    maxIter = 100; 
end

if isfield(opts,'tol')  
    tol = opts.tol; 
else
    tol = 1e-6; 
end

if isfield(opts,'beta')
    beta = opts.beta; 
else
    beta = ones(5,1)*1e-2/mean(abs(y)); 
    beta(4) = beta(4)*10;
end

if isfield(opts,'isCont')
    isCont = opts.isCont; 
    contpar = opts.contpar; 
    contfactor = opts.contfactor;
 else
    isCont = true; 
    contpar = 0.95; 
    contfactor = 1.15;
end
disp(opts);
%% initial variables 
h    = sizeD(1);
w    = sizeD(2);
s    = sizeD(3);
%% initial x1 x,and u,v
Aty        = A'*y;
Aty_mat    = reshape(Aty,h*w,s);
[U,Sig,V]  = MySVD(Aty_mat);
rk         = rank;
x1         = U(:,1:rk)*Sig(1:rk,1:rk)*V(:,1:rk)';
x1         = x1(:);
normD      = norm(Aty,'fro');
if isfield(opts,'init_x')
    x = opts.init_x;
else
    x = x1;
end
% L  = reshape(x1,sizeD);
% figure;
% imshow(L(:,:,1),[]);title('Clean image');
[mp,sm,er] = msqia(reshape(x1,[h,w,s])/255,opts.trX/255);
fprintf('mp = %f, sm = %f, er =%f.\n',mp,sm,er);

e   = zeros(size(x));
tv_x           = diff_x(x1,sizeD);
tv_x           = reshape(tv_x,[h*w,s]);
[U_x,S_x,V_x]  = svd(tv_x,'econ');
U_x            = U_x(:,1:rk)*S_x(1:rk,1:rk);
V_x            = V_x(:,1:rk);

tv_y           = diff_y(x1,sizeD);
tv_y           = reshape(tv_y,[h*w,s]);
[U_y,S_y,V_y]  = svd(tv_y,'econ');
U_y            = U_y(:,1:rk)*S_y(1:rk,1:rk);
V_y            = V_y(:,1:rk);

tv_z           = diff_z(x1,sizeD);
tv_z           = reshape(tv_z,[h*w,s]);
[U_z,S_z,V_z]  = svd(tv_z,'econ');
U_z            = U_z(:,1:rk)*S_z(1:rk,1:rk);
V_z            = V_z(:,1:rk);
%% admm parameter
psnrpath  = zeros(maxIter,1);
lambda    = 1;
M1        = zeros([h*w,s]); 
M2        = M1;
M3        = M2;             %\Delta_i(x1)=u*v
M4        = zeros(size(y)); % A*x-y
M5        = zeros(h*w*s,1);             % the multiplier for x-x1-e
%% main loop
for iter = 1 : maxIter
    x_pre = x;
    fprintf('\n**********************iter: %d ***********************\n', iter');     
    %% e subproblem
    e = (beta(5)*(x-x1)+M5)/(lambda+beta(5));
    %% x subproblem
    x = x(:);
    T = beta(4)*(A'*y) + (A'*M4) - M5;
    T = T(:) +beta(5)*(x1+e) ;
    x = myPCG_CS(x,A,T,beta);   
    %% U_x,U_y,U_z subproblem
    tmp_x         = reshape(diff_x(x1,sizeD),[h*w,s]);
    tmp_x         = tmp_x+M1/beta(1);
    U_x           = softthre(tmp_x*V_x, weight(1)/beta(1));
    tmp_y         = reshape(diff_y(x1,sizeD),[h*w,s]);
    tmp_y         = tmp_y+M2/beta(2);
    U_y           = softthre(tmp_y*V_y, weight(2)/beta(2));
    tmp_z         = reshape(diff_z(x1,sizeD),[h*w,s]);
    tmp_z         = tmp_z+M3/beta(3);
    U_z           = softthre(tmp_z*V_z, weight(3)/beta(3));
    %% V_x,V_y,V_z subproblem
    [u,~,v]       = svd(tmp_x'*U_x,'econ');
    V_x           = u*v';
    [u,~,v]       = svd(tmp_y'*U_y,'econ');
    V_y           = u*v';
    [u,~,v]       = svd(tmp_z'*U_z,'econ');
    V_z           = u*v';
    %% x1 subproblem
    x1   = x1(:);
    T    = beta(1)*diff_xT(U_x*V_x',sizeD)+beta(2)*diff_yT(U_y*V_y',sizeD)+beta(3)*diff_zT(U_z*V_z',sizeD);
    T    = T(:)-diff_xT(M1,sizeD)-diff_yT(M2,sizeD)-diff_zT(M3,sizeD) +beta(5)*(x-e)+M5(:);
    x1   = myPCG_TV(x1,T,beta,sizeD);
    %% - updating multipliers
    leq1 = reshape(diff_x(x1,sizeD),[h*w,s])- U_x*V_x';
    leq2 = reshape(diff_y(x1,sizeD),[h*w,s])- U_y*V_y';
    leq3 = reshape(diff_z(x1,sizeD),[h*w,s])- U_z*V_z';
    leq4 = y-A*x;
    leq5 = x - x1 -e;
    M1   = M1 + beta(1)*leq1;
    M2   = M2 + beta(2)*leq2;
    M3   = M3 + beta(3)*leq3;
	M4   = M4 + beta(4)*leq4;
	M5   = M5 + beta(5)*leq5;
    
    %- terminating the algorithm and saving some ralted information 
    relChgX = norm(x-x_pre,'fro')/max(1,norm(x_pre,'fro'));
    fprintf('stopFormulaVal: %4.4e \n', relChgX);
    fprintf('stopFormulaVal2: %4.4e \n', norm(leq4));
    if (iter> 200) ||  (relChgX < tol ) 
          disp(' !!!stopped by termination rule!!! ');  
          break;       
    end
    
    
    if isCont
        nr1 = max(abs(leq1(:)));
        nr2 = max(abs(leq2(:)));
        nr3 = max(abs(leq3(:)));
        nr4 = norm(leq4,'fro')/normD;
        nr5 = norm(leq5,'fro')/normD;
        if iter >1 && nr1 > contpar * nr1_pre
            beta(1) = min(contfactor*beta(1),1e6);
        end
        if iter>1 && nr2 > contpar * nr2_pre
            beta(2) = min(contfactor*beta(2),1e6);
        end
        if iter>1 && nr3 > contpar * nr3_pre
            beta(3) = min(contfactor*beta(3),1e6);
        end
		if iter>1 && nr4 > contpar * nr4_pre
            beta(4) = min(contfactor*beta(4),1e6);
        end
        if iter>1 && nr5 > contpar * nr5_pre
            beta(5) = min(contfactor*beta(5),1e6);
        end  
        nr1_pre = nr1; nr2_pre = nr2; nr3_pre = nr3;
        nr4_pre = nr4; nr5_pre = nr5;
    end
    fprintf('beta value is: %f, %f, %f, %f\n',beta(1),beta(2),beta(3),beta(4));
    x_ten = reshape(x,sizeD);
    L     = reshape(x1,sizeD);
    [mp,sm,er] = msqia(L/255,opts.trX/255);
    fprintf('The clean image: mp = %f, sm = %f, er =%f.\n',mp,sm,er);
    [mp,sm,er] = msqia(x_ten/255,opts.trX/255);
    fprintf('The full image: mp = %f, sm = %f, er =%f.\n',mp,sm,er);
%     figure;subplot 121;imshow(x_ten(:,:,1),[]);subplot 122;imshow(x_ten(:,:,110),[])
%     figure;
%     subplot 121;imshow(L(:,:,1),[]);title('Clean image');
%     subplot 122;imshow(x_ten(:,:,1),[]);title('Original image');
end

end