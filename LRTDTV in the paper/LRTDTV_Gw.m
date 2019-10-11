%% =========================== Frist part notes ===========================
% ADMM algorithm: tensor denosing
%    solve the following problem 
%           (SSTV regularized low rank tucker decomposition problem)
%    
%      The approximate model:  
%                         min tau*||F||_1 + lambda*||E||_1 + beta*||N||_1
%                             s.t. D = X+E, X = Z,DwZ = F
%                                  X = Core*U1*U2*U3,rank(Core_i)=ri 
%      where,D is weight SSTV difference operator in the above model
%
%      The lagrange function is :
%       tau*||F||_1 + lambda*||E||_1 +beta*||N||_1+ <M1, D-X-E-N> + <M2,X-Z> 
%   +<Gamma,DZ - F> + beta/2*( ||D-X-E-N||_F^2 + ||X-Z||_F^2 + ||DZ - F||_F^2 )        
%
%% =========================== Second part notes===========================

% Reference paper: Hyperspectral Image Restoration via Total Variation 
%                      Regularized Low-rank Tensor Decomposition
% Author: Yao Wang, Jiangjun Peng
% E-mail addresses: andrew.pengjj@gmail.com
% -------------------------------------------------------------------------

%% =========================== Thrid part notes =========================== 
% INPUT:
%   Noi:     noisy 3-D image of size M*N*p normalized to [0,1] band by band
%   tau:     the trade-off parameter (recommended value 1)              
%   lambda:  sparse noise coefficient
%   beta :   The Gaussian noise cofficient
%   rank:    rank constraint,[0.8*M,0.8*N,r]
% OUTPUT:
%  clean_iamge:   3-D denoised image
%  S:             The noise term
%  out_value:     MPSNR and MSSIM and ERGAS valuses of each iteration 
%  ========================================================================

function [clean_image,S,G, errList,time] = LRTDTV_Gw(Noi, tau,lambda,beta,rank,weight)
tic
sizeD           = size(Noi);
normD           = norm(Noi(:)); 
n               = prod(sizeD);
maxIter         = 40;
epsilon         = 1e-6;  
mu              = 0.01;             % The ascending multiplier value

out_value       = [];
out_value.SSIM  = [];
out_value.PSNR  = [];
out_value.ERGAS = [];
%%  Initialization 
X               = rand(sizeD);      % X : The clean image
Z               = X;                % Z : auxiliary variable for X
G               = zeros(sizeD);     % G : Gaussian noise
S               = zeros(sizeD);     % S : sparse noise 
F               = zeros(3*n,1);     % F : auxiliary variable for tv
Gamma           = F;                % The multiplier for DZ-F
M1              = zeros(size(Noi)); % The multiplier for          
M2              = M1;
%% main loop

for iter = 1: maxIter
    preX       = X;
    %% - update X
    temp       =  1/2*(Noi-S-G+Z+(M1-M2)/mu);
    X          = tucker_hooi(temp,rank);
    %% - update Z
    z          = Z(:);
    z          = myPCG1(z,X,M2,F,Gamma,mu,sizeD);
    Z          = reshape(z,sizeD);
    %% - update F
    diff_Z     = diff3_weight(Z(:), sizeD,weight); 
    F          = softthre( diff_Z+ Gamma/mu, tau/mu );   
    %% - update S
    S          = softthre(Noi-X-G+M1/mu,lambda/mu);
    %% - update G
    G          = (mu*(Noi -X -S) + M1)/(mu+2*beta);    
    %% - update M
    M1         = M1+ mu*(Noi-X-S-G);
    M2         = M2 + mu*(X-Z); 
    %% update Gamma
    Gamma      =  Gamma+mu*(diff_Z-F);            
    mu         =  mu * 1.5; 
    %% compute the error
    errList    = norm(X(:)-preX(:)) / normD;
    fprintf('LRTDTV: iterations = %d   difference=%f\n', iter, errList);
    if errList < epsilon
        break;  
    end 
    %% output SSIM and PSNR values of each step
%     load simu_indian
%     OriData3 = simu_indian;
%     [out_value.PSNR(iter),out_value.SSIM(iter),out_value.ERGAS(iter)]=msqia(OriData3,X);
end
clean_image = X;
fprintf('HALRTC ends: total iterations = %d   difference=%f\n\n', iter, errList);
toc
time=toc;
end