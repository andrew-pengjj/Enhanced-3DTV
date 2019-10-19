function [Z,Znorm,Zrank,out]=Tensor_IHT(X,Omega,nRank,opts)
% DESCRIPTION:
% Problem: 0.5*||(Z-X).*Omega||_{F}^{2}  s.t. rank(Z_i}<=r_i
%
% Algorithm: (IHT)
%            Y_k<---Z_k- tao*Omega'*(Omega(Z_k)-Omega(X));
%            for mode =1: nDims
%                M_mode=Shrink(Unfold(Y_k,mode),r_mode);
%            end
%            Z_{k+1}<----\sum_mode w_mode Fold(M_mode,mode);
% Ref.: Min Zhang, Lei Yang and Zheng Hai Huang. Mimimum n-Rank
% Approximation via Iterative Hard Thresholding, manuscirpt 2014.
%
% INPUTS:
%        X ----------- Observed tensor
%        Omega---------Observatin mask
%        nRank -------- Tucker rank (Note it is a vector)
%
%        opts.
%        Z_init: the initilization of Z
%        tao : iteration step
%        epsilon: error bound
%        Zstar:  the true Z  
%        maxiter: the number of the maximum iteration
%        kesi: the parameter for predicting n-rank
%          
% OUTPUTS:
%       Z------------- The iterative result
%       Znorm---------- The nuclear norm for each mode 
%       Zrank---------- The Tucker rank of the recovered Z
%
%       out:
%       iter 
%       time
%       relErr
%       all_relchg
%
%  E-mail: caowenf2006@163.com  @ Xi'an Jiaotong University
%
%% Check the inputs and Initilization
if nargin<3
    error('At least three essential input variables');
end
sizeD=size(X);
nDim =ndims(X);

if isempty(nRank)
    Rank_Predict=true;
else
    Rank_Predict=false;
    if length(nRank)~=nDim
        error('nRank: note nRank is one vector');
    end
end

if ~exist('opts','var'); opts=[]; end;
if isfield(opts,'Z_init'); Z_init=opts.Z_init;else Z_init=X;end;
if isfield(opts,'tao'); tao=opts.tao; else tao=1; end;
if isfield(opts,'epsilon');epsilon=opts.epsilon; else epsilon=1e-5;end;
if isfield(opts,'Zstar'); Zstar=opts.Zstar; isTrueZ=1;else isTrueZ=0; end;
if isfield(opts,'maxiter'); maxiter=opts.maxiter;else maxiter=300; end;
if isfield(opts,'kesi'); kesi=opts.kesi; else kesi=1e-2;end
if isfield(opts,'w'); w=opts.w; else w=ones(nDim,1)*(1/nDim);end

if Rank_Predict
    nRank_init=floor(sqrt(sizeD));%floor(sizeD/2);
    nRank =nRank_init;
end


%% Main loop
tic;
Znorm=zeros(nDim,1);
Zrank=zeros(nDim,1);
Z = Z_init;
all_relchg=zeros(maxiter,1);
ZXerr =zeros(maxiter,1);
sigvec =cell(nDim,1);
if isTrueZ
    relErr=zeros(maxiter,1);
end
for iter=1:maxiter
    
    Z_old=Z;
    
    % Gradient descent
    Y=Z.*(1-tao*Omega)+tao*Omega.*X;
    
    % Hard Thresholding 
    Z=zeros(sizeD);
    for mode =nDim
        tempY=Unfold(Y,sizeD,mode);
        [ThresY,sigvec{mode}]=Hard_Thres(tempY,nRank(mode));
        Znorm(mode)=sum(sigvec{mode}(1:nRank(mode)));
        Zrank(mode)=nRank(mode);
        Z=Z+w(mode)*Fold(ThresY,sizeD,mode);
    end
    
    % Rank Prediction
    TempZX=(Z-X).*Omega;
    ZXerr(iter)=norm(TempZX(:),'fro');
    if Rank_Predict
        for mode=1:nDim
           nRank(mode)= length(find(sigvec{mode}>kesi*sigvec{mode}(1)));
        end
        
        if iter>2 && ZXerr(iter)>ZXerr(iter-1) && ~any(nRank+1>sizeD)%% Rank shouldn't larger than the sizeD
           nRank = nRank +1;
        end
    end
    
    % Check the termination condition
    relchg=norm(Z(:)-Z_old(:),'fro')/max(1,norm(Z_old(:),'fro'));
    all_relchg(iter)=relchg;
    if iter>4 && all( all_relchg(iter-3 : iter)<epsilon )
        break;
    end
    
    % Compute the recovery bound if the true Z0 exists
    if isTrueZ
        relErr(iter)=norm(Z(:)-Zstar(:),'fro');
    end
    
end%%END FOR

%% Output 
if isTrueZ
    out.relErr=relErr(1:iter);
else
    out.relErr=[];
end
out.time=toc;
out.all_relchg=all_relchg(1:iter);
out.iter=iter;

return;
end


function [ThresY,sigvec_mode]=Hard_Thres(Y,Rank_mode)%
%
% min_X 0.5*||X-A||^2 s.t. rank(X)<= r;
%
%

[U,D,V]=MySVD(Y);
sigvec_mode=diag(D);
Vt=V';
ThresY =U(:,1:Rank_mode)*diag(sigvec_mode(1:Rank_mode))*Vt(1:Rank_mode,:);

end

    
    
    
    






