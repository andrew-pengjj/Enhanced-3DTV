function [rX,rgX,alpha,Dict,out] = decKron3D(X,G,opts)
%*************************************************************************
% min 0.5*||X - \sum_g \sum_{p \in I_g} alpha_p \kron Dict_p  ||^2
% with alpha_p and Dict_p satisfying some constraint
% such as. nonnegative, orthogonal, alpha_p is low rank in the third mode.
%--------------------------------------------------------------------------
% INPUTs:
% X - input tensor
% G - size for pattern Dict
% opts:
%  - maxIter
%  - tol
%  - display
% OUTPUTs:
% rX  - rec. X
% rgX - rec. X in each group
% alpha, Dict -                   root
%                      /     /      |    \     \
%                     g1    g2      g3    g4    g5   <-----alpha{g}
%
%Note:  alpha{g} is one tensor of 4-D order
%-------------------------------------------------------------------------
% caowenf2006@163.com
%*************************************************************************
%
%
%% check and initilization
if ~exist('opts','var'); opts=[];end
if isfield(opts,'maxIter'); maxIter = opts.maxIter; else maxIter = 100; end
if isfield(opts,'tol');   tol = opts.tol; else tol = 1e-6; end;
if isfield(opts,'init_method'); init_method = opts.init_method; else init_method ='adpative';end
nDim  = ndims(X);
if nDim ~= 3; error('the input argument should be 3-order tensor'); end
con_type = 'orth';
%% pad the input argument X
psizeD = size(X);
X      = preProcX(X,G); % pad each frame in X in symmetric manner
sizeD  = size(X);
row    = sizeD(1); col = sizeD(2); dep = sizeD(3);
nG     = size(G,1); % e.g. G=[ 3,4,1; 8 8,1; 4 8 1; 8 4 1]

%% initilizating alpha and Dict
alpha = cell(nG,1);
Dict  = cell(nG,1);
for g = 1 : nG
    
    dict_szg  = G(g,:);
    alpha_szg = sizeD./dict_szg;
    matX      = mexkronUnFold3D(X,alpha_szg,dict_szg);
    rk        = rank(matX);
    if strcmp('random',init_method)
        alpha{g}   = rand([alpha_szg,rk]);
        Dict{g}    = rand([dict_szg,rk]);
        
    else %adaptive 
        [U1,S1,V1] = svd(matX,'econ');
%         dgD     = max(diag(S1));
%         rk      = numel(find(diag(S1)>dgD*0.005));
        
        U1 = U1(:,1:rk); S1 = S1(1:rk,1:rk); V1 = V1(:,1:rk);
        alpha{g}   = reshape(U1*S1,[alpha_szg,rk]);
        Dict{g}    = reshape(V1,[dict_szg,rk]);
    end
    
end      
%% main loop
tic;
cost_hist = zeros(maxIter,1);
for iter = 1 : maxIter
    
    fprintf('--->iter: %d \n', iter);
    
    % updatin alpha_p and Dict_p

    resX = X - allDecKronAppr(alpha,Dict,sizeD);
    for g = 1 : nG
        gAppr1      = gDecKronAppr(alpha{g},Dict{g});
        resX        = resX + gAppr1;
        
        dict_szg    = G(g,:);
        alpha_szg   = sizeD./dict_szg;
        matResX     = mexkronUnFold3D(resX,alpha_szg,dict_szg);
        [alpha{g},Dict{g}] = gDecKronSolver(matResX,alpha_szg,dict_szg,con_type);
        
        gAppr2      = gDecKronAppr(alpha{g},Dict{g});
        resX        = resX - gAppr2;
    end
    
    % stopping rule
    cost_hist(iter) = norm(resX(:),'fro');
    fprintf('cost_hist: %4.4e  \n',cost_hist(iter));
    if iter > 2  && all( abs( cost_hist(iter-1:iter) - cost_hist(iter-2:iter-1) ) / max(1,max(cost_hist(iter-2:iter-1))) < tol )
        disp('stopped by the termination rule.');
        break;
    end
    
    
end

rgX           = cell(nG,1);
rX            = zeros(psizeD);
for g = 1 : nG
    rgX{g}    = postProcX(gDecKronAppr(alpha{g},Dict{g}), psizeD);
    rX        = rX + rgX{g};
end
out.time     = toc;
out.cost_his = cost_hist(1:iter);
out.iter     = iter;

return;
end
        

function Y = preProcX(X,G)
% G = [size of row,size of col, size of dep ]  for each group in Dict
% 
%
[row,col,dep] = size(X);
rSz = G(:,1); cSz = G(:,2); dSz = G(:,3);
if numel(rSz) ==1
    rCM = rSz;
end
if numel(cSz) == 1
    cCM = cSz;
end
if numel(dSz) == 1
    dCM = dSz;
end
for i = 2 : numel(rSz);
    if i == 2
        rCM  = lcm(rSz(i),rSz(i-1));
        cCM  = lcm(cSz(i),cSz(i-1));
        dCM  = lcm(dSz(i),dSz(i-1));
    else
        rCM = lcm(rCM,rSz(i));
        cCM = lcm(cCM,rSz(i));
        dCM = lcm(dCM,dSz(i));
    end
end

rInc    = mod(row,rCM);
cInc    = mod(col,cCM);
dInc    = mod(dep,dCM);

rAdd    = (rInc~=0)*(rCM-rInc);
cAdd    = (cInc~=0)*(cCM-cInc);
dAdd    = (dInc~=0)*(dCM-dInc);
if rAdd==0 && cAdd==0 && dAdd ==0
    Y   = X;
else
    Y   = padarray(X,[rAdd,cAdd,dAdd],'replicate','post');
end

end
    
function Y     = postProcX(X,psizeD)

row  = psizeD(1);
col  = psizeD(2);
nfrm = psizeD(3);
Y   =  X(1:row,1:col,1:nfrm);
end
  

function  rgX   = gDecKronAppr(Ag,Bg)

nDim1   = ndims(Ag);
nDim2   = ndims(Bg);
if ~isequal(nDim1,nDim2)
    error('the dim. is not equal.');
end
sizeD1 = size(Ag);
sizeD2 = size(Bg);

rgX   = zeros(sizeD1(1)*sizeD2(1),sizeD1(2)*sizeD2(2),sizeD1(3)*sizeD2(3));
for k = 1 : sizeD1(4)
    rgX = rgX + mexKron3D(Ag(:,:,:,k),Bg(:,:,:,k));
end

end

function C = allDecKronAppr(A,B,sizeD)

nG = numel(A);
C  = zeros(sizeD);
for g = 1 : nG
    C = C + gDecKronAppr(A{g},B{g});
end

end

function [Ag,Bg] = gDecKronSolver(matResX,a_szg,d_szg,con_type)
%
%  0.5*|| X - \sum_p vec(A_p) * vec(B_p)' ||^2, 
%   s.t. (1) orth
%        (2) nonnegative
%        (3) Laplacian smoothness
%        (4) low rank
%        ????
%

if strcmp('orth',con_type)
    rk      = rank(matResX);
    [U,D,V] = svd(matResX,'econ');
%     dgD     = max(diag(D));
%     rk      = numel(find(diag(D)>dgD*0.005))
%     rk = 15;
    U       = U(:,1:rk); D = D(1:rk,1:rk); V = V(:,1:rk);
    Ag      = reshape(U*D,[a_szg,rk]);
    Bg      = reshape(V,[d_szg,rk]);
elseif strcmp('nonneg',con_type)
    error(' to be updated later...');
elseif strcmp('low_rank',con_type)
    error(' to be updated later....');
end

end
            
        

