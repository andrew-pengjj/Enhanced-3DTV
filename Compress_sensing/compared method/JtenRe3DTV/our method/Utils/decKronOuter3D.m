function [Y,gY,A,B,C,out] = decKronOuter3D(X,nTerm,G,opts)
% **********************************************************************
% min || X - sum_p ( sum_{g_p} A_{g_p} \kron B_{g_p} ) \outprod c_p ||^2
% with A_{g_p} and B_{g_p} satisfying orthogonal constraint.
% 
%----------------------------------------------------------------------
%INPUTs:
% X:     input 3D tensor
% nTerm: the number of decomposed terms
% G:     pattern for each term
% opts:
% - maxIter
% - tol
%
%OUTPUTs:
% Y: output 3D tensor (Y = sum_p (sum_{g_p} A_{g_p} \kron B_{g_p} )\out c_p )
% A,B,C. A: intersity, B: pattern, C: factor for mode 3
% A:  A{ii}{g}     root
%                 /    \
%          Term1        Term2      <--- ii
%       /   /   /        \  \  \ 
%      g1  g2  g3        g1 g2 g3  <---- g
%-------------------------------------------------------------------------
% caowenf2006@163.com
% author: Yao Wang and W.F. Cao 
% xi'an jiaotong university
%*************************************************************************
%
%% process 'opts'
if ~exist('opts','var'); opts=[];end
if isfield(opts,'maxIter'); maxIter = opts.maxIter; else maxIter = 100; end
if isfield(opts,'tol'); tol = opts.tol; else tol = 1e-6; end;
nDim  = ndims(X);
if nDim ~= 3; error('the input argument should be 3-order tensor'); end
%% pad the input argument X
psizeD= size(X);
X     = preProcX(X,G); % pad each frame in X in symmetric manner
sizeD = size(X);
row   = sizeD(1); col = sizeD(2); dep = sizeD(3);

%% initilizing C and A B: A{ii}{g}
% X_dep    = Unfold(X,sizeD,3);
% [U,D,V]  = svd(X_dep,'econ');
% C        = U(:,1:nTerm);
% matX     = reshape(X,sizeD(1)*sizeD(2),sizeD(3));
% [U,D,V]  = svd(matX,'econ');
% C        = V(:,1:nTerm);
% tmpX     = reshape(U(:,1:nTerm)*D(1:nTerm,1:nTerm),sizeD(1),sizeD(2),nTerm);
% idz      = randsample(dep,nTerm);
% tmpX     = X(:,:,idz);
nG       = size(G,1); % e.g. G=[ 4 8; 1 4; 4 1; 8  8] 
for ii = 1 : nTerm
% %     tmpMat = tmpX(:,:,ii);
    tmpMat = rand(row,col); %randomized initialization
    for g = 1 :1: nG
        szgB     = G(g,:);
        szgA     = [row, col]./szgB; 
        tmpFlat  = kronUnFold(tmpMat,szgA,szgB);
        [U,D,V]  = svd(tmpFlat,'econ');
        ng       = ceil(size(U,2)*1);
        U        = U(:,1:ng); V = V(:,1:ng); D = D(1:ng,1:ng);
        A{ii}{g} = reshape(U*D,[szgA,ng]);
        B{ii}{g} = reshape(V,  [szgB,ng]);
    end
end
    
%% BCD method is used to slove this optimizaiton problem
cost_hist = zeros(maxIter,1);
tic;
for  iter = 1 : maxIter
    
    fprintf('--->iter: %d  \n',iter);
    
    %- undating C
    AB   = termSum(A,B,row,col);    % AB: (row*col) x nTerm
    matAB= reshape(AB,row*col,nTerm);
    matX = reshape(X,row*col,dep);
    C    = ( pinv(matAB)*matX )';
%     [C,~]= qr(C,0);
    C    = C*diag( 1./sqrt(sum(C.^2,1)) );
    
    %- updating A B
%     Y = allKronDecAppr(A,B,C,sizeD);
%     for ii = 1 : nTerm
%         for g = 1 :1 : nG
%             gAppr1   = gKronDecAppr(A{ii}{g},B{ii}{g},C(:,ii));
%             resX     = X - Y + gAppr1;
%             szgB     = G(g,:);  % the size of B in g-th group
%             szgA     = [row, col]./szgB; % the size of A in g-th group
%             [A{ii}{g},B{ii}{g}] = gKronDec(resX,C(:,ii),szgA,szgB); 
% %             size(A{ii}{g})
% %             size(B{ii}{g})
%             gAppr    = gKronDecAppr(A{ii}{g},B{ii}{g},C(:,ii));
%             Y        = Y - gAppr1 + gAppr;
%         end
%     end
    resX = X-allKronDecAppr(A,B,C,sizeD);
    for ii = 1 : nTerm
        for g = 1 :1: nG
             gAppr_g1 = gKronDecAppr(A{ii}{g},B{ii}{g},C(:,ii));
             resX     = resX + gAppr_g1;
             szgB     = G(g,:);
             szgA     = [row,col]./szgB;
             [A{ii}{g},B{ii}{g}] = gKronDec(resX,C(:,ii),szgA,szgB);
             gAppr_g2 = gKronDecAppr(A{ii}{g},B{ii}{g},C(:,ii));
             resX     = resX - gAppr_g2;
        end
    end
      
    
    %- stopping rule
    cost_hist(iter) = norm(resX(:),'fro');
    fprintf('cost_hist: %4.4e  \n',cost_hist(iter));
    if iter > 5  && all( abs( cost_hist(iter-3:iter) - cost_hist(iter-4:iter-1) ) / max(1,max(cost_hist(iter-4:iter-1))) < tol )
        disp('stopped by the termination rule.');
        break;
    end

end

%% remove the padded elements
gY       = cell(2,nG);
for g = 1 : nG
    for ii = 1 : nTerm
        gY{ii,g} = gKronDecAppr(A{ii}{g},B{ii}{g},C(:,ii));
        gY{ii,g} = postProcX(gY{ii,g},psizeD);
    end
end

% Y = gY{1};
% for g = 2 : nG
%     Y = gY{g} + Y;
% end
Y = postProcX(allKronDecAppr(A,B,C,sizeD),psizeD);
out.iter = iter;
out.time = toc;
out.cost_hist = cost_hist(1:iter);
return;

end

function  AB   = termSum(A,B,row,col)    
% AB: (row*col) x nTerm
% sum_p (A_p \kron B_p)
% A{ii}{g} B{ii}{g}
% 
nTerm1 = numel(A);
nTerm2 = numel(B);
if ~isequal(nTerm1,nTerm2)
    error('the number of terms in A must be equal to that in B!');
end
nG    = numel(A{1});
AB    = zeros(row,col,nTerm1);
for ii = 1 : nTerm1
    for g = 1 : nG
        AB(:,:,ii) = AB(:,:,ii) + tenKron(A{ii}{g},B{ii}{g});
    end
end

end

function C = tenKron(tenA,tenB)
nfrm1 = size(tenA,3);
nfrm2 = size(tenB,3);
C     = zeros(size(tenA,1)*size(tenB,1), size(tenA,2)*size(tenB,2));
if ~isequal(nfrm1,nfrm2)
    error('the size of tenA must be equal to that of tenB!');
end
for k = 1 : nfrm1
    C = C + kron(tenA(:,:,k),tenB(:,:,k));
end

end
function gAppr = gKronDecAppr(tenA,tenB,vecC)
%
% ( sum_g(A_p \kron B_p) ) \outprod (vecC)
% 
sizeD1 = size(tenA);
sizeD2 = size(tenB);
dep   = numel(vecC);
tmp   = zeros(sizeD1(1)*sizeD2(1),sizeD1(2)*sizeD2(2));
tmp   = tmp + tenKron(tenA,tenB);
gAppr = reshape(tmp(:)*vecC', sizeD1(1)*sizeD2(1),sizeD1(2)*sizeD2(2),dep);
end

function  Y = allKronDecAppr(A,B,C,sizeD)
%
% Y = sum_term (sum_g gAppr)
% 
nTerm = size(C,2);
Y     = zeros(sizeD);
nG    = numel(A{1});
for ii = 1 : nTerm
    for g = 1 : nG
        Y = Y + gKronDecAppr(A{ii}{g},B{ii}{g},C(:,ii));
    end
end

end

function  [A_ig,B_ig] = gKronDec(resX,vecC,szgA,szgB)
% || resX - (\sum_p A_p \kron B_p) \outprod vecC ||^2 
% s.t. A_P is orthgonal, B_p is orthgonal
%
[row,col,dep] = size(resX);
tmp           = reshape(resX,row*col,dep);
tmp_mat       = reshape(tmp*vecC,row,col); % note ||vecC||_2 =1
kron_mat      = kronUnFold(tmp_mat,szgA,szgB);
[U,D,V]       = svd(kron_mat,'econ');
ng            = ceil(size(U,2)*1); %%%????
U             = U(:,1:ng); V = V(:,1:ng); D = D(1:ng,1:ng);
A_ig          = reshape(U*D,[szgA,ng]);
B_ig          = reshape(V, [szgB,ng]);
end

function Y = preProcX(X,G)
% G = [size of row,size of col ]  for each group in B
% 
%
[row,col,nfrm] = size(X);
rSz = G(:,1); cSz = G(:,2);
if numel(rSz) ==1
    rCM = rSz;
end
if numel(cSz) == 1
    cCM = cSz;
end
for i = 2 : numel(rSz);
    if i == 2
        rCM  = lcm(rSz(i),rSz(i-1));
        cCM  = lcm(cSz(i),cSz(i-1));
    else
        rCM = lcm(rCM,rSz(i));
        cCM = lcm(cCM,rSz(i));
    end
end

rInc    = mod(row,rCM);
if rInc ~= 0
    row_new = row + (rCM-rInc);
else
    rInc    = 0;
    row_new = row;
end

cInc    = mod(col,cCM);
if cInc ~= 0
    col_new = col + (cCM-cInc);
else
    cInc    = 0;
    col_new = col;
end

Y       = zeros(row_new,col_new,nfrm);
if cInc ~=0 && rInc ~=0
    for k = 1 : nfrm
        Y(:,:,k) = padarray(X(:,:,k),[rCM-rInc,cCM-cInc],'replicate','post');
    end
elseif cInc~=0 && rInc ==0
    for k = 1 : nfrm
        Y(:,:,k) = padarray(X(:,:,k),[0,cCM-cInc],'replicate','post');
    end
elseif cInc ==0 && rInc~=0
    for k = 1 : nfrm
        Y(:,:,k) = padarray(X(:,:,k),[rCM-rInc,0],'replicate','post');
    end
else
    Y = X;
end;

end
    
function Y     = postProcX(Y,psizeD)

row = psizeD(1);
col = psizeD(2);
Y   = Y(1:row,1:col,:);
end
  
