function [Z]=Ten_ShrinkNuclear_PG(X,lambda)
% minimize_Z 0.5||Z-X||_{F}^{2} + \sum_i lambda_i ||Z_(i)||_{*}
%
% % tic;
% % sizeD=size(X);
% % nDim =ndims(X);
% % Zm=cell(nDim,1);
% % for mode = 1:nDim
% %     tempX=Unfold(X,sizeD,mode);
% %     [U,D,V]=MySVD(tempX);
% %     d_shrink=max(diag(D)-lambda(mode),0);
% %     idx=find(d_shrink>0);
% %     tempZ=U(:,idx)*diag(d_shrink(idx))*V(:,idx)';
% %     Zm{mode}=Fold(tempZ,sizeD,mode);
% % end
% % 
% % Z=zeros(sizeD);
% % for mode=1:nDim
% %     Z=Z+Zm{mode};
% % end
% % Z=Z/nDim;
% % toc
% % end
%
% tic;
mu=1e10;
sizeD=size(X);
nDim =ndims(X);
Z=X;%zeros(sizeD);
M=cell(nDim,1);
tol=1e-6;
maxiter=20;
for iter =1:maxiter
    
    Zold=Z;
    
    for mode = 1:nDim
        tempX=Unfold(X,sizeD,mode);
        tempZ=Unfold(Z,sizeD,mode);
        XZ =((1/nDim)*tempX + mu*tempZ )/(1/nDim+mu);

        [U,D,V]=MySVD(XZ);
        d_shrink=max(diag(D)-lambda(mode)/(1/nDim+mu),0);
%         idx=find(d_shrink>0);
%         tempM=U(:,idx)*diag(d_shrink(idx))*V(:,idx)';
        tempM=U*diag(d_shrink)*V';
        rank(tempM)
        M{mode}=Fold(tempM,sizeD,mode);
    end

    Z=M{1};
    for mode=2:nDim
        Z=Z+M{mode};
    end
    Z=Z/nDim;
    
    relchg=norm(Z(:)-Zold(:))/norm(Zold(:));
    if  relchg<tol; break; end;
    
end
% iter
% toc
end
