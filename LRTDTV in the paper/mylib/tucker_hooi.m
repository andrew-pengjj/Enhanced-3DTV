function X =tucker_hooi(X_all,rank)
[Core, U]  =  tensorSVD(X_all,rank);% initial with HOSVD
ilastX     =  X_all;
ndim       =  length(size(Core));
sizeD      =  size(X_all);
normD      =  norm(X_all(:));
for j=1:ndim
    Ut{j}  =  U{j}';
end
for Initer =1:6
    for j=1:ndim
        RRank =rank;
        tempC    = my_ttm(X_all,Ut,[1:j-1,j+1:ndim],sizeD,rank,ndim);
        RRank(j) =sizeD(j);
        UnfoldC  = Unfold(tempC,RRank,j);
        [V1,~,~]= svd(UnfoldC,'econ');
        U{j}    = V1(:,1:rank(j));
        Ut{j}       = U{j}';
    end
    Core = my_ttm(X_all,Ut,1:ndim,sizeD,rank,ndim);
    X    = my_ttm(Core,U,1:ndim,rank,sizeD,ndim);
    % stop criterion
    errT  = ilastX-X;
    err   = norm(errT(:))/normD;
    if err<=1e-3
%         fprintf('  ++++ Innerloop iterations = %d+++  \n',Initer);
        break; 
    else
        ilastX=X;
    end
end