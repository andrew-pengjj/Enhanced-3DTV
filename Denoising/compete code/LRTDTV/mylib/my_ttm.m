function X = my_ttm(X,V,dims,sizeOld,sizeNew,Ndim)
% tensor times matrix,V is cell
for n = dims(1:end) 
    order = [n,1:n-1,n+1:Ndim];
    temp  = reshape(permute(X,order), sizeOld(n), []);
    temp  = V{n}*temp;
    sizeOld(n) = sizeNew(n);    
    X = ipermute(reshape(temp, [sizeOld(n),sizeOld(1:n-1),sizeOld(n+1:Ndim)]),order);
end