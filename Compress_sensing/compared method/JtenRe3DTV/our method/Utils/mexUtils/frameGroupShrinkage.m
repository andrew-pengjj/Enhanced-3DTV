function  [z] = frameGroupShrinkage(X,G,lambda,beta)


[~, n ]= size(X);
z = cell(1,n);
for k = 1: n %% please speed up!!!
    
    tempZ = diag(sparse(X(:,k)))*G + lambda{k}/beta;
    z{k} = shrinkageMex(tempZ, [1,beta]);
    
end

end