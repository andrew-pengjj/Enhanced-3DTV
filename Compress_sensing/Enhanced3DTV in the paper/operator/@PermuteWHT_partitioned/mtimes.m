function res=mtimes(A,x)
if ~isa(A,'PermuteWHT_partitioned')
    error('In A*x, A should be a object of name PermuteWHT_partitioned');
end

if ~A.adjoint
    
    x      = reshape(x, A.m, A.L);
    X      = [x; zeros(A.M - A.m, A.L)];
    FX     = fdWHtrans(X(A.perm,:));
%     FX=fwht(X(A.perm,:));
    res     = FX(A.picks);
else
    
    FY = zeros(A.M, A.L);
    FY(A.picks) = x;
    res  = zeros(A.M, A.L);
    res(A.perm,:) = fdWHtrans(FY);
%     res(A.perm,:)=fwht(FY);
    res = res(1:A.m,:);
    res = res(:);
end
end