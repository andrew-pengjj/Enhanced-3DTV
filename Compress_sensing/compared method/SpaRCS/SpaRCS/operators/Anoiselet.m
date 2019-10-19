function b = Anoiselet(X, indx, indx2)
X = X(indx2);
y = realnoiselet(X(:))/sqrt(length(X(:)));
b = y(indx);


%%%%%%%
%indx = randperm(prod(size(X)); indx = indx(1:M);
