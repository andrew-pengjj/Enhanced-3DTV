function X = Atnoiselet(b, indx, indx2, siz)

y = zeros(prod(siz), 1);
y(indx) = b;

X1 = realnoiselet(y)/sqrt(length(y));
X = 0*y;
X(indx2) = X1;
X = reshape(X, siz);