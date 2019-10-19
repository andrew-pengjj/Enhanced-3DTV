function y = A_I(A,x,I,N);
% assume consistency all around
% this function takes as input 
% a) function handle A \in R^{MxN}
% b) vector x \in R^{length(I)x1}
% and computes y = A(:,I)*x
z = zeros(N,1); z(I) = x;
y = A(z);