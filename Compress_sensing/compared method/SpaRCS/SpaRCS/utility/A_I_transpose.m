function y = A_I_transpose(A,x,I);
% assume consistency all around
% this function takes as input a function handle A and a vector x
% and computes y = A(:,I)'*x
z = A(x);
y = z(I);