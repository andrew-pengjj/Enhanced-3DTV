function res=ctranspose(A)
res=A;
res.adjoint=xor(res.adjoint,1);
end