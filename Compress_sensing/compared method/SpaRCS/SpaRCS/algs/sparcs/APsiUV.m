function b = APsiUV(alpha,A,PsiU, PsiV)

z = zeros(size(PsiU, 1), size(PsiV, 1));
%for kk=1:length(alpha)
%    z= z+alpha(kk)*PsiU(:, kk)*PsiV(:, kk)';
%end
z = PsiU*diag(alpha)*PsiV';

b = A(z(:));