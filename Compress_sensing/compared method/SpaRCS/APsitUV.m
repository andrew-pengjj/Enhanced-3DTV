function alpha = APsitUV(b,At,PsiU, PsiV)

t = At(b);

alpha = PsiU'*t*PsiV;
alpha = diag(alpha);