function Clean_HSI=RPCA_WSNM_ite(noised_HSI,nSig,p,C,lamb,delta)
N_Img = noised_HSI;
E_Img = N_Img;
for iter = 1:3
    ti=clock;
    par = ParWsnmSet(nSig);
    E_Img = E_Img + delta*(N_Img - E_Img);
    E_Img = RPCA_WSNM(E_Img, N_Img, nSig, p, C, lamb, iter, par);
    Clean_HSI = E_Img;
    tf=clock;
    fprintf('run time1 = %g secs \n',etime(tf,ti));
end
