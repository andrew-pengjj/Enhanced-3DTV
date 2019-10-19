% clear;
% clc;
addpath(genpath(pwd));

save_path = 'HSI_results\';

% load('real_data\simulated_australia1.mat');
% 
% noised_HSI = hsi_normalized;
noised_HSI   = D;
nSig = 20/255;
p = 0.7;
C = 0.007;
lamb = 1.2;
N_Img = noised_HSI;
E_Img = N_Img;
delta = 0.1;
%% Iterative Restoration Method for HSI
for iter = 1:3 
    ti=clock;
    par = ParSet(nSig);
    E_Img = E_Img + delta*(N_Img - E_Img);
    E_Img = RPCA_WSNM(E_Img, N_Img, nSig, p, C, lamb, iter, par);
    Clean_HSI = E_Img;
    save([save_path 'Clean_HSI_iter=' num2str(iter) '.mat'],'Clean_HSI');
    tf=clock;
    fprintf('run time1 = %g secs \n',etime(tf,ti));
end
