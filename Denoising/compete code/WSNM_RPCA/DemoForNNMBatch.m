%%仿真数值实验，将WNNM的RPCA与过去的RPCA对比
%%
clc
clear
addpath exact_alm_rpca inexact_alm_rpca
m=400;

%num=5;
Accur_NNM=0;rank_NNM=0;E_NNM_0=0;

pr = [0.01:0.01:0.5];
ps = [0.01:0.01:0.5];
error_mat = zeros(length(pr), length(ps));
rank_mat  = zeros(length(pr), length(ps));
sparse_num = zeros(length(pr), length(ps));
for i=20:length(pr)
    for j=30:length(ps)
        r=round(pr(i)*m);
        EL0=round(m*m*ps(j));
        
        %construct low rank matrix for recovering
        U=normrnd(0,1,m,r);V=normrnd(0,1,m,r);
        A0=U*V';
        E=zeros(m,m);
        Ind = randperm(m*m);
        % 在低秩矩阵上加入稀疏噪声
        E(Ind(1:EL0))=1000*rand(1,EL0)-500;
        D=A0+E;
        
        %perform RPCA
        fprintf('----pr=%3f--------ps=%3f------\n', pr(i), ps(j));
        [A_NNM E_NNM ]=exact_alm_rpca(D,1 / sqrt(m));
        Accur_NNM=(sum(sum((A_NNM-A0).^2))).^0.5/(sum(sum(A0.^2))).^0.5;
        rank_NNM=rank(A_NNM);
        E_NNM_0=length(find(abs(E_NNM)>0));
        disp(['NNM精度:  ' num2str(Accur_NNM ) ',  秩：'  num2str(rank_NNM)  ',  |E|_0;'  num2str(E_NNM_0)]);
        disp('--------------------------------------');
        
        error_mat(i,j) = Accur_NNM;
        rank_mat(i,j) = rank_NNM;
        sparse_num(i,j) = E_NNM_0;
    end
end

save('D:\SR\WSNM-new-exp\WSNM_RPCA\batch_plot\NNM_error_matrix.mat','error_mat', 'rank_mat', 'sparse_num');
% disp(['inexactWNNM精度:   ' num2str(Accur_iWNNM ) '.     秩：'  num2str(rank_iWNNM)  '.    |E|_0;'  num2str(E_iWNNM_0)  '.    时间：' num2str(time_iWNNM) '.' ])
% disp(['iNNM精度:          ' num2str(Accur_iNNM ) '.     秩：'  num2str(rank_iNNM)  '.     |E|_0;'  num2str(E_iNNM_0)  '.   时间：' num2str(time_iNNM) '.' ])
% disp(['秩'  num2str(r) '|E|_0' num2str(EL0)])
% disp(['inexactNNM精度:    ' num2str(Accur_iNNM ) '.     秩：'  num2str(rank_iNNM)  '.     时间：' num2str(time_iNNM) '.' ])
% %
