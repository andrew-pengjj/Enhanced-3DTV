function [mpsnr,mssim,ergas,Time]=RunAllMethod(options,Omsi,Nmsi)
% 1.1 ========== METHOD SET==========
METHODSET = {'Nosiy','NNM','WNNM','LRMR','BM4D','TDL','WSNM','LRTV','LRTDTV','Enhanced3DTV'};
% 1.2 ========== METHOD SWITCH==========
STATUS_NNM   = options.STATUS_NNM;  % set to 0 for turning off;
STATUS_WNNM  = options.STATUS_WNNM;
STATUS_LRMR  = options.STATUS_LRMR;%adjusted
STATUS_BM4D  = options.STATUS_BM4D;
STATUS_TDL   = options.STATUS_TDL;
STATUS_WSNM  = options.STATUS_WSNM;
STATUS_LRTV  = options.STATUS_LRTV;
STATUS_LRTDTV = options.STATUS_LRTDTV;
STATUS_Enhanced3DTV = options.STATUS_Enhanced3DTV;
gausssigma   = options.gausssigma;
% spsigma      = options.spsigma;
msi_sz    =    size(Omsi);
[M,N,p]   =    size(Omsi);

index  = 1;
[mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi, Nmsi);
Time(index) = 0;

index = index+1;
%% Use NNM
if STATUS_NNM
    lambda = 1/sqrt(msi_sz(1)*msi_sz(2));
    tol=1e-6;
    maxIter=100;
    fprintf('\n');
    disp(['performing ',METHODSET{index}, ' ... ']);
    D = zeros(M*N,p);
    for i=1:p
        bandp  = Nmsi(:,:,i);
        D(:,i) = bandp(:); 
    end
    tic;
    [A,~,~] = inexact_alm_rpca(D, lambda, tol, maxIter);
    Img = reshape(A,msi_sz);
    Time(index) = toc;
    %   filename = strcat('Result','\',pre_,METHODSET{index},'.mat');
    %     save(filename,'Img');
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
end

%% WNNM
index  = index +1;
if STATUS_WNNM
    C = 0.02;
    tol=1e-6;
    maxIter=100;
    fprintf('\n');
    disp(['performing ',METHODSET{index}, ' ... ']);
    D = zeros(M*N,p);
    for i=1:p
        bandp  = Nmsi(:,:,i);
        D(:,i) = bandp(:); 
    end
    tic;
    [A_hat,~,~] = inexact_alm_WNNMrpca(D,C,tol, maxIter);
    Img     = reshape(A_hat,[M,N,p]);
    Time(index) = toc;
%     save(filename,'Img');
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
end

%% LRMR
index  = index +1;
if STATUS_LRMR
    r = 8;
    slide =20;
    s = 0.05;
    stepsize = 8;
    fprintf('\n');
    disp(['performing ',METHODSET{index}, ' ... ']);
    tic;
    Img    = LRMR_HSI_denoise( Nmsi,r,slide,s,stepsize );
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
end

%% BM4D
index  = index +1;
if STATUS_BM4D
%     sigma=0.1;
    fprintf('\n');
    disp(['performing ',METHODSET{index}, ' ... ']);
    tic
    if gausssigma~=0  %有Gaussian
        [~, Img,sigma_est] = bm4d(Omsi, Nmsi, gausssigma);
    else
        [~, Img,sigma_est] = bm4d(Omsi, Nmsi, 0);
    end
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
end

%% Use TDL
index = index+1;
if STATUS_TDL
    disp(['performing ',METHODSET{index}, ' ... ']);
    tic;
    if gausssigma~=0  %有Gaussian 
        vstbmtf_params.peak_value = 0;
        vstbmtf_params.nsigma = gausssigma;
        vstbmtf_params.peak_value = 0;
        vstbmtf_params.nsigma = mean(sigma_est(:));
    else
        vstbmtf_params.peak_value = 0;
        vstbmtf_params.nsigma = mean(sigma_est(:));
    end
    tic;
    Img = TensorDL(Nmsi, vstbmtf_params);
    Time(index) = toc;
    [n1,n2,n3]=size(Img);
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi(1:n1,1:n2,1:n3),Img);
end

%% WSNM_RPCA
index  = index +1;
if STATUS_WSNM
%     nSig = 0.1;
    if gausssigma~=0  %有Gaussian
        nSig = gausssigma;
    else
        nSig = mean(sigma_est(:));
    end
%     nSig = 20/255;
    pp = 0.7;
    C = 0.007;
    lamb = 1.2;
    delta = 0.1;
    fprintf('\n');
    disp(['performing ',METHODSET{index}, ' ... ']);
    tic;
    Img= RPCA_WSNM_ite(Nmsi,nSig,pp,C,lamb,delta); 
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi(1:end-4,1:end-4,:),Img(1:end-4,1:end-4,:));
%     [n1,n2,n3]=size(Img);
%     [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi(1:n1,1:n2,1:n3),Img);
end

%% LRTV
index  = index +1;
if STATUS_LRTV
    tau = 0.015;% keep tau = 0.01
    lambda =20/sqrt(M*N);% fine tune lambda parameter
    rank = 10;% fine tune rank parameter
    fprintf('\n');
    disp(['performing ',METHODSET{index}, ' ... ']);
    tic
    [Img, ~] = LRTV(Nmsi, tau,lambda, rank);
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
end

%% LRTDTV
index  = index +1;
if STATUS_LRTDTV
    disp(['performing ',METHODSET{index}, ' ... ']);
    if gausssigma~=0  %有Gaussian
        tau    = 1;
        lambda = 100000;beta=100;
        Rank   = [round(M*0.8),round(N*0.8),10];
        fprintf('\n');
        tic
        Img                = LRTDTV_G(Nmsi, tau,lambda,beta,Rank);
        Time(index) = toc;
        [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
    else
        tau    = 1;
        lambda = 10;
        Rank   = [round(M*0.8),round(N*0.8),10];
        fprintf('\n');
        tic
        Img                = LRTDTV(Nmsi, tau,lambda,Rank);
        Time(index) = toc;
        [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
    end
end

%% Enhanced3DTV
index  = index +1;
if STATUS_Enhanced3DTV
    rank   = [13,13,13];
    tau    = 0.004 *sqrt(M*N);
    fprintf('\n');
    disp(['performing ',METHODSET{index}, ' ... ']);
    tic
    Img = EnhancedTV(Nmsi,tau,rank);
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,Img);
end