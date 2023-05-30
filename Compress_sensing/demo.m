clc,clear; close all
clear all;clc;
addpath(genpath('compared method'));
addpath(genpath('quality assess'));
addpath(genpath('Enhanced3DTV in the paper'));
seed = 2015; 
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

%% read data
load x_dc
% x_dc = x_dc(1:100,1:100,1:30);
x=trans255(double(x_dc));
[w,h,s] = size(x);
%% Generate measurements
ratio  = 0.05; 
N      = h*w;
% A      = PermuteWHT2(N,s,ratio);
A      = PermuteWHT_partitioned(N,s,ratio);
AT     = A';
% A      = PermuteWHT2(N,s,ratio);
b      = A*x(:);

% GCS
% fprintf('=========== Enhanced3DTV ============\n');
% clear opts;
% opts.maxIter = 100;
% opts.tol     = 1e-6;
% opts.trX     = x;
% opts.gauss_frag = 1;
% Rk           = 7;
% weight       = [0.015,0.015,0.015];
% [x_gcs,x1,e,psnrpath] = EnhancedTV_CS(A, b, size(x), Rk, weight, opts);
% xrec_ttv=reshape(x_gcs,[w,h,s]);   
% i=1;[mp(i),sm(i),er(i)] = msqia(xrec_ttv/255,opts.trX/255);
% % lrtdtv
% fprintf('=========== lrtdtv ============\n');
% clear opts;
% opts.maxIter = 100;
% opts.tol     = 1e-9;
% opts.trX     = x;
% opts.gauss_frag = 1;
% rk           = [ceil(h*0.6),ceil(w*0.6),6];
% lam=0.001;
% x_rec_re= LrApprReTV(A,b,size(x),rk,lam,opts);
% xrec_rettv=reshape(x_rec_re,size(x)); 
% i=2;[mp(i),sm(i),er(i)] = msqia(xrec_rettv/255,opts.trX/255);
% % lrtv
% fprintf('=========== lrtv ============\n');
% lambda=0.1; mu=[0.01 150 0.2]; iter=100;                                    
% xrec=funHSI(A,b,[w,h,s],lambda,mu,iter);
% xrec_lrtv=reshape(xrec,w,h,s);
% i=3;[mp(i),sm(i),er(i)] = msqia(xrec_lrtv/255,x/255);
% %% sparcs
fprintf('=========== sparcs ============\n');
r=15;%15
K = ceil(0.05*N*s); %0.02 
[L,S,err]=sparcs(b, r, K, A, [w*h,s],'svds', 5e-4, 20, 1);
xrec_sparcs=reshape((L+S),w,h,s);
i=4;[mp(i),sm(i),er(i)] = msqia(xrec_sparcs/255,x/255);
% 
%% KCS
% regularization parameter
% tau = ratio*max(abs(A'*b));%0.05
% first_tau_factor = 0.5*(max(abs(A'*b))/tau);
% steps = 5;%5
% debias = 0;
% 
% stopCri=5;
% tolA=1.e-5;
% fprintf('=========== KCS ============\n');
% [x_rec,x_debias,objective,times,debias_start,mses,taus]= ...
%     GPSR_BB(b,A,AT,tau,...
%          'Debias',debias,...
%          'Monotone',0,...
%          'Initialization',2,...
%          'MaxiterA',100,...
%          'True_x',x(:),...
%          'StopCriterion',stopCri,...
%        	 'ToleranceA',tolA,...
%          'Verbose',0);
% xrec_kcs = reshape(x_rec,[w,h,s]);
% i=5;[mp(i),sm(i),er(i)] = msqia(xrec_kcs/255,x/255);

%% NCS
% clear opts
% opts.mu    = 2^12;
% opts.beta  = 2^7;%2^5
% opts.tol   = 1E-3;
% opts.maxit = 300;
% opts.TVnorm = 1;
% opts.nonneg = false;
% t = cputime;
% fprintf('=========== NCS ============\n');
% [U, out] = TVAL3(A,b,N,s,opts);
% t = cputime - t;
% xrec_ncs = reshape(U,[w,h,s]);
% i=6;[mp(i),sm(i),er(i)] = msqia(xrec_ncs/255,x/255);
