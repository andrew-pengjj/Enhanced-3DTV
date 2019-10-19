clc,clear; close all
seed = 2015; 
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

%% read data
load x_dc
x=trans255(x_dc);
[w,h,s] = size(x);

%% Generate measurements
ratio  = 0.05; 
N      = h*w;
% m      = ceil(N*s*ratio);
% Nn_ext = 2^(ceil(log2(N*s)));
% p      = randperm(Nn_ext);
% picks  = sort(p(1:m),'ascend');  picks(1) = 1;
% p      = randperm(Nn_ext);
% A      = partialRPermuteWHT(h*w*s,picks,p); 
A      = PermuteWHT_partitioned(N,s,ratio);
y=A*x(:);
 
r=5;
K = ceil(0.02*N*s);  
[L,S,err]=sparcs(y, r, K, A, [w*h,s],'svds', 5e-4, 20, 1);
xrec_sparcs=reshape((L+S),w,h,s);
% save xrec_sparcs.mat xrec_sparcs;
