clc,clear; close all
seed = 2015; 
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

%% read data
load x_dc
x_dc = x_dc(1:100,1:100,1:30);
x=trans255(x_dc);
[w,h,s] = size(x);
%% Generate measurements
ratio  = 0.05; 
N      = h*w;
A      = PermuteWHT_partitioned(N,s,ratio);
% A      = PermuteWHT2(N,s,ratio);
b      = A*x(:);
%% tenspec
clear opts;
opts.maxIter = 3;
opts.tol     = 1e-6;
opts.trX     = x;
opts.gauss_frag = 1;
% Rk           = 4;
% weight       = [0.015,0.015,0.015];
Rk           = 7;
weight       = [0.015,0.015,0.015];
for i=1:1
    fprintf('================== %d/4 ===================\n',i)
    [x_rec,x1,e,psnrpath] = EnhancedTV_CS(A, b, size(x), Rk(i), weight, opts);
    xrec_ttv=reshape(x_rec,[w,h,s]);   
    [mp(i),sm(i),er(i)] = msqia(xrec_ttv/255,opts.trX/255);
end
% weight       = [1,1,0.015];
% [x,x1,e,psnrpath] = TVS_CS(A, b, size(x), rk, weight, opts);
% xrec_ttv=reshape(x,[w,h,s]);   
% [mp,sm,er] = msqia(xrec_ttv/255,opts.trX/255);

% save xrec_ttv.mat xrec_ttv;