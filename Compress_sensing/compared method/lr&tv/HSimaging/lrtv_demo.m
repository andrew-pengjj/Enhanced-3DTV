clc,clear; close all
seed = 2015; 
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

%% read data
load x_dc
x= trans255(x_dc);
[w,h,s] = size(x);
%% Generate measurements
ratio=0.05;
N      = h*w;
% m      = ceil(N*s*ratio);
% Nn_ext = 2^(ceil(log2(N*s)));
% p      = randperm(Nn_ext);
% picks  = sort(p(1:m),'ascend');  picks(1) = 1;
% p      = randperm(Nn_ext);
% A      = partialRPermuteWHT(h*w*s,picks,p);
A      = PermuteWHT_partitioned(N,s,ratio);
%% Take measurements  
y=A*x(:); 
%% Set parameters and run the algorithm
lambda=0.1; mu=[0.01 150 0.2]; iter=100;                                    % Iteration number
                                   % Iteration number
xrec=funHSI(A,y,[w,h,s],lambda,mu,iter);
xrec_lrtv=reshape(xrec,w,h,s);
% imshow(xrec_lrtv(:,:,20),[]);
i=1;[mp(i),sm(i),er(i)] = msqia(xrec_lrtv/255,x/255)