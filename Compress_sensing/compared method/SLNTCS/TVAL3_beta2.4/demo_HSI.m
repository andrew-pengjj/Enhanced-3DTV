clc,clear; %close all
seed = 2015; 
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

%% read data
load x_urban
x=trans255(x_urban);
[w,h,s] = size(x);
% x_mat   = zeros(w*s,h);
% for col=1:h
%     for band = 1:s
%         tmp = x(:,col,band);
%         x_mat((band-1)*w+1:band*w,col)=tmp(:);
%     end
% end
%% Generate measurements
ratio  = 0.05; 
N      = w*h;
A      = PermuteWHT_partitioned(N,s,ratio);
y      = A*x(:);
%% Run TVAL3
clear opts
opts.mu    = 2^12;
opts.beta  = 2^5;
opts.tol   = 1E-3;
opts.maxit = 300;
opts.TVnorm = 1;
opts.nonneg = false;
t = cputime;
[U, out] = TVAL3(A,y,N,s,opts);
t = cputime - t;
xrec_ncs = reshape(U,[w,h,s]);
% xrec_ncs=zeros(w,h,s);
% for col=1:h
%     tmp  = U(:,col);
%     tmpM = reshape(tmp,[w,s]);
%     for band = 1:s
%         xrec_ncs(:,col,band)= tmpM(:,band);
%     end
% end
i=1;[mp(i),sm(i),er(i)] = msqia1(x/255,(xrec_ncs)/255);