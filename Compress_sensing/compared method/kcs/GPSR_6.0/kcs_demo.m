clc,clear; %close all
seed = 2015; 
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));

%% read data
load x_urban
x=trans255(x_urban);
[w,h,s] = size(x);
%% Generate measurements
ratio  = 0.1; 
N      = w*h;
A      = PermuteWHT_partitioned(N,s,ratio);
AT     = A';
y      = A*x(:);

%% Run GPRS
% regularization parameter
tau = 0.05*max(abs(A'*y));
first_tau_factor = 0.5*(max(abs(A'*y))/tau);
steps = 5;
debias = 0;

stopCri=3;
tolA=1.e-5;
t = cputime;
[x_rec,x_debias,objective,times,debias_start,mses,taus]= ...
    GPSR_BB(y,A,AT,tau,...
         'Debias',debias,...
         'Monotone',0,...
         'Initialization',2,...
         'MaxiterA',10000,...
         'True_x',x(:),...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'Verbose',0);
t = cputime - t;
xrec_kcs = reshape(x_rec,[w,h,s]);
i=1;[mp(i),sm(i),er(i)] = msqia1(x/255,(xrec_kcs)/255);
figure;imshow(xrec_kcs(:,:,103),[])