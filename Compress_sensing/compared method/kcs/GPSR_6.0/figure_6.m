
close all
clear all

% signal length
n = 2^13;
% number of spikes to put down
% n_spikes = floor(.01*n);
n_spikes = 100;
% number of observations to make
k = 2^9;

% random +/- 1 signal
f = zeros(n,1);
q = randperm(n);
f(q(1:n_spikes)) = sign(randn(n_spikes,1));
%f(q(1:n_spikes)) = randn(n_spikes,1);

% measurement matrix
disp('Creating measurement matrix...');
R = randn(k,n);
R = orth(R')';
disp('Done.');

% noisy observations
sigma = 0.01;
y = R*f + sigma*randn(k,1);

taus = [0.05:0.025:0.275]*max(abs(R'*y));
trials = 10;

iters = [];
warm_durations = [];
nz = [];
cold_iters = [];
cold_durations = [];

debias = 0;
stopCri = 3;
tolA = 0.001;


for i=1:length(taus)
    
    cold_durations(i) = 0;
    warm_durations(i) = 0;

    for tr = 1:trials

        if i==1
           init = 0;
        else
           init = theta;
        end
        [theta,theta_debias,obj,times,debias_start,mses]= ...
                 GPSR_BB(y,R,taus(i),...
                 'Debias',debias,...
                 'AT',R',... %'TrueTheta',f,...
                 'Monotone',1,...
                 'Initialization',init,...
                 'StopCriterion',stopCri,...
                 'ToleranceA',tolA,...
                 'ToleranceD',0.0001);
        iters(i) = length(times);     
        warm_durations(i) =  warm_durations(i) + times(end);

        [theta,theta_debias,obj,times2,debias_start,mses]= ...
                 GPSR_BB(y,R,taus(i),...
                 'Debias',debias,...
                 'AT',R',... %'TrueTheta',f,...
                 'Monotone',1,...
                 'Initialization',0,...
                 'StopCriterion',stopCri,...
                 'ToleranceA',tolA,...
                 'ToleranceD',0.0001);
        if i==1
           cold_iters(i) = iters(i);     
           cold_durations(i) =  cold_durations(i) + times(end);
        else
           cold_iters(i) = length(times);     
           cold_durations(i) =  cold_durations(i) + times2(end);
        end
    end
end

cold_durations = cold_durations / trials;
warm_durations = warm_durations / trials;


% ================= Plotting results ==========

figure(1)
plot(taus,cold_durations,'bs')
hold on
plot(taus,warm_durations,'r*')
hold off
legend('Cold start','Warm start')
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title('Warm versus cold start')
xlabel('Value of \tau')
ylabel('CPU time')

