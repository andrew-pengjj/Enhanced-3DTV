
close all
% n is the original signal length
n = 2^12;

% k is number of observations to make
k = 2^10;

% number of spikes to put down
% n_spikes = floor(.01*n);
n_spikes = 160;


% random +/- 1 signal
f = zeros(n,1);
q = randperm(n);
f(q(1:n_spikes)) = sign(randn(n_spikes,1));
%f(q(1:n_spikes)) = randn(n_spikes,1);

% measurement matrix
disp('Creating measurement matrix...');
R = randn(k,n);

% orthonormalize rows
R = orth(R')';

if n == 8192  
   % in this case, we load a precomputed
   % matrix to save some time
   load Rmatrix_2048_8192.mat
end
%
disp('Finished creating matrix');

hR = @(x) R*x;
hRt = @(x) R'*x;

% noisy observations
sigma = 0; 
y = hR(f) + sigma*randn(k,1);

% regularization parameter
tau = 0.005*max(abs(R'*y));

first_tau_factor = 0.8*(max(abs(R'*y))/tau);
steps = 5;

debias = 0;

stopCri=3;
tolA=1.e-5;

disp('Starting GPSR BB monotonic')
[x_BB_mono,x_debias_BB_mono,obj_BB_mono,...
    times_BB_mono,debias_start_BB_mono,mse_BB_mono]= ...
         GPSR_BB(y,R,R',tau,...
         'Debias',debias,...
         'Monotone',1,...
         'Initialization',0,...
         'MaxiterA',10000,...
         'True_x',f,...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'Verbose',0);
t_BB_mono = times_BB_mono(end);


disp('Starting GPSR BB monotonic with continuation')
[x_BB_mono_cont,x_debias_BB_mono_cont,obj_BB_mono_cont,...
    times_BB_mono_cont,debias_start_BB_mono,mse_BB_mono_cont]= ...
         GPSR_BB(y,R,tau,...
         'Debias',debias,...
         'Continuation',1,...
         'ContinuationSteps',steps,...
         'FirstTauFactor',first_tau_factor,...
         'Monotone',1,...
         'Initialization',0,...
         'True_x',f,...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'Verbose',0);
t_BB_mono_cont = times_BB_mono_cont(end);

disp('Starting GPSR Basic')
[x_Basic,x_debias_Basic,obj_Basic,...
    times_Basic,debias_start_Basic,mse_Basic]= ...
         GPSR_Basic(y,R,tau,...
         'Debias',debias,...
         'Initialization',0,...
         'MaxiterA',10000,...
         'True_x',f,...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'Verbose',0);
t_Basic = times_Basic(end);

disp('Starting GPSR Basic with continuation')
[x_Basic_cont,x_debias_Basic_cont,obj_Basic_cont,...
    times_Basic_cont,debias_start_Basic,mse_Basic_cont]= ...
         GPSR_Basic(y,R,tau,...
         'Debias',debias,...
         'Continuation',1,...
         'ContinuationSteps',steps,...
         'FirstTauFactor',first_tau_factor,...
         'Initialization',0,...
         'True_x',f,...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'Verbose',0);
t_Basic_cont = times_Basic_cont(end);



fprintf(1,'\n\n-------------------------------------------------\n')   
fprintf(1,'-------------------------------------------------\n')   
fprintf(1,'Problem: n = %g,  k = %g, number of spikes = %g\n',n,k,n_spikes)
fprintf(1,'Parameters: sigma = %g, tau = %g, debiasing = %g\n',sigma,tau,debias)
fprintf(1,'All BB algorithms initialized with zeros\n')
fprintf(1,'-------------------------------------------------\n')


fprintf(1,'\nGPSR-BB monotone; cpu: %6.2f secs (%d iterations)\n',...
        t_BB_mono,length(obj_BB_mono))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          obj_BB_mono(end),(1/n)*norm(x_BB_mono-f)^2)      


fprintf(1,'\nGPSR-BB monotone continuation; cpu: %6.2f secs (%d iterations)\n',...
        t_BB_mono_cont,length(obj_BB_mono_cont))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          obj_BB_mono_cont(end),(1/n)*norm(x_BB_mono_cont-f)^2)      

fprintf(1,'\nGPSR-Basic; cpu: %6.2f secs (%d iterations)\n',...
        t_Basic,length(obj_Basic))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          obj_Basic(end),(1/n)*norm(x_Basic-f)^2)      

fprintf(1,'\nGPSR-Basic-continuation; cpu: %6.2f secs (%d iterations)\n',...
        t_Basic_cont,length(obj_Basic_cont))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          obj_Basic_cont(end),(1/n)*norm(x_Basic_cont-f)^2)      

fprintf(1,'-------------------------------------------------\n')
fprintf(1,'-------------------------------------------------\n')


% ================= Plotting results ==========

figure(1)
plot(mse_BB_mono,'LineWidth',2)
hold on
plot(mse_BB_mono_cont,'k:','LineWidth',2)
%plot(history(2,:)/2,'g-.','LineWidth',2)
%legend('GPSR-BB monotone','GPSR-BB monotone with continuation','l1-ls')
legend('GPSR-BB monotone','GPSR-BB monotone with continuation')
set(gca,'FontName','Times','FontSize',16)
xlabel('Iterations')
ylabel('MSE')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off

figure(2)
plot(times_BB_mono,mse_BB_mono,'LineWidth',2)
hold on
plot(times_BB_mono_cont,mse_BB_mono_cont,'k:','LineWidth',2)
%plot(history(7,:),history(2,:)/2,'g-.','LineWidth',2)
%legend('GPSR-BB monotone','GPSR-BB monotone with continuation','l1-ls')
legend('GPSR-BB monotone','GPSR-BB monotone with continuation')
set(gca,'FontName','Times','FontSize',16)
xlabel('CPU time (seconds)')
ylabel('MSE')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off

figure(3)
plot(mse_Basic,'LineWidth',2)
hold on
plot(mse_Basic_cont,'k:','LineWidth',2)
%plot(history(2,:)/2,'g-.','LineWidth',2)
%legend('GPSR-Basic','GPSR-Basic with continuation','l1-ls')
legend('GPSR-Basic','GPSR-Basic with continuation')
set(gca,'FontName','Times','FontSize',16)
xlabel('Iterations')
ylabel('MSE')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off

% this is a plot, not in the paper, which also includes l1_ls
figure(4)
plot(times_Basic,mse_Basic,'LineWidth',2)
hold on
plot(times_Basic_cont,mse_Basic_cont,'k:','LineWidth',2)
%plot(history(7,:),history(2,:)/2,'g-.','LineWidth',2)
%legend('GPSR-Basic','GPSR-Basic with continuation','l1-ls')
legend('GPSR-Basic','GPSR-Basic with continuation')
set(gca,'FontName','Times','FontSize',16)
xlabel('CPU time (seconds)')
ylabel('MSE')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off
% 


figure(7)
scrsz = get(0,'ScreenSize');
set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(3,1,1)
plot(f,'LineWidth',1.1)
top = max(f(:));
bottom = min(f(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Original (n = %g, number of nonzeros = %g)',n,n_spikes))
axis(v)

subplot(3,1,2)
plot(x_BB_mono,'LineWidth',1.1)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
axis(v)
title(sprintf('GPSR-BB (%g iterations, k = %g, tau = %5.3g, MSE = %5.3g)',...
    length(obj_BB_mono),k,tau,(1/n)*norm(x_BB_mono-f)^2))

subplot(3,1,3)
plot(x_BB_mono_cont,'LineWidth',1.1)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('GPSR-BB with continuation (%g iterations, MSE = %0.4g)',...
      length(obj_BB_mono_cont),(1/n)*norm(x_BB_mono_cont-f)^2))
axis(v)




