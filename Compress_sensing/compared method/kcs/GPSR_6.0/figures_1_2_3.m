
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
sigma = 0.01;
y = hR(f) + sigma*randn(k,1);

% regularization parameter
tau = 0.1*max(abs(R'*y));

debias = 0;

%
% Need to call l1_ls with tau/2 because it assumes
% a different objective:
%  || y - R*x||_2^2 + tau ||x||_1
% instead of the one assumed by GPSR which is 
%  (1/2)*|| y - R*x||_2^2 + tau ||x||_1
%
[x_l1_ls,status,history] = l1_ls(R,y,2*tau,0.01);
t_l1_ls = history(7,end);

stopCri = 4;

%
% for the reason explained above the objective
% function of l1_ls is twice that of GPSR, so we
% need the following correction.
tolA = history(2,end)/2;

[x_BB_mono,x_debias_BB_mono,obj_BB_mono,...
    times_BB_mono,debias_start_BB_mono,mse]= ...
         GPSR_BB(y,hR,tau,...
         'Debias',debias,...
         'AT',hRt,... 
         'Monotone',1,...
         'Initialization',0,...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'ToleranceD',0.001);
t_BB_mono = times_BB_mono(end);

[x_BB_notmono,x_debias_BB_notmono,obj_BB_notmono,...
    times_BB_notmono,debias_start_BB_notmono,mse]= ...
         GPSR_BB(y,hR,tau,...
         'Debias',debias,...
         'AT',hRt,... 
         'Monotone',0,...
         'Initialization',0,...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'ToleranceD',0.0001);
t_BB_notmono = times_BB_notmono(end);

[x_GPSR_Basic,x_debias_GPSR_Basic,obj_GPSR_Basic,...
    times_GPSR_Basic,debias_start_Basic,mse]= ...
	 GPSR_Basic(y,hR,tau,...
         'Debias',debias,...
         'AT',hRt,... 
         'Initialization',0,...
    	 'StopCriterion',stopCri,...
	     'ToleranceA',tolA,...
         'ToleranceD',0.0001);
t_GPSR_Basic = times_GPSR_Basic(end);


fprintf(1,'\n\n-------------------------------------------------\n')   
fprintf(1,'-------------------------------------------------\n')   
fprintf(1,'Problem: n = %g,  k = %g, number of spikes = %g\n',n,k,n_spikes)
fprintf(1,'Parameters: sigma = %g, tau = %g, debiasing = %g\n',sigma,tau,debias)
fprintf(1,'All BB algorithms initialized with zeros\n')
fprintf(1,'-------------------------------------------------\n')


fprintf(1,'\nQP-BB-monotone; cpu: %6.2f secs (%d iterations)\n',...
        t_BB_mono,length(obj_BB_mono))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          obj_BB_mono(end),(1/n)*norm(x_BB_mono-f)^2)      
fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_BB_mono~=0))


fprintf(1,'QP-BB-non-monotone; cpu: %6.2f secs (%d iterations)\n',...
        t_BB_notmono,length(obj_BB_notmono))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          obj_BB_notmono(end),(1/n)*norm(x_BB_notmono-f)^2)      
fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_BB_notmono~=0))

fprintf(1,'QP GPSR-Basic; cpu: %6.2f secs (%d iterations)\n',...
        t_GPSR_Basic,length(obj_GPSR_Basic))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          obj_GPSR_Basic(end),(1/n)*norm(x_GPSR_Basic-f)^2)
fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_GPSR_Basic~=0))


fprintf(1,'l1_ls; cpu: %6.2f secs (%d iterations)\n',t_l1_ls,length(history(2,:)))
fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
          history(2,end)/2,(1/n)*norm(x_l1_ls-f)^2)
fprintf(1,'number of non-zero estimates = %g\n\n',sum(x_l1_ls~=0))

fprintf(1,'-------------------------------------------------\n')
fprintf(1,'-------------------------------------------------\n')


% ================= Plotting results ==========

% This is the first plot in figure 2 of the paper
figure(1)
plot(obj_BB_mono,'LineWidth',2)
hold on
plot(obj_BB_notmono,'r--','LineWidth',2)
plot(obj_GPSR_Basic,'k:','LineWidth',2)
legend('GPSR-BB monotone','GPSR-BB non-monotone','GPSR-Basic')
set(gca,'FontName','Times','FontSize',16)
xlabel('Iterations')
ylabel('Objective function')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off

% This is the second plot in figure 2 of the paper
figure(2)
plot(times_BB_mono,obj_BB_mono,'LineWidth',2)
hold on
plot(times_BB_notmono,obj_BB_notmono,'r--','LineWidth',2)
plot(times_GPSR_Basic,obj_GPSR_Basic,'k:','LineWidth',2)
legend('GPSR-BB monotone','GPSR-BB non-monotone','GPSR-Basic')
set(gca,'FontName','Times','FontSize',16)
xlabel('CPU time (seconds)')
ylabel('Objective function')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off

% this is a plot, not in the paper, which also includes l1_ls
figure(6)
plot(times_BB_mono,obj_BB_mono,'LineWidth',2)
hold on
plot(times_BB_notmono,obj_BB_notmono,'r--','LineWidth',2)
plot(times_GPSR_Basic,obj_GPSR_Basic,'k:','LineWidth',2)
plot(history(7,:),history(2,:)/2,'g-.','LineWidth',2)
legend('GPSR-BB monotone','GPSR-BB non-monotone','GPSR-Basic','l1-ls')
set(gca,'FontName','Times','FontSize',16)
xlabel('CPU time (seconds)')
ylabel('Objective function')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off



% Now we run with debiasing, do figure 3 of the paper
[x_BB_mono,x_debias_BB_mono,obj_BB_mono,...
    times_BB_mono,debias_start_BB_mono,mse_BB_mono]= ...
         GPSR_BB(y,hR,tau,...
         'Debias',1,...
         'AT',hRt,... 
         'True_x',f,...
         'Monotone',1,...
         'Initialization',0,...
         'StopCriterion',stopCri,...
       	 'ToleranceA',tolA,...
         'ToleranceD',0.001);

% This is the first plot of figure 3 of the paper
figure(3)
plot(times_BB_mono,obj_BB_mono,'LineWidth',2)
v = axis;
line([times_BB_mono(debias_start_BB_mono),...
      times_BB_mono(debias_start_BB_mono)],...
      [v(3),v(4)],'LineStyle','--','Color',[1 0 0])
set(gca,'FontName','Times','FontSize',16)
text(times_BB_mono(debias_start_BB_mono)+0.01*(v(2)-v(1)),...
     v(3)+0.8*(v(4)-v(3)),'Debiasing','FontName','Times','FontSize',16)
xlabel('CPU time')
ylabel('Objective function')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off

% This is the second plot of figure 3 of the paper
figure(4)
plot(times_BB_mono,mse_BB_mono,'LineWidth',2)
v = axis;
line([times_BB_mono(debias_start_BB_mono),...
      times_BB_mono(debias_start_BB_mono)],...
      [v(3),v(4)],'LineStyle','--','Color',[1 0 0])
text(times_BB_mono(debias_start_BB_mono)+0.01*(v(2)-v(1)),...
     v(3)+0.8*(v(4)-v(3)),'Debiasing','FontName','Times','FontSize',16)
set(gca,'FontName','Times','FontSize',16)
xlabel('CPU time (seconds)')
ylabel('MSE')
title(sprintf('n=%d, k=%d, tau=%g',n,k,tau))
hold off


% This is figure 1 of the paper.
debias = 1

figure(5)
scrsz = get(0,'ScreenSize');
set(5,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
if debias
    subplot(4,1,1)
else
    subplot(3,1,1)
end
plot(f,'LineWidth',1.1)
top = max(f(:));
bottom = min(f(:));
v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Original (n = %g, number of nonzeros = %g)',n,n_spikes))
axis(v)

if debias
    subplot(4,1,2)
else
    subplot(3,1,2)
end
plot(x_BB_mono,'LineWidth',1.1)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
axis(v)
title(sprintf('GPSR reconstruction (k = %g, tau = %5.3g, MSE = %5.3g)',...
    k,tau,(1/n)*norm(x_BB_mono-f)^2))

if debias
    subplot(4,1,3)
    plot(x_debias_BB_mono,'LineWidth',1.1)
    set(gca,'FontName','Times')
    set(gca,'FontSize',14)
    top = max(f(:));
    bottom = min(f(:));
    v = [0 n+1 bottom-0.15*(top-bottom)  top+0.15*((top-bottom))];
    axis(v)
    title(sprintf(...
     'Debiased (MSE = %0.4g)',(1/n)*norm(x_debias_BB_mono-f)^2))
end

if debias
    subplot(4,1,4)
else
    subplot(3,1,3)
end
pseudo = hRt(y);
plot(pseudo,'LineWidth',1.1)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Minimum norm solution (MSE = %0.4g)',(1/n)*norm(pseudo-f)^2))
top = max(pseudo(:));
bottom = min(pseudo(:));
v = [0 n+1 bottom-0.15*(top-bottom)  top+0.15*((top-bottom))];
axis(v)




