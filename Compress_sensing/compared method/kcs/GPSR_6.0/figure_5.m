% code to generate Figure 5.

% clear everything
close all
clear all
pack

% number of trials for each size (set to 10 for the final run)
trials = 10;

% signal lengths to consider
ns = [1e4 3e4 6e4 1e5 3e5 6e5 1e6];

% initialize arrays for average time computations
ts_l1_ls = zeros(1,length(ns));
ts_QP_BB_mono = zeros(1,length(ns));
ts_QP_BB_notmono = zeros(1,length(ns));
ts_GPSR_Basic = zeros(1,length(ns));
ts_IST = zeros(1,length(ns));

for len = 1:length(ns)

    % set signal length
    n = ns(len);

    % set number of observations (matrix rows)
    k = n/10;

    % set number of spikes
    n_spikes = 0.25*n;


    for tr = 1:trials
        
    % clear the previous matrix, to save memory
    clear R

    % random +/- 1 signal
    f = zeros(n,1);
    q = randperm(n);
    f(q(1:n_spikes)) = sign(randn(n_spikes,1));

    %  create measurement matrix
    disp('Creating measurement matrix...');
    R = sprandn(k,n,30/n);
    disp('Done.');

    % noisy observations
    sigma = 0.01;
    y = R*f + sigma*randn(k,1);

    % Same selection of tau as in the l1_ls paper
    tau = 0.1*max(abs(R'*y));

    t0 = cputime;
    [theta_l1_ls,status,history] = l1_ls(R,y,2*tau,0.001);
    t_l1_ls = cputime - t0;


    t0 = cputime;
    [theta_QP_BB_mono,theta_debias_QP_BB_mono,obj_QP_BB_mono,...
        times_QP_BB_mono,debias_start,mses]= ...
             GPSR_BB(y,R,tau,...
             'Debias',0,...
             'AT',R', ... 
             'Monotone',1,...
             'Initialization',0,...
             'StopCriterion',4,...
             'ToleranceA',0.5*history(2,end),...
             'ToleranceD',0.0001);
    t_QP_BB_mono = cputime - t0;


    t0 = cputime;
    [theta_QP_BB_notmono,theta_debias_QP_BB_notmono,obj_QP_BB_notmono,...
        times_QP_BB_notmono,debias_start,mses]= ...
             GPSR_BB(y,R,tau,...
             'Debias',0,...
             'AT',R',... 
             'Monotone',0,...
             'Initialization',0,...
             'StopCriterion',4,...
             'ToleranceA',0.5*history(2,end),...
             'ToleranceD',0.0001);
    t_QP_BB_notmono = cputime - t0;


    t0 = cputime;      
    [theta_GPSR_Basic,theta_debias_GPSR_Basic,obj_GPSR_Basic,...
        times_GPSR_Basic,debias_start,mses]= ...
         GPSR_Basic(y,R,tau,...
             'Debias',0,...
             'AT',R',...
             'Initialization',0,...
             'StopCriterion',4,...
             'ToleranceA',0.5*history(2,end),...
            'ToleranceD',0.0001);
    t_GPSR_Basic = cputime - t0;

    % IST assumes that matrix R has unit norm; 
    % because of that, we rescale the problem accordingly
    %
    gamma = svds(R,1);
    t0 = cputime;
    [theta_IST,theta_debias_IST,obj_IST,...
        times_IST,debias_start,mses]= ...
             IST(y/gamma,R/gamma,tau/gamma^2,...
             'Debias',0,...
             'AT',R', ...
             'Initialization',0,...
             'StopCriterion',4,...
             'ToleranceA',obj_QP_BB_mono(end)/gamma^2,...
             'ToleranceD',0.0001);
    t_IST = cputime - t0;


    ts_l1_ls(len) = ts_l1_ls(len) + t_l1_ls;
    ts_QP_BB_mono(len) = ts_QP_BB_mono(len) + t_QP_BB_mono;
    ts_QP_BB_notmono(len) = ts_QP_BB_notmono(len) + t_QP_BB_notmono;
    ts_GPSR_Basic(len) = ts_GPSR_Basic(len) + t_GPSR_Basic;
    ts_IST(len) = ts_IST(len) + t_IST;

    end

end

ts_l1_ls = ts_l1_ls / trials;
ts_QP_BB_mono = ts_QP_BB_mono / trials;
ts_QP_BB_notmono = ts_QP_BB_notmono / trials;
ts_GPSR_Basic = ts_GPSR_Basic / trials;
ts_IST = ts_IST / trials;

% find the empirical asymptotical exponents
pp_l1_ls = polyfit(log(ns),log(ts_l1_ls),1);
pp_QP_BB_mono = polyfit(log(ns),log(ts_QP_BB_mono),1);
pp_QP_BB_notmono = polyfit(log(ns),log(ts_QP_BB_notmono),1);
pp_GPSR_Basic = polyfit(log(ns),log(ts_GPSR_Basic),1);
pp_IST = polyfit(log(ns),log(ts_IST),1);

% array of points for plotting the fitted lines
nnss = exp([min(log(ns)):.01:max(log(ns))]);

% compute fitted lines
fit_l1_ls = exp(pp_l1_ls(2))*nnss.^pp_l1_ls(1);
fit_QP_BB_mono = exp(pp_QP_BB_mono(2))*nnss.^pp_QP_BB_mono(1);
fit_QP_BB_notmono = exp(pp_QP_BB_notmono(2))*nnss.^pp_QP_BB_notmono(1);
fit_GPSR_Basic = exp(pp_GPSR_Basic(2))*nnss.^pp_GPSR_Basic(1);
fit_IST = exp(pp_IST(2))*nnss.^pp_IST(1);

% do the plots
loglog(ns,ts_GPSR_Basic,'r+')
hold on
loglog(ns,ts_QP_BB_mono,'go')
loglog(ns,ts_QP_BB_notmono,'b*')
loglog(ns,ts_l1_ls,'mx')
loglog(ns,ts_IST,'cs')
loglog(nnss,fit_QP_BB_mono,'g','LineWidth',1.8)
loglog(nnss,fit_QP_BB_notmono,'b','LineWidth',1.8)
loglog(nnss,fit_GPSR_Basic,'r','LineWidth',1.8)
loglog(nnss,fit_l1_ls,'m','LineWidth',1.8)
loglog(nnss,fit_IST,'c','LineWidth',1.8)


leg = legend(strcat('GPSR-BB monotone (\alpha = ',sprintf(' %0.3g)',pp_QP_BB_mono(1))),...
             strcat('GPSR-BB non-monotone (\alpha = ',sprintf(' %0.3g)',pp_QP_BB_notmono(1))),...
             strcat('GPSR-Basic (\alpha = ',sprintf(' %0.3g)',pp_GPSR_Basic(1))),...
             strcat('l1-ls (\alpha = ',sprintf(' %0.3g)',pp_l1_ls(1))),...
             strcat('IST (\alpha = ',sprintf(' %0.3g)',pp_IST(1))) )
set(gca,'FontName','Times')
set(gca,'FontSize',14)
set(leg,'FontName','Times')
set(leg,'FontSize',12)
xlabel('Problem size (n)')
ylabel('Average CPU time')
title('Empirical asymptotic exponents O(n^\alpha)')
hold off





