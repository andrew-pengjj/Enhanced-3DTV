
close all
clear all

% signal length
n = 2^12;
% number of observations
k = 2^10;

spikes = [5:15:245];

trials = 1;

for spi = 1:length(spikes)
    
    timing_GPSR(spi) = 0;
    timing_OMP(spi) = 0;
    mse_GPSR(spi) = 0;
    mse_OMP(spi) = 0;
    timing_OMP_sparselab(spi) = 0;
    mse_OMP_sparselab(spi) = 0;

    n_spikes = spikes(spi);
    fprintf(1,'\n Number of spikes = %d\n',n_spikes);

    for tr = 1:trials
        disp('Creating measurement matrix...');
        R = randn(k,n);
        % normalize the columns of R, as required by OMP
        for ii=1:n
            %ii
            R(:,ii) = R(:,ii)/norm(R(:,ii));
        end
        disp('Finished building matrix');

        % random spikes signal
        f = zeros(n,1);
        q = randperm(n);
        f(q(1:n_spikes)) = sign(randn(n_spikes,1));

        % noisy observations
        sigma = 0.01;  
        y = R*f + sigma*randn(k,1);

        % regularization parameter
        tau = 0.1*max(abs(R'*y));

        debias = 1;
        stopCri = 3;
        tolA = 0.01;

        t0 = cputime;
        [theta_QP_BB_mono,theta_debias_QP_BB_mono,obj_QP_BB_mono,...
            times_QP_BB_mono,debias_start,mses]= ...
                 GPSR_BB(y,R,tau,...
                 'Debias',debias,...
                 'Monotone',1,...
                 'Initialization',0,...
                 'StopCriterion',stopCri,...
                 'ToleranceA',tolA,...
                 'ToleranceD',0.00001,...
                 'Verbose',0);
        t_QP_BB_mono = cputime - t0;
        if debias 
           theta_GPSR = theta_debias_QP_BB_mono;
        else
           theta_GPSR = theta_QP_BB_mono;
        end

        goal_mse = (1/k)*norm(y-R*theta_GPSR)^2;

        timing_GPSR(spi) = timing_GPSR(spi) + t_QP_BB_mono;
        mse_GPSR(spi) = mse_GPSR(spi) + (1/n)*norm(theta_GPSR-f)^2;


        disp('Starting Greedlab OMP')
        t0 = cputime;
        [theta_OMP, err_mse, times_OMP] = greed_omp_qr(y,R,n,...
              'stopCrit','mse','stopTol',goal_mse,'verbose',true); % ,'original',f);
        t_OMP_sparsify = cputime - t0;
        timing_OMP(spi) = timing_OMP(spi) + t_OMP_sparsify;
        mse_OMP(spi) = mse_OMP(spi) + (1/n)*norm(theta_OMP-f)^2;
        disp('Finished Greedlab OMP')


        disp('Starting SparseLab OMP')
        t0 = cputime;
        [sols, iters_sparselab, activationHist] = SolveOMP(R, y, n, n, 0, 0, 0, ...
                                        norm(y-R*theta_GPSR)/norm(y));
        t_OMP_sparselab = cputime - t0;
        timing_OMP_sparselab(spi) = timing_OMP_sparselab(spi) + t_OMP_sparselab;
        mse_OMP_sparselab(spi) = mse_OMP_sparselab(spi) + (1/n)*norm(sols-f)^2;
        disp('Finished SparseLab OMP')


        fprintf(1,'\n-------------------------------------------------\n')   
        fprintf(1,'Problem: n = %g,  k = %g, number of spikes = %g\n',n,k,n_spikes)
        fprintf(1,'Parameters: sigma = %g, tau = %g, debiasing = %g\n',sigma,tau,debias)
        fprintf(1,'All BB algorithms initialized with zeros\n')
        fprintf(1,'-------------------------------------------------\n')


        fprintf(1,'\nQP-BB-monotone; cpu: %6.2f secs (%d iterations)\n',...
                t_QP_BB_mono,length(obj_QP_BB_mono))
        fprintf(1,'final value of the objective function = %6.3e, \nMSE of the solution = %6.3e\n',...
                  obj_QP_BB_mono(end),(1/n)*norm(theta_GPSR-f)^2)      
        fprintf(1,'number of non-zero estimates = %g\n\n',sum(theta_GPSR~=0))

        fprintf(1,'\nOMP (Sparsify); cpu: %6.2f secs (%d iterations)\n',...
                t_OMP_sparsify,length(times_OMP))
        fprintf(1,'MSE of the solution = %6.3e\n',(1/n)*norm(theta_OMP-f)^2)      
        fprintf(1,'number of non-zero estimates = %g\n\n',sum(theta_OMP~=0))

        fprintf(1,'\nOMP (SparseLab); cpu: %6.2f secs (%d iterations)\n',...
                t_OMP_sparselab,iters_sparselab)
        fprintf(1,'MSE of the solution = %6.3e\n',(1/n)*norm(sols-f)^2)      
        fprintf(1,'number of non-zero estimates = %g\n\n',sum(sols~=0))

        fprintf(1,'-------------------------------------------------\n')
        fprintf(1,'-------------------------------------------------\n')
    end  % the trials loop
end % end the number of spikes loop

timing_OMP = timing_OMP / trials;
timing_GPSR = timing_GPSR / trials;

mse_GPSR = mse_GPSR / trials;
mse_OMP = mse_OMP / trials;

mse_OMP_sparselab = mse_OMP_sparselab / trials;
timing_OMP_sparselab = timing_OMP_sparselab / trials;


figure(1)
plot(spikes,mse_GPSR,'r','LineWidth',2)
hold on
plot(spikes,mse_OMP,'b--','LineWidth',2)
plot(spikes,mse_OMP_sparselab,'g-.','LineWidth',2)
hold off
set(gca,'FontName','Times')
set(gca,'FontSize',14)
xlabel('Number of non-zero components')
ylabel('MSE')
title('MSE versus sparseness degree')
legend('GPSR','Sparsify OMP','SparseLab OMP')


figure(2)
plot(spikes,timing_GPSR,'r','LineWidth',2)
hold on
plot(spikes,timing_OMP,'b--','LineWidth',2)
plot(spikes,timing_OMP_sparselab,'g-.','LineWidth',2)
hold off
set(gca,'FontName','Times')
set(gca,'FontSize',14)
xlabel('Number of non-zero components')
ylabel('CPU time (seconds)')
title('CPU time versus sparseness degree')
legend('GPSR','Sparsify OMP','SparseLab OMP')



