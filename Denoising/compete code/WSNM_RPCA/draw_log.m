load('D:\SR\WSNM-new-exp\WSNM_RPCA_Batch\batch_plot\WNNM_error_matrix.mat');
error_plot = log(error_mat);
%error_plot = flipud(error_plot);

ps = [0.01:0.01:0.5];
pr = [0.01:0.01:0.5];

imagesc(ps, pr, error_plot, [-12,0]);
figure(gcf);colorbar;colormap(jet); 

set(gca, 'XTick', [0.1:0.1:0.5]);
set(gca, 'YTick', [0.1:0.1:0.5]);
axis xy;

% set(gca,'XMinorGrid','on');
% set(gca,'YMinorGrid','on');
% grid on; 