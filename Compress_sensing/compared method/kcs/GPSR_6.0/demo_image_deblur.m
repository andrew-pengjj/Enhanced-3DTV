% This demo shows the reconstruction of a sparse image
% of randomly placed plus an minus ones
% 
close all
clear
clf

f = double(imread('Camera.tif'));
[m n] = size(f);

scrsz = get(0,'ScreenSize');
figure(1)
set(1,'Position',[10 scrsz(4)*0.05 scrsz(3)/4 0.85*scrsz(4)])
subplot(3,1,1)
imagesc(f)
colormap(gray(255))
axis off
axis equal
title('Original image','FontName','Times','FontSize',14)



% create observation operator; in this case 
% it will be a blur function composed with an
% inverse weavelet transform
disp('Creating observation operator...');

middle = n/2 + 1;

% uncomment the following lines for Experiment 1 (see paper)
% sigma = 0.56;
% h = zeros(size(f));
% for i=-4:4
%    for j=-4:4
%       h(i+middle,j+middle)= 1; 
%    end
% end


% uncomment the following lines for Experiment 2 (see paper)
sigma = sqrt(2);
h = zeros(size(f));
for i=-4:4
   for j=-4:4
      h(i+middle,j+middle)= (1/(1+i*i+j*j));
   end
end

% uncomment the following lines for Experiment 3 (see paper)
% sigma = sqrt(8);
% h = zeros(size(f));
% for i=-4:4
%    for j=-4:4
%       h(i+middle,j+middle)= (1/(1+i*i+j*j));
%    end
% end

% % center and normalize the blur
h = fftshift(h);   
h = h/sum(h(:));

% definde the function handles that compute 
% the blur and the conjugate blur.
R = @(x) real(ifft2(fft2(h).*fft2(x)));
RT = @(x) real(ifft2(conj(fft2(h)).*fft2(x)));

% define the function handles that compute 
% the products by W (inverse DWT) and W' (DWT)
wav = daubcqf(2);
W = @(x) midwt(x,wav,3);
WT = @(x) mdwt(x,wav,3);

%Finally define the function handles that compute 
% the products by A = RW  and A' =W'*R' 
A = @(x) R(W(x));
AT = @(x) WT(RT(x));

% generate noisy blurred observations
y = R(f) + sigma*randn(size(f));
figure(1)
subplot(3,1,2)
imagesc(y)
colormap(gray(255))
axis off
axis equal
title('Blurred image','FontName','Times','FontSize',14)

% regularization parameter
tau = .35;

% set tolA
tolA = 1.e-3;

% Run IST until the relative change in objective function is no
% larger than tolA
[theta_ist,theta_debias,obj_IST,times_IST,debias_s,mses_IST]= ...
	 IST(y,A,tau,...
	'Debias',0,...
	'AT',AT,... 
    'True_x',WT(f),...
	'Initialization',AT(y),...
	'StopCriterion',1,...
	'ToleranceA',tolA);

% Now, run the GPSR functions, until they reach the same value
% of objective function reached by IST.
[theta,theta_debias,obj_GPSR_Basic,times_GPSR_Basic,debias_s,mses_GPSR_Basic]= ...
	GPSR_Basic(y,A,tau,...
	'Debias',0,...
	'AT',AT,... 
    'True_x',WT(f),...
	'Initialization',AT(y),...
	'StopCriterion',4,...
	'ToleranceA',obj_IST(end));

% [theta,theta_debias,obj_QP_BB_notmono,times_QP_BB_notmono,debias_s,mses_QP_BB_notmono]= ...
% 	GPSR_BB(y,A,tau,...
% 	'AT', AT,...
% 	'Debias',0,...
% 	'Initialization',AT(y),...
%     'True_x',WT(f),...
% 	'Monotone',0,...
% 	'StopCriterion',4,...
% 	'ToleranceA',obj_IST(end));

% [theta,theta_debias,obj_QP_BB_mono,times_QP_BB_mono,debias_start,mses_QP_BB_mono]= ...
% 	GPSR_BB(y,A,tau,...
% 	'AT', AT,...
% 	'Debias',0,...
% 	'Initialization', AT(y),...
%     'True_x',WT(f),...
% 	'Monotone',1,...
% 	'StopCriterion',4,...
% 	'ToleranceA',obj_IST(end));
% 
% 
% % ================= Plotting results ==========
figure(1)
subplot(3,1,3)
if prod(size(theta_debias))~=0
   imagesc(W(theta_debias))
else
   imagesc(W(theta_ist))
end
colormap(gray)
axis off
axis equal
title('Deblurred image','FontName','Times','FontSize',14)
%   
% figure(2)
% scrsz = get(0,'ScreenSize');
% lft = 0.55*scrsz(3)-10;
% btm = 0.525*scrsz(4);
% wdt = 0.45*scrsz(3);
% hgt = 0.375*scrsz(4);
% 
% set(2,'Position',[lft btm wdt hgt])
% plot(obj_QP_BB_mono,'b','LineWidth',1.8);
% hold on
% plot(obj_QP_BB_notmono,'r--','LineWidth',1.8);
% plot(obj_GPSR_Basic,'g:','LineWidth',1.8);
% plot(obj_IST,'m-.','LineWidth',1.8);
% hold off
% leg = legend('GPSR-BB monotone','GPSR-BB non-monotone','GPSR-Basic','IST')
% v = axis;
% if debias_start ~= 0
%    line([debias_start,debias_start],[v(3),v(4)],'LineStyle',':')
%    text(debias_start+0.01*(v(2)-v(1)),...
%    v(3)+0.8*(v(4)-v(3)),'Debiasing')
% end
% ylabel('Objective function','FontName','Times','FontSize',16)
% xlabel('Iterations','FontName','Times','FontSize',16)
% 
% 
% set(leg,'FontName','Times')
% set(leg,'FontSize',16)
% set(gca,'FontName','Times')
% set(gca,'FontSize',16)
%     
% 
% figure(3)
% scrsz = get(0,'ScreenSize');
% lft = 0.55*scrsz(3)-10;
% btm = 0.025*scrsz(4);
% wdt = 0.45*scrsz(3);
% hgt = 0.375*scrsz(4);
% 
% set(3,'Position',[lft btm wdt hgt])
% plot(times_QP_BB_mono,obj_QP_BB_mono,'b','LineWidth',1.8);
% hold on
% plot(times_QP_BB_notmono,obj_QP_BB_notmono,'r--','LineWidth',1.8);
% plot(times_GPSR_Basic,obj_GPSR_Basic,'g:','LineWidth',1.8);
% plot(times_IST,obj_IST,'m-.','LineWidth',1.8);
% hold off
% leg = legend('GPSR-BB monotone','GPSR-BB non-monotone','GPSR-Basic','IST')
% v = axis;
% if debias_start ~= 0
%    line([debias_start,debias_start],[v(3),v(4)],'LineStyle',':')
%    text(debias_start+0.01*(v(2)-v(1)),...
%    v(3)+0.8*(v(4)-v(3)),'Debiasing')
% end
% ylabel('Objective function','FontName','Times','FontSize',16)
% xlabel('CPU time (seconds)','FontName','Times','FontSize',16)
% 
% 
% set(leg,'FontName','Times')
% set(leg,'FontSize',16)
% set(gca,'FontName','Times')
% set(gca,'FontSize',16)
%     
% 
% figure(4)
% 
% plot(times_QP_BB_mono,mses_QP_BB_mono,'b','LineWidth',1.8);
% hold on
% plot(times_QP_BB_notmono,mses_QP_BB_notmono,'r--','LineWidth',1.8);
% plot(times_GPSR_Basic,mses_GPSR_Basic,'g:','LineWidth',1.8);
% plot(times_IST,mses_IST,'m-.','LineWidth',1.8);
% hold off
% leg = legend('GPSR-BB monotone','GPSR-BB non-monotone','GPSR-Basic','IST')
% v = axis;
% if debias_start ~= 0
%    line([debias_start,debias_start],[v(3),v(4)],'LineStyle',':')
%    text(debias_start+0.01*(v(2)-v(1)),...
%    v(3)+0.8*(v(4)-v(3)),'Debiasing')
% end
% ylabel('Deconvolution MSE','FontName','Times','FontSize',16)
% xlabel('CPU time (seconds)','FontName','Times','FontSize',16)
% 
% 
% set(leg,'FontName','Times')
% set(leg,'FontSize',16)
% set(gca,'FontName','Times')
% set(gca,'FontSize',16)
%     
% 
% 
% 
% 
