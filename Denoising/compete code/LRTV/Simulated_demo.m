% clear all;
% close all;
% clc;	

addpath ./prox_operators
addpath ./PROPACK
%% real data experiment
% % load 'Indian_145_220';
% % load 'hypercube166_200';
% % load 'URBAN_200_307';
%% simulated experiment
%----------------------------------------------load image---------------------------------------------------------
load simu_indian
OriData3 = simu_indian;
ratio = 0.15*ones(1,224);            % for case 1
noiselevel = 0.075*ones(1,224);      % for case 1

% load Simu_ratio                       % for case 2
% load Simu_noiselevel                  % for case 2
%----------------------------------------------noise simulated---------------------------------------------------------
oriData3_noise = OriData3;
[M N p] = size(OriData3);
% Gaussian noise
for i =1:p
     oriData3_noise(:,:,i)=OriData3(:,:,i)  + noiselevel(i)*randn(M,N);
%        mixed(vv,:)=addnoise(mixed(vv,:),'mg',ratio);
end
%S&P noise
for i =1:p
     oriData3_noise(:,:,i)=imnoise(oriData3_noise(:,:,i),'salt & pepper',ratio(i));
end
%% LRTV denoising
%
tau = 0.015;
lambda =20/sqrt(M*N);
rank = 10;
[ output_image out_value] = LRTV(oriData3_noise, tau,lambda, rank);
%% ALM + SSAHTV
% lambda =1/sqrt(M*N);
% tol = 1e-8;
% maxIter = 500;
% [output_image] = inexact_alm_rpca(oriData3_noise, lambda, tol, maxIter); 
% 
% lambda = 2;
% sa=1;
% output_image = tvdenoise(output_image,lambda,sa,'l1');
%% SSAHTV  need tvrestore.m tvrestoresa.m shrink2.m uupdategs.m getWeights.m
% lambda = 6;
% sa=1;
% output_image = tvdenoise(oriData3_noise,lambda,sa,'l1');
%% LRMR
% r = 7;
% slide =20;
% % s= 3000;
% s = 0.1;
% stepsize = 8;
% [ output_image ] = LRMR_HSI_denoise( oriData3_noise,r,slide,s,stepsize );
%%
PSNRvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);

    I=255*output_image(:,:,i);

      PSNRvector(1,i)=PSNR(J,I,M,N);
end
% dlmwrite('PSNRvector.txt',PSNRvector,'delimiter','\t','newline','pc');
% 
PSNRvec = mean(PSNRvector);
SSIMvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);
%     Jnoise=oriData3_noise(:,:,i);
    I=255*output_image(:,:,i); 
%      [ SSIMnoisevector(1,i),ssim_map] = ssim(J,Jnoise);
      [ SSIMvector(1,i),ssim_map] = ssim(J,I);
end
mean(SSIMvector);
% dlmwrite('SSIMvector.txt',SSIMvector,'delimiter','\t','newline','pc');
