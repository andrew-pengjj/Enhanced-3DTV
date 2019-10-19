clear all;
clc;	
%% simulated experiment 3
% ----------------------------load image-----------------------------------
load simu_indian
Ohsi       = simu_indian;
ratio      = 0.15*ones(1,224);           
noiselevel = 0.075*ones(1,224); 
% ------------------------ Simulation experiment --------------------------
Nhsi      = Ohsi;
[M,N,p]   = size(Ohsi);
%% Gaussian noise
for i = 1:p
     Nhsi(:,:,i)=Ohsi(:,:,i)  + noiselevel(i)*randn(M,N);
end
%% S&P noise
for i = 1:p
     Nhsi(:,:,i)=imnoise(Nhsi(:,:,i),'salt & pepper',ratio(i));
end
%% TV sparsity denoising
rank   = [13,13,13];
tau    = 0.004 *sqrt(M*N);
output_image = EnhancedTV(Nhsi,tau,rank);
[mpsnr,mssim,ergas]=msqia(Ohsi,output_image);