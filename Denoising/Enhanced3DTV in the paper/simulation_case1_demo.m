clear all;
clc;	
%% simulated experiment 1
% ----------------------------load image-----------------------------------
addpath('TV_operator')
addpath('quality assess')
load simu_indian
Ohsi       = simu_indian;
noiselevel = 0.1*ones(1,224); 
% ------------------------ Simulation experiment --------------------------
Nhsi      = Ohsi;
[M,N,p]   = size(Ohsi);
%% Gaussian noise
for i = 1:p
     Nhsi(:,:,i)=Ohsi(:,:,i)+noiselevel(i)*randn(M,N);
end
%% TV sparsity denoising
rank   = [13,13,13];
tau    = 0.004 *sqrt(M*N);
output_image = EnhancedTV(Nhsi,tau,rank);
[mpsnr,mssim,ergas]=msqia(Ohsi,output_image);

