clear all;
clc;	

addpath ./prox_operators
addpath ./mylib

%% simulated experiment 1
% ----------------------------load image-----------------------------------
load phantomo;
Omsi       = phantomo;%(61:200,31:210,:);
ratio      = 0.5;
noiselevel = ratio*ones(1,224); 

% ------------------------ Simulation experiment --------------------------
Nmsi      = Omsi;
[M,N,p]   = size(Omsi);
%% Gaussian noise
for i = 1:p
     Nmsi(:,:,i)=Omsi(:,:,i)  + noiselevel(i)*randn(M,N);
end
for i = 1:p
    tmp = Nmsi(:,:,i);
    tmp(Ind)=0;
    Nmsi(:,:,i)=tmp;
end
%% LRTD denoising
tau    = 1/4;
lambda = (1/ratio)^2;
Rank   = [round(M*0.6),round(N*0.6),3];r=2;
it     = 1;
weight = [1,1,1];
clean_image                    = LRTDTV_w(Nmsi, tau,lambda,Rank,weight);
[mpsnr(it),mssim(it),ergas(it)]= msqia(Omsi, clean_image);
for i = 1:p
    tmp = clean_image(:,:,i);
    tmp(Ind)=0;
    output_image(:,:,i)=tmp;
end
it =2;
[mpsnr(it),mssim(it),ergas(it)]= msqia(Omsi, output_image);

