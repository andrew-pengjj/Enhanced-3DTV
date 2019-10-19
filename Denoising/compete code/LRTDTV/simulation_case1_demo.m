clear all;
clc;	

addpath ./prox_operators
addpath ./mylib

%% simulated experiment 1
% ----------------------------load image-----------------------------------
load simu_indian
Omsi       = simu_indian;
noiselevel = 0.1*ones(1,224); 

% ------------------------ Simulation experiment --------------------------
Nmsi      = Omsi;
[M,N,p]   = size(Omsi);
%% Gaussian noise
for i = 1:p
     Nmsi(:,:,i)=Omsi(:,:,i)  + noiselevel(i)*randn(M,N);
end
% for i=91:130
%     indp=randperm(10,1)+2;
%     ind=randperm(N-1,indp);
%     an=funrand(2,length(ind));
%     % searching the location of an which value is 1,2,3
%     loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
%     Nmsi(:,ind(loc1),i)=0; 
%     Nmsi(:,ind(loc2):ind(loc2)+1,i)=0;
%     Nmsi(:,ind(loc3)-1:ind(loc3)+1,i)=0;
% end  
%% LRTD denoising
% lambda =100;
% Rank   = [120,120,10];
% [clean_image,S,out_value,time] = LRTD(Nmsi,lambda,Rank);
% it     = 1;
% [mpsnr(it),mssim(it),ergas(it)]= msqia(Omsi, clean_image);
tau    = 1;
lambda = 100000;beta=100;
Rank   = [120,120,10];
it     = 1;
output_image                   = LRTDTV_G(Nmsi, tau,lambda,beta,Rank);
[mpsnr(it),mssim(it),ergas(it)]= msqia(Omsi, output_image);
