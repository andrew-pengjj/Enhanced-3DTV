clear all;
clc;	

addpath ./prox_operators
addpath ./mylib

%% simulated experiment 3
% ----------------------------load image-----------------------------------
load DC_pure
Omsi       = simu_indian;           
noiselevel = 0.1*ones(1,224); 

% ------------------------ Simulation experiment --------------------------
Nmsi      = Omsi;
[M,N,p]   = size(Omsi);
%% Gaussian noise
for i = 1:p
     Nmsi(:,:,i)=Omsi(:,:,i)  + noiselevel(i)*randn(M,N);
end
% %% S&P noise
% for i = 1:p
%      Nmsi(:,:,i)=imnoise(Nmsi(:,:,i),'salt & pepper',ratio(i));
% end
%% dead line 
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
tau    = 1;
Lambda = 50;
Rank   = [160,160,6];
for it=1:length(Lambda)
    fprintf('===============%d/%d==============\n',it,length(Lambda));
clean_image                    = LRTDTV(Nmsi, tau,Lambda(it),Rank);
[mpsnr(it),mssim(it),ergas(it)]= msqia(Omsi, clean_image);
end
% it=2;
% weight=[1,1,0.5];
% clean_image                    = LRTDTV_w(Nmsi, tau,lambda,Rank,weight);
% [mpsnr(it),mssim(it),ergas(it)]= msqia(Omsi, clean_image);

