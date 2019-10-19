clear all;
clc;	
%% simulated experiment 6
% ----------------------------load image-----------------------------------
load simu_indian
Ohsi = simu_indian;
load Simu_ratio                       
load Simu_noiselevel                 

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
%% dead line 
for i=91:130
    indp=randperm(10,1)+2;
    ind=randperm(N-1,indp);
    an=funrand(2,length(ind));
    % searching the location of an which value is 1,2,3
    loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
    Nhsi(:,ind(loc1),i)=0; 
    Nhsi(:,ind(loc2):ind(loc2)+1,i)=0;
    Nhsi(:,ind(loc3)-1:ind(loc3)+1,i)=0;
end  
%% stripe line
for band=161:190
    num = 19+randperm(21,1);
    loc = ceil(N*rand(1,num));
    t = rand(1,length(loc))*0.5-0.25;
    Nhsi(:,loc,band) = bsxfun(@minus,Nhsi(:,loc,band),t);
end
%% TV sparsity denoising
rank   = [13,13,13];
tau    = 0.004 *sqrt(M*N);
output_image = EnhancedTV(Nhsi,tau,rank);
[mpsnr,mssim,ergas]=msqia(Ohsi,output_image);