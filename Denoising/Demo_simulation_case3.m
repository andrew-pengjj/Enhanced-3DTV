%==========================================================================
% This script compares multi-spectral imagery (MSI) noise removal methods
% listed as follows:
%   1. NNM
%   2. WNNM
%   3. LRMR
%   4. BM4D
%   5. TDL
%   6. WSNM
%   7. LRTV
%   8. LRTDTV
%   9. Ours
%==========================================================================

%% 1: initiation
% 1.1 ========== DATA SET==========
clear all;clc;
addpath(genpath('compete code'));
addpath(genpath('quality assess'));
addpath(genpath('Enhanced3DTV in the paper'));
load simu_indian
Ohsi       = simu_indian;
[M,N,p]    = size(Ohsi);

% 1.3 ========== NOISE PATTERN SET==========
noiselevel = 0.075*ones(1,p); 
gausssigma = mean(noiselevel);
ratio      = 0.15*ones(1,p);
spsigma    = mean(ratio);
% ------------------------ Simulation experiment --------------------------
Nhsi      = Ohsi;
%% Gaussian noise
for i = 1:p
     Nhsi(:,:,i)=Ohsi(:,:,i)+noiselevel(i)*randn(M,N);
end  
%% S&P noise
for i = 1:p
     Nhsi(:,:,i)=imnoise(Nhsi(:,:,i),'salt & pepper',ratio(i));
end
% 1.4 ========== GCP ACCELARATE SET==========
%     if isempty(gcp)
%         parpool('local',4,'IdleTimeout', 6000); % If your computer's memory is less than 8G, do not use more than 4 workers.
%     end 
%% 2 START RUNNING THE METHOD
options.STATUS_NNM  = 0;
options.STATUS_WNNM = 0;
options.STATUS_LRMR = 0;
options.STATUS_BM4D = 0;
options.STATUS_TDL  = 0;
options.STATUS_WSNM = 0;
options.STATUS_LRTV = 1;
options.STATUS_LRTDTV = 0;
options.STATUS_Enhanced3DTV  = 0;
options.spsigma = spsigma;
options.gausssigma = gausssigma;
if spsigma ~=0
    options.gausssigma = 0;
end
[mpsnr,mssim,ergas,Time]=RunAllMethod(options,Ohsi,Nhsi);