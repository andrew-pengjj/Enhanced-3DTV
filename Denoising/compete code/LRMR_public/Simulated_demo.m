% clear all;
% close all;
% clc;	

load 'simu_indian'
OriData3=simu_indian;
[M N p] = size(OriData3);
oriData3_noise = OriData3;
%% simulated experiment
ratio = 0.15*ones(1,224);            % for case 1
noiselevel = 0.075*ones(1,224);      % for case 1
%noise simulated
for i =1:p
     oriData3_noise(:,:,i)=OriData3(:,:,i)  + noiselevel(i)*randn(M,N);
%        mixed(vv,:)=addnoise(mixed(vv,:),'mg',ratio);
end
% oriData3_noise = OriData3 + 0.075*randn(M,N,p);
for i =1:p
     oriData3_noise(:,:,i)=imnoise(oriData3_noise(:,:,i),'salt & pepper',ratio(i));
%        mixed(vv,:)=addnoise(mixed(vv,:),'mg',ratio);
end
 
%% LRMR
r = 7;
blocksize =20;
% s= 3000;  % for LRMR implemented with GoDec 
s = 0.15;   % for LRMR implemented with ssGoDec 
stepsize = 8;
[ output_image ] = LRMR_HSI_denoise( oriData3_noise,r,blocksize,s,stepsize );
%%
% PSNR_in=10*log10((mean((OriData3(:)).^2))/mean((oriData3_noise(:)-OriData3(:)).^2));
% PSNR_out=10*log10((mean((OriData3(:)).^2))/mean((output_image(:)-OriData3(:)).^2));

PSNRvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);
    I=255*output_image(:,:,i);
    PSNRvector(1,i)=PSNR(J,I,M,N);
end
% dlmwrite('PSNRvector.txt',PSNRvector,'delimiter','\t','newline','pc'); 
MPSNR = mean(PSNRvector);
