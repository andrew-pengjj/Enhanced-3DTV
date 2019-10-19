load simu_indian
OriData3 = simu_indian;
ratio = 0.15*ones(1,224);            % for case 1 ,   0.15
noiselevel = 0.075*ones(1,224);      % for case 1,      0.075
oriData3_noise = OriData3;
[M,N,p] = size(OriData3);
num=6;
MPSNR=zeros(num,p);
MSSIM=MPSNR;
for iter=1:num
    % Gaussian noise
    for i =1:p
        oriData3_noise(:,:,i)=OriData3(:,:,i)  + noiselevel(i)*randn(M,N);
    end
    %S&P noise
    for i =1:p
        oriData3_noise(:,:,i)=imnoise(oriData3_noise(:,:,i),'salt & pepper',ratio(i));
    end
    %% LRMR
    fprintf('=============================   %d ´ÎÑ­»·  =============================\n',iter);
    r = 7;
    blocksize =20;

    s = 0.15;   
    stepsize = 8;
    [ output_image ] = LRMR_HSI_denoise( oriData3_noise,r,blocksize,s,stepsize );

    %%  LRMR index
    PSNRvector=zeros(1,p);
    for i=1:1:p
        J=255*OriData3(:,:,i);
        I=255*output_image(:,:,i);
        PSNRvector(i)=PSNR(J,I,M,N);
    end
    MPSNR(iter,:)=PSNRvector;
    
    SSIMvector=zeros(1,p);
    for i=1:1:p
        J=255*OriData3(:,:,i);
        I=255*output_image(:,:,i); 
        [ SSIMvector(i),ssim_map] = ssim(J,I);
    end
    MSSIM(iter,:)=SSIMvector;
end
