    close all;
    clear all;
    addpath(genpath('Utilize\'));
    addpath(genpath('Data\'));    
    nSig = 20;    

    load('chart_and_stuffed_toy_msUint8.mat')
    O_Img = O_ImgUint8;
%     O_Img = O_ImgUint8(200:300,200:300,:);
    randn('seed', 0);
    N_Img = O_Img + nSig * randn(size(O_Img));                                   %Generate noisy image   
    PSNR  =  csnr( N_Img, O_Img, 0, 0 );
    SSIM  = cal_ssim( N_Img, O_Img, 0, 0 ); 
    fprintf( 'Noisy Image: nSig = %2.3f, PSNR = %2.2f, SSIM = %2.3f,\n', nSig, PSNR,SSIM );
    
    Par   = ParSet(nSig); 
    tic
    [E_Img]= LLRT_DeNoising( N_Img, O_Img, Par);  %LLRT denoisng function
    toc
    PSNR  = csnr( O_Img, E_Img, 0, 0 );
    fprintf( 'Estimated Image: nSig = %2.3f, PSNR = %2.2f \n', nSig, PSNR );
