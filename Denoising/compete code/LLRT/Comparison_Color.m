% clc;
clear all;
close all;
addpath(genpath('Utilize\'));
addpath(genpath('Data\'));

noise = [10,20,30,40];
randn('seed', 0);

for j =1:length(noise)
nSig =noise(j);

Original_image_dir    =    'Data\BSD';
Output_dir    =   ['Results\Color_Denoising_results\nsig_',num2str(nSig)];
mkdir(Output_dir); 
pre           =   'ALL_Denoising_';
fpath         =   fullfile(Original_image_dir, '*.jpg');
im_dir        =   dir(fpath);
im_num        =   length(im_dir);
sum_psnr      =   0;
sum_ssim      =   0;
time0         =   clock;
fn_txt        =   strcat( pre, 'PSNR_SSIM.txt' ); 
fd_txt        =   fopen( fullfile(Output_dir, fn_txt), 'wt');
PSNR_Value    =   zeros(3,im_num);
SSIM_Value    =   zeros(3,im_num);
for i = 1:im_num
    O_Img      =   double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
    N_Img      =   O_Img + nSig * randn(size(O_Img)); 
 
%% LSCD restore results
    param      =  SetParameter(nSig);
    LSCD_Img   =  imclsmooth( N_Img/256, param );
    PSNR_LSCD  =  csnr( O_Img, LSCD_Img*256, 0, 0 );
    SSIM_LSCD  =  cal_ssim( O_Img, LSCD_Img*256, 0, 0 ); 
    fprintf( sprintf('%s: PSNR_LSCD = %3.2f  SSIM_LSCD = %f \n', im_dir(i).name, PSNR_LSCD, SSIM_LSCD) );
    fprintf(fd_txt, '%s : PSNR_LSCD = %2.2f  SSIM_LSCD = %2.4f \n', im_dir(i).name, PSNR_LSCD, SSIM_LSCD);   

    
%% CBM3D restore results
    [~, CBM3D_Img] = CBM3D(1, N_Img, nSig);
    PSNR_CBM3D =  csnr( O_Img, CBM3D_Img*256, 0, 0 );
    SSIM_CBM3D  =  cal_ssim( O_Img, CBM3D_Img*256, 0, 0 ); 
    fprintf( sprintf('%s: PSNR_CBM3D = %3.2f  SSIM_CBM3D = %f \n', im_dir(i).name, PSNR_CBM3D, SSIM_CBM3D) );
    fprintf(fd_txt, '%s : PSNR_CBM3D = %2.2f  SSIM_CBM3D = %2.4f \n', im_dir(i).name, PSNR_CBM3D, SSIM_CBM3D);    

%% LLRT restore results
    Par   = ParSet_Color(nSig);  
    [E_Img]= LLRT_DeNoising( N_Img, O_Img, Par);       
    PSNR_LLRT  =  csnr( O_Img, E_Img, 0, 0 );
    SSIM_LLRT  =  cal_ssim( O_Img, E_Img, 0, 0 ); 
    fprintf( sprintf('%s: PSNR_LLRT = %3.2f  SSIM_LLRT = %f \n', im_dir(i).name, PSNR_LLRT, SSIM_LLRT) );
    fprintf(fd_txt, '%s : PSNR_LLRT = %2.2f  SSIM_LLRT = %2.4f \n', im_dir(i).name, PSNR_LLRT, SSIM_LLRT);   

%% Quantitative Assessment
    PSNR_Value(:,i) = [PSNR_LSCD;PSNR_CBM3D;PSNR_LLRT];
    SSIM_Value(:,i) = [SSIM_LSCD;SSIM_CBM3D;SSIM_LLRT];
end
fclose(fd_txt);
end
fprintf(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));