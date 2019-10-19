clc;
clear all;
close all;
addpath(genpath('Utilize\'));
addpath(genpath('Data\'));
noise = [10,30,50,100];
randn('seed', 0);

for j =1:length(noise)
nSig =noise(j);

Original_image_dir    =    'Data\CAVE';
Output_dir    =   ['Results\Denoising_results\nsig_',num2str(nSig)];
mkdir(Output_dir); 
pre           =   'All_Denoising_';
fpath         =   fullfile(Original_image_dir, '*.mat');
im_dir        =   dir(fpath);
im_num        =   length(im_dir);
sum_psnr      =   0;
sum_ssim      =   0;
time0         =   clock;
fn_txt        =   strcat( pre, 'PSNR_SSIM.txt' ); 
fd_txt        =   fopen( fullfile(Output_dir, fn_txt), 'wt');
PSNR_Value    =   zeros(12,im_num);
SSIM_Value    =   zeros(12,im_num);
SAM_Value     =   zeros(12,im_num);
ERGAS_Value   =   zeros(12,im_num);
for i = 1:im_num
    load(fullfile(Original_image_dir, im_dir(i).name))
    O_Img      =   O_ImgUint8(200:300,200:300,:);
    N_Img      =   O_Img + nSig * randn(size(O_Img)); 

%% BM4D restore results
    [y_est, ~] = bm4d(N_Img, 'Gauss', (~0)*nSig, 'mp', 1, 1);
    PSNR_BM4D  =  csnr( O_Img, y_est, 0, 0 );
    SSIM_BM4D  =  cal_ssim( O_Img, y_est, 0, 0 ); 
    SAM_BM4D   = SpectAngMapper(O_Img, y_est);
    ERGAS_BM4D = ErrRelGlobAdimSyn(O_Img, y_est);
    fprintf( sprintf('%s: PSNR_BM4D = %3.2f  SSIM_BM4D = %f SAM_BM4D =%f  ERGAS_BM4D =%f \n', im_dir(i).name, PSNR_BM4D, SSIM_BM4D, SAM_BM4D, ERGAS_BM4D) );
    fprintf(fd_txt, '%s : PSNR_BM4D = %2.2f  SSIM_BM4D = %2.4f SAM_BM4D =%f  ERGAS_BM4D =%f \n', im_dir(i).name, PSNR_BM4D, SSIM_BM4D, SAM_BM4D, ERGAS_BM4D); 
    
%% PARAFAC restore results
    [PARAFAC_Img, k] = PARAFAC(N_Img);
    PSNR_PARAFAC  =  csnr( O_Img, PARAFAC_Img, 0, 0 );
    SSIM_PARAFAC  =  cal_ssim( O_Img, PARAFAC_Img, 0, 0 );
    SAM_PARAFAC = SpectAngMapper(O_Img, PARAFAC_Img);
    ERGAS_PARAFAC = ErrRelGlobAdimSyn(O_Img, PARAFAC_Img);
    fprintf( sprintf('%s: PSNR_PARAFAC = %3.2f  SSIM_PARAFAC = %f SAM_PARAFAC = %f  ERGAS_PARAFAC = %f \n', im_dir(i).name, PSNR_PARAFAC, SSIM_PARAFAC, SAM_PARAFAC, ERGAS_PARAFAC) );
    fprintf(fd_txt, '%s : PSNR_PARAFAC = %2.2f  SSIM_PARAFAC = %2.4f SAM_PARAFAC = %f  ERGAS_PARAFAC = %f \n', im_dir(i).name, PSNR_PARAFAC, SSIM_PARAFAC, SAM_PARAFAC, ERGAS_PARAFAC);     

%% LRTA restore results
    LRTA_Img   =  double(LRTA(tensor(N_Img)));
    PSNR_LRTA  =  csnr( O_Img, LRTA_Img, 0, 0 );
    SSIM_LRTA  =  cal_ssim( O_Img, LRTA_Img, 0, 0 ); 
    SAM_LRTA   = SpectAngMapper(O_Img, LRTA_Img);
    ERGAS_LRTA = ErrRelGlobAdimSyn(O_Img, LRTA_Img);
    fprintf( sprintf('%s: PSNR_LRTA = %3.2f  SSIM_LRTA = %f SAM_LRTA = %f  ERGAS_LRTA = %f \n', im_dir(i).name, PSNR_LRTA, SSIM_LRTA, SAM_LRTA , ERGAS_LRTA ) );
    fprintf(fd_txt, '%s : PSNR_LRTA = %2.2f  SSIM_LRTA = %2.4f SAM_LRTA = %f  ERGAS_LRTA = %f \n', im_dir(i).name, PSNR_LRTA, SSIM_LRTA, SAM_LRTA, ERGAS_LRTA);     

%% SDS restore results
    [SDS_Img,~] = denoise(shiftdim(N_Img,2));
    PSNR_SDS  =  csnr( O_Img, shiftdim(SDS_Img,1), 0, 0 );
    SSIM_SDS  =  cal_ssim( O_Img, shiftdim(SDS_Img,1), 0, 0 ); 
    SAM_SDS   = SpectAngMapper(O_Img, shiftdim(SDS_Img,1));
    ERGAS_SDS = ErrRelGlobAdimSyn(O_Img, shiftdim(SDS_Img,1));
    fprintf( sprintf('%s: PSNR_SDS = %3.2f  SSIM_SDS = %f SAM_SDS = %f  ERGAS_SDS = %f \n', im_dir(i).name, PSNR_SDS, SSIM_SDS, SAM_SDS , ERGAS_SDS ) );
    fprintf(fd_txt, '%s : PSNR_SDS = %2.2f  SSIM_SDS = %2.4f SAM_SDS = %f  ERGAS_SDS = %f \n', im_dir(i).name, PSNR_SDS, SSIM_SDS, SAM_SDS, ERGAS_SDS);           

%% SSTV restore results
    SSTV_Img  =   funSSTV(N_Img,30,0.1,3,0.2);
    PSNR_SSTV  =  csnr( O_Img, SSTV_Img, 0, 0 );
    SSIM_SSTV  =  cal_ssim( O_Img, SSTV_Img, 0, 0 ); 
    SAM_SSTV   = SpectAngMapper(O_Img,SSTV_Img);
    ERGAS_SSTV = ErrRelGlobAdimSyn(O_Img, SSTV_Img);
    fprintf( sprintf('%s: PSNR_SSTV = %3.2f  SSIM_SSTV = %f SAM_SSTV = %f  ERGAS_SSTV = %f \n', im_dir(i).name, PSNR_SSTV, SSIM_SSTV, SAM_SSTV , ERGAS_SSTV ) );
    fprintf(fd_txt, '%s : PSNR_SSTV = %2.2f  SSIM_SSTV = %2.4f SAM_SSTV = %f  ERGAS_SSTV = %f \n', im_dir(i).name, PSNR_SSTV, SSIM_SSTV, SAM_SSTV, ERGAS_SSTV);               
    
%% ANLM restore results
    min_val = min(N_Img(:));
    TN_Img  = N_Img - min_val;
    fimau = MABONLM3D(TN_Img, 3, 1, 0);
    fimao = MABONLM3D(TN_Img, 3, 2, 0);
    ANLM_Img = mixingsubband(fimau, fimao) + min_val;
    PSNR_ANLM  =  csnr( O_Img, ANLM_Img, 0, 0 );
    SSIM_ANLM =  cal_ssim( O_Img, ANLM_Img, 0, 0 ); 
    SAM_ANLM   = SpectAngMapper(O_Img,ANLM_Img);
    ERGAS_ANLM = ErrRelGlobAdimSyn(O_Img, ANLM_Img);
    fprintf( sprintf('%s: PSNR_ANLM = %3.2f  SSIM_ANLM = %f SAM_ANLM = %f  ERGAS_ANLM = %f \n', im_dir(i).name, PSNR_ANLM, SSIM_ANLM, SAM_ANLM , ERGAS_ANLM ) );
    fprintf(fd_txt, '%s : PSNR_ANLM = %2.2f  SSIM_ANLM = %2.4f SAM_ANLM = %f  ERGAS_ANLM = %f \n', im_dir(i).name, PSNR_ANLM, SSIM_ANLM, SAM_ANLM, ERGAS_ANLM);    

%% NMF restore results
    [NMF_Img, A, S] = nmf_denoising(N_Img, [7,7], 2, 98, nSig*3);
    PSNR_NMF  =  csnr( O_Img, NMF_Img, 0, 0 );
    SSIM_NMF  =  cal_ssim( O_Img, NMF_Img, 0, 0 ); 
    SAM_NMF   = SpectAngMapper(O_Img,NMF_Img);
    ERGAS_NMF = ErrRelGlobAdimSyn(O_Img, NMF_Img);
    fprintf( sprintf('%s: PSNR_NMF = %3.2f  SSIM_NMF = %f SAM_NMF = %f  ERGAS_NMF = %f \n', im_dir(i).name, PSNR_NMF, SSIM_NMF, SAM_NMF , ERGAS_NMF ) );
    fprintf(fd_txt, '%s : PSNR_NMF = %2.2f  SSIM_NMF = %2.4f SAM_NMF = %f  ERGAS_NMF = %f \n', im_dir(i).name, PSNR_NMF, SSIM_NMF, SAM_NMF, ERGAS_NMF);    
    
%% TDL restore results
    vstbmtf_params.peak_value = 1;
    vstbmtf_params.nsigma = nSig/256;
    TDL_Img = TensorDL(N_Img/256, vstbmtf_params);
    PSNR_TDL  =  csnr( O_Img, 256*TDL_Img, 0, 0 );
    SSIM_TDL  =  cal_ssim( O_Img, 256*TDL_Img, 0, 0 ); 
    SAM_TDL   = SpectAngMapper(O_Img,256*TDL_Img);
    ERGAS_TDL = ErrRelGlobAdimSyn(O_Img, 256*TDL_Img);
    fprintf( sprintf('%s: PSNR_TDL = %3.2f  SSIM_TDL = %f SAM_TDL = %f  ERGAS_TDL = %f \n', im_dir(i).name, PSNR_TDL, SSIM_TDL, SAM_TDL , ERGAS_TDL ) );
    fprintf(fd_txt, '%s : PSNR_TDL = %2.2f  SSIM_TDL = %2.4f SAM_TDL = %f  ERGAS_TDL = %f \n', im_dir(i).name, PSNR_TDL, SSIM_TDL, SAM_TDL, ERGAS_TDL);   
    
%% BM3D restore results
    BM3D_Img   =  N_Img;
    for k = 1:size(N_Img,3)
        [~, BM3D_Img(:,:,k)] = BM3D(1,N_Img(:,:,k), nSig);
    end
    PSNR_BM3D  =  csnr( O_Img, 256*BM3D_Img, 0, 0 );
    SSIM_BM3D  =  cal_ssim( O_Img, 256*BM3D_Img, 0, 0 ); 
    SAM_BM3D   = SpectAngMapper(O_Img,256*BM3D_Img);
    ERGAS_BM3D = ErrRelGlobAdimSyn(O_Img, 256*BM3D_Img);
    fprintf( sprintf('%s: PSNR_BM3D = %3.2f  SSIM_BM3D = %f SAM_BM3D = %f  ERGAS_BM3D = %f \n', im_dir(i).name, PSNR_BM3D, SSIM_BM3D, SAM_BM3D , ERGAS_BM3D ) );
    fprintf(fd_txt, '%s : PSNR_BM3D = %2.2f  SSIM_BM3D = %2.4f SAM_BM3D = %f  ERGAS_BM3D = %f \n', im_dir(i).name, PSNR_BM3D, SSIM_BM3D, SAM_BM3D, ERGAS_BM3D);   

%% LRMR restore results
    Par        =  ParamSet(nSig);
    LRMR_Img   =  LRMR_DeNoising( N_Img, Par );
    PSNR_LRMR  =  csnr( O_Img, LRMR_Img, 0, 0 );
    SSIM_LRMR  =  cal_ssim( O_Img, LRMR_Img, 0, 0 ); 
    SAM_LRMR   = SpectAngMapper(O_Img,LRMR_Img);
    ERGAS_LRMR = ErrRelGlobAdimSyn(O_Img, LRMR_Img);
    fprintf( sprintf('%s: PSNR_LRMR = %3.2f  SSIM_LRMR = %f SAM_LRMR = %f  ERGAS_LRMR = %f \n', im_dir(i).name, PSNR_LRMR, SSIM_LRMR, SAM_LRMR , ERGAS_LRMR ) );
    fprintf(fd_txt, '%s : PSNR_LRMR = %2.2f  SSIM_LRMR = %2.4f SAM_LRMR = %f  ERGAS_LRMR = %f \n', im_dir(i).name, PSNR_LRMR, SSIM_LRMR, SAM_LRMR, ERGAS_LRMR);    


%% ISTreg restore results
    ISTReg_Img = ITS_DeNoising(N_Img/256,nSig/256,O_Img/256);
    PSNR_ISTReg  =  csnr( O_Img, 256*ISTReg_Img, 0, 0 );
    SSIM_ISTReg  =  cal_ssim( O_Img, 256*ISTReg_Img, 0, 0 ); 
    SAM_ISTReg   = SpectAngMapper(O_Img,256*ISTReg_Img);
    ERGAS_ISTReg = ErrRelGlobAdimSyn(O_Img, 256*ISTReg_Img);
    fprintf( sprintf('%s: PSNR_ISTReg = %3.2f  SSIM_ISTReg = %f SAM_ISTReg = %f  ERGAS_ISTReg = %f \n', im_dir(i).name, PSNR_ISTReg, SSIM_ISTReg, SAM_ISTReg , ERGAS_ISTReg ) );
    fprintf(fd_txt, '%s : PSNR_ISTReg = %2.2f  SSIM_ISTReg = %2.4f SAM_ISTReg = %f  ERGAS_ISTReg = %f \n', im_dir(i).name, PSNR_ISTReg, SSIM_ISTReg, SAM_ISTReg, ERGAS_ISTReg);
    
%% LLRT restore results
    Par   = ParSet(nSig);  
    [E_Img]= LLRT_DeNoising( N_Img, O_Img, Par);       
    PSNR_LLRT  =  csnr( O_Img, E_Img, 0, 0 );
    SSIM_LLRT  =  cal_ssim( O_Img, E_Img, 0, 0 ); 
    SAM_LLRT   = SpectAngMapper(O_Img,E_Img);
    ERGAS_LLRT = ErrRelGlobAdimSyn(O_Img, E_Img);
    fprintf( sprintf('%s: PSNR_LLRT = %3.2f  SSIM_LLRT = %f SAM_LLRT = %f  ERGAS_LLRT = %f \n', im_dir(i).name, PSNR_LLRT, SSIM_LLRT, SAM_LLRT , ERGAS_LLRT ) );
    fprintf(fd_txt, '%s : PSNR_LLRT = %2.2f  SSIM_LLRT = %2.4f SAM_LLRT = %f  ERGAS_LLRT = %f \n', im_dir(i).name, PSNR_LLRT, SSIM_LLRT, SAM_LLRT, ERGAS_LLRT);      
%% Quantitative Assessment
    PSNR_Value(:,i) = [PSNR_BM4D;PSNR_PARAFAC;PSNR_LRTA;PSNR_SDS;PSNR_SSTV;PSNR_ANLM;PSNR_NMF;PSNR_TDL;PSNR_BM3D;PSNR_LRMR;PSNR_ISTReg;PSNR_LLRT];
    SSIM_Value(:,i) = [SSIM_BM4D;SSIM_PARAFAC;SSIM_LRTA;SSIM_SDS;SSIM_SSTV;SSIM_ANLM;SSIM_NMF;SSIM_TDL;SSIM_BM3D;SSIM_LRMR;SSIM_ISTReg;SSIM_LLRT];
    SAM_Value(:,i) = [SAM_BM4D;SAM_PARAFAC;SAM_LRTA;SAM_SDS;SAM_SSTV;SAM_ANLM;SAM_NMF;SAM_TDL;SAM_BM3D;SAM_LRMR;SAM_ISTReg;SAM_LLRT];
    ERGAS_Value(:,i) = [ERGAS_BM4D;ERGAS_PARAFAC;ERGAS_LRTA;ERGAS_SDS;ERGAS_SSTV;ERGAS_ANLM;ERGAS_NMF;ERGAS_TDL;ERGAS_BM3D;ERGAS_LRMR;ERGAS_ISTReg;ERGAS_LLRT];
end
fclose(fd_txt);
end
fprintf(sprintf('Total elapsed time = %f min\n', (etime(clock,time0)/60) ));