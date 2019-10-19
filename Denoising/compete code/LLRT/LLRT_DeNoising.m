function [E_Img]   =  LLRT_DeNoising( N_Img, O_Img, Par )

E_Img            = N_Img;                                                         % Estimated Image
[Height, Width, Band]  = size(E_Img);   
TotalPatNum      = (Height-Par.patsize+1)*(Width-Par.patsize+1);                  % Total Patch Number in the image
Average          = mean(N_Img,3);                                                 % Calculate the average band for fast spatial non-local searching
[Neighbor_arr, Num_arr, Self_arr] =	NeighborIndex(Average, Par);                  % PreCompute all the patch index in the searching window 

for iter = 1 : Par.Iter        
    Average             =   mean(E_Img,3);
    [CurPat, Mat, Sigma_arr]	=	Cub2Patch( E_Img, N_Img, Average, Par );      % Cubic to patch and estimate local noise variance
    
    if (mod(iter-1,2)==0)
        Par.patnum = Par.patnum - 10;                                             % Lower Noise level, less NL patches
        NL_mat  =  Block_matching(Mat, Par, Neighbor_arr, Num_arr, Self_arr);     % Caculate Non-local similar patches for each
        if(iter==1)
            Sigma_arr = Par.nSig * ones(size(Sigma_arr));                         % First Iteration use the input noise parameter
        end
    end
    
    [Spa_EPat, Spa_W]    =  NLPatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par);      % Estimate the spatial patches
    [Spa_Img, Spa_Wei]   =  Patch2Cub( Spa_EPat, Spa_W, Par.patsize, Height, Width, Band );       % Patch to Cubic
%          E_Img = Spa_Img./Spa_Wei;
    [E_Img]         =  ImUpdata(N_Img, Spa_Img, Spa_Wei, sum(Sigma_arr)/TotalPatNum, Par );       % Spectral sparsity constraint. Comment out this line means no spectral sparsity constraint. 
    % For HSI with large spectral bands, the spectral constraint could 
    % improve the denoising result; while for color image denoising, the 
    % spectral gradient sparsity assumption may be invalid.
    % The specific rgb correlationship constraint is expected (out of the
    % scope of this work). So when you test the color image denoising,
    % please comment out this line.
    
    PSNR  = csnr( O_Img, E_Img, 0, 0 );
    SSIM  = cal_ssim( O_Img, E_Img, 0, 0 );
    fprintf( 'Iter = %2.3f, PSNR = %2.2f, SSIM = %2.3f, NoiseLevel = %2.3f \n', iter, PSNR, SSIM, sum(Sigma_arr)/TotalPatNum);
end
return;


