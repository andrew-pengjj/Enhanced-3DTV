% Main function
% Usage :
%        main_dtdwt
% INPUT :
%        Raw Lena image
% OUTPUT :
%        PSNR value of the denoised image
%
% Load clean image
fid = fopen('boat','r');
s  = fread(fid,[512 512],'unsigned char');
fclose(fid);
N = 512;

% Noise variance
sigma_n = 30;
n = sigma_n*randn(N);

% Add noise 
x = s + n;

% Run local adaptive image denoising algorithm using dual-tree DWT. 
y = denoising_dtdwt(x);

% Calculate the error 
err = s - y;

% Calculate the PSNR value
PSNR = 20*log10(256/std(err(:)))