% Main function
% Usage :
%        main
% INPUT :
%        Raw Lena image
% OUTPUT :
%        PSNR value of the denoised image
%
% Load clean image
fid = fopen('barbara','r');
s  = fread(fid,[512 512],'unsigned char');
fclose(fid)
%s = imread('zoe.bmp','bmp');
%s = double(s);
N = 512;

% Noise variance
sigma_n = 25;
n = sigma_n*randn(N);

% Add noise 
x = s + n;

% Run local adaptive image denoising algorithm
y = denoising_dwt(x);

% Calculate the error 
err = s - y;

% Calculate the PSNR value
PSNR = 20*log10(256/std(err(:)))