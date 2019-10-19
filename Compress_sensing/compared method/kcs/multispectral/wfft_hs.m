% wfft_hs.m
%
% Fourier/wavelet tensor transform for hyperspectral imaging
%
% Usage: w = wvlt_hs(x, h)
%
% x - input signal of size YxMxN, where MxN is the image size and Y is the
% number of spectral bands
%
% h - Wavelet filter name
%
% w - 3D coefficient matrix for x
%
% Written by: Marco F. Duarte, Rice University
% Created: February 25 2008

function w = wfft_hs(x, h)

N = size(x,1);
Y = size(x,3);
if Y > 1,
    wvlts = zeros(size(x));
    for i=1:N,
        for j=1:N,
            % Perform 1D FFT (spectra) on pixel (i,j)
            wvlts(i,j,:) = dct(x(i,j,:));
        end
    end
else
    wvlts = x;
end
w = zeros(size(wvlts));
for i=1:Y,
    % Perform 2D DWT on snapshot i (space)
    tmp = reshape(mdwt(squeeze(wvlts(:,:,i)),h,log2(N)),[ N N 1]);
    w(:,:,i) = tmp(:,:,1);
end
w = w(:);
