% iwfft_hs.m
%
% Fourier/wavelet tensor transform for hyperspectral images
%
% Usage: x = iwfft_hs(w, h, S)
%
% w - input coefficient matrix
%
% h - Wavelet filter
%
% S - number of spectral bands
%
% x - output data cube
%
% Written by: Marco F. Duarte, Rice University
% Created: February 25 2008
%
% Uses Rice Wavelet Toolbox, http://dsp.rice.edu/rwt

function x = iwfft_hs(w, h,S)
M = sqrt(length(w)/S);
w = reshape(w,[M M S]);
wvltx = zeros(size(w));
for i=1:S,
    % Perform 2D IWT on snapshot i (space)
    tmp = reshape(midwt(squeeze(w(:,:,i)),h,log2(M)),[M M 1]);
    wvltx(:,:,i) = tmp(:,:,1);
end

if S == 1,
    x = wvltx;
else
    x = zeros(size(wvltx));
    for i=1:M,
        for j=1:M,
            % Perform 1D FFT (spectra) on pixel (i,j)
            x(i,j,:) = idct(wvltx(i,j,:));
        end
    end
end
