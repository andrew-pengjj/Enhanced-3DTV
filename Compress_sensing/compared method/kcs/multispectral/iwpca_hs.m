% iwpca_hs.m
%
% Fourier/wavelet tensor transform for hyperspectral images
%
% Usage: x = iwpca_hs(w, h, S)
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

function x = iwpca_hs(w, h, v)

S = size(v,1);
M = sqrt(length(w)/S);
w = reshape(w,[S M M]);
wvltx = zeros(size(w));
for i=1:S,
    % Perform 2D IWT on snapshot i (space)
    tmp = reshape(midwt(squeeze(w(i,:,:)),h,log2(M)),[1 M M]);
    wvltx(i,:,:) = tmp(1,:,:);
end

if S == 1,
    x = wvltx;
else
    x = zeros(size(wvltx));
    for i=1:M,
        for j=1:M,
            % Perform 1D PCA synthesis (spectra) on pixel (i,j)
            x(:,i,j) = v*wvltx(:,i,j);
        end
    end
end
