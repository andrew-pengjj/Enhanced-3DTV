% iwvlt_hs.m
%
% 3D inverse wavelet transform for hyperspectral images
%
% Usage: x = iwvlt_hs(w, h, S)
%
% w - input wavelet coefficient vector
%
% h - Wavelet filter
%
% S - number of spectral bands
%
% x - output data cube
%
% Written by: Marco F. Duarte, Rice University
% Created: November 5 2006
%
% Uses Rice Wavelet Toolbox, http://dsp.rice.edu/rwt

function x = iwvlt_hs(w, h, S)

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
            % Perform 1D DWT (spectra) on pixel (i,j)
            x(:,i,j) = midwt(wvltx(:,i,j),h,log2(S));
        end
    end
end
