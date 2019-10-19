% wvlt_hs.m
%
% 3D wavelet transform for hyperspectral imaging
%
% Usage: w = wvlt_hs(x, h)
%
% x - input signal of size YxMxN, where MxN is the image size and Y is the
% number of spectral bands
%
% h - Wavelet filter name
%
% w - 3D wavelet coefficient vector for x
%
% Written by: Marco F. Duarte, Rice University
% Created: November 13 2006

function w = wvlt_hs(x, h)

N = size(x,2);
Y = size(x,1);
if Y > 1,
    wvlts = zeros(size(x));
    for i=1:N,
        for j=1:N,
            % Perform 1D DWT (spectra) on pixel (i,j)
            wvlts(:,i,j) = mdwt(x(:,i,j),h,log2(Y));
        end
    end
else
    wvlts = x;
end
w = zeros(size(wvlts));
for i=1:Y,
    % Perform 2D DWT on snapshot i (space)
    tmp = reshape(mdwt(squeeze(wvlts(i,:,:)),h,log2(N)),[1 N N]);
    w(i,:,:) = tmp(1,:,:);
end
w = w(:);
