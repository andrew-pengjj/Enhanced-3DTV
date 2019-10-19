% wpca_hs.m
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

function w = wpca_hs(x, h, v)

N = size(x,2);
Y = size(x,1);
if Y > 1,
    wvlts = zeros(size(x));
    for i=1:N,
        for j=1:N,
            % Perform 1D PCA analysis (spectra) on pixel (i,j)
            wvlts(:,i,j) = v'*x(:,i,j);
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
