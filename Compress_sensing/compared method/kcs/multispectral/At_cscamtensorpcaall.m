% At_cscamtensorpca.m
%
% Adjoint for Mersenne measurements for hyperspectral 
% signals with wavelet/PCA tensor product sparsity
%
% Usage: w = At_cscamtensorpca(b, OMEGA, idx, P, h)
%
% b - K*S vector = K-length measurement vector for each of S spectral 
% bands, concatenated
%
% OMEGA - K vector denoting which Walsh coefficients to use
%
% idx - Permutation vector for Walsh transform
%
% P - Permutation to apply to the input vector.
%
% h - Wavelet filter
%
% w - N^2*S tensor coefficient vector for x(t)
%
% Written by: Marco F. Duarte, Rice University
% Created: June 26 2006
%
% Based on code by Justin Romberg, Georgia Institute of Technology
% Uses Rice Wavelet Toolbox, http://dsp.rice.edu/rwt

function [w,x,tmp] = At_cscamtensorpcaall(b, OMEGA, idx, P, h, v)

S = size(v,1);
N = sqrt(length(P)/S);

% Apply transpose of measurement matrix
fx = zeros(N^2*S,1);
fx(OMEGA) = b/(N*sqrt(S));
tmp = zeros(size(fx));
% tmp(P) = ifastwalsh(fx,idx);
tmp(P) = fasterwalsh(fx,idx);

x = reshape(tmp,[S N N]);
% Apply wavelet transform
w = wpca_hs(x,h,v);

% Vectorize
w = w(:);
