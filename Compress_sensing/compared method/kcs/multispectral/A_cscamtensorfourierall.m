% A_cscamtensorfourierall.m
%
% Takes Mersenne global measurements of hyperspectral signals with 
% foirier/wavelet tensor product sparsity.
%
% Usage: b = A_cscamtensorfourierall(w, OMEGA, idx, P, h)
%
% w -   SxN^2 vector of tensor coefficients
%
% b -   K*S vector = K-length measurement for each of S spectral 
%       bands, concatenated
%
% OMEGA - K vector denoting which Walsh coefficients to use
%
% idx - Permutation vector for Walsh transform
%
% P - Permutation to apply to the input vector.
%
% h -   Wavelet filter
%
% Written by: Marco F. Duarte, Rice University
% Created: June 25 2007
%
% Based on code by Justin Romberg, Georgia Institute of Technology

function [b,bm,tmp,x] = A_cscamtensorfourierall(w, OMEGA, idx, P, h, S)

N = length(P);
M = length(OMEGA);

% Reconstruct from wavelet coefficients
x = iwfft_hs(w,h,S);
x = x(:);
% Apply measurement matrix
% fx = sqrt(N)*fastwalsh(x(P),idx);
fx = fasterwalsh(x(P),idx)/sqrt(N);
b = fx(OMEGA);
