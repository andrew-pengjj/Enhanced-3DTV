% A_cscamtensorfourier.m
%
% Takes Mersenne measurements of hyperspectral signals with 
% foirier/wavelet tensor product sparsity.
%
% Usage: b = A_cscamtensorfourier(w, OMEGA, idx, P, h)
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

function [b,bm,tmp,x] = A_cscamtensorfourier(w, OMEGA, idx, P, h)

N = length(P);
S = length(w)/N;
M = length(OMEGA);

% Reconstruct from wavelet coefficients
x = iwfft_hs(w,h,S);
% Apply measurement matrix to each slice
tmp = reshape(x,S,N)';
bm = zeros(M,S);
for i=1:S,
    % fx = sqrt(N)*fastwalsh(tmp(P,i),idx);
    fx = fasterwalsh(tmp(P,i),idx)/sqrt(N);
    bm(:,i) = fx(OMEGA);
end
b = bm(:);
