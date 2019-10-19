% At_cscamtensorfourier.m
%
% Adjoint for Mersenne measurements for hyperspectral 
% signals with wavelet/fourier tensor product sparsity
%
% Usage: w = At_cscamtensorfourier(b, OMEGA, idx, P, h)
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

function [w,x,tmp,bm] = At_cscamtensorfourier(b, OMEGA, idx, P, h)

N = sqrt(length(P));
M = length(OMEGA);
S = length(b)/M;
bm = reshape(b,M,S);

% Apply transpose of measurement matrix
tmp = zeros(N^2,S);
for i=1:S,
    fx = zeros(N^2,1);
    fx(OMEGA) = bm(:,i);
    % tmp(P,i) = ifastwalsh(fx,idx)/N;
    tmp(P,i) = fasterwalsh(fx,idx)/N;
end

x = reshape(tmp',[S N N]);
% Apply wavelet transform
w = wfft_hs(x,h);

% Vectorize
w = w(:);
