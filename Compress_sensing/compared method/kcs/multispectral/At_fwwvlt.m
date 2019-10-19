% At_fwwvlt.m
%
% Adjoint for "scrambled fast Walsh" measurements for wavelet sparse signals.
%
% Usage: w = At_fwwvlt(b, N, OMEGA, P, h)
%
% b - K length measurement vector
%
% OMEGA - K vector denoting which Walsh coefficients to use
%
% idx - Permutation vector for Walsh transform
%
% P - Permutation to apply to the input vector.  
%
% h - Wavelet filter
%
% w - N wavelet coefficient vector for x(t)
%
% Written by: Marco F. Duarte, Rice University
% Created: October 2006
%
% Based on code by Justin Romberg, Georgia Institute of Technology
% Uses Rice Wavelet Toolbox, http://dsp.rice.edu/rwt


function w = At_fwwvlt(b, OMEGA, idx, P, h)

N = length(P);

K = length(b);
fx = zeros(N,1);
fx(OMEGA) = b;
x = zeros(N,1);
% x(P) = ifastwalsh(fx,idx)/sqrt(N);
x(P) = fasterwalsh(fx,idx)/sqrt(N);
w = mdwt(reshape(x,[sqrt(N) sqrt(N)]),h,log2(N)/2);
w = w(:);
