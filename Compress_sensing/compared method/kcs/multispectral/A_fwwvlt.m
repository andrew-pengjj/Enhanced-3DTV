% A_fwwvlt.m
%
% Takes "scrambled fast Walsh" measurements for wavelet-sparse signals.
%
% Usage: b = A_fwwvlt(w, OMEGA, idx, P, h)
%
% w - N vector of wavelet coefficients
%
% b - K vector
%
% OMEGA - K vector denoting which Walsh coefficients to use
%
% idx - Permutation vector for Walsh transform
%
% P - Permutation to apply to the input vector.
%
% h - Wavelet filter
%
% Written by: Marco F. Duarte, Rice University
% Created: December 2006
%
% Based on code by Justin Romberg, Georgia Institute of Technology
% Uses Rice Wavelet Toolbox, http://dsp.rice.edu/rwt

function b = A_fwwvlt(w, OMEGA, idx, P, h)

N = length(w);

x = midwt(reshape(w,[sqrt(N) sqrt(N)]),h,log2(N)/2);
x = x(:);
% fx = sqrt(N)*fastwalsh(x(P),idx);
fx = fasterwalsh(x(P),idx)/sqrt(N);
b = fx(OMEGA);
