% At_fw.m
%
% Adjoint for "scrambled Fast Walsh" measurements for wavelet sparse signals.
%
% Usage: x = At_fw(b, OMEGA, idx, P)
%
% b - K length measurement vector
%
% OMEGA - K vector denoting which Walsh coefficients to use
%
% idx - Permutation vector for Walsh transform
%
% P - Permutation to apply to the input vector.  
%
% x - N length signal vector
%
% Written by: Marco F. Duarte, Rice University
% Created: October 2006
%
% Based on code by Justin Romberg, Georgia Institute of Technology

function x = At_fw(b, OMEGA, idx, P)

N = length(P);

K = length(b);
fx = zeros(N,1);
fx(OMEGA) = b/sqrt(N);
x = zeros(N,1);
% x(P) = ifastwalsh(fx,idx);
x(P) = fasterwalsh(fx,idx);