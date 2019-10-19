% A_fw.m
%
% Takes "scrambled Fast Walsh" measurements
%
% Usage: b = A_fw(w, OMEGA, idx, P)
%
% w - N length signal vector
%
% b - K vector
%
% OMEGA - K vector denoting which Walsh coefficients to use
%
% idx - Permutation vector for Walsh transform
%
% P - Permutation to apply to the input vector.
%
% Written by: Marco F. Duarte, Rice University
% Created: December 2006
%
% Based on code by Justin Romberg, Georgia Institute of Technology

function b = A_fw(x, OMEGA, idx, P)

N = length(x);

x = x(:);
% fx = fastwalsh(x(P),idx)*sqrt(N);
fx = fasterwalsh(x(P),idx)/sqrt(N);
b = fx(OMEGA);
