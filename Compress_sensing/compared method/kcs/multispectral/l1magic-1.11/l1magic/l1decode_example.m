% l1decode_example.m
%
% Test out l1decode code.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

path(path, './Optimization');

% source length
N = 256;

% codeword length
M = 4*N;

% number of perturbations
T = round(.2*M);

% coding matrix
G = randn(M,N);

% source word
x = randn(N,1);

% code word
c = G*x;

% channel: perturb T randomly chosen entries
q = randperm(M);
y = c;
y(q(1:T)) = randn(T,1);

% recover
x0 = inv(G'*G)*G'*y;
xp = l1decode_pd(x0, G, [], y, 1e-4, 30);

% large scale
% gfun = @(z) G*z;
% gtfun = @(z) G'*z;
% xp = l1decode_pd(x0, gfun, gtfun, y, 1e-3, 25, 1e-8, 200);
