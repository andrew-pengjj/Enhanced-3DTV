% tvqc_quantization_example.m
%
% Takes random measurements of the 'boats' image, quantizes, and
% reconstructs using quadractically constrained TV.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

path(path, './Optimization');
path(path, './Measurements');
path(path, './Data');

% boats
n = 256;
N = n*n;
load boats
img = boats;
I = img/norm(img(:));
I = I - mean(I(:));
x = reshape(I,N,1);

% num obs
K = 25000;

% permutation P and observation set OMEGA
P = randperm(N)';
q = randperm(N/2-1)+1;
OMEGA = q(1:K/2)';

% measurement implicit matrices
A = @(z) A_f(z, OMEGA, P);
At = @(z) At_f(z, N, OMEGA, P);

% take measurements
Ax = A(x);
% quantization codebook
mxax = max(abs(Ax));
bins = 10;
C = 2*mxax*(linspace(1/(2*bins),1-1/(2*bins),bins)-1/2);
% quantize
b = zeros(K,1);
symbols = zeros(K,1);
for kk = 1:K
  [mn,mni] = min(abs(Ax(kk)-C));
  b(kk) = C(mni);
  symbols(kk) = mni;
end
% error (just for future reference)
e = b - Ax;  

% the error should be roughly this size
epsilon = (2*mxax/bins)/sqrt(12)*sqrt(K);


% initial point  = min l2
x0 = At(b);


tvI = sum(sum(sqrt([diff(I,1,2) zeros(n,1)].^2 + [diff(I,1,1); zeros(1,n)].^2 )));
disp(sprintf('original TV = %8.3f', tvI));

time0 = clock;
xp = tvqc_logbarrier(x0, A, At, b, epsilon, 1e-1, 5);
Ip = reshape(xp, n, n);
disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));
