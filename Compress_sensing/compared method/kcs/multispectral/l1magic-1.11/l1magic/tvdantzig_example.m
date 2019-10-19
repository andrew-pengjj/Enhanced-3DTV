% tvdantzig_example.m
%
% Test out tvdantzig code (TV minimization with bounded residual correlation).
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

% use implicit, matrix-free algorithms ?  
largescale = 0;

path(path, './Optimization');
path(path, './Measurements');
path(path, './Data');


% test image = 32x32 piece of cameraman's arm
load camera
I = camera(81:112,37:68);
n = 32;
N = n*n;
I = I/norm(I(:));
I = I - mean(I(:));
x = reshape(I,N,1);


% num obs
K = 300;

% permutation P and observation set OMEGA
P = randperm(N)';
q = randperm(N/2-1)+1;
OMEGA = q(1:K/2)';



% measurement matrix
if (largescale)
  A = @(z) A_f(z, OMEGA, P);
  At = @(z) At_f(z, N, OMEGA, P);
  % obsevations
  b = A(x);
  % initial point
  x0 = At(b);
else
  FT = 1/sqrt(N)*fft(eye(N));
  A = sqrt(2)*[real(FT(OMEGA,:)); imag(FT(OMEGA,:))];
  A = [1/sqrt(N)*ones(1,N); A];
  At = [];
  % observations
  b = A*x;
  % initial point
  x0 = A'*b;
end


epsilon = 5e-3;



tvI = sum(sum(sqrt([diff(I,1,2) zeros(n,1)].^2 + [diff(I,1,1); zeros(1,n)].^2 )));
disp(sprintf('Original TV = %.3f', tvI));

time0 = clock;
xp =  tvdantzig_logbarrier(x0, A, At, b, epsilon, 1e-3, 5, 1e-8, 1500);
Ip = reshape(xp, n, n);
disp(sprintf('Total elapsed time = %f secs\n', etime(clock,time0)));



                                                                                                                                     
