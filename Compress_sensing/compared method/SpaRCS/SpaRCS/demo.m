%Some simple test code for SpaRCS
clear all
close all

addpath(genpath('./algs/'))
addpath(genpath('./operators/'))
addpath(genpath('./PROPACK/'))
addpath(genpath('./utility/'))

%User params
N = 256;            %The size of the matrix
r = [2];            %The rank of L
K = ceil(.02*N^2);  %The sparsity level of S
p = floor(.2*N^2);  %The number of measurements to acquire

%Create the data matrix X = L+S
x = randn(N,r); y = randn(N,r);
L = x*y';

S = zeros(size(L)); 
indices = randperm(N^2); indices = indices(1:K);
S(indices) = randn(size(indices));
S = S*norm(L,'fro')/norm(S,'fro');

prob = 'noiselets';

switch prob
    case 'identity'
        %This solves the Robust PCA problem
        A = @(z) A_id(z);
        At = @(z) At_id(z,N,N);
	case 'noiselets'
        %This solves the SpaRCS problem
		idx = randperm(prod(size(L))); idx = idx(1:p);
		idx2 = randperm(N^2);
		A = @(z) Anoiselet(z,idx, idx2);
		At = @(z) Atnoiselet(z,idx, idx2, size(L));
	case 'robustMC'
        %Robust matrix completion
		idx = randperm(prod(size(L))); idx = idx(1:p);
		[ix, iy] = ind2sub(size(L), idx);
		idx = idx(:);
		A = @(z) z(idx);
		At = @(z) full(sparse(ix, iy, z, N, N));
        
        %Need to redefine S so that all K terms are supported by the
        %sampling set idx
        indices = randperm(length(idx)); indices = idx(indices); 
        indices = indices(1:K);
        S = zeros(size(L));
        S(indices) = randn(size(indices));
        S = S*norm(L,'fro')/norm(S,'fro');
end

b = A(L + S);

tic
[L1 S1 err] = sparcs(b, r, K, A, At, 'svds',5e-4, 20, 1); 
toc

disp(sprintf('RSNR for L: %3.3f dB', -20*log10(norm(L1 - L,'fro')/norm(L, 'fro'))))
disp(sprintf('RSNR for S: %3.3f dB', -20*log10(norm(S1 - S,'fro')/norm(S, 'fro'))))
disp(sprintf('RSNR for M_hat: %3.3f dB', -20*log10(norm(L+S-S1-L1,'fro')/norm(L+S, 'fro'))))
