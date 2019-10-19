function S = DimDetectAIC( D )

%==========================================================================
% Detect the dimensionality of inherent subspace of data matrix D via
%   Akaike Information Criteria (AIC[2]) following [1]. Each row of D is a
%   P-dimensional vector (a sample of a P-dimensional signal / an
%   observation of a P-dimensional random vector) and D has N rows in total
%   (N samples / N observations).
%
% S = DimDetectAIC( D ) returns the AIC scores of Q-dimensional model where
%   Q varies from 1 to P-1. Thus S is a (P-1)-dimensional score vector.
%
% [1] M. Wax, T. Kailath, "Detection of Signals by Information Theoretic
%     Criteria", IEEE International Conference on Acoustics Speech and
%     Signal Processing, vol. 33, pp. 387-392, 1985.
% [2] H. Akaike, "A new look at the statistical model identification", IEEE
%     Trans. Automat. Contr., vol. AC-19, pp. 716-723, 1974.
%
% See also DimDetectMDL
%
% by Yi Peng
%==========================================================================

[n, p] = size(D);

% remove the average from each random vector
Dmean = mean(D);
D = D - repmat(Dmean, n, 1);

% only D's singular value is used, so choose the smaller covariance matrix
if n >= p
    [~, l] = eig(D'*D/n);
else
    [~, l] = eig(D*D'/n);
    tmp = n; n = p; p = tmp;
end
q = 1:p-1;
l = diag(l)';               % eigenvalues of data covariance, ascending order
ari_avg = cumsum(l)./(1:p);
ari_avg = ari_avg(p-q);     % arithmetic average of smallest eigenvalues
log_sum = cumsum(log(l));
log_sum = log_sum(p-q);     % summation of log of smallest eigenvalues
fv = 2*p*q - q.^2 + 1;      % number of free adjusted variables
S = -2*n*log_sum(q) + 2*n*(p-q).*log(ari_avg(q)) + 2*fv;