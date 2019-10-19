function S = DimDetectMDL( D )

%==========================================================================
% Detect the dimensionality of inherent subspace of data matrix D via
%   Minimum Description Length (MDL[2][3]) following [1]. Each row of D is
%   a P-dimensional vector (a sample of a P-dimensional signal / an
%   observation of a P-dimensional random vector) and D has N rows in total
%   (N samples / N observations).
%
% S = DimDetectMDL( D ) returns the MDL scores of Q-dimensional model where
%   Q varies from 1 to P-1. Thus S is a (P-1)-dimensional score vector.
%
% [1] M. Wax, T. Kailath, "Detection of Signals by Information Theoretic
%     Criteria", IEEE International Conference on Acoustics Speech and
%     Signal Processing, vol. 33, pp. 387-392, 1985.
% [2] G. Schwartz, "Estimating the dimension of a model", Ann. Stat., vol.
%     6, pp. 461-464, 1978.
% [3] J. Rissanen, "Modeling by shortest data description", Automatica,
%     vol. 14, pp. 465-471, 1978.
%
% See also DimDetectAIC
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
if size(l, 2) ~= p
    printf('Error: dimension do not agree.\n');
end
ari_avg = cumsum(l)./(1:p);
ari_avg = ari_avg(p-q);     % arithmetic average of smallest eigenvalues
log_sum = cumsum(log(l));
log_sum = log_sum(p-q);     % summation of log of smallest eigenvalues
fv = 2*p*q - q.^2 + 1;      % number of free adjusted variables
S = -n*log_sum(q) + n*(p-q).*log(ari_avg(q)) + 0.5*fv*log(n);