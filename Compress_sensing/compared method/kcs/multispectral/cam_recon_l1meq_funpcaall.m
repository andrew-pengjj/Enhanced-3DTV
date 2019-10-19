%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cam_recon_l1meq_funpcaall
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pkeep = percent to keep
% pkeep must be between 0 and 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Written by: Marco F. Duarte, Rice University
% Created: June 2006

function estIm = cam_recon_l1meq_funpcaall(coeffs, pkeep, sidelength, S);

%Setup Input
y = coeffs;

% Image size
N = (sidelength^2)*S;

% Number of measurements
M = (round(N*pkeep));
% Get reasonable number of samples
y = y(1:M);

% Load permutation vectors
eval(['load WalshParams_' num2str(sidelength) '_' num2str(S)]);
OMEGA = permy2(1:M);
eval(['load aviris' num2str(S) '_pca']);

h = daubcqf(8);
Af = @(z) A_cscamtensorpcaall(z, OMEGA, idx2, permx2, h, v);
Ab = @(z) At_cscamtensorpcaall(z, OMEGA, idx2, permx2, h, v);
S
%Result Values
x0 = Ab(y);
estCoeffs = l1eq_pd(x0, Af, Ab, y, norm(y)/1000);

%Result Image
estIm = iwpca_hs(estCoeffs,h,v);
