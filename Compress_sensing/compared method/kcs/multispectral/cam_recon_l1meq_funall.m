%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cam_recon_l1meq_funall
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

function estIm = cam_recon_l1meq_funall(coeffs, pkeep, sidelength, S);

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

h = daubcqf(8);
Af = @(z) A_cscamtensorall(z, OMEGA, idx2, permx2, h, S);
Ab = @(z) At_cscamtensorall(z, OMEGA, idx2, permx2, h, S);

%Result Values
x0 = Ab(y);
estCoeffs = l1eq_pd(x0, Af, Ab, y, norm(y)/1000);

%Result Image
estIm = iwvlt_hs(estCoeffs,h,S);
