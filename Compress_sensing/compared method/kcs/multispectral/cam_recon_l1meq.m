%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cam_recon_l1meq
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

function estIm = cam_recon_l1meq(coeffs, pkeep, sidelength);

%Setup Input
y = coeffs;

% Image size
N = sidelength^2;
S = size(y,2);

% Number of measurements
M = (round(N*pkeep));
% Get reasonable number of samples
y = y(1:M,:);

% Build Measurement Matrix
eval(['load WalshParams_' num2str(sidelength)]);
OMEGA = permy(1:M);

h = daubcqf(8);
L = log2(sidelength);
Af = @(z) A_fwwvlt(z, OMEGA, idx, permx, h);
Ab = @(z) At_fwwvlt(z, OMEGA, idx, permx, h);

%Result Values
estIm = zeros(S,sidelength,sidelength);
for i=1:S,
    x0 = Ab(y(:,i));
    tmp = l1eq_pd(x0, Af, Ab, y(:,i), norm(y(:,i))/1000);
    %Result Image
    estIm(i,:,:) = reshape(midwt(reshape(tmp,sqrt(N),sqrt(N)),h,L),[1 sidelength sidelength]);
end
