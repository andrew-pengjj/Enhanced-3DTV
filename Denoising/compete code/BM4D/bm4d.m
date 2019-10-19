function [PSNR, y_est, sigma_est] = bm4d(y, z, sigma, is_rician, alpha, do_wiener, use_mod_profile, print_to_screen)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM4D is an algorithm for attenuation of additive white Gaussian noise in 
%  volumetric data. This algorithm reproduces the results from the article:
%
%  [1] M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, "A Nonlocal 
%      Transform-Domain Filter for Volumetric Data Denoising and 
%      Reconstruction", submitted to IEEE TIP, 2011.
%      preprint at http://www.cs.tut.fi/~foi/GCF-BM3D.
%
%
%  FUNCTION INTERFACE:
% 
%  INPUTS (only z and sigma are required, see examples below):
%     1) y               (3D array) : noise-free volume (only for computing PSNR),
%                                     replace with the scalar 1 if not available.
%     2) z               (3D array) : noisy volume (intensities in the range [0,1])
%     3) sigma             (double) : standard deviation of the noise
%                                     (in the range [0,1])
%                                     If unknown put 0 to enable noise estimation
%     4) is_rician        (logical) : 0 --> z is Gaussian distributed 
%                                   : 1 --> z is Rician distributed
%     5) alpha             (double) : alpha-rooting parameter for sharpening
%                                     (default is 1, which means no sharpening)
%     6) do_wiener        (logical) : perform collaborative wiener filtering
%                                     (default is 1)
%     7) use_mod_profile  (logical) : 0 --> 'normal' Profile 
%                                   : 1 --> 'modified' Profile
%                                     (default is 1)
%     8) print_to_screen  (logical) : 0 --> do not print output information
%                                     1 --> print information to screen
%                                     (default is 0)
%
%  OUTPUTS:
%     1) PSNR              (double) : Output PSNR (dB), only if the original 
%                                     volume y is available, otherwise PSNR = 0
%     2) y_est           (3D array) : Final estimate (in the range [0,1])
%     3) sigma_est       (3D array) : Voxel-wise standard deviation estimate
%
%
%  TYPICAL USAGE EXAMPLES:
%
%  Case: Original volume y is available
%   % load volume and scale to range [0,1]
%      y = fread(fopen(phantom));
%   % Add the AWGN with zero mean and standard deviation 'sigma' in [0,1]
%      z = y + sigma*randn(size(y));
%   % Denoise 'z'. The denoised volume is 'y_est', and 'PSNR' is 
%   % equal to 10*log10(1/mean((y(:)-y_est(:)).^2))
%      [PSNR, y_est] = bm4d(y, z, sigma);
%
%  Case: Only noisy volume z is available
%   % Denoise 'z'. The denoised volume is 'y_est', and 'dummy = 1' because 
%   % the noise-free volume was not provided.
%      [dummy, y_est] = bm4d(1, z, sigma, is_rician);
%
%
%  NOTES:
%
%  This function interface is based upon the BM3D.m function by K. Dabov,
%  with the main notable difference in the scaling of the sigma input,
%  which here is consistent with that of the signal (whereas in BM3D.m the
%  sigma input is always scaled with respect to [0,255] range, even when z
%  is given on the [0,1] range).
%
%  The external VST package downloadable at http://www.cs.tut.fi/~foi/RiceOptVST/
%  is required to process Rician-distributed data
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2010-2011 Tampere University of Technology.
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% AUTHOR:
%     Matteo Maggioni, email: matteo.maggioni _at_ tut.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Specify the std. dev. of the corrupting noise
if (exist('sigma','var') ~= 1),
    sigma = 0; %% default std of the AWGN (enables noise estimation)
end
if (exist('is_rician','var') ~= 1),
    is_rician = 0; %% default noise distribution
end
if (exist('alpha','var') ~= 1),
    alpha = 1.0; %% default sharpening parameter
end
if (exist('do_wiener','var') ~= 1),
    do_wiener = 1; %% 
end
if (exist('use_mod_profile','var') ~= 1),
    use_mod_profile = 1; %%
end
dump_output_information = 1;
if (exist('print_to_screen','var') == 1) && (print_to_screen == 0),
    dump_output_information = 0; %% 
end
background = -1; %% intensities below this threshold will not be included in PSNR calculation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Following are the parameters for the Normal Profile.
%%%%

%%%% Select transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
transform_3D_HT_name     = 'bior1.5'; %% transform used for the HT filtering
transform_3D_Wiener_name = 'dct';     %% transform used for the Wiener filtering
transform_4D_name        = 'haar';    %% transform used in the grouping dim, 
                                      %  the same for HT and Wiener filt.

%%%% Hard-thresholding (HT) parameters:
N1                  = 4;    %% cube has size (N1 x N1 x N3)
N3                  = N1;   %% cube has size (N1 x N1 x N3)
Nstep               = 3;    %% sliding step to process every next reference cube
N2                  = 16;   %% maximum number of similar cubes
Ns                  = 11;   %% length of the side of the search neighborhood for full-search cube-matching
tau_match           = 3000; %% threshold for the cube-distance (d-distance)
lambda_thr4D        = 2.7;  %% threshold parameter for the hard-thresholding in 4D transform domain

%%%% Wiener filtering parameters:
N1_wiener           = N1;
N3_wiener           = N1_wiener;
Nstep_wiener        = Nstep;
N2_wiener           = 32;
Ns_wiener           = Ns;
tau_match_wiener    = 400;

%%%% Cube-matching parameters:
synchronous = 0;  %% if 1, the grouped cubes have coordinates lying in the same slice
decLevel    = 0;  %% decimation levels of the dyadic wavelet 2D transform 
                  %  0 means full decimation, higher values decrease the dec. number

if use_mod_profile
    N2                  = 32;
    N1_wiener           = 5;
    N3_wiener           = N1_wiener;
    lambda_thr4D        = 2.8;
    tau_match_wiener    = 3500;
    tau_match           = 25000;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create transform matrices, etc.
%%%%
[Tfor, Tinv]     = getTransfMatrix(N1, transform_3D_HT_name, decLevel);
[TforW, TinvW]   = getTransfMatrix(N1_wiener, transform_3D_Wiener_name, 0);
[Tfor3, Tinv3]   = getTransfMatrix(N3, transform_3D_HT_name, decLevel);
[Tfor3W, Tinv3W] = getTransfMatrix(N3_wiener, transform_3D_Wiener_name, 0);

Tfor4 = cell(1,max(N2,N2_wiener));
Tinv4 = cell(1,max(N2,N2_wiener));
for hpow = 0:ceil(log2(max(N2,N2_wiener))),
    h = 2^hpow;
    [Tfor4rd, Tinv4rd] = getTransfMatrix(h, transform_4D_name, 0);
    Tfor4{h} = double(Tfor4rd);
    Tinv4{h} = double(Tinv4rd');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% If needed, read images, generate noise, etc.
%%%%
if (exist('y','var') && ndims(y)==3 && ~exist('z','var'))
    if sigma~=0
        randn('seed', 0);
        if is_rician
            z = sqrt( (y+sigma.*randn(size(y))).^2 + (sigma.*randn(size(y))).^2 );
        else
            z = y + sigma*randn(size(y));
        end
    else
        error('Standard deviation sigma must be greater than 0.')
    end
else
    % convert z and y to double precision if needed
    z = double(z);
    y = double(y);
end
if (size(z,4) ~= 1) || (size(z,3) == 1),
    error('BM4D accepts only grayscale volumetric images.');
end
% Check if the true image y is a valid one; if not, then we cannot compute PSNR, etc.
y_is_invalid_image = (length(size(z)) ~= length(size(y))) | ...
    (size(z,1) ~= size(y,1)) | (size(z,2) ~= size(y,2)) | (size(z,3) ~= size(y,3));
if (y_is_invalid_image),
    dump_output_information = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Print image information to the screen
%%%%
if dump_output_information
    if is_rician
        fprintf('Rician noise distribution \n')
    else
        fprintf('Gaussian noise distribution. \n')
    end
    if sigma==0
        fprintf('Volumetric image: (%dx%dx%d) \nAdaptive noise variance estimation enabled \n', ...
            size(z,1), size(z,2), size(z,3));
    else
        fprintf('Volumetric image: (%dx%dx%d), sigma: %.1f%% \n', ...
            size(z,1), size(z,2), size(z,3), sigma*100);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Rician data forward stabilization
%%%%
if is_rician && ((sigma~=0 && ~exist('riceVST.m','file')) || ...
        (sigma==0 && ~exist('ricePairInversion.m','file')))
    error(['VST framework for Rician-distributed data not found. ',...
        'Please download it from http://www.cs.tut.fi/~foi/RiceOptVST/.'])
end
if is_rician && sigma~=0
    smoothVST = 'A';
    sigma_rice = sigma;
    sigma = 1;
    z = riceVST(z,sigma_rice,smoothVST);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalizing input
%%%%
maxtransformed = max(z(:));
mintransformed = min(z(:));
scale_range = 0.7;
scale_shift = (1-scale_range)/2;
z = (z-mintransformed) / (maxtransformed-mintransformed);
z = z*scale_range+scale_shift;
sigma = sigma/(maxtransformed-mintransformed);
sigma = sigma*scale_range;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 1. Produce the basic estimate by HT filtering
%%%%
tic;
if is_rician && sigma==0
    [y_hat, sigma_hat] = bm4d_thr_rice_mex(z, Nstep, N1, N2, N3,...
        lambda_thr4D, tau_match*N1*N1*N1/(255*255), (Ns-1)/2, synchronous, ...
        0, alpha, Tfor, Tinv', Tfor3, Tinv3', Tfor4, Tinv4 );
else
    [y_hat, sigma_hat] = bm4d_thr_mex(z, Nstep, N1, N2, N3,...
        lambda_thr4D, tau_match*N1*N1*N1/(255*255), (Ns-1)/2, synchronous, ...
        sigma, alpha, Tfor, Tinv', Tfor3, Tinv3', Tfor4, Tinv4 );
end
estimate_elapsed_time = toc;

if dump_output_information
    ind = y>background;
    fprintf('Basic estimate');
    if ~is_rician || (is_rician && sigma==0)
        PSNR_INITIAL_ESTIMATE = 10*log10(1 / mean(( y(ind)-...
            ((y_hat(ind)-scale_shift)/scale_range*...
            (maxtransformed-mintransformed)+mintransformed)).^2));
        fprintf(', PSNR: %.2f dB', PSNR_INITIAL_ESTIMATE);
    else
        PSNR_INITIAL_ESTIMATE = 10*log10(1 / mean(( y(ind)-...
            riceVST_EUI(((y_hat(ind)-scale_shift)/scale_range*...
            (maxtransformed-mintransformed)+mintransformed),sigma_rice,smoothVST) ).^2));
        fprintf(', PSNR: %.2f dB', PSNR_INITIAL_ESTIMATE);
    end
    fprintf(' \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2. Produce the final estimate by Wiener filtering (using the
%%%%  hard-thresholding initial estimate)
%%%%
tic;
if alpha~=1.0 || ~do_wiener
    y_est     = y_hat;
    sigma_est = sigma_hat;
else
    if is_rician && sigma==0
        [y_est, sigma_est] = bm4d_wie_rice_mex(z, y_hat, Nstep_wiener, N1_wiener, N2_wiener, N3_wiener, ...
            tau_match_wiener*N1_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, ...
            synchronous, 0, TforW, TinvW', Tfor3W, Tinv3W', Tfor4, Tinv4 );
    else
        [y_est, sigma_est] = bm4d_wie_mex(z, y_hat, Nstep_wiener, N1_wiener, N2_wiener, N3_wiener, ...
            tau_match_wiener*N1_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, ...
            synchronous, sigma, TforW, TinvW', Tfor3W, Tinv3W', Tfor4, Tinv4 );
    end
end
wiener_elapsed_time = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Restoring original range
%%%%
y_est     = (y_est-scale_shift)/scale_range;
y_est     = y_est*(maxtransformed-mintransformed)+mintransformed;
sigma_est = sigma_est/scale_range;
sigma_est = sigma_est*(maxtransformed-mintransformed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Rician data inverse stabilization
%%%%
if is_rician && sigma~=0
    y_est = riceVST_EUI(y_est,sigma_rice,smoothVST);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Calculate the final estimate's PSNR, print it, and show the
%%%% denoised image next to the noisy one
%%%%
PSNR = 0; %% Remains 0 if the true image y is not available
if ~y_is_invalid_image
    ind  = y>background;
    PSNR = 10*log10(1/mean((y(ind)-y_est(ind)).^2)); % y is valid
end

if dump_output_information
    fprintf('Final estimate, PSNR: %.2f dB \nTotal execution time: %.1f sec \n', ...
        PSNR, wiener_elapsed_time + estimate_elapsed_time);
end

return;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [Transf_Matrix, invTransf_Matrix] = getTransfMatrix (N, transform_type, Nden)
%
% Create forward and inverse transform matrices, which allow for perfect
% reconstruction. The forward transform matrix is normalized so that the
% l2-norm of each basis element is 1.
%
% [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
%  INPUTS:
%
%   N           --> Size of the transform (for wavelets, must be 2^K)
%
%   transform_type  --> 'dct', 'dst', 'hadamard', or anything that is
%                       listed by 'help wfilters' (bi-orthogonal wavelets)
%                       'DCrand' -- an orthonormal transform with a DC and all
%                       the other basis elements of random nature
%
%   dec_levels      --> If a wavelet transform is generated, this is the
%                       desired decomposition level. Must be in the
%                       range [0, log2(N)-1], where "0" implies
%                       full decomposition.
%
%  OUTPUTS:
%
%   Tforward        --> (N x N) Forward transform matrix
%
%   Tforward        --> (N x N) Inverse transform matrix
%

if ~exist('Nden','var')
    Nden = 0;
end

if N == 1,
    Transf_Matrix = 1;
elseif strcmp(transform_type, 'dct') == 1,
    Transf_Matrix    = dct(eye(N));
elseif strcmp(transform_type, 'dst') == 1,
    Transf_Matrix    = dst(eye(N));
elseif strcmp(transform_type, 'DCrand') == 1,
    x = randn(N); x(1:end,1) = 1; [Q,R] = qr(x);
    if (Q(1) < 0),
        Q = -Q;
    end;
    Transf_Matrix = Q';
elseif strcmp(transform_type, 'hadamard') == 1,
    Transf_Matrix    = hadamard(N);
else %% wavelet transform

    %%% Set periodic boundary conditions, to preserve bi-orthogonality
    dwtmode('per','nodisp');

    [LO_D, HI_D, LO_R, HI_R] = wfilters(transform_type);
    for i = 1:N
        Transf_Matrix(i,:)=waverec(circshift([1 zeros(1,N-1)],[Nden i-1]), ...
            2.^[Nden Nden:log2(N)], LO_D, -HI_D);  %% construct transform matrix
    end
end

%%% Normalize the basis elements
Transf_Matrix = (Transf_Matrix' * diag(sqrt(1./sum(Transf_Matrix.^2,2))))';

%%% Compute the inverse transform matrix
invTransf_Matrix = inv(Transf_Matrix);

return;

