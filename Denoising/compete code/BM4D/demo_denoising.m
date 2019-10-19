% --------------------------------------------------------------------------------------------
%
%     Demo software for BM4D volumetric data denoising 
%               Release ver. 2.0  (17 April 2012)
%
% --------------------------------------------------------------------------------------------
%
% The software implements the BM4D denoising algorithm described in the papers:
%
% M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, "A Nonlocal Transform-Domain Filter 
% for Volumetric Data Denoising and Reconstruction", submitted to IEEE TIP, 2011.
%
% M. Maggioni, A. Foi, "Nonlocal Transform-Domain Denoising of Volumetric Data With 
% Groupwise Adaptive Variance Estimation", in Proc. of SPIE Electronic Imaging (EI), 
% Jan. 2012, San Francisco, CA, USA.
%
% --------------------------------------------------------------------------------------------
%
% authors:               Matteo Maggioni
%                        Alessandro Foi
%
% web page:              http://www.cs.tut.fi/~foi/GCF-BM3D
%
% contact:               firstname.lastname@tut.fi
%
% --------------------------------------------------------------------------------------------
% Copyright (c) 2010-2012 Tampere University of Technology.
% All rights reserved.
% This work should be used for nonprofit purposes only.
% --------------------------------------------------------------------------------------------
%
% Disclaimer
% ----------
%
% Any unauthorized use of these routines for industrial or profit-oriented activities is
% expressively prohibited. By downloading and/or using any of these files, you implicitly
% agree to all the terms of the TUT limited license (included in the file Legal_Notice.txt).
% --------------------------------------------------------------------------------------------
%

clear all;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modifiable parameters
phantom           = 't1_icbm_normal_1mm_pn0_rf0.rawb'; % name of the phantom raw data
sigma             = 11;      % noise standard deviation (%)
is_rician         = 0;       % noise distribution (Gaussian --> 0, Rician --> 1)
sharpAlphaRoot    = 1.0;     % alpha-rooting parameter for sharpening (1 means no sharpening)
do_wiener         = 1;       % perform collaborative wiener filtering
use_mod_profile   = 1;       % always use the modified profile (refer to [1] in the README)
verbose           = 1;       % verbose mode

crop_phantom      = 1;       % experiment on smaller phantom
estimate_sigma    = 0;       % enable sigma estimation
save_mat          = 0;       % save result to matlab .mat file
variable_noise    = 0;       % enable spatially varying noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MODIFY BELOW THIS POINT ONLY IF YOU KNOW WHAT YOU ARE DOING       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate original and noisy phantom
if ~exist(phantom,'file')
    error(['Could not read BrainWeb phantom "',phantom,'" file. ',...
        'Please download it from http://mouldy.bic.mni.mcgill.ca/'...
        'brainweb/selection_normal.html ',...
        '(Modality T1, Slice thickness 1mm, ',...
        'Noise 0%, Intensity non-uniformity 0%)']);
end
y = reshape(fread(fopen(phantom),181*217*181),[181 217 181])/255;
if crop_phantom
    y = y(51:100,51:100,51:100);
end
if variable_noise
    disp(' ')
    disp('Spatially-varying noise ')
    estimate_sigma = 1;
    % noise modulation field
    s = size(y);
    map = ones(3,3,3);
    map(2,2,2)=3;
    [x1,y1,z1] = meshgrid(1:3,1:3,1:3);
    [x2,y2,z2] = meshgrid(1:2/(s(2)-1):3,1:2/(s(1)-1):3,1:2/(s(3)-1):3);
    map = interp3(x1,y1,z1,map,x2,y2,z2,'cubic'); 
    clear s x2 y2 z2;
else
    disp('Uniform noise ')
    map = 1;
end
randn('seed',0);
rand('seed',0);
if is_rician
    z = sqrt( (y+sigma/100*map.*randn(size(y))).^2 + (sigma/100*map.*randn(size(y))).^2 );
else
    z = y + sigma/100*map.*randn(size(y));
end

% perform filtering
disp('Denoising started (might take a while...)')
[PSNR, y_est, sigma_est] = bm4d(y, z, (~estimate_sigma)*sigma/100, ...
    is_rician, sharpAlphaRoot, do_wiener, use_mod_profile, verbose);

% plot historgram of the estimated standard deviation
if estimate_sigma && ~variable_noise
    figure,
    hold all
    hist(sigma_est(y>10/255)*100,50);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','w','EdgeColor','k')
    line([sigma sigma],get(gca,'ylim'), 'Color','k','LineWidth',2);
    box
    xlabel('$\sigma (\%)$','Interpreter','Latex')
    ylabel('Voxels','Interpreter','Latex')
    hold off
end

% show 3-D cross-sections
volumes = {y,z,y_est};
fig_title = {'Original phantom',...
    sprintf('Noisy phantom - sigma %.2f%%',sigma),...
    sprintf('Final estimate - %.2fdB',PSNR)};
visualizeXsect( volumes, fig_title );

% save experiments workspace
if save_mat
    save([phantom,'_sigma',num2str(sigma),'.mat'],...
        'y','z','y_est','sigma_est','sigma','PSNR')
end
