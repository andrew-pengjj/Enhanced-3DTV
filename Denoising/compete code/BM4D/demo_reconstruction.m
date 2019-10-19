% --------------------------------------------------------------------------------------------
%
%     Demo software for BM4D volumetric data reconstruction
%               Release ver. 2.0  (17 April 2012)
%
% --------------------------------------------------------------------------------------------
%
% The software implements the iterative reconstruction algorithm described in the paper:
%
% M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, "A Nonlocal Transform-Domain Filter
% for Volumetric Data Denoising and Reconstruction", submitted to IEEE TIP, 2011.
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
% Copyright (c) 2011-2012 Tampere University of Technology.
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

clear all;
close all;

% load constants
constantsSparseTraj3D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modifiable parameters
phantom    = SHEPPLOGAN;   % BRAINWEB SHEPPLOGAN
trajectory = RADIAL;       % RADIAL SPIRAL LOG_SPIRAL LIM_ANGLE SPHERICAL
data       = COMPLEX;      % REAL COMPLEX

n          = 128;          % size of the 3d phantom (power of 2 <=256)
data_std   = 5.0;          % AWGN standard deviation in the initial observed data (%)
coverage   = 30.0;         % percentage of sampled pixels (%)
low_pass   = 9;            % number of retained phase coefficients (per dimension)
excursion  = 4;            % phase excursion (>=1)
min_norm   = eps;          % minimum normalized p-norm required to continue
pnorm      = 2;            % norm type used in early-stop condition

iter_nbr   = 1e3;          % max number of iterations
alpha      = 1.010;        % alpha noise-excitation
beta       = 5.0e2;        % beta noise-excitation

tol        = 1.0;          % tolerance error of coverage (%)
rot_deg    = 1*1e0;        % rotation degrees between consecutive trajectories
line_nbr   = 1;            % number of subsampling trajectories per slice
line_std   = 0.0;          % noise in subsampling trajectories (%)

lapse      = 10;           % lapse between saved slices during reconstruction
verbose    = IMAGE;        % NONE TEXT IMAGE
save_mat   = 1;            % save result to matlab .mat file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       MODIFY BELOW THIS POINT ONLY IF YOU KNOW WHAT YOU ARE DOING       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if phantom==SHEPPLOGAN || ~exist('y','var')
    % load Shepp-Logan phantom
    load SheppLogan3D.mat
    y = double(SheppLogan3D{log2(n)});
    clear SheppLogan3D
end
if phantom==BRAINWEB
    % padding phantom
    brainweb = 't1_icbm_normal_1mm_pn0_rf0.rawb';
    if ~exist(brainweb,'file')
        error(['Could not read BrainWeb phantom "',brainweb,'" file. ',...
            'Please download it from http://mouldy.bic.mni.mcgill.ca/'...
            'brainweb/selection_normal.html ',...
            '(Modality T1, Slice thickness 1mm, ',...
            'Noise 0%, Intensity non-uniformity 0%)']);
    end
    % load BrainWeb phantom and scale it to size n
    y = double(reshape(fread(fopen(brainweb),181*217*181),[181 217 181]));
    w = zeros([217,217,217]);
    w(19:181+18,:,19:181+18) = y;
    % scale x-y
    s = n/217;
    w = imresize(w,s,'bilinear');
    % scale z
    y = zeros(n,n,n);
    for i=1:n
        for j=1:n
            y(i,j,:) = imresize(squeeze(w(i,j,:)),s,'bilinear');
        end
    end
    clear w s
end
% normalize phantom
y = (y-min(y(:))) / (max(y(:))-min(y(:)));
[size1 size2 size3] = size(y);

% BM4D parameters
sharpAlphaRoot    = 1.0;     % alpha-rooting parameter for sharpening (1 means no sharpening)
do_wiener         = 0;       % perform collaborative wiener filtering
use_mod_profile   = 0;       % modified profile

% setting transform functions
if trajectory==SPHERICAL
    % 3-D FFT
    transform  = @fftn;
    itransform = @ifftn;
else
    % multi slice 2-D FFT
    transform  = @msfft2;
    itransform = @imsfft2;
end

% initial seeds
randn('seed',0);
rand('seed',0);

% synthetic phase
if data==COMPLEX
    phi = 2*pi*0.25*randn(size(y));
    mask = zeros(size(phi));
    mask(1:low_pass,1:low_pass,1:low_pass) = 1;
    phi = idct3(dct3(phi).*mask);
    R = max(phi(:))-min(phi(:));
    phi = mod(2*pi*phi/R*max(1,excursion)+pi,2*pi)-pi;
    z = y.*exp(1i*phi);
else
    z = y;
end
% adding noise (if any) to initial data
data_std = data_std/100;
z = z + data_std/sqrt(2)*(randn(size(y)) + 1i*randn(size(y)));

% subsampling trajectory
S = sampling( trajectory, n, rot_deg, line_nbr, line_std/100, coverage/100, tol/100 );

% observed data and initial back-projections
theta_hat_0 = transform(z) .* S;
if data==COMPLEX
    y_hat_0     = abs(itransform(theta_hat_0));
    phi_hat_0   = angle(itransform(theta_hat_0));
else
    y_hat_0     = real(itransform(theta_hat_0));
    phi_hat_0   = zeros(size(y_hat_0));
end

% getting sigma to be used in the first iteration
if data_std~=0
    sigma0 = data_std;
else
    sigma0=1000/eps;
end

% excitation noise
min_std_excite  = 10^(-255/20)/sqrt(1-sum(S(:))/numel(y_hat_0));
min_std_excite  = (min_std_excite~=Inf) * min_std_excite;
std_excite      = max(min_std_excite,sqrt(alpha.^(-(1:iter_nbr)-beta))) + data_std;
std_excite(end) = 0;

% parameters' initialization
y_hat_k         = y_hat_0;
y_tilde_k       = y_hat_0;
y_hat_excite_k  = y_hat_0;
phi_hat_k       = phi_hat_0;
phi_tilde_k     = phi_hat_0;
sigma           = sigma0;
realcoverage    = sum(S(:))/numel(S)*100;
psnr_tilde      = zeros(1,iter_nbr);
psnr_hat        = zeros(1,iter_nbr);
ssim_tilde      = zeros(1,iter_nbr);
ssim_hat        = zeros(1,iter_nbr);
psnr_tilde_phi  = zeros(1,iter_nbr);
psnr_hat_phi    = zeros(1,iter_nbr);
ssim_tilde_phi  = zeros(1,iter_nbr);
ssim_hat_phi    = zeros(1,iter_nbr);
xSection        = ceil(size3/2);
psnr_ind        = y>10/255;
sw              = [1 1 1];

buffer_cplx     = (y_hat_0.*exp(1i*phi_hat_0)) ./ sigma0^2;
weight          = 1/sigma0^2;

prev_y_tilde_k  = abs(buffer_cplx) ./ weight;

% initializing progression variable
if data==COMPLEX
    progress        = zeros([size1 3*size2 iter_nbr/lapse+1]);
    progress(:,:,1) = [y_hat_0(:,:,xSection) y_hat_0(:,:,xSection) ...
        phi_hat_0(:,:,xSection)/2/pi+0.5];
else
    progress        = zeros([size1 2*size2 iter_nbr/lapse+1]);
    progress(:,:,1) = [y_hat_0(:,:,xSection) y_hat_0(:,:,xSection)];
end
idx = 2;

% print initial information
if verbose~=NONE
    fprintf('\n%s %s phantom of size %dx%dx%d \n', OBS{data},DATA{phantom},size1,size2,size3);
    fprintf('%s trajectory with %.2f%% subsampling \n', TRAJ{trajectory},realcoverage);
    if data_std>0
        fprintf('AWGN with sigma %.2f%% \n\n',data_std*100);
    end
end

start = tic;
k = 1;
early_stop = 0;
while k<=iter_nbr && ~early_stop
    % start counting
    iter = tic;
    
    if k>1
        % getting right sigma
        sigma = std_excite(k-1);
        
        % magnitude regularization
        [PSNR, y_hat_k] = bm4d(1, y_hat_excite_k, sigma, (data==COMPLEX),...
            sharpAlphaRoot, do_wiener, use_mod_profile, 0);
        
        % phase regularization
        if data==COMPLEX
            rand_phase_shift = 2*rand*pi-pi;
            [PSNR, phi_hat_k] = bm4d(1, mod(phi_hat_k+pi+rand_phase_shift,2*pi)/2/pi, sigma, ...
                0, sharpAlphaRoot, do_wiener, use_mod_profile, 0);
            phi_hat_k    = phi_hat_k*2*pi-pi-rand_phase_shift;
            phi_hat_k    = mod(phi_hat_k+pi,2*pi)-pi;
        end
        
        % buffer update
        weight = weight + 1/sigma^2;
        buffer_cplx = buffer_cplx + (y_hat_k.*exp(1i*phi_hat_k))/sigma^2;
        
        if data==COMPLEX
            y_tilde_k   = abs(buffer_cplx./ weight);
            phi_tilde_k = angle(buffer_cplx./ weight);
        else
            y_tilde_k   = real(buffer_cplx./ weight);
            phi_tilde_k = zeros(size(y_tilde_k));
        end
        
        % early stop condition
        early_stop      = (sum(abs(prev_y_tilde_k(:)-y_tilde_k(:)).^pnorm)).^(1/pnorm) / ...
            (numel(y_tilde_k)).^(1/pnorm) < min_norm;
        prev_y_tilde_k  = y_tilde_k;
        
        % excitation (prepare)
        if data==COMPLEX
            y_hat_excite_k = y_hat_k.*exp(1i*phi_tilde_k);
        else
            y_hat_excite_k = y_hat_k;
        end
    end
    
    % excitation (excite)
    y_hat_excite_k    = y_hat_excite_k  + std_excite(k)*randn(size(y_hat_excite_k)) + ...
        1i*std_excite(k)*randn(size(y_hat_excite_k));
    theta_hat_k       = transform(y_hat_excite_k);
    theta_hat_k(S==1) = theta_hat_0(S==1);
    i_t_theta_hat_k   = itransform(theta_hat_k);
    y_hat_excite_k    = abs(i_t_theta_hat_k);
    if data==COMPLEX
        phi_hat_k = angle(i_t_theta_hat_k).*(y_hat_excite_k~=0) + ...
            (y_hat_excite_k==0).*phi_hat_k;
    end
    
    % performances
    psnr_tilde(k)     = 10*log10(1/mean((y(psnr_ind)-y_tilde_k(psnr_ind)).^2));
    psnr_hat(k)       = 10*log10(1/mean((y(psnr_ind)-y_hat_k(psnr_ind)).^2));
    ssim_tilde(k)     = ssim_index3d(y_tilde_k*255,y*255,sw,psnr_ind);
    ssim_hat(k)       = ssim_index3d(y_hat_k*255,y*255,sw,psnr_ind);
    if data==COMPLEX
        psnr_tilde_phi(k) = 10*log10(1/mean(((phi(psnr_ind)-phi_tilde_k(psnr_ind))/2/pi).^2));
        psnr_hat_phi(k)   = 10*log10(1/mean(((phi(psnr_ind)-phi_hat_k(psnr_ind))/2/pi).^2));
        ssim_tilde_phi(k) = ssim_index3d((phi_tilde_k/2/pi+0.5)*255,(phi/2/pi+0.5)*255,sw,psnr_ind);
        ssim_hat_phi(k)   = ssim_index3d((phi_hat_k/2/pi+0.5)*255,(phi/2/pi+0.5)*255,sw,psnr_ind);
    end
    
    % stop counting
    iter_time = toc(iter);
    
    % saving cross section of the progression
    if mod(k,lapse)==0
        if data==COMPLEX
            progress(:,:,idx) = [y_hat_excite_k(:,:,xSection) ...
                y_tilde_k(:,:,xSection) phi_tilde_k(:,:,xSection)/2/pi+0.5];
        else
            progress(:,:,idx) = [y_hat_excite_k(:,:,xSection) y_tilde_k(:,:,xSection)];
        end
        idx = idx+1;
    end
    
    % storing data to disk
    if save_mat
        if data==COMPLEX
            save([DATA{phantom},'_',TRAJ{trajectory},'_cov',...
                num2str(coverage),'_sigma',num2str(data_std*100),'_COMPLEX.mat'],...
                'S','y','progress','psnr_hat','psnr_tilde','ssim_tilde','ssim_hat',...
                'psnr_hat_phi','psnr_tilde_phi','ssim_tilde_phi','ssim_hat_phi',...
                'iter_nbr','alpha','beta','std_excite','data_std','y_hat_k','y_tilde_k','y_hat_0',...
                'phi','phi_hat_k','phi_hat_0','phi_tilde_k','low_pass','excursion')
        else
            save([DATA{phantom},'_',TRAJ{trajectory},'_cov',...
                num2str(coverage),'_sigma',num2str(data_std*100),'_REAL.mat'],...
                'S','y','progress','psnr_hat','psnr_tilde','ssim_tilde','ssim_hat',...
                'iter_nbr','alpha','beta','std_excite','data_std','y_hat_k','y_tilde_k','y_hat_0')
        end
    end
    
    %%% OUTPUT CODE
    if verbose==IMAGE
        if ~ishandle(1)
            figure(1),
        end
        if data==COMPLEX
            imshow([ y(:,:,xSection) ...
                phi(:,:,xSection)/(2*pi)+0.5...
                y_hat_0(:,:,xSection) ...
                phi_hat_0(:,:,xSection)/(2*pi)+0.5;
                y_hat_k(:,:,xSection)...
                phi_hat_k(:,:,xSection)/(2*pi)+0.5...
                y_tilde_k(:,:,xSection) ...
                phi_tilde_k(:,:,xSection)/(2*pi)+0.5], ...
                'InitialMagnification','fit');
            title(sprintf(['Iteration #%04d (%.1fsec) - Excitation sigma %.2f%% \n',...
                'Magnitude PSNR %.2fdB SSIM %.2f - Phase PSNR %.2fdB SSIM %.2f \n',...
                '%s subsampling %.2f%% - Initial sigma %.2f%% \n'],...
                k, iter_time, sigma*100, psnr_tilde(k), ssim_tilde(k), ...
                psnr_tilde_phi(k), ssim_tilde_phi(k), TRAJ{trajectory}, ...
                realcoverage, data_std*100));
            text(size2*0.5,0*size1-(1.5),'$y$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','bottom');
            text(size2*1.5,0*size1-(1.5),'$\phi$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','bottom');
            text(size2*2.5,0*size1-(1.5),'$\hat{y}^{(0)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','bottom');
            text(size2*3.5,0*size1-(1.5),'$\hat{\phi}^{(0)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','bottom');
            text(size2*0.5,2*size1+(1.5),'$\hat{y}^{(k)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            text(size2*1.5,2*size1+(1.5),'$\hat{\phi}^{(k)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            text(size2*2.5,2*size1+(1.5),'$\tilde{y}^{(k)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            text(size2*3.5,2*size1+(1.5),'$\tilde{\phi}^{(k)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            pos_caption_y = 2;
            pos_caption_x = 2;
        else
            imshow([ y(:,:,xSection) ...
                y_hat_0(:,:,xSection) ...
                y_hat_excite_k(:,:,xSection) ...
                y_hat_k(:,:,xSection)...
                y_tilde_k(:,:,xSection) ], ...
                'InitialMagnification','fit');
            title(sprintf(['Iteration #%04d (%.1fsec) - ',...
                'Excitation sigma %.2f%% - PSNR %.2fdB \n',...
                '%s of subsampling %.2f%% - Initial sigma %.2f%% \n'],...
                k, iter_time, sigma*100, psnr_tilde(k), ...
                TRAJ{trajectory}, realcoverage, data_std*100));
            text(size2*0.5,size1+(1.5),'$y$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            text(size2*1.5,size1+(1.5),'$\hat{y}^{(0)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            text(size2*2.5,size1+(1.5),'$\hat{y}_{\mathrm{excite}}^{(k)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            text(size2*3.5,size1+(1.5),'$\hat{y}^{(k)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            text(size2*4.5,size1+(1.5),'$\tilde{y}^{(k)}$','Interpreter','Latex',...
                'HorizontalAlignment','center','VerticalAlignment','top');
            pos_caption_y = 1;
            pos_caption_x = 2.5;
        end
        text(pos_caption_x*size2,pos_caption_y*size1+4,...
            {' ',['Phantom size ',num2str(size1),'x',num2str(size2),'x',...
            num2str(size3),',  displaying cross-section ',num2str(xSection),' of ',num2str(size3)]},...
            'HorizontalAlignment','center','VerticalAlignment','top');
        pause(eps)
        drawnow
    end
    if verbose==TEXT
        if data==COMPLEX
            fprintf('#%04d | %.2fsec - sigma %.2f%% - Magnitude %.2fdB - Phase %.2fdB \n',...
                k, iter_time, std_excite(k)*100, psnr_tilde(k), psnr_tilde_phi(k));
        else
            fprintf('#%04d | %.2fsec - sigma %.2f%% - %.2fdB \n',...
                k, iter_time, std_excite(k)*100, psnr_tilde(k));
        end
    end
    k = k + 1;
end
exec_time = toc(start);
fprintf('Reconstruction from trajectory %s (%.2f%%) done (%.2fsec - %.2fdB) \n', ...
    TRAJ{trajectory}, realcoverage, exec_time, psnr_tilde(end));

if verbose==IMAGE
    % show PSNR progression
    figure,
    hold all
    plot(1:iter_nbr,psnr_tilde)
    plot(1:iter_nbr,psnr_hat)
    box
    grid
    h = legend('$\tilde{y}^{(k)}$','$\hat{y}^{(k)}$', 'Location','Best');
    set(h,'Interpreter','Latex')
    xlabel('Iteration','Interpreter','Latex','FontSize',11)
    ylabel('PSNR (dB)','Interpreter','Latex','FontSize',11)
    set(gca,'xtick',0:iter_nbr/10:iter_nbr)
    hold off
    
    % compare y_tilde and y_hat PSNRs
    figure,
    hold all
    plot(1:iter_nbr,psnr_tilde./psnr_hat,'-k')
    box
    grid
    xlabel('Iteration','Interpreter','Latex','FontSize',11)
    ylabel('PSNR $\tilde{y}^{(k)}$ / PSNR $\hat{y}^{(k)}$',...
        'Interpreter','Latex','FontSize',11)
    set(gca,'xtick',0:iter_nbr/10:iter_nbr)
    hold off
    
    if data==COMPLEX
        % same as before, but relative to phase reconstruction
        figure,
        hold all
        plot(1:iter_nbr,psnr_tilde_phi)
        plot(1:iter_nbr,psnr_hat_phi)
        box
        grid
        h = legend('$\tilde{\phi}^{(k)}$','$\hat{\phi}^{(k)}$', 'Location','Best');
        set(h,'Interpreter','Latex')
        xlabel('Iteration','Interpreter','Latex','FontSize',11)
        ylabel('PSNR (dB)','Interpreter','Latex','FontSize',11)
        set(gca,'xtick',0:iter_nbr/10:iter_nbr)
        hold off
        
        figure,
        hold all
        plot(1:iter_nbr,psnr_tilde_phi./psnr_hat_phi,'-k')
        box
        grid
        xlabel('Iteration','Interpreter','Latex','FontSize',11)
        ylabel('PSNR $\tilde{\phi}^{(k)}$ / PSNR $\hat{\phi}^{(k)}$',...
            'Interpreter','Latex','FontSize',11)
        set(gca,'xtick',0:iter_nbr/10:iter_nbr)
        hold off
    end
    
    % show 3-D cross-sections
    if data==COMPLEX
        volumes = {y,phi/2/pi+0.5,y_tilde_k,phi_tilde_k/2/pi+0.5};
        fig_title = {'Original magnitude','Original phase',...
            'Magnitude estimate','Phase estimate'};
        visualizeXsect( volumes, fig_title );
    else
        volumes = {y,y_tilde_k};
        fig_title = {'Original phantom','Final estimate'};
        visualizeXsect( volumes, fig_title );
    end
end
