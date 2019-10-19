function testfig6(numbands)
% This script generates data for Figure 6 of 
% "Kronecker Compressive Sensing" 
% by Marco F. Duarte and Richard G. Baraniuk
%
% To generate the files fig6results_XX_YY.mat, 
% XX = 5, 10, ... 50, run 'testfig8(YY)', 
% setting YY = 16, 32, 64, 128.
%
% Written by: Marco F. Duarte, Duke University
% Created: March 9 2011


load xhs_others
xhs = xhs2;
sidelength = 128;
%numbands = 128;
if ((128/numbands)-round(128/numbands)),
    error('Number of bands must be a power of 2 between 1 and 128')
end
eval(['load WalshParams_' num2str(sidelength)])
eval(['load WalshParams_' num2str(sidelength) '_' num2str(numbands)])
spoff = [permx(1)-(sidelength*(ceil(permx(1)/sidelength)-1)) ceil(permx(1)/sidelength)]

if numbands == 128,
    xh = xhs;
else
    xh = zeros(size(xhs,1),numbands);
    for i=1:numbands,
        xh(:,i) = sum(xhs(:,((i-1)*(128/numbands)+1):(i*(128/numbands))),2);
    end
end

h = daubcqf(8);
L = log2(sidelength);

pkeep=1/4;
    % Calculate measurements
%     M = round(pkeep*sidelength*sidelength);
%     OMEGA = permy(1:M);
%     
%     Af = @(z) A_fw(z, OMEGA, idx, permx);
%     yhyper = zeros(M,size(xh,2));
%     for i=1:size(xh,2),
%         yhyper(:,i) = Af(xh(:,i));
%     end

    % Calculate global measurements
    M2 = round(pkeep*sidelength*sidelength*numbands);
    OMEGA2 = permy2(1:M2);
    Af2 = @(z) A_fw(z, OMEGA2, idx2, permx2);
    xh2 = xh';
    yhyper2 = Af2(xh2(:));

%     clear xl1all xl1ind xl1glo xl1fall xl1fglo tl1fall tl1all tl1fglo tl1ind tl1glo
%     disp('All, wavelet/wavelet')
%     tic
%     xl1all = cam_recon_l1meq_fun(yhyper, pkeep, sidelength);
%     xl1all(:,spoff(1),spoff(2)) = xl1all(:,spoff(1),spoff(2)+1);
%     tl1all = toc;
%     eval(['save fig6results_' num2str(100*pkeep) '_' num2str(numbands) ' xl1all tl1all pkeep'])
%     disp('Independent')
%     tic
%     xl1ind = cam_recon_l1meq(yhyper, pkeep, sidelength);
%     xl1ind(:,spoff(1),spoff(2)) = xl1ind(:,spoff(1),spoff(2)+1);
%     tl1ind = toc;
%     eval(['save fig6results_' num2str(100*pkeep) '_' num2str(numbands) ' xl1all xl1ind tl1all tl1ind pkeep'])
%     disp('Global, wavelet/wavelet')
%     tic
%     xl1glo = cam_recon_l1meq_funall(yhyper2, pkeep, sidelength, numbands);
%     xl1glo(:,spoff(1),spoff(2)) = xl1glo(:,spoff(1),spoff(2)+1);
%     tl1glo = toc;
%     eval(['save fig6results_' num2str(100*pkeep) '_' num2str(numbands) ' xl1all tl1all xl1ind tl1ind xl1glo tl1glo pkeep'])
%     disp('All, wavelet/Fourier')
%     tic
%     xl1fall = cam_recon_l1meq_funfourier(yhyper, pkeep, sidelength);
%     xl1fall(:,spoff(1),spoff(2)) = xl1fall(:,spoff(1),spoff(2)+1);
%     tl1fall = toc;
%     eval(['save fig6results_' num2str(100*pkeep) '_' num2str(numbands) ' xl1all tl1all xl1ind tl1ind xl1glo tl1glo xl1fall tl1fall pkeep'])
%     disp('Global, wavelet/Fourier')
%     tic
    xl1fglo = cam_recon_l1meq_funfourierall(yhyper2, pkeep, sidelength, numbands);
    xl1fglo(:,spoff(1),spoff(2)) = xl1fglo(:,spoff(1),spoff(2)+1);
%     tl1fglo = toc;
%     eval(['save fig6results_' num2str(100*pkeep) '_' num2str(numbands) ' xl1all tl1all xl1ind tl1ind xl1glo tl1glo xl1fall tl1fall xl1fglo tl1fglo pkeep'])
%     disp('All, wavelet/KLT')
%     tic
%     xl1pall = cam_recon_l1meq_funpca(yhyper, pkeep, sidelength);
%     xl1pall(:,spoff(1),spoff(2)) = xl1pall(:,spoff(1),spoff(2)+1);
%     tl1pall = toc;
%     eval(['save fig6results_' num2str(100*pkeep) '_' num2str(numbands) ' xl1all tl1all xl1ind tl1ind xl1glo tl1glo xl1fall tl1fall xl1fglo tl1fglo xl1pall tl1pall pkeep'])
%     disp('Global, wavelet/KLT')
%     tic
%     xl1pglo = cam_recon_l1meq_funpcaall(yhyper2, pkeep, sidelength, numbands);
%     xl1pglo(:,spoff(1),spoff(2)) = xl1pglo(:,spoff(1),spoff(2)+1);
%     tl1pglo = toc;
%     eval(['save fig6results_' num2str(100*pkeep) '_' num2str(numbands) ' xl1all tl1all xl1ind tl1ind xl1glo tl1glo xl1fall tl1fall xl1fglo tl1fglo xl1pall tl1pall xl1pglo tl1pglo pkeep'])

