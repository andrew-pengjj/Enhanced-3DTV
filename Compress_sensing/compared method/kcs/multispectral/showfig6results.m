function showfig6results(numbands,verbose)
% This script generates Figures 5 and 6 of 
% "Kronecker Compressive Sensing" by Marco F. Duarte and Richard G. Baraniuk
% Running the script figure6(XX) needs the files 
% fig6results_YY_1_XXX.mat, YY = 5, 10, ... 50,
% which are generated using testfig8.m
%
% Written by: Marco F. Duarte, Duke University
% Created: March 9 2011


if nargin < 1, verbose = 0; end
load xhs_others
eval(['load aviris' num2str(numbands) '_pca'])
xhs = xhs2;
bands=7;
sidelength = 128;
N = sidelength;
S = numbands;
if ((128/numbands)-round(128/numbands)),
    error('Number of bands must be a power of 2 between 1 and 128')
end

if numbands == 128,
    xh = xhs;
else
    xh = zeros(size(xhs,1),numbands);
    for i=1:numbands,
        xh(:,i) = sum(xhs(:,((i-1)*(128/numbands)+1):(i*(128/numbands))),2);
    end
end
x = reshape(xh',[numbands sidelength sidelength]);
h = daubcqf(8);
w = wvlt_hs(x,h);
ws = sort(abs(w(:)),'descend');
werr = (sum(ws.^2)-cumsum(ws.^2)).^0.5;

f = wfft_hs(x,h);
fs = sort(abs(f(:)),'descend');
ferr = (sum(fs.^2)-cumsum(fs.^2)).^0.5;

p = wpca_hs(x,h,v);
ps = sort(abs(p(:)),'descend');
perr = (sum(ps.^2)-cumsum(ps.^2)).^0.5;

% Get 3D wavelet filters (pick one)
%%% FARRAS
[af, sf] = farras;
J = log2(N)-3;
%%% HAAR
% af = [1 1;1 -1]/sqrt(2);sf = af(end:-1:1, :);
% J = log2(sidelength)-1;
if N == S,
    w3 = dwt3D_vectorized(x,J,af);
    w3s = sort(abs(w3(:)),'descend');
    w3err = (sum(w3s.^2)-cumsum(w3s.^2)).^0.5;
end
wsep = zeros(N^2,S);
for i=1:S,
    tmp = squeeze(x(i,:,:));
    tmpw = mdwt(tmp,h,log2(N));
    wsep(:,i) = sort(abs(tmpw(:)),'descend');
end
wseperr = (sum(wsep(:).^2)-sum(cumsum(wsep.^2),2)).^0.5;
wfreq = zeros(S,N,N);
for i=1:N,
    for j=1:N,
        tmpw = mdwt(x(:,i,j),h,log2(S));
        wfreq(:,i,j) = sort(abs(tmpw(:)),'descend');
    end
end
wfreqerr = (sum(wfreq(:).^2)-sum(sum(cumsum(wfreq.^2),2),3)).^0.5;
wfpca = zeros(S,N,N);
for i=1:N,
    for j=1:N,
        tmpw = v'*x(:,i,j);
        wfpca(:,i,j) = sort(abs(tmpw(:)),'descend');
    end
end
wfpcaerr = (sum(wfpca(:).^2)-sum(sum(cumsum(wfpca.^2),2),3)).^0.5;

% wss = sort(abs(ws(:)),'descend');
kmax = 100000;
figure(1),clf,plot((1:5000:kmax),20*log10(norm(ws(:))./werr(1:5000:kmax)),'o','LineWidth',2)
hold on%,plot((1:kmax),20*log10(norm(ws(:))./ferr(1:kmax)),'--or','LineWidth',2)
plot((1:5000:kmax),20*log10(norm(ws(:))./perr(1:5000:kmax)),'go','LineWidth',2)
if N == S,
    plot((1:kmax),20*log10(norm(w3s(:))./w3err(1:kmax)),':k','LineWidth',2)
end
plot(S*(1:ceil(kmax/S)),20*log10(norm(wsep(:))./wseperr(1:ceil(kmax/(S)))),'--k','LineWidth',2)
plot((N^2)*(1:ceil(kmax/(N^2))),20*log10(norm(wfreq(:))./wfreqerr(1:ceil(kmax/(N^2)))),'-','LineWidth',2)
%plot((N^2)*(1:ceil(kmax/(N^2))),20*log10(norm(wffour(:))./wffourerr(1:ceil(kmax/(N^2)))),'-.g','LineWidth',2)
plot((N^2)*(1:ceil(kmax/(N^2))),20*log10(norm(wfpca(:))./wfpcaerr(1:ceil(kmax/(N^2)))),'-.g','LineWidth',2)
plot((1:kmax),20*log10(norm(ws(:))./werr(1:kmax)),'-','LineWidth',2)
plot((1:kmax),20*log10(norm(ws(:))./perr(1:kmax)),'-.g','LineWidth',2)
axisfortex('','Number of coefficients, K','Transform coding compression SNR, dB')
axis([0 1e5 0 43])
if N == S,
    legend('Hyperbolic Wavelet','Wavelet/KLT','Isotropic Wavelet','Space Wavelet','Frequency Wavelet','Frequency KLT','Location','SouthEast')
else
    legend('Hyperbolic Wavelet','Wavelet/KLT','Space Wavelet','Frequency Wavelet','Frequency KLT','Location','SouthEast')
end
figure(4),clf,loglog((1:5000:kmax),werr(1:5000:kmax)./norm(ws(:)),'-o','LineWidth',2)
hold on%,loglog((1:kmax),ferr(1:kmax)./norm(ws(:)),'--','LineWidth',2)
loglog((1:5000:kmax),perr(1:5000:kmax)./norm(ws(:)),'go-.','LineWidth',2)
if N == S,
    loglog((1:kmax),w3err(1:kmax)./norm(w3s(:)),':k','LineWidth',2)
end
loglog(S*(1:ceil(kmax/S)),wseperr(1:ceil(kmax/(S)))./norm(wsep(:)),'--k','LineWidth',2)
loglog((N^2)*(1:ceil(kmax/(N^2))),wfreqerr(1:ceil(kmax/(N^2)))./norm(wfreq(:)),'-','LineWidth',2)
loglog((N^2)*(1:ceil(kmax/(N^2))),wfpcaerr(1:ceil(kmax/(N^2)))./norm(wfpca(:)),'-.g','LineWidth',2)
loglog((1:kmax),werr(1:kmax)./norm(ws(:)),'-','LineWidth',2)
loglog((1:kmax),perr(1:kmax)./norm(ws(:)),'g-.','LineWidth',2)
axisfortex('','Number of coefficients, K (log scale)','Normalized error magnitude (log scale), dB')
axis tight
if N == S,
    legend('Hyperbolic Wavelet','Wavelet/KLT','Isotropic Wavelet','Space Wavelet','Frequency Wavelet','Frequency KLT','Location','SouthWest')
else
    legend('Hyperbolic Wavelet','Wavelet/KLT','Space Wavelet','Frequency Wavelet','Frequency KLT','Location','SouthWest')
end
errormat = zeros(10,7);
timemat = zeros(10,7);
for pkeep = 5:5:50,
    pkeep2 = pkeep;
    eval(['load fig6results_' num2str(pkeep) '_' num2str(numbands)]) % xl1all xl1ind xsompall tsompall tl1all tl1ind'])
    errormat(pkeep2/5,1) = 20*log10(norm(x(:))/norm(xl1ind(:)-x(:)));
    errormat(pkeep2/5,2) = 20*log10(norm(x(:))/norm(xl1all(:)-x(:)));
    errormat(pkeep2/5,3) = 20*log10(norm(x(:))/norm(xl1glo(:)-x(:)));
    errormat(pkeep2/5,4) = 20*log10(norm(x(:))/norm(xl1fall(:)-x(:)));
    errormat(pkeep2/5,5) = 20*log10(norm(x(:))/norm(xl1fglo(:)-x(:)));
    errormat(pkeep2/5,6) = 20*log10(norm(x(:))/norm(xl1pall(:)-x(:)));
    errormat(pkeep2/5,7) = 20*log10(norm(x(:))/norm(xl1pglo(:)-x(:)));
    timemat(pkeep2/5,1) = tl1ind;
    timemat(pkeep2/5,2) = tl1all;
    timemat(pkeep2/5,3) = tl1glo;
    timemat(pkeep2/5,4) = tl1fall;
    timemat(pkeep2/5,5) = tl1fglo;
    timemat(pkeep2/5,6) = tl1pall;
    timemat(pkeep2/5,7) = tl1pglo;
end

if verbose,
    figure(11),clf,imagesc(makehstiles(x)),colormap gray,axis image,title('Original'),rmaxis(11)
    figure(12),clf,imagesc(makehstiles(xl1all),[min(x(:)) max(x(:))]),colormap gray,axis image,title('Distributed W/W'),rmaxis(12)
    figure(13),clf,imagesc(makehstiles(xl1glo),[min(x(:)) max(x(:))]),colormap gray,axis image,title('Global W/W'),rmaxis(13)
    figure(14),clf,imagesc(makehstiles(xl1fall),[min(x(:)) max(x(:))]),colormap gray,axis image,title('Distributed W/F'),rmaxis(14)
    figure(15),clf,imagesc(makehstiles(xl1fglo),[min(x(:)) max(x(:))]),colormap gray,axis image,title('Global W/F'),rmaxis(15)
    figure(16),clf,imagesc(makehstiles(xl1pall),[min(x(:)) max(x(:))]),colormap gray,axis image,title('Distributed W/P'),rmaxis(16)
    figure(17),clf,imagesc(makehstiles(xl1pglo),[min(x(:)) max(x(:))]),colormap gray,axis image,title('Global W/P'),rmaxis(17)
    figure(18),clf,imagesc(makehstiles(xl1ind),[min(x(:)) max(x(:))]),colormap gray,axis image,title('Separate'),rmaxis(18)
end

% K = (1:25)*(128^3)/200;
figure(2),clf,
plot((0.05:0.05:0.5),errormat(:,7),'-.og','LineWidth',2)
hold on
plot((0.05:0.05:0.5),errormat(:,6),'-.g','LineWidth',2)
plot((0.05:0.05:0.5),errormat(:,3),'-o','LineWidth',2)
% plot((0.05:0.05:0.5),errormat(:,5),'--or','LineWidth',2)
% plot((0.05:0.05:0.5),errormat(:,4),'--r','LineWidth',2)
plot((0.05:0.05:0.5),errormat(:,2),'LineWidth',2)
plot((0.05:0.05:0.5),errormat(:,1),'-.k','LineWidth',2)
axisfortex('','Normalized number of measurements, M/N','SNR, dB')
switch numbands
    case 16,
        axis([0.05 0.5 5 21])
    case 32,
        axis([0.05 0.5 5 21])
    case 64,
        axis([0.05 0.5 6 21])
    case 128,
        axis([0.05 0.5 7 21])
end
%legend('Independent Recovery','Kronecker KLT/Global','KCS-KLT','Kronecker Wavelet/Global','Kronecker Fourier/Global','KCS Fourier','KCS Wavelet','Location','SouthEast')
legend('Kronecker KLT/Global','KCS-KLT','Kronecker Wavelet/Global','KCS-Wavelet','Independent Recovery','Location','SouthEast')
% figure(3),clf,plot((0.05:0.05:0.5),timemat(:,2),'LineWidth',2)
% hold on,plot((0.05:0.05:0.5),timemat(:,3),'b-o','LineWidth',2)
% plot((0.05:0.05:0.5),timemat(:,6),'r--','LineWidth',2)
% plot((0.05:0.05:0.5),timemat(:,7),'r--o','LineWidth',2)
% plot((0.05:0.05:0.5),timemat(:,1),'k-.','LineWidth',2)
% axisfortex('','Normalized number of measurements, M/N','Time, seconds')
%axis([0.05 0.5 0 3600])
mean(timemat,1)
%axis tight
%legend('Independent CS','KCS Wavelet','Kronecker Wavelet/Global','KCS Fourier','Kronecker Fourier/Global','KCS KLT','Kronecker KLT/Global','Location','Best')
%legend('KCS Wavelet','Kronecker Wavelet/Global','KCS KLT','Kronecker KLT/Global','Independent CS','Location','Best')
eval(['print -depsc2 -f1 hspec_transformcoding_' num2str(numbands) '.eps'])
eval(['print -dtiff -f1 hspec_transformcoding_' num2str(numbands) '.tif'])
eval(['print -dtiff -f2 hspec_cssnr_pca_' num2str(numbands) '.tif'])
eval(['print -depsc2 -f2 hspec_cssnr_pca_' num2str(numbands) '.eps'])
eval(['print -dtiff -f4 hspec_transformcodingerate_' num2str(numbands) '.tif'])
eval(['print -depsc2 -f4 hspec_transformcodingerate_' num2str(numbands) '.eps'])
% eval(['print -dtiff -f3 hspec_cstime_pca_' num2str(numbands) '.tif'])
% eval(['print -depsc2 -f3 hspec_cstime_pca_' num2str(numbands) '.eps'])
