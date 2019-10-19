% This script generates Figure 4 of 
% "Kronecker Compressive Sensing" 
% by Marco F. Duarte and Richard G. Baraniuk
%
% Written by: Marco F. Duarte, Rice University
% Created: February 25 2008

clear all
load dc
dc=double(dc);
D=dc(500:627,100:227,1:128);
[width,height,band] = size(D);
sizeD               = size(D);
for i=1:band
    D(:,:,i)= D(:,:,i)/max(max( D(:,:,i)));
end

xhs =reshape( D,128*128,128);
N = 128;
S = 128;
if S == 128,
    xh = xhs;
else
    xh = zeros(size(xhs,1),S);
    for i=1:S,
        xh(:,i) = sum(xhs(:,((i-1)*(128/S)+1):(i*(128/S))),2);
    end
end
x = reshape(xh',[S N N]);
xtiles = makehstiles(x);
figure(1),clf,imagesc(xtiles),colormap(gray),axis image,rmaxis(1)
m =1;
for i=1:7,
    m = [m 2^i*ones(2^(i-1)); 2^i*ones(2^(i-1),2^i)];
end
mk = kron(ones(size(xtiles)/N),m);
h = daubcqf(8);
% Approach 1: Wavelet tensor product
wt = wvlt_hs(x,h);
figure(2),clf,imagesc(makehstiles(reshape(abs(wt),[S N N])).*mk),axis image,rmaxis(2)
% Approach 2: Wavelets in space
ws = zeros(size(x));
for i=1:S,
    tmp = reshape(mdwt(squeeze(x(i,:,:)),h,log2(N)),[1 N N]);
    ws(i,:,:) = tmp(1,:,:);
end
figure(3),clf,imagesc(makehstiles(reshape(abs(ws),[S N N])).*mk),axis image,rmaxis(3)
wss = sort(abs(ws(:)),'descend');
% Approach 3: Wavelets in spectra
wp = zeros(size(x));
for i=1:N,
    for j=1:N,
        wp(:,i,j) = mdwt(x(:,i,j),h,log2(S));
    end
end
figure(4),clf,imagesc(makehstiles(reshape(abs(wp),[S N N]))),axis image,rmaxis(4)
wps = sort(abs(wp(:)),'descend');

