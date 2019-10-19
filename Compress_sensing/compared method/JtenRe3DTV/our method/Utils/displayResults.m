function displayResults(Orig,xs1,xs2,height,width,nfrm)
% INPUTS:
% Orig : original video.
% xs1 : background.
% xs2 : foreground.
% height, width, nfrm.
%
Orig = reshape(Orig,height,width,nfrm);
a1=reshape(xs1,height,width,nfrm);
a2=reshape(xs2,height,width,nfrm);
% a1=wcodemat(a1,255);
% a2=wcodemat(a2,255);
% figure,
% for i=1:nfrm
%    subplot(131), imagesc(Orig(:,:,i));% imshow(Orig(:,:,i),[]);
%    subplot(132), imagesc(a1(:,:,i));% imshow(a1(:,:,i),[]);
%    subplot(133), imagesc(a2(:,:,i));% imshow(a2(:,:,i),[]);
%    pause(0.05);
% end
for i=1:nfrm
   subplot(131), imshow(Orig(:,:,i),[]);
   subplot(132), imshow(a1(:,:,i),[]);
   subplot(133), imshow(a2(:,:,i),[]);
   pause(0.05);
end
end