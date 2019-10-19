function  [ ssim ] = SSIM(video0,video1,sizeD)

V0 = reshape(video0,sizeD);
V1 = reshape(video1,sizeD);

nfrm = sizeD(3);
ssim = zeros(nfrm,1);
for k = 1 : nfrm
    
     [mssim1 ssim_map1] = ssim_index(V0(:,:,k),V1(:,:,k));
     ssim(k)            = mssim1;    
end

ssim = mean(ssim);
end
