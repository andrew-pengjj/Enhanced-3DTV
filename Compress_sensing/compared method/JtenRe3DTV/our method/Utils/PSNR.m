function psnr=PSNR(sig0,sig1)

mse=mean((sig0(:)-sig1(:)).^2);

psnr=10*log10(255^2/mse);

end