function cpsnr=myPSNR(org,recon,skip)

org=org(skip+1:end-skip,skip+1:end-skip,:);
recon=recon(skip+1:end-skip,skip+1:end-skip,:);
[m, n,~]=size(org);
    %this is the sum of square error for each band
    sse=sum(sum((org-recon).^2));   
    mse=sse./(m*n);  %mean square error of each band.
    rmse=sqrt(sum(mse)/numel(mse)); 
    maxval=max( max(abs(org(:))),max(abs(recon(:))));
    
    cpsnr=20*log10(maxval/rmse);
end