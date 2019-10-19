function [snr,bgErr,fScore,Mask,Thres] = evelResults(tenD,tenB,tenS,xs1,xs2,sizeD)

snr   = SNR(tenD(:),xs1+xs2);
bgErr = norm(xs1-tenB(:),'fro')/norm(tenB(:),'fro');

S1    = reshape(xs2,sizeD);

nfrm     = sizeD(3);
fScore   = zeros(nfrm,1);
Thres    = zeros(nfrm,1);
Mask     = zeros(sizeD);
for i = 1 : nfrm
   [fScore(i), Mask(:,:,i), Thres(i)] = findFMeasure(S1(:,:,i), tenS(:,:,i));
end

end


