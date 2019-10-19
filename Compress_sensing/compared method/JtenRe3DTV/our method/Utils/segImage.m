function I_s = segImage(I,S)

    S = single(S);
    [cx,cy] = gradient(S);
    ccc = (abs(cx)+abs(cy))~=0;
    ccc = uint8(ccc)*255;
%     I_s = I;
    [m,n] = size(I);
    I_s = uint8(zeros(m,n,3));
    I_s(:,:,1) = I;
    I_s(:,:,2) = I;
    I_s(:,:,3) = I;
    I_s(:,:,1) = max(I_s(:,:,1),ccc);
    I_s(:,:,2) = min(I_s(:,:,2),255-ccc);
    I_s(:,:,3) = min(I_s(:,:,3),255-ccc);
