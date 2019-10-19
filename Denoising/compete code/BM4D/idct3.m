function [ out ] = idct3( in )
    % assuming input to be a 3-D volume of size NxMxP
    
    M1 = dctmtx(size(in,1));
    M2 = dctmtx(size(in,2));
    M3 = dctmtx(size(in,3));
    
    out = zeros(size(in));
    for i=1:size(in,3)
        out(:,:,i) = M1' * in(:,:,i) * M2;
    end
    
    for i=1:size(in,1)
        for j=1:size(in,2)
            out(i,j,:) = squeeze(out(i,j,:))' * M3;
        end
    end
    
end
