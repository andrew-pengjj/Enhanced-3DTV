function [ phantom ] = imsfft2( theta )
    phantom = zeros(size(theta));
    for i=1:size(theta,3)
        phantom(:,:,i) = ifft2(theta(:,:,i));
    end
end
