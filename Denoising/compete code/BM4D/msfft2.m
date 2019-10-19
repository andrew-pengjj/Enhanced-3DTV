function [ theta ] = msfft2( phantom )
    theta = zeros(size(phantom));
    for i=1:size(phantom,3)
        theta(:,:,i) = fft2(phantom(:,:,i));
    end
end

