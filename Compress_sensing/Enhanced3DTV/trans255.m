function [X_out ] = trans255( X_in)
band=size(X_in,3);
X_out=zeros(size(X_in));
for i=1:band
    X_out(:,:,i)=(X_in(:,:,i)-min(min(X_in(:,:,i))))/(max(max(X_in(:,:,i)))-min(min(X_in(:,:,i))))*255;
end



end

