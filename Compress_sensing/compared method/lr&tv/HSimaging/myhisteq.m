function img=myhisteq(img)
[~,~,dim]=size(img);

for i=1:dim
    img(:,:,i)=histeq(img(:,:,i));
end
