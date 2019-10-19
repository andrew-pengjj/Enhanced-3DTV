function  writeGif_real(name,D,Bg,Fg,imSize)
D  = reshape(D,imSize);
Bg = reshape(Bg,imSize);
Fg = reshape(Fg,imSize);
% Fg = wcodemat(Fg,255);

nfrm = imSize(3);
% frame =struct;
for k = 1:nfrm

    figure(1),
%     imshow(Fg(:,:,k),[]); 
%     subplot(131),imshow(D(:,:,k),[]); title('Orig. image');
%     subplot(132),imshow(Bg(:,:,k),[]); title('Rec. Bg');
%     subplot(133),imshow(Fg(:,:,k),[]); title('Rec. Fg');
    subplot(2,2,1);imshow(Bg(:,:,k)/255);title('背景');subplot(2,2,2);imshow(abs(Fg(:,:,k))/255);title('前景');
    subplot(2,2,3);imshow(Bg(:,:,k),[]);title('放大对比度背景');subplot(2,2,4);imshow(Fg(:,:,k),[]);title('放大对比度前景');


    
    frame(k)=getframe(gcf); % get the frame
end

dt = 0.1;
writegif2(name,frame,dt);
close all;
end


function writegif2(name,frames,dt)
    
nframe=length(frames);
for i=1:nframe
    [image,map]=frame2im(frames(i));
    [im,map2]=rgb2ind(image,128);
    if i==1
        imwrite(im,map2,name,'gif','writeMode','overwrite','delaytime',dt,'loopcount',inf);
    else
        imwrite(im,map2,name,'writeMode','append','delaytime',dt); %,'loopcount',inf);
    end
end
end