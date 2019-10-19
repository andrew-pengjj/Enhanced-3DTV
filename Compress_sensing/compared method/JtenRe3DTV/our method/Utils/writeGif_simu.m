function writeGif_simu(origIm,origFg,recFg,recBg,name,imSize)
origIm = reshape(origIm, imSize);
origFg = reshape(origFg, imSize);
recFg  = reshape(recFg, imSize);
recBg  = reshape(recBg, imSize);
recFg = wcodemat(recFg,255);


for i = 1:imSize(3)
    
    figure(1); clf;
    
    subplot(2,2,1);
    imagesc(origIm(:,:,i))%, axis off, colormap gray; 
    title('Orig. Bg','fontsize',12);
    
    subplot(2,2,2);
    imagesc(origFg(:,:,i))%, axis off,colormap gray; 
    title('Orig. Fg','fontsize',12);
    
    subplot(2,2,3);
    imagesc(recBg(:,:,i))%, axis off,colormap gray;
    title('Rec. Bg','fontsize',12);
    
    subplot(2,2,4);
    imagesc((recFg(:,:,i)))%, axis off, colormap gray; 
    title('Rec. Fg','fontsize',12);
    
    pause(0.05);
    
    frame(i)=getframe(gcf); % get the frame
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