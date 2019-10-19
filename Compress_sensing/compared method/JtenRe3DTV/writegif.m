function  writegif(name,D,imSize)
D  = reshape(D,imSize);
nfrm = imSize(3);
% frame =struct;
for k = 1:nfrm
    figure(1)
   imshow(D(:,:,k),[]);title('kcs');
   set(gcf,'color','w');
    set(gca,'units','pixels','Visible','off');
    q=get(gca,'position');
    q(1)=0;
    q(2)=0;
    set(gca,'position',q);  
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
    [im,map2]=rgb2ind(image,256);
    if i==1
        imwrite(im,map2,name,'gif','writeMode','overwrite','delaytime',dt,'loopcount',inf);
    else
        imwrite(im,map2,name,'writeMode','append','delaytime',dt); %,'loopcount',inf);
    end
end
end