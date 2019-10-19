function displayVideo(video)

nDim = ndims(video);
% figure,
%  set(gcf,'Position',[5 0 1000, 1000]);
% set(gcf,'position',[200,200,1000,1000]);

if nDim < 4
    nfrm = size(video,3);
    
    for i = 1:nfrm
        imshow(video(:,:,i),[]);
        pause(0.01);
    end
else
    
    nfrm = size(video,4);
    
    for i =1 : nfrm
        imshow(video(:,:,:,i),[]);
        pause(0.01);
    end
end

end
       