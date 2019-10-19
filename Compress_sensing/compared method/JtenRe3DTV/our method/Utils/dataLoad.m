function videoOut = dataLoad(inputFolder,param,varargin)
%frag = 1 indicates reading avi vedio
%frag = 0 indicates reading image sequences.

if nargin ==3 
% switching variable for video or image seqences.
    frag       = 1; 
    fileName   = varargin{1}; 
else
    frag       = 0;
    partFName  = param.partFName; 
    fileFormat = param.fileFormat;
end; 

scale    = param.scale;
grayFrag = param.grayFrag;
startIdx = param.startIdx;
endIdx   = param.endIdx;
nfrm     = endIdx - startIdx+1;

if frag % read vedio

    fullName = [inputFolder,fileName];
    
    if grayFrag
    
        videoOut =aviReadGray(fullName,startIdx,endIdx,nfrm,scale);
    else
        
        videoOut = aviReadColor(fullName,startIdx,endIdx,nfrm,scale);
    end
    
else % read image sequensces
    
    imgPath=strcat([inputFolder,partFName,num2str(startIdx),'.',fileFormat]);
    im = imread(imgPath);
    
    [height,width, channel] = size(imresize(im,scale));
    
    if channel == 1
        
        videoOut = zeros(height,width,nfrm);
        idx = 1;
        for i = startIdx : endIdx
            
            imgPath=strcat([inputFolder,partFName,num2str(i),'.',fileFormat]);
            im = imread(imgPath);
            
            videoOut(:,:,idx) = imresize(im, scale);
            
            idx = idx+1;
        end
        
    else
      %%% processing multiple channels
     if grayFrag

        videoOut = zeros(height,width,nfrm);
        idx =1;
        for i = startIdx : endIdx
            
            imgPath=strcat([inputFolder,partFName,num2str(i),'.',fileFormat]);
            im = imread(imgPath);

            videoOut(:,:,idx) = imresize( rgb2gray(im) , scale);
            
            idx = idx+1;
        end
    else
        
        videoOut = zeros(height,width,channel,nfrm);
        idx =1;
        for i = startIdx : endIdx
            
            imgPath = strcat([inputFolder,partFName,num2str(i),'.',fileFormat]); 
            im = imread(imgPath);
            
            videoOut(:,:,:,idx) = imresize(im, scale);
            
            idx = idx +1;
        end  
     end
        
    end% END IF
    
end% END IF

return;
       
end



function Gray=aviReadGray(vfile,startIdx,endIdx,nfrm,scale)
% 
%
Mobj = VideoReader(vfile);
nFrames = Mobj.NumberofFrames;
height = Mobj.Height;
width = Mobj.Width;

if nfrm < nFrames
    
    temp  =read(Mobj,1);
    [height,width] = size(imresize( rgb2gray(temp), scale));
    Gray=zeros(height,width,nfrm);
    idx = 1;
    for k=startIdx : endIdx
        temp=read(Mobj,k);
        Gray(:,:,idx)=imresize( rgb2gray(temp), scale);
        idx = idx +1;
    end
    
else
    error('Excced the number of frames in this video.');
end

end

function Color=aviReadColor(vfile,startIdx,endIdx,nfrm,scale)
% 
%
Mobj=VideoReader(vfile);
nFrames=Mobj.NumberofFrames;
height=Mobj.Height;
width=Mobj.Width;

if nfrm < nFrames
   
    temp  =read(Mobj,1);
    [height,width,channel] = size(imresize(temp, scale));
    Color=zeros(height,width,channel,nfrm);
    idx =1;
    for k=startIdx:endIdx
        
        Color(:,:,:,idx)= imresize( read(Mobj,k), scale );
        idx = idx +1;
    end
    
else
    error('Excced the number of frames in this video.');
end

end


    