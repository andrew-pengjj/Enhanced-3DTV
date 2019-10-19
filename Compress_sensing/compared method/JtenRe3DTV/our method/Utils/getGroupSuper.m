function G = getGroupSuper(I, slicParam)
    

    N = size(I,1) * size(I, 2);

    G = sparse(N,N);
 
    %imlab = vl_xyz2lab(vl_rgb2xyz(I)) ;
    imlab = single(I);
    %slicParam = param.superpixel.slicParam;
     
    groupCount = 0;
    for i=1:size(slicParam,1)
                  
        segments = vl_slic(imlab, slicParam(i,1),  slicParam(i,2)) ;
       
        %I_sp = segImageRegion(I,segments);
        %I_s = segImage(I,segments);
        %figure,
        %subplot(121),imshow(I_sp);
        %subplot(122),imshow(I_s);
        
        segments = segments(:);
       
        maxLabel = max(segments);
        for j = 0:maxLabel
            G(:, groupCount+j+1) = segments==j;
        end
        groupCount  = groupCount+ maxLabel+1;
        
        
    end
    
    % croping G
    G =G(:, 1:groupCount);
     
%     G = G';
end