function G = getGroupSuperColor(I, slicParam)
    
    
    
    N = size(I,1) * size(I, 2);

    g = sparse(N, N);
  
    imlab = single(I);
     
    groupCount = 0;
    for i=1:size(slicParam,1)
                  
        segments = vl_slic(imlab, slicParam(i,1),  slicParam(i,2)) ;
%         
%         I_sp = segImageRegion(uint8(I),segments);
%         I_s = segImage(uint8(I),segments);
%         figure,
%         subplot(121),imshow(I_sp);
%         subplot(122),imshow(I_s);
        
        segments = segments(:);
       
        maxLabel = max(segments);
        for j = 0:maxLabel
            g(:, groupCount+j+1) = segments==j;
        end
        groupCount  = groupCount+ maxLabel+1;
        
        
    end
    
    % croping G
    g =g(:, 1:groupCount);
    G = sparse(3*N, 3.0*double(groupCount)); 
    G(1:N, 1:groupCount) = g;
    G(N+1:2*N, groupCount+1:2*groupCount) = g;
    G(2*N+1:3*N, 2*groupCount+1:3*groupCount) = g;
    
end