function Idx = vecGrid(row,dep)

if log2(row) < dep
    error('the depth is set too large!\n');
end
        
Idx=[1,fix(row/2),row]; % for dep = 1
q=1;
while q<dep 
      
    tempIdx = zeros(1,(numel(Idx)-1)*2 +1);
    tempIdx(1) = 1;
    
    tempCnt = 2;
    for p = 2 : numel(Idx)
        
        tempIdx(tempCnt) = fix((Idx(p) + Idx(p-1))/2);
        tempIdx(tempCnt+1) = Idx(p);
        tempCnt = tempCnt + 2;
        
    end
    
    Idx = tempIdx;
    
    q = q +1;
end
        
end