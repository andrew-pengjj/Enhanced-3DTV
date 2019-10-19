function g =  getGroupOverlapWithHierarchical(row, col,rdepth,cdepth)
%************************************************************************
% Decription: this script aims to generate a overlap-group of tree
% struture.
%-------------------------------------------------------------------------
%INPUTS:
% - row, col: the size of row and col.
% - depth : the depth of tree. 
%OUTPUTS:
% - g : one matrix extrating the overlapped group of tree structure.
%wenfeic2006@163.com
%**************************************************************************

    if floor(log2(row))<rdepth || floor(log2(col)) < cdepth
        error('The depth is set too large!\n');
    end
    
    N = row*col;
    
    g = sparse(zeros(N,1));
    g = diag(g);
    
    % the root in the tree (depth =0)
    g(1,:) = ones(1,N);
    
    % the first level in the tree (depth = 1,2,3,...)
   if rdepth < cdepth;

        cnt = 1;
        for k = 1 : rdepth
            
          curr_cnt = cnt;
          [tempg,nGp] = extractGp(row,col,k,k);
          cnt = cnt + nGp;
          g(curr_cnt+1:cnt,:) = tempg;
          
        end
        
        for k = rdepth+1:cdepth
            
            curr_cnt = cnt;
            [tempg, nGp] =extractGp(row,col,rdepth,k);
            cnt = cnt + nGp;
            g(curr_cnt+1:cnt,:) = tempg;
            
        end
   else
        cnt = 1;
        for k = 1 : cdepth
            
          curr_cnt = cnt;
          [tempg,nGp] = extractGp(row,col,k,k);
          cnt = cnt + nGp;
          g(curr_cnt+1:cnt,:) = tempg;
          
        end
        
        if cdepth < rdepth
            for k = cdepth+1:rdepth

                curr_cnt = cnt;
                [tempg, nGp] =extractGp(row,col,k,cdepth);
                cnt = cnt + nGp;
                g(curr_cnt+1:cnt,:) = tempg;

            end
        end
        
   end
         
    
    if cnt < N
        g = g(1:cnt,:);
    end
   
    g = g';
    return;
    
end

function [tempg, nGp]= extractGp(row, col, kr, kc)

[rIdx] = vecGrid(row, kr);
[cIdx] = vecGrid(col, kc);

N = row*col;
nGp = (numel(rIdx)-1) * (numel(cIdx)-1);
tempg = sparse(nGp,N,0);
allIdx = reshape([1:N],row,col);
vecZero = zeros(1,N);
cnt = 1;
for ii = 1 : numel(rIdx) - 1
   
    for jj = 1 : numel(cIdx) - 1
        
        if ii == numel(rIdx) -1 && jj == numel(cIdx) - 1
            
            cellIdx = allIdx(rIdx(ii) : rIdx(ii+1), cIdx(jj) : cIdx(jj+1) );
        elseif ii == numel(rIdx) -1 && jj ~= numel(cIdx) -1
            
            cellIdx = allIdx(rIdx(ii) : rIdx(ii+1), cIdx(jj) : (cIdx(jj+1)-1) );
        elseif ii ~= numel(rIdx) -1 && jj == numel(cIdx) -1
            
            cellIdx = allIdx(rIdx(ii) : (rIdx(ii+1)-1), cIdx(jj) : cIdx(jj+1));
        else
            
            cellIdx = allIdx( rIdx(ii) : (rIdx(ii+1)-1), cIdx(jj) : (cIdx(jj+1) -1) );
        end
        
        vecZero(cellIdx) = 1;
        figure(1),imshow(reshape(vecZero,row,col),[]);
             
        tempg(cnt,:) = sparse(vecZero);
        
        vecZero = zeros(1,N);
        
        cnt = cnt +1;
    end
    
end

function Idx = vecGrid(row,dep)
        
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

end


        
        
        
        
        





    
    
    
    