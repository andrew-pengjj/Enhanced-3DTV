function tenG = genTenGroupOverlap(height,width,nfrm)
% tenG : N x number of groups
%
%
sizeD = [height,width,nfrm];
N     = height*width*nfrm;
tenG  = sparse(zeros(N,1));
tenG  = diag(tenG);

idx  = 1:N;
tenIdx = reshape(idx, sizeD);
cnt  = 0;
tempIdx = zeros(10,1);
for k = 1 : nfrm -1
    
	%-------------------- The first case ----------------%
	%top left corner
	tempIdx(1) = tenIdx(1,1,k);
	tempIdx(2) = tenIdx(2,1,k);
	tempIdx(3) = tenIdx(1,2,k);
	tempIdx(4) = tenIdx(2,2,k);
	tempIdx(5) = tenIdx(1,1,k+1);

	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%bottom left corner
	tempIdx(1) = tenIdx(height,1,k);
	tempIdx(2) = tenIdx(height-1,1,k);
	tempIdx(3) = tenIdx(height,2,k);
	tempIdx(4) = tenIdx(height-1,2,k);
	tempIdx(5) = tenIdx(height,1,k+1);
	
	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%top right corner
	tempIdx(1) = tenIdx(1,width,k);
	tempIdx(2) = tenIdx(1,width-1,k);
	tempIdx(3) = tenIdx(2,width,k);
	tempIdx(4) = tenIdx(2,width-1,k);
	tempIdx(5) = tenIdx(1,width,k+1);
	
	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%bottom right corner
	tempIdx(1) = tenIdx(height,width,k);
	tempIdx(2) = tenIdx(height,width-1,k);
	tempIdx(3) = tenIdx(height-1,width,k);
	tempIdx(4) = tenIdx(height-1,width-1,k);
	tempIdx(5) = tenIdx(height,width,k+1);
	
	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%-------------The second case -----------------------%
	%top row 
	for j = 2:width-1
	 tempIdx(1) = tenIdx(1,j,k);
	 tempIdx(2)  = tenIdx(1,j-1,k);
	 tempIdx(3) = tenIdx(1,j+1,k);
	 tempIdx(4) = tenIdx(2,j-1,k);
	 tempIdx(5) = tenIdx(2,j,k);
	 tempIdx(6) = tenIdx(2,j+1,k);
	 tempIdx(7) = tenIdx(1,j,k+1);
	
	 cnt = cnt + 1;
	 tenG(cnt,tempIdx(1:7)) = 1;
	end
	
	%bottom row
	for j = 2:width-1
	 tempIdx(1) = tenIdx(height,j,k);
	 tempIdx(2) = tenIdx(height,j-1,k);
	 tempIdx(3) = tenIdx(height,j+1,k);
	 tempIdx(4) = tenIdx(height-1,j-1,k);
	 tempIdx(5) = tenIdx(height-1,j,k);
	 tempIdx(6) = tenIdx(height-1,j+1,k);
	 tempIdx(7) = tenIdx(height,j,k+1);
	
	 cnt = cnt + 1;
	 tenG(cnt,tempIdx(1:7)) = 1;
	end
	
    %left col
	for i = 2:height-1
	 tempIdx(1) = tenIdx(i,1,k);
	 tempIdx(2) = tenIdx(i-1,1,k);
	 tempIdx(3) = tenIdx(i+1,1,k);
	 tempIdx(4) = tenIdx(i-1,2,k);
	 tempIdx(5) = tenIdx(i,2,k);
	 tempIdx(6) = tenIdx(i+1,2,k);
	 tempIdx(7) = tenIdx(i,1,k+1);
	
	 cnt = cnt + 1;
	 tenG(cnt,tempIdx(1:7)) = 1;
	end
	
	%right col
	for i = 2:height -1
	  tempIdx(1) = tenIdx(i,width,k);
	  tempIdx(2) = tenIdx(i-1,width,k);
	  tempIdx(3) = tenIdx(i+1,width,k);
	  tempIdx(4) = tenIdx(i,width-1,k);
	  tempIdx(5) = tenIdx(i-1,width-1,k);
	  tempIdx(6) = tenIdx(i+1,width-1,k);
	  tempIdx(7) = tenIdx(i,width,k+1);
	  
	  cnt = cnt +1;
	  tenG(cnt,tempIdx(1:7)) = 1;
	 end
	 
	%-----------The third case ----------------------%
	for i = 2 : height -1
        for j = 2 : width -1
            
			tempIdx(1) = tenIdx(i,j,k);
			tempIdx(2) = tenIdx(i,j-1,k);
			tempIdx(3) = tenIdx(i,j+1,k);
			tempIdx(4) = tenIdx(i-1,j,k);
			tempIdx(5) = tenIdx(i-1,j-1,k);
			tempIdx(6) = tenIdx(i-1,j+1,k);
			tempIdx(7) = tenIdx(i+1,j,k);
			tempIdx(8) = tenIdx(i+1,j-1,k);
			tempIdx(9) = tenIdx(i+1,j+1,k);
			tempIdx(10)= tenIdx(i,j,k+1);
			
			cnt = cnt +1;
			tenG(cnt,tempIdx(1:10)) = 1;
        end
   end		
	
end


%%% for the last frame
k = nfrm;
%-------------------- The first case ----------------%
	%top left corner
	tempIdx(1) = tenIdx(1,1,k);
	tempIdx(2) = tenIdx(2,1,k);
	tempIdx(3) = tenIdx(2,2,k);
	tempIdx(4) = tenIdx(1,2,k);
	tempIdx(5) = tenIdx(1,1,k-1);

	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%bottom left corner
	tempIdx(1) = tenIdx(height,1,k);
	tempIdx(2) = tenIdx(height-1,1,k);
	tempIdx(3) = tenIdx(height,2,k);
	tempIdx(4) = tenIdx(height-1,2,k);
	tempIdx(5) = tenIdx(height,1,k-1);
	
	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%top right corner
	tempIdx(1) = tenIdx(1,width,k);
	tempIdx(2) = tenIdx(1,width-1,k);
	tempIdx(3) = tenIdx(2,width,k);
	tempIdx(4) = tenIdx(2,width-1,k);
	tempIdx(5) = tenIdx(1,width,k-1);
	
	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%bottom right corner
	tempIdx(1) = tenIdx(height,width,k);
	tempIdx(2) = tenIdx(height,width-1,k);
	tempIdx(3) = tenIdx(height-1,width,k);
	tempIdx(4) = tenIdx(height-1,width-1,k);
	tempIdx(5) = tenIdx(height,width,k-1);
	
	cnt = cnt +1;
	tenG(cnt,tempIdx(1:5)) = 1;
	
	%-------------The second case -----------------------%
	%top row 
	for j = 2:width-1
	 tempIdx(1) = tenIdx(1,j,k);
	 tempIdx(2) = tenIdx(1,j-1,k);
	 tempIdx(3) = tenIdx(1,j+1,k);
	 tempIdx(4) = tenIdx(2,j-1,k);
	 tempIdx(5) = tenIdx(2,j,k);
	 tempIdx(6) = tenIdx(2,j+1,k);
	 tempIdx(7) = tenIdx(1,j,k-1);
	
	 cnt = cnt + 1;
	 tenG(cnt,tempIdx(1:7)) = 1;
	end
	
	%bottom row
	for j = 2:width-1
	 tempIdx(1) = tenIdx(height,j,k);
	 tempIdx(2) = tenIdx(height,j-1,k);
	 tempIdx(3) = tenIdx(height,j+1,k);
	 tempIdx(4) = tenIdx(height-1,j-1,k);
	 tempIdx(5) = tenIdx(height-1,j,k);
	 tempIdx(6) = tenIdx(height-1,j+1,k);
	 tempIdx(7) = tenIdx(height,j,k-1);
	
	 cnt = cnt + 1;
	 tenG(cnt,tempIdx(1:7)) = 1;
	end
	
    %left col
	for i = 2:height-1
	 tempIdx(1) = tenIdx(i,1,k);
	 tempIdx(2) = tenIdx(i-1,1,k);
	 tempIdx(3) = tenIdx(i+1,1,k);
	 tempIdx(4) = tenIdx(i-1,2,k);
	 tempIdx(5) = tenIdx(i,2,k);
	 tempIdx(6) = tenIdx(i+1,2,k);
	 tempIdx(7) = tenIdx(i,1,k-1);
	
	 cnt = cnt + 1;
	 tenG(cnt,tempIdx(1:7)) = 1;
	end
	
	%right col
	for i = 2:height -1
	  tempIdx(1) = tenIdx(i,width,k);
	  tempIdx(2) = tenIdx(i-1,width,k);
	  tempIdx(3) = tenIdx(i+1,width,k);
	  tempIdx(4) = tenIdx(i,width-1,k);
	  tempIdx(5) = tenIdx(i-1,width-1,k);
	  tempIdx(6) = tenIdx(i+1,width-1,k);
	  tempIdx(7) = tenIdx(i,width,k-1);
	  
	  cnt = cnt +1;
	  tenG(cnt,tempIdx(1:7)) = 1;
	 end
	 
	%-----------The third case ----------------------%
	 for i = 2 : height -1
        for j = 2 : width -1
            
			tempIdx(1) = tenIdx(i,j,k);
			tempIdx(2) = tenIdx(i,j-1,k);
			tempIdx(3) = tenIdx(i,j+1,k);
			tempIdx(4) = tenIdx(i-1,j,k);
			tempIdx(5) = tenIdx(i-1,j-1,k);
			tempIdx(6) = tenIdx(i-1,j+1,k);
			tempIdx(7) = tenIdx(i+1,j,k);
			tempIdx(8) = tenIdx(i+1,j-1,k);
			tempIdx(9) = tenIdx(i+1,j+1,k);
			tempIdx(10)= tenIdx(i,j,k-1);
			
			cnt = cnt +1;
			tenG(cnt,tempIdx(1:10)) = 1;
        end
   end		
	
	
tenG(1:cnt,:) = tenG;
	
tenG = tenG';	%tenG : N x cnt
	
return;
end
