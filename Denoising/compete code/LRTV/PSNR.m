function result = PSNR(ReferBuffer,UnReferBuffer,lHeight,lWidth)

    result = 0;
	for j=1:lWidth*lHeight
		temp = ReferBuffer(j)-UnReferBuffer(j);
		result = result + double(temp*temp);
    end;
    
	if(result==0)
		result =100;
	else 
		result = 10*log10(255*255/result*lWidth*lHeight);
    end;