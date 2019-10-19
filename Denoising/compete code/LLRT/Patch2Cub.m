function  [E_Img,W_Img]  =  Patch2Cub( ImPat, WPat, PatSize, ImageH, ImageW,ImageB )
    TempR        =   ImageH-PatSize+1;
    TempC        =   ImageW-PatSize+1;
    TempOffsetR  =   1:TempR;
    TempOffsetC  =   1:TempC;    

    E_Img  	=  zeros(ImageH,ImageW,ImageB);
    W_Img 	=  zeros(ImageH,ImageW,ImageB);
    k        =   0;
for o = 1:ImageB
    for i  = 1:PatSize
        for j  = 1:PatSize
            k    =  k+1;
            E_Img(TempOffsetR-1+i,TempOffsetC-1+j,o)  =  E_Img(TempOffsetR-1+i,TempOffsetC-1+j,o) + reshape( ImPat(k,:)', [TempR TempC]);
            W_Img(TempOffsetR-1+i,TempOffsetC-1+j,o)  =  W_Img(TempOffsetR-1+i,TempOffsetC-1+j,o) + reshape( WPat(k,:)',  [TempR TempC]);
        end
    end
end

    