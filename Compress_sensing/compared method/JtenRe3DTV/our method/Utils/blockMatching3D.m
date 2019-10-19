function  pos_arr   = blockMatching3D(im, win,step,nblk)

if ndims(im) ~= 3
    error('3-D image required!\n');
end

[h, w, d] = size(im);

S         =   16;
f         =   win;
f2        =   f^2;
s         =   step;

N         =   h-f+1;
M         =   w-f+1;
r         =   [1:s:N];
r         =   [r r(end)+1:N];
c         =   [1:s:M];
c         =   [c c(end)+1:M];
L         =   N*M;
X         =   zeros(f*f,L,d, 'single'); %% f*f x L x d

k    =  0;
for i  = 1:f
    for j  = 1:f
        k    =  k+1;
        blk  =  im(i:end-f+i,j:end-f+j,:);
        blk  = reshape(blk,size(blk,1)*size(blk,2),size(blk,3));
        X(k,:,:) =  blk;
    end
end

% Index image
I     =   (1:L);
I     =   reshape(I, N, M);
N1    =   length(r);
M1    =   length(c);
pos_arr   =  zeros(nblk, N1*M1 );

for  i  =  1 : N1
    for  j  =  1 : M1  
        
        row     =   r(i);
        col     =   c(j);
        off     =  (col-1)*N + row;
        off1    =  (j-1)*N1 + i;
                
        rmin    =   max( row-S, 1 );
        rmax    =   min( row+S, N );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, M );
         
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        B       =   X(:,idx, :); %extracting all pathces in this local window      
        v       =   X(:,off, :); %examplar patch. 
        
        B      =    permute(B,[1,3,2]);
        B      =    reshape(B,size(B,1)*size(B,2),size(B,3));
        v      =    permute(v,[1,3,2]);
        v      =    v(:);
        
        
        dis     =   (B(1,:) - v(1)).^2;
        for k = 2 : f2*d
            dis   =  dis + (B(k,:) - v(k)).^2;
        end
        dis   =  dis./f2;
        [~,ind]   =  sort(dis);        
        pos_arr(:,off1)  =  idx( ind(1:nblk) );        
    end
end