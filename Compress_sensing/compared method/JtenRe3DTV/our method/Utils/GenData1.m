function [D,S]=GenData(D,W)

Range =max(D(:)) - min(D(:));
Template = rand(W,W)-0.5;
Template = Template*Range;
S=(zeros(size(D)));

[height,width,nfrm]=size(D);
d0=1;
step=10;
i=1;
o=3;
for k=1:nfrm
    if d0+(i-1)*step+W-1<=height
       
%         if o+W-1>width
%             o=1;
%             step=step+4;
%         end     
        dlw=min(o+W-1,width);
        D(d0+(i-1)*step : d0+(i-1)*step+W-1,o : dlw,k)=Template(:,1:dlw-o+1);
        S(d0+(i-1)*step : d0+(i-1)*step+W-1,o : dlw,k)=1;
        i = i+1;
        
    else
        i=1;
        o=o+4;
%         if o+W-1>width
%             o=1;
%             step=step+4;
%         end
        dlw=min(o+W-1,width);
        D(d0+(i-1)*step : d0+(i-1)*step+W-1,o : dlw,k)=Template(:,1:dlw-o+1);
        S(d0+(i-1)*step : d0+(i-1)*step+W-1,o : dlw,k)=1;
        i=i+1;
    end
end
end