function [D,S] = genData(B,W)
% parameters:
% cStp, rStp
% d0 and o
% W
% D    = B;

rStp = W - 6; %speed
cStp = floor(W/5); %3;  %3 4 5

maxD = max(B(:));
minD = min(B(:));
szGrid = 4;
r = 1: szGrid : W;
c = r;

Template = zeros(W,W);
for ii = 1 : numel(r)
    for jj = 1:numel(c)
        
         ridx = r(ii); cidx = c(jj);
         val  = ( (maxD-minD)*rand-minD );
         
         rend = min(ridx+szGrid-1,W);
         cend = min(cidx+szGrid-1,W);
         
         Template(ridx: rend, cidx :cend) = val;
         
    end
end
S=(zeros(size(B)));

[height,width,nfrm]=size(B);

%starting point for col
d0    = 1;
j     = 1;
% cStp  = 4; %moving speed

%starting point for row
o     =  15;
% rStp  = 4;

frag = 1;
for k = 1 : nfrm
    if frag
        minR = min(o+W-1,height);
        B(o:minR,d0+(j-1)*cStp : d0+(j-1)*cStp + W-1,k) = Template(1:minR-o+1,:);
        S(o:minR,d0+(j-1)*cStp : d0+(j-1)*cStp + W-1,k) = 1;
        j = j+1;
        if d0 + (j-1)*cStp + W -1 > width
            frag = 0;
            o = o + rStp;
        end
        
    else
        j = j -1;
        minR = min(o+W-1,height);
        B(o:minR,d0+(j-1)*cStp : d0+(j-1)*cStp + W-1,k) = Template(1:minR-o+1,:);
        S(o:minR,d0+(j-1)*cStp : d0+(j-1)*cStp + W-1,k) = 1;
        if d0 + (j-1)*cStp  < d0 + cStp
            frag = 1;
            o = o + rStp;
            j = 1;
        end
    end
end

D = B;
return;

end
        

% %%% starting point for row
% d0     = 1;
% i      = 1;
% step_r = 3;%moving speed
% 
% %%% starting point for col
% o      = fix(width/3);
% step_c = 3;%moving speed
% frag   = 1;
% 
% for k=1:nfrm
%     
%         
%     if frag
%         
%         dlw=min(o+W-1,width);
%         D(d0+(i-1)*step_r : d0+(i-1)*step_r+W-1,o : dlw,k)=Template(:,1:dlw-o+1);
%         S(d0+(i-1)*step_r : d0+(i-1)*step_r+W-1,o : dlw,k)=1;
%         i = i+1;
%         if d0+(i-1)*step_r+W-1 > height
%             frag = 0; o=o + step_c;
%         end
%         
%     else
%         
%         i=i-1;
%         dlw=min(o+W-1,width);
%         D(d0+(i-1)*step_r : d0+(i-1)*step_r+W-1,o : dlw,k)=Template(:,1:dlw-o+1);
%         S(d0+(i-1)*step_r : d0+(i-1)*step_r+W-1,o : dlw,k)=1;
%         if d0 + (i-1)*step_r < d0 + step_r
%             frag =1; o=o + step_c;
%             i=1;
%         end
%     end
% end
