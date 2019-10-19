function x = funHSI(A,y,sizex,lambda,mu,maxiter)
%% Initialization

 mu1=mu(1); mu2=mu(2) ; mu3=mu(3);
 
 h=sizex(1);w=sizex(2);n=h*w;s=sizex(3);


B1=zeros(n,s); B2=zeros(n,s); B3=zeros(n,s);   % Bregman Variables

% Create Total variation matrices
[Dh,Dv]=opTV2(h,w);                                     % Vertical and Horizontal Variation matrices

x=zeros(n,s);

%% Main iteration
for i=1:maxiter
    
    Q=MySoftTh(Dh*x+B1,lambda/mu1); 
    R=MySoftTh(Dv*x+B2,lambda/mu1);
    S=NucTh(x+B3,mu2/mu3);
   
    

    temp=A'*y;
    tempr=reshape(temp,n,s);
    
    bigY=tempr+Dh'*(mu1*(Q-B1))+Dv'*(mu1*(R-B2))+mu3*(S-B3);  
    

 
    [x(:),~]=lsqr(@afun,bigY(:),1e-6,5,[],[],x(:));
   

    B1=B1+Dh*x-Q;
    B2=B2+Dv*x-R;
    B3=B3+x-S;
    
       if rem(i,10)==0    
           fprintf(' %d iteration done of %d \n',i, maxiter);
       end     
end

function y = afun(z,str)
% 
zmat=reshape(z,n,s);
temp1=reshape(mu1*(Dh'*(Dh*zmat)),n*s,1);
temp2=reshape(mu1*(Dv'*(Dv*zmat)),n*s,1);

tt= (A'*(A*z))+temp1+ temp2 +mu3*z;
        switch str
            case 'transp'
                y = tt;
            case 'notransp'
                y = tt;
        end
end



end


 
%% This is soft thresholding operation
function X= MySoftTh(B,lambda)

   X=sign(B).*max(0,abs(B)-(lambda/2));
end
%% This is nuclear norm thresholding
function X=NucTh(X,lambda)
if isnan(lambda)
    lambda=0;
end
[u,s,v]= svd(X,0);
s1=MySoftTh(diag(s),lambda);
X=u*diag(s1)*v';
end
%% 
% function newZ = SolveLeastSquare(A,Y,Z)
%     [mn,dim]=size(Z);
%     newZ=zeros(mn,dim);
%     for i=1:dim
%         [newZ(:,i),~]=lsqr(A,Y(:,i),1e-6,5,[],[],Z(:,i));                
%                   
%     end   
% end

%% This is a function to make total variation matrix
function [Dh, Dv]=opTV2(h,w)

Dh = spdiags([-ones(w,1) ones(w,1)],[0 1],w,w);
Dh(w,:) = 0;
Dh = kron(Dh,speye(h));
   
Dv = spdiags([-ones(h,1) ones(h,1)],[0 1],h,h);
Dv(h,:) = 0;
Dv = kron(speye(w),Dv);
end
