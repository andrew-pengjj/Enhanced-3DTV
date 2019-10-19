 function [tempx2] = RelatedGradient2(X2,G,lambda_x2,zx2,beta2)
 
  [m,n]= size(X2);
  tempx2 = zeros(m,n);
  for nfrm = 1:n %% please speed up!!!
%      temp = G.*(lambda_x2{nfrm} - beta2*zx2{nfrm}) + beta2*diag(sparse(X2(:,nfrm)))*G
     tempx2(:,nfrm) = sum(G.*(lambda_x2{nfrm} - beta2*zx2{nfrm}),2) +...
                                beta2*sum(diag(sparse(X2(:,nfrm)))*G, 2);
                      
  end;