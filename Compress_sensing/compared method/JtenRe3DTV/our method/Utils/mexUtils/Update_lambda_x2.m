function  [lambda_x2,norm_diff2] = Update_lambda_x2(X2,lambda_x2,zx2,...
                                               G,beta2,gamma,isCont)
% speed up !!!
  

   nfrm =size(X2,2);
   norm_diff2 = 0;
   
   for k1 =1 : nfrm
       
       diagX2 = diag(sparse(X2(:,k1)));
       
       lambda_x2{k1} = lambda_x2{k1} - gamma*beta2*(zx2{k1} -diagX2*G);
       
       if isCont
           norm_diff2 = norm_diff2 + (norm(zx2{k1} - diagX2*G,'fro'))^2;
       end
   end
end