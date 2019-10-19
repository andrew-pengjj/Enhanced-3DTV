function x= myPCG_ss(x,X,M2,M3,F,beta,dim)
  temp=beta*X(:)+beta*diffT3(F,dim)+M2(:)-diffT3(M3,dim);                               
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
         y=beta*x+beta*diffT3(diff3(x,dim),dim);
    end
end
