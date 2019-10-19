function x= myPCG1(x,Ysum,M2,F,Gamma,beta,dim)
  temp=beta*Ysum(:)+M2(:)+beta*diffT3(F,dim)-diffT3(Gamma,dim);                                                
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
         y=beta*x+beta*diffT3(diff3(x,dim),dim);
    end
end
