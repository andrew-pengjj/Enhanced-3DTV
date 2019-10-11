function x= myPCG1_w(x,Ysum,M2,F,Gamma,weight,beta,dim)
  temp=beta*Ysum(:)+M2(:)-diffT3_weight(Gamma,dim,weight)+beta*diffT3_weight(F,dim,weight);                                                
  [x, ~] = pcg(@(x) Fun(x),temp,1e-4,1000,[],[],x);   
    function y = Fun(x)
         y=beta*x+beta*diffT3_weight(diff3_weight(x,dim,weight),dim,weight);
    end
end
