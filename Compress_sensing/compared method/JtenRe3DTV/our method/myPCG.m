function x= myPCG(A,y,L,x,p,v1,v2,v3,e,beta,sizex )
  temp=beta(3)*(e+L)+v3-diffT3(v1,sizex)+beta(1)*diffT3(p,sizex)+beta(2)*(A'*y)-A'*(v2);                                                
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,4000,[],[],x);   
    function y = Fun(x)
         y=beta(3)*x+beta(1)*diffT3(diff3(x,sizex),sizex)+beta(2)*(A'*(A*x));
    end
end