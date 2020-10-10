function x= myPCG_CS(x,A,T,beta)
  temp =T(:);                                              
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
         y = beta(4)*(A'*(A*x));
         y = y(:) + beta(5)*x;
    end
end