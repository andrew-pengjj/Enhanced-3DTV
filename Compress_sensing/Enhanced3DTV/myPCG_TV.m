function x   = myPCG_TV(x,T,beta,dim)
  temp=T(:);                                               
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
         y = beta(5)*x+beta(1)*diff_xT(diff_x(x,dim),dim)+beta(2)*diff_yT(diff_y(x,dim),dim);
         y = y + beta(3)*diff_zT(diff_z(x,dim),dim);
    end
end