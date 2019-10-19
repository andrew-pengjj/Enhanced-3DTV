function x= myPCG_sstv(x,Cha,Mu_x,Mu_y,Mu_z,M1,M2,M3,M4,beta,dim)
  temp=beta*Cha(:)+beta*(diff_xT(Mu_x,dim)+diff_yT(Mu_y,dim)+diff_zT(Mu_z,dim))+...
  M1(:)-diff_xT(M2,dim)-diff_yT(M3,dim)-diff_zT(M4,dim);                                                
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
         y=beta*x+beta*(diff_xT(diff_x(x,dim),dim)+diff_yT(diff_y(x,dim),dim)+diff_zT(diff_z(x,dim),dim));
    end
end
