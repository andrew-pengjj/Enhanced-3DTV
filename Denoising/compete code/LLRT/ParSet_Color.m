function  [par]=ParSet_Color(nSig)

%% Patch-based Iteration Parameters
par.nSig        =   nSig;                               % Variance of the noise image
par.SearchWin   =   30;                                 % Non-local patch searching window
par.c1          =   1.5*sqrt(2);                        % Constant num for color image
%% Image-based Iteration Parameters
par.Innerloop_X   =   20;                               % InnerLoop Num of estimating clear image
par.kappa    = 0.05;
par.alpha    = 10;    par.rho     = 1;
par.belta    = 1;     par.lambda = 0.8;
%% Patch and noise Parameters
if nSig<=20
    par.patsize       =   6;
    par.patnum        =   200;                          
    par.Iter          =   6;
    par.lamada        =   0.54;     
else
    par.patsize       =   7;
    par.patnum        =   200;
    par.Iter          =   6;
    par.lamada        =   0.56; 
end
par.step      =   floor((par.patsize)/2-1); 
