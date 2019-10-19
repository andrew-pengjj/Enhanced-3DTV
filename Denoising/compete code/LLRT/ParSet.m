function  [par]=ParSet(nSig)

%% Patch-based Iteration Parameters
par.nSig        =   nSig;                               % Variance of the noise image
par.SearchWin   =   30;                                 % Non-local patch searching window
par.c1          =   10*sqrt(2);                         % Constant num for HSI
%% Image-based Iteration Parameters
par.Innerloop_X   =   20;                               % InnerLoop Num of estimating clear image
par.kappa    = 0.05;
par.alpha    = 10;    par.rho     = 1;
par.belta    = 1;     par.lambda = 0.8;
par.step     = 5;
%% Patch and noise Parameters
if nSig<=10/255
    par.patsize       =   6;
    par.patnum        =   300;                          % Increase the patch number and iterations could further improve the performance, at the cost of running time.
    par.Iter          =   4;
    par.lamada        =   0.54;     
elseif nSig <= 30/255
    par.patsize       =   7;
    par.patnum        =   300;
    par.Iter          =   4;
    par.lamada        =   0.56; 
else
    par.patsize       =   8;
    par.patnum        =   300;
    par.Iter          =   4;
    par.lamada        =   0.58; 
end

