function [X]         =  ImUpdata(N_Img, Spa_Img, Spa_Wei, nSig, opts )
%% initialization
[m,n,frame] = size(N_Img);
template_time = zeros(2,2,2);    % construct the template
m1 = [1 0;0 0 ];
m2 = [-1 0 ;0 0 ];
template_time(:,:,1) = m1;
template_time(:,:,2) = m2;
otfDz = psf2otf(template_time,[m,n, frame]);

Dz    = zeros(m,n,frame);
Bz    = zeros(m,n,frame);
M     = zeros(m,n,frame);
J1    = zeros(m,n,frame);

%% innerloop iteration
for iter = 1:opts.Innerloop_X
    
    %% X - subproblem
    Denom  =   opts.belta * abs(otfDz).^2 + opts.alpha + opts.kappa;  
    Fx     =   ( opts.kappa * fftn( N_Img ) + conj(otfDz).*fftn(opts.belta * Dz - Bz)+ fftn(opts.alpha * M - J1))./Denom;     % nFFT domain computing
    X      =   real( ifftn( Fx ) );  

    %% M - subproblem
    M      =   (opts.rho*Spa_Img + opts.alpha * X +  J1)./(opts.alpha + opts.rho*Spa_Wei);              % Spatial reconstruction
    J1     =    J1 + opts.alpha*(X - M);
    opts.alpha = opts.alpha*1.02;
 
    %% Dz- subproblem

    Du_z   =   real(ifftn(otfDz.*fftn(X)));                                         % spectral sparsity
    Dz     =   soft_L1_shrink(Bz/opts.belta + Du_z, nSig^2*opts.lambda/opts.belta); % At first, we thought Lp is proper for modeling the spectral gradient. 
    % While considering the Fig. 3(e) and experiment results, we find the L1 norm is slightly better than Lp norm. 
    
    %% updata Bz

    Bz     =    Bz + opts.belta*(Du_z - Dz);
    opts.belta = opts.belta*1.02;

end
end