function clean_HSI = RPCA_WSNM(data, N_Img, nSig, p, C, lamb, it, par)

    %% Rearrange subcube to 2D matrix D;
    [X num] = Cube2Matrices(data, par);
    [N_Y numy] = Cube2Matrices(N_Img, par);
    [w h nump] = size(X);
    Y = reshape(X, [w*h nump]);
    N_Y = reshape(N_Y, [w*h nump]);
    %% Update \eta
    SigmaArr = par.gama*sqrt(abs(repmat(nSig^2,1,size(Y,2))-mean((N_Y-Y).^2)));
    if it==1
        SigmaArr = nSig*ones(1, num);
    end
    Low_Rank = zeros(size(X));
    lambda = lamb / sqrt(size(Low_Rank,1));
    tol = 1e-7;
    maxIter = 1000;
    %% Estimate A
    for i = 1:num
        [Low_Rank(:,:,i),S,iter] = inexact_alm_rpca_wsnm(X(:,:,i), C, lambda, SigmaArr(i), p, tol, maxIter);
    end
    clean_HSI = Matrices2Cube(Low_Rank, data, par);
        
end