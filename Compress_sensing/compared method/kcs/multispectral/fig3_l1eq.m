% This script generates data for Figure 3 of 
% "Kronecker Compressive Sensing" 
% by Marco F. Duarte and Richard G. Baraniuk
%
% Written by: Marco F. Duarte, Rice University
% Created: February 25 2008

% Size of dim 1
S = 8;
% Size of dims 2-3
N = 8;
% Wavelet filter
h = daubcqf(8);
% Sparsity
K = 10;
% Number of iterations
numiter = 100;

% Psi1 - Tensor product of wavelets
Psi1 = zeros(N^2*S);
for i=1:(N^2*S),
    tmp = zeros(S,N,N);
    tmp(i) = 1;
    Psi1(:,i) = wvlt_hs(tmp,h);
end
Psi1 = Psi1';

% Psi2 - Tensor product of wavelets and Fourier
Psi2 = zeros(N^2*S);
for i=1:(N^2*S),
    tmp = zeros(S,N,N);
    tmp(i) = 1;
    Psi2(:,i) = wfft_hs(tmp,h);
end
Psi2 = Psi2';

% Psi3 - 2D Wavelet sparsity
Psi3 = zeros(N^2);
for i=1:(N^2),
    tmp = zeros(N);
    tmp(i) = 1;
    Psi3(:,i) = reshape(midwt(tmp,h,log2(N)),N^2,1);
end


% For pseudorandom measurements
[a,idxN2] = bitrevorder(ones(N^2,1));
permxN2 = randperm(N^2);
permyN2 = randperm((N^2)-1);
permyN2 = [1 1+permyN2(1:(((N^2))-1))];
[a,idxS] = bitrevorder(ones(S,1));
permxS = randperm(S);
permyS = randperm(S);
[a,idxN2S] = bitrevorder(ones((N^2*S),1));
permxN2S = randperm(N^2*S);
permyN2S = randperm(N^2*S-1);
permyN2S = [1 1+permyN2S(1:(((S*N^2))-1))];


% Phi1 - Band-wise PWT measurements
Phi1 = zeros((N^2*S),N^2*S);
for i=1:((N^2)),
    tmp = zeros((N^2),1);
    tmp(i) = 1;
    tmp2 = At_fw(tmp,permyN2,idxN2,permxN2);
    for j=1:S,
        Phitmp = zeros(S,N,N);
        Phitmp(j,:,:) = reshape(tmp2,[1 N N]);
        Phi1(((j-1)*(N^2))+i,:) = reshape(Phitmp,N^2*S,[])';
    end
end

% Phi2 - Tensor product PWT measurements
Phi2 = zeros((N^2*S),N^2*S);
for i=1:((N^2)),
    tmp = zeros((N^2),1);
    tmp(i) = 1;
    tmp2 = At_fw(tmp,permyN2,idxN2,permxN2);
    for j=1:S,
        Phitmp = zeros(S,N,N);
        Phitmp2 = zeros(S,N,N);
        Phitmp(j,:,:) = reshape(tmp2,[1 N N]);
        for k=1:8,
            for l=1:8,
                Phitmp2(:,k,l) = A_fw(Phitmp(:,k,l),permyS,idxS,permxS);
            end
        end
        Phi2(((j-1)*(N^2))+i,:) = reshape(Phitmp2,N^2*S,[])';
    end
end

% Phi3 - Global measurements
Phi3 = zeros((N^2*S),N^2*S);
for i=1:(N^2*S),
    tmp = zeros((N^2*S),1);
    tmp(i) = 1;
    Phi3(i,:) = At_fw(tmp,permyN2S,idxN2S,permxN2S)';
end

% Phi4 - Random measurements
Phi4 = randn(N^2*S)/(N*sqrt(S));

% Phi5 - 2D measurements
Phi5 = zeros(N^2,N^2);
for i=1:(N^2),
    tmp = zeros(N^2,1);
    tmp(i) = 1;
    Phi5(i,:) = At_fw(tmp,permyN2,idxN2,permxN2)';
end

% Test reconstruction

Mmax = floor((N^2*S)/10);
snrout = zeros(10,Mmax);
success = zeros(10,Mmax);

for iter=1:numiter,
    iter

    idx = randperm(N^2*S);
    x = zeros(N^2*S,1);
    x(idx(1:K)) = randn(K,1);
    y11 = Phi1*Psi1*x;
    y12 = Phi2*Psi1*x;
    y13 = Phi3*Psi1*x;
    y14 = Phi4*Psi1*x;
    y21 = Phi1*Psi2*x;
    y22 = Phi2*Psi2*x;
    y23 = Phi3*Psi2*x;
    y24 = Phi4*Psi2*x;
    yref = Phi4*x;

    for i=Mmax:-1:1,
        M = round(i*10/S);
        idx = 1:M;
        for j=2:S,
            idx = [idx ((j-1)*(N^2))+(1:M)];
        end
        x11 = l1eq_pd(Psi1'*Phi1(idx,:)'*y11(idx), Phi1(idx,:)*Psi1, [], y11(idx), norm(y11(idx))/1000);
        snrout(1,i) = snrout(1,i)+20*log10(norm(x)/norm(x-x11));
        success(1,i) = success(1,i)+(norm(x-x11) <= 1e-3*norm(x));
        x12 = l1eq_pd(Psi1'*Phi2(idx,:)'*y12(idx), Phi2(idx,:)*Psi1, [], y12(idx), norm(y12(idx))/1000);
        snrout(2,i) = snrout(2,i)+20*log10(norm(x)/norm(x-x12));
        success(2,i) = success(2,i)+(norm(x-x12) <= 1e-3*norm(x));
        x13 = l1eq_pd(Psi1'*Phi3(idx,:)'*y13(idx), Phi3(idx,:)*Psi1, [], y13(idx), norm(y13(idx))/1000);
        snrout(3,i) = snrout(3,i)+20*log10(norm(x)/norm(x-x13));
        success(3,i) = success(3,i)+(norm(x-x13) <= 1e-3*norm(x));
        x14 = l1eq_pd(Psi1'*Phi4(idx,:)'*y14(idx), Phi4(idx,:)*Psi1, [], y14(idx), norm(y14(idx))/1000);
        snrout(4,i) = snrout(4,i)+20*log10(norm(x)/norm(x-x14));
        success(4,i) = success(4,i)+(norm(x-x14) <= 1e-3*norm(x));
        xtmp = zeros(S,N,N);
        for j=1:S,
            ysm = y11(idx(((j-1)*M)+(1:M)));
            tmp = l1eq_pd(Psi3'*Phi5(1:M,:)'*ysm, Phi5(1:M,:)*Psi3, [], ysm, norm(ysm)/1000);
            xtmp(j,:,:) = reshape(midwt(reshape(tmp,N,N),h,log2(N)),[1 N N]);
        end
        x15 = wvlt_hs(xtmp,h);
        snrout(5,i) = snrout(5,i)+20*log10(norm(x)/norm(x-x15));
        success(5,i) = success(5,i)+(norm(x-x15) <= 1e-3*norm(x));

        x21 = l1eq_pd(Psi2'*Phi1(idx,:)'*y21(idx), Phi1(idx,:)*Psi2, [], y21(idx), norm(y21(idx))/1000);
        snrout(6,i) = snrout(6,i)+20*log10(norm(x)/norm(x-x21));
        success(6,i) = success(6,i)+(norm(x-x21) <= 1e-3*norm(x));
        x22 = l1eq_pd(Psi2'*Phi2(idx,:)'*y22(idx), Phi2(idx,:)*Psi2, [], y22(idx), norm(y22(idx))/1000);
        snrout(7,i) = snrout(7,i)+20*log10(norm(x)/norm(x-x22));
        success(7,i) = success(7,i)+(norm(x-x22) <= 1e-3*norm(x));
        x23 = l1eq_pd(Psi2'*Phi3(idx,:)'*y23(idx), Phi3(idx,:)*Psi2, [], y23(idx), norm(y23(idx))/1000);
        snrout(8,i) = snrout(8,i)+20*log10(norm(x)/norm(x-x23));
        success(8,i) = success(8,i)+(norm(x-x23) <= 1e-3*norm(x));
        x24 = l1eq_pd(Psi2'*Phi4(idx,:)'*y24(idx), Phi4(idx,:)*Psi2, [], y24(idx), norm(y24(idx))/1000);
        snrout(9,i) = snrout(9,i)+20*log10(norm(x)/norm(x-x24));
        success(9,i) = success(9,i)+(norm(x-x24) <= 1e-3*norm(x));
        tau = max(norm(yref(idx))/sqrt(S*M)/4,sqrt(3)*std(yref(idx)))/4;
        xtmp = zeros(S,N,N);
        for j=1:S,
            ysm = y21(idx(((j-1)*M)+(1:M)));
            tmp = l1eq_pd(Psi3'*Phi5(1:M,:)'*ysm, Phi5(1:M,:)*Psi3, [], ysm, norm(ysm)/1000);
            xtmp(j,:,:) = reshape(midwt(reshape(tmp,N,N),h,log2(N)),[1 N N]);
        end
        x25 = wfft_hs(xtmp,h);
        snrout(10,i) = snrout(10,i)+20*log10(norm(x)/norm(x-x25));
        success(10,i) = success(10,i)+(norm(x-x25) <= 1e-3*norm(x));

        eval(['save fig3_l1eqresults_N' num2str(N) '_S' num2str(S) ' x11 x12 x13 x14 x15 x21 x22 x23 x24 x25 x snrout success i '])
    end
end
snrout = snrout/numiter;
success = success*100/numiter;
eval(['save fig3_l1eqresults_N' num2str(N) '_S' num2str(S) ' x11 x12 x13 x14 x15 x21 x22 x23 x24 x25 x snrout success'])