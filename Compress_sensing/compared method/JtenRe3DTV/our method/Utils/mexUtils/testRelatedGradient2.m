clc,clear;

t1=zeros(10,1);
t2 =zeros(10,1);

nfrm = 20;
m = 5000;
g = 5;
for k = 1:1
    X = randn(m,nfrm);
    temp =rand(m,g);
    G=sparse(temp<0.8);

    lambda = cell(1,nfrm);
    z = cell(1,nfrm);
    for k = 1: nfrm
        lambda{k} =rand(m,g).*G; % sparse(m,g);%
        z{k} = rand(m,g).*G; %sparse(m,g);%
    end

    beta = 1e-3;

    tic;
    [out1] = RelatedGradient2(X,G,lambda,z,beta);
    t1(k)=toc;

    lambda1 =cell2mat(lambda);
    z1 = cell2mat(z);
    tic;
    [out2] = mexRelatedGradient2(X,G,lambda1,z1,beta);
    t2(k)=toc;

    disp(norm(out1(:)-out2(:)));

end
disp(mean(t1)/mean(t2));