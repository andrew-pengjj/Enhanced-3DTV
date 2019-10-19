clc,clear;

t1=zeros(10,1);
t2 =zeros(10,1);

nfrm = 200;
for k =1:1
X = randn(5000,nfrm);
temp = rand(5000,300);
G = sparse(temp<0.5);
lambda = cell(1,nfrm);
for i =1:nfrm
    
  lambda{i} = sparse(5000,300);
%   lambda{i} = rand(5000,300).*G;
  
end
beta = 1e-1;
tic;
[z1] = frameGroupShrinkage(X,G,lambda,beta);
t1(k)=toc;
Z1 = cell2mat(z1);

lambda1= cell2mat(lambda);
tic;
Z2 = mexFrameGroupShrinkage(X,G,lambda1,beta);
t2(2)=toc;

disp(norm(Z1-Z2,'fro'));
end

disp(mean(t1)/mean(t2));