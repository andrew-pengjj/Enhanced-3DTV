clc,clear;
t1 = zeros(10,1);
t2 = zeros(10,1);

m = 10000;
n = 100;
g =500;
for iter = 1:1
    
X2 =rand(m,n);
temp=rand(m,g);
G=sparse((temp<0.3));
lambda_x2=cell(1,n);
for k =1:n
%     lambda_x2{k}=randn(size(G)).*G;
    lambda_x2{k}= sparse(size(G,1),size(G,2));
end
zx2=lambda_x2;
beta2=1e1;
gamma =1e1;
isCont =1;

tic;
[lambda_x1,norm_diff1] = Update_lambda_x2(X2,lambda_x2,zx2,...
                                               G,beta2,gamma,isCont);
t1(iter)=toc;
                                           
tic;
[lambda_x2,norm_diff2] = mexUpdate_lambda_x2(X2,lambda_x2,zx2,...
                                                 G,beta2,gamma,isCont); %
t2(iter)=toc;

a1=cell2mat(lambda_x1);
a2=cell2mat(lambda_x2);

disp(norm(a1(:)-a2(:)));
disp(norm_diff1-norm_diff2);
end
disp(mean(t1)/mean(t2));