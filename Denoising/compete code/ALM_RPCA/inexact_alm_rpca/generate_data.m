clear all
clc
load('urbanpart.mat')
Data=urbanpart(51:250,51:250,1:60);
ind=randperm(60,10);
Data(:,70:10:120,ind)=0;
for i=1:60
    O(:,i)=reshape(Data(:,:,i),200*200,1)/max(max(Data(:,:,i)));%干净数据
end
ind=randperm(60,10);
Y=O+normrnd(0,0.04,size(O));% 加Guassian噪音
