%% 提出图片中的主体部分
function P_Image=remove_edge(D)
%  input: 
%         D：原始图像，带有边框
%  output
%         P_Image:经过处理不带边框的图片

ind1=mean(D,1);
column=find(ind1~=255);
tmp=D(:,column);
ind2=mean(D,2);
row=find(ind2~=255);
P_Image=tmp(row,:);