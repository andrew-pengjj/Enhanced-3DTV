function [SigmaX,svp]=ClosedWNNM(SigmaY,C,oureps)
%加权核范数最小化问题的解析解
%这里的目标函数如下
%         sum(w*SigmaY)+1/2*||Y-X||_F^2
% 其中w_i =C/(sigmaX_i+oureps),oureps是一个足够小的常数
%
temp=(SigmaY-oureps).^2-4*(C-oureps*SigmaY);
ind=find (temp>0);
svp=length(ind);
SigmaX=max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;

end