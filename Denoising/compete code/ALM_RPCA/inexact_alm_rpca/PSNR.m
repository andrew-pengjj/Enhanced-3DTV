%% PSNR 
function [mse,psnr]=PSNR(O,C)
%    input 
%          O:原图
%          C：待比较的图片
%    output
%          mse:均方差
%          psnr:峰值信噪比
mse=sum(sum(abs(O-C).*abs(O-C)))/(size(O,1)*size(O,2));
psnr=20*log10(1/sqrt(mse));