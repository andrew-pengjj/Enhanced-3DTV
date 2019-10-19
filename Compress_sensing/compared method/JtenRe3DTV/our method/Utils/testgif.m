%Ê¾Àý£º
clear all;
clc;
set(gca,'nextplot','replacechildren','box','off','color','b','xgrid','on');
title('Get 20 frames of current window');
%%
for j=1:20
    
    plot(fft(eye(j+16)));
    axis([-1. 1. -1. 1.]);
    frame(j)=getframe(gcf); % get the frame
end
writegif('test.gif',frame,0.1);