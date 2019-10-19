function rmaxis(n)

% rmaxis(n)
% removes the axis for figure n

figure(n)
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'ZTick',[])