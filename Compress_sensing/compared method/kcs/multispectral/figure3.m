% This script generates Figures 5 and 6 of 
% "Kronecker Compressive Sensing" 
% by Marco F. Duarte and Richard G. Baraniuk
%
% This script needs the file fig3_l1eqresults_N8_S8.mat, 
% which is generated using fig3_l1eq.m
%
% Written by: Marco F. Duarte, Rice University
% Created: February 25 2008

close all
load fig3_l1eqresults_N8_S8
Mmax = size(snrout,2);
figure(1),clf,plot(10*(1:Mmax),snrout(4,:)','k--','LineWidth',2)
hold on,plot(10*(1:Mmax),snrout(1,:)','LineWidth',2)
hold on,plot(10*(1:Mmax),snrout(5,:)','r-.','LineWidth',2)
axisfortex('','Number of measurements M','SNR, dB')
legend('Global Measurements','KCS Recovery','Independent Recovery','Location','SouthEast')
axis([0 512 0 125])
figure(2),clf,plot(10*(1:Mmax),success(4,:)','k--','LineWidth',2)
hold on,plot(10*(1:Mmax),success(1,:)','LineWidth',2)
hold on,plot(10*(1:Mmax),success(5,:)','r-.','LineWidth',2)
axisfortex('','Number of measurements M',['Probability of successful recovery, ' 37])
legend('Global Measurements','KCS Recovery','Independent Recovery','Location','SouthEast')
axis([0 512 0 101])
print -depsc2 -f2 tensorprob.eps
print -dtiff -f2 tensorprob.tif
