clear all 
close all

[af, sf] = farras;              % analysis and synthesis filter
x=zeros(1,64);                  % create generic signal whose values are all 0
w = dwt(x,3,af);                % analysis filter banks (3 stages)
w{3}(4)=1;                      % setting one of the wavelet coefficients to be 1
y = idwt(w,3,sf);               % synthesis filter banks (3 stages)
plot(y);            % Plot the wave
title('Standard 1-D wavelet transform') 
xlabel('t');                    
ylabel('\psi(t)');



