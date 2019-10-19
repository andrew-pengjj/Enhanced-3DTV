function x=ifastwalsh(data,idx)
% The function implement the 1D sequency(Walsh) ordered fast Walsh-Hadamard transform,
% wchich can be used in signal processing, pattern recongnition and Genetic alogorithms.
% This algorithm uses a Cooley-Tukey type signal flow graph and is implemented in N log2 N
% additions and subtractions. Data sequence length must be an integer power of 2.
%
% The inverse transform is the same as the forward transform except for the multiplication factor N.
% Inversre tranform form can be easily achieved by deleting the last line
% i.e  x=inv(N)*x;
%
% Example:
% x=[1 2 1 1];
% W=FWHT(x);
%
% Author: Gylson Thomas
% e-mail: gylson_thomas@yahoo.com
% Asst. Professor, Electrical and Electronics Engineering Dept.
% MES College of Engineering Kuttippuram,
% Kerala, India, February 2005.
% copyright 2005.
% Reference: N.Ahmed, K.R. Rao, "Orthogonal Transformations for
% Digital Signal Processing" Spring Verlag, New York 1975. page-111.
N = length(data); 
x = data(idx); 
L=log2(N); 
k1=N; 
k2=1;
k3=N/2;
for i1=1:L  %Iteration stage
    L1=1;
    for i2=1:k2
        for i3=1:k3
            i=i3+L1-1; j=i+k3;
            temp1= x(i); temp2 = x(j);
            if(mod(i2,2) == 0)
              x(i) = temp1 - temp2;
              x(j) = temp1 + temp2;
            else
              x(i) = temp1 + temp2;
              x(j) = temp1 - temp2;
            end
        end
            L1=L1+k1;
    end
        k1 = k1/2;  k2 = k2*2;  k3 = k3/2;
end
