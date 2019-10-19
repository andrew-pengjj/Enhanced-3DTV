function [ SigmaX, svp ] = IterativeWSNM( SigmaY, C, p )
% weighted schatten p-norm minimization
% Objective function:
% min |A|^{p}_w,p + |E|_1, s.t. A + E = D
% w_i = C*sqrt(m*n)/(SigmaX_i + eps);

% input: SigmaY is the singular value of Y (maybe partial)
% C = C*sqrt(m*n)/mu;

% here, we use GISA to iteratively solve the following lp-based
% minimization:
% X_{k+1} = argmin_X |X|^{p}_w,p + \mu/2 ||X - (Y - E_{k+1} + \mu^{-1} L_k)||_F
 
p1    =   p;    
Temp  =   SigmaY;
s     =   SigmaY;
s1    =   zeros(size(s));

for i=1:3
   W_Vec    =   C./( (Temp).^(1/p1) + eps );               % Weight vector
   [s1, svp]=   solve_Lp_w(s, W_Vec, p1);
   Temp     =   s1;
end
SigmaX = s1;

end

