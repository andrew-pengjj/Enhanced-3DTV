function [dx] = soft_L1_shrink(dx,Thresh)
%SHRINK2   Vectorial shrinkage (soft-threholding)
%   [dxnew,dynew] = SHRINK2(dx,dy,Thresh)

s = abs( dx );
% s = max(s - Thresh,0)./max(1e-12,s);
% dx = s.*dx;
dx = max(s - Thresh,0).*sign(dx);
