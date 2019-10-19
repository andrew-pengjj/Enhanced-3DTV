function w = dwt3D_vectorized(x, J, af)

% Modified by Mike Wakin from Ivan Selesnick's original function

LLL = x;

sidelength = size(x,1);

w = zeros(sidelength^3,1);

for k = 1:J
  sidelength = sidelength/2;
  [LLL,LLH,LHL,LHH,HLL,HLH,HHL,HHH] = afb3D_vectorized(LLL, af, af, af);
  
  w(sidelength^3+1:8*sidelength^3) = [LLH(:);LHL(:);LHH(:);HLL(:);HLH(:);HHL(:);HHH(:)];
  
end
w(1:sidelength^3) = LLL(:);

