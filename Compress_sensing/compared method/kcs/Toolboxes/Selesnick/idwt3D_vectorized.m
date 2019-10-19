function y = idwt3D_vectorized(w, J, sf)

% Modified by Mike Wakin from Ivan Selesnick's original function

sidelength = round((length(w)^(1/3))/(2^J));
LLL = reshape(w(1:sidelength^3),[sidelength sidelength sidelength]);

for k = J:-1:1
  s = sidelength^3;
  LLH = reshape(w(  s+1:2*s),[sidelength sidelength sidelength]);
  LHL = reshape(w(2*s+1:3*s),[sidelength sidelength sidelength]);
  LHH = reshape(w(3*s+1:4*s),[sidelength sidelength sidelength]);
  HLL = reshape(w(4*s+1:5*s),[sidelength sidelength sidelength]);
  HLH = reshape(w(5*s+1:6*s),[sidelength sidelength sidelength]);
  HHL = reshape(w(6*s+1:7*s),[sidelength sidelength sidelength]);
  HHH = reshape(w(7*s+1:8*s),[sidelength sidelength sidelength]);
  
  LLL = sfb3D_vectorized(LLL,LLH,LHL,LHH,HLL,HLH,HHL,HHH, sf, sf, sf);
  sidelength = sidelength*2;
end

y = LLL;