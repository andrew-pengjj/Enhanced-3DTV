function res=PermuteWHT_partitioned(m,L,ratio)
%INPUTS:
%      m ---------------the length of one signal
%      L -------------- the number of signals
%     samp--------- -sampling index
%     permute-------one random permutation
%OUTPUTS:
%      res---------- one object
%

n      = ceil(m * ratio);
M     = 2^(ceil(log2(m)));
picks = false(M, L);
for i = 1 : L
     p       =  randsample(M, M); 
     idx     =  sort(p(1: n), 'ascend');   idx(1) = 1;
     picks(idx, i)  = 1;   
end
permu = randsample(M,  M);

res.M        = M;
res.m        = m;
res.L        = L;
res.perm     = permu;
res.adjoint  = 0;
res.picks     = picks;

res = class(res, 'PermuteWHT_partitioned');
end