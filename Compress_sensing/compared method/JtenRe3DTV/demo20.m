clear;clc;close all;
seed = 2015; 
fprintf('Seed = %d\n',seed);
RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
load x_dc
x=trans255(x_dc);
[w,h,s]=size(x);

ratio  = 0.2; 
N      = h*w;
A      = PermuteWHT2(N,s,ratio);

b      = A*x(:);
clear opts;
opts.maxIter = 200;
opts.tol     = 1e-9;
opts.trX     = x;
opts.gauss_frag = 1;
rk           = [ceil(h*0.6),ceil(w*0.6),3];
lam=0.001;
x_rec= LrApprTVFast(A,b,size(x),rk,lam,opts);
xrec_ttv=reshape(x_rec,size(x));  
[everyindex_ttv,outindex_ttv,globalpsnr_ttv ,ergas_ttv, msam_ttv] = MSIQA(x, xrec_ttv);
save ('F:\CS of HSI partitioned\results\dc\0.2\xrec_ttv', 'xrec_ttv');
save ('F:\CS of HSI partitioned\results\dc\0.2\index_ttv',' everyindex_ttv',' outindex_ttv',' globalpsnr_ttv', 'ergas_ttv' ,'msam_ttv');
x_rec_re= LrApprReTV(A,b,size(x),rk,lam,opts);
xrec_rettv=reshape(x_rec_re,size(x)); 
[everyindex_rettv,outindex_rettv,globalpsnr_rettv ,ergas_rettv, msam_rettv] = MSIQA(x, xrec_rettv);
save ('F:\CS of HSI partitioned\results\dc\0.2\xrec_rettv', 'xrec_rettv');
save ('F:\CS of HSI partitioned\results\dc\0.2\index_rettv',' everyindex_rettv',' outindex_rettv',' globalpsnr_rettv', 'ergas_rettv' ,'msam_rettv');

