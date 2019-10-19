% This package contains a matlab implementation of the LRMR  algorithm [1].
%
% Pleae run the "Simulated_demo.m" file,
% which will run our LRMR algorithm on the simulated indian data [2].
%
%
%[1]H. Zhang, W. He, L. Zhang, H. Shen, and  Q. Yuan, ¡°Hyperspectral Image 
%   Restoration Using Low-Rank Matrix Recovery,¡±
%   IEEE Trans. Geosci. Remote Sens. , vol. 52, pp. 4729-4743, Aug. 2014.
%[2] W. He, H. Zhang, L. Zhang, and  H. Shen, ¡°Total-Variation-Regularized 
%    Low-Rank Matrix Factorization for Hyperspectral Image Restoration,¡± 
%    IEEE Trans. Geosci. Remote Sens. , vol. 54, pp. 178-188, Jan. 2016.
%[3] W. He, H. Zhang, L. Zhang, H. Shen, "Hyperspectral Image Denoising 
%    via Noise-Adjusted Iterative Low-Rank Matrix Approximation", IEEE 
%    Journal of Selected Topics in Applied Earth Observations and Remote
%    Sensing, vol. 8, no. 6, pp. 3050 - 3061, 2015.
%
% NOTE:  LRMR is implemented with ssGoDec, which is different from the one 
%        mentioned in the original paper [1]. The selection of parameter s is 
%        changed correspondingly.

%
% Author: Wei He (November, 2014)
% E-mail addresses:(weihe1990@whu.edu.cn)
