
GPSR 5.0  

December 5, 2007

---------------------------------------------------------------------------
Copyright (2007): Mario Figueiredo, Robert Nowak, Stephen Wright

GPSR is distributed under the terms of the GNU General Public License 2.0.

Permission to use, copy, modify, and distribute this software for
any purpose without fee is hereby granted, provided that this entire
notice is included in all copies of any software which is or includes
a copy or modification of this software and in all copies of the
supporting documentation for such software.
This software is being provided "as is", without any express or
implied warranty.  In particular, the authors do not make any
representation or warranty of any kind concerning the merchantability
of this software or its fitness for any particular purpose."
---------------------------------------------------------------------------


This set of MATLAB files contain an implementation of the
gradient projection algorithms described in the paper
"Gradient Projection for Sparse Reconstruction: Application
to Compressed Sensing and Other Inverse Problems"
by Mario A. T. Figueiredo, Robert D. Nowak, Stephen J. Wright,
Journal of Selected Topics in Signal Processing: SPecial Issue
on COnvex Optimization for Signal Processing, December 2007.
 
BOth the paper and the code are available at
http://www.lx.it.pt/~mtf/GPSR/


There are two main files (GPSR_BB.m and GPSR_BAsic.m)  which 
contain two versions of the algorithm, for solving the convex problem
min_theta = 0.5*|| y - A theta ||_2^2 + tau ||theta||_1
as described in the paper.

For usage details, type  "help GSPR_BB" or "help GSPR_Basic"
at the MATLAB prompt.


New in version 3.0: this new version includes a continuation option,
which makes the algorithm work much faster when using small values
of tau. Although it's not obvious to defined what we mean by small,
a simple rule of thumb can be to consider tau small if it is less
tham 0.05 || A^T y  ||_infinity

There are 5 demos included, four of which can be used to 
reproduce most of the figures in the paper. 

The demo "demo_image_deblur.m" requires the presence of 
the Rice wavelet toolbox somewhere on the MATLAB  path. 
This toolbox can be freely downloaded from 
http://www-dsp.rice.edu/software/rwt.shtml

This package also includes an implementation of the IST
algorithm, described in the paper
M. Figueiredo, and R. Nowak, "A Bound Optimization 
Approach to Wavelet-Based Image Deconvolution", IEEE 
International Conference on Image Processing - ICIP'2005, 
Genoa, Italy, September  2005,
which is available at
http://www.lx.it.pt/~mtf/Figueiredo_Nowak_ICIP05.pdf

Attention: the IST algorithm assumes that (and only works if)
the larges singular value of matrix A is no larger than 1.


This package also includes the l1_ls algorithm, which is
described and can be downloaded at
http://www.stanford.edu/~boyd/l1_ls/



Contacts: mario.figueiredo@lx.it.pt
          nowak@ece.wisc.edu
          swright@cs.wisc.edu



This code is in development stage; any comments or bug reports are 
very welcome.
