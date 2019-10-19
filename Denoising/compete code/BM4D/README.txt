-------------------------------------------------------------------

  BM4D software for volumetric data denoising and reconstruction
            Public release ver. 2.0 (14 May 2012) 

-------------------------------------------------------------------

Copyright (c) 2011 Tampere University of Technology. 
All rights reserved.
This work should be used for nonprofit purposes only.

Authors:                     Matteo Maggioni
                             Alessandro Foi


BM4D web page:               http://www.cs.tut.fi/~foi/GCF-BM3D


-------------------------------------------------------------------
 Contents
-------------------------------------------------------------------

The package contains these files

*) demo_denoising.m         : denoising demo script
*) demo_reconstruction.m    : reconstruction demo script
*) bm4d.m                   : BM4D volumetric denoising filter [1]
*) sampling.m               : 3-D sampling trajectories generator
*) (i)msfft2.m              : multi-slice 2-D FFT (inverse) transform
*) (i)dct3.m                : 3-D DCT (inverse) transform
*) visualizeXsect.m         : displays phantom cross-sections
*) constantsSparseTraj3D.m  : useful constants used by master scripts
*) ssim_index3d.m           : 3-D SSIM index [3,4]
*) phantom3d.mat            : 3-D Shepp-Logan phantom
*) t1_icbm_normal_1mm_pn0_rf0.rawb   : BrainWeb T1 phantom [2]

-------------------------------------------------------------------
 Installation & Usage
-------------------------------------------------------------------

Unzip BM4D.zip (contains codes) in a folder that is in the MATLAB 
path. Execute the script "demo_reconstruction.m" to run the
reconstruction demo, and execute the script "demo_denoising.m" to
run a volumetric denoising demo. You can freely change the 
parameters involved in the filtering by modifying their initial 
value at the beginning of the master scripts. Comments will help 
you to understand their meaning.

-------------------------------------------------------------------
 Requirements
-------------------------------------------------------------------

*) MS Windows 64 bit, Linux 64 bit or Mac OS X 64 bit
*) Matlab v.7.1 or later with installed:
   -- Image Processing Toolbox (for visualization with "imshow")
   -- Wavelet Toolbox
*) VST framework for Rician-distributed data. Downloadable 
   from http://www.cs.tut.fi/~foi/RiceOptVST/. Required for the
   denoising of Rician data, and for the reconstruction of phantom
   data with non-zero phase.

-------------------------------------------------------------------
 Change log
-------------------------------------------------------------------
v2.0   (17 April 2012)
 + reconstruction of volumetric phantom data with non-zero phase 
   from noisy and incomplete Fourier-domain (k-space) measurements
 + adaptive denoising for data corrupted by spatially varying noise

v1.01  (18 July 2011)
 + fixed few typos, corrected lambda_thr4D in modified profile

v1.0   (17 July 2011)
 + initial version

-------------------------------------------------------------------
 References
-------------------------------------------------------------------

[1] M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, "A Nonlocal 
    Transform-Domain Filter for Volumetric Data Denoising and 
    Reconstruction", submitted to IEEE Trans. Image Process., 2011.

[2] R. Vincent, "Brainweb:  Simulated  brain  database", online at
    http://mouldy.bic.mni.mcgill.ca/brainweb/, 2006.

[3] Z. Wang, A. Bovik, H. Sheikh, E. Simoncelli, "Image quality 
    assessment: from error visibility to structural similarity",
    IEEE Transactions on Image Processing, vol. 13, no. 4, 2004.

[4] J. V. Manjon, P. Coupe, A. Buades, D. L. Collins, M. Robles, 
    "New methods for MRI denoising based on sparseness and 
     self-similarity,” Medical Image Analysis (in press), 2011.
 
-------------------------------------------------------------------
 Disclaimer
-------------------------------------------------------------------

Any unauthorized use of these routines for industrial or profit-
oriented activities is expressively prohibited. By downloading 
and/or using any of these files, you implicitly agree to all the 
terms of the TUT limited license, as specified in the document
Legal_Notice.txt (included in this package) and online at
http://www.cs.tut.fi/~foi/GCF-BM3D/legal_notice.html

-------------------------------------------------------------------
 Feedback
-------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact    Matteo Maggioni   at  matteo.maggioni<at>tut.fi


