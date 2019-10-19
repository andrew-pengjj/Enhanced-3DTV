The following files contain a MATLAB implementation for the SpaRCS algorithm developed in:

A. E. Waters, A. C. Sankaranarayanan, and R. G. Baraniuk. SpaRCS: recovering low rank and sparse matrices from compressive measurements. In Neural Information Processing Symposium (NIPS), Granada, Spain, Dec. 2011.


Some notes on the contents of the directory:

1) demo.m is a simple file that can be used to test SpaRCS.

2) setup.m is a simple file to install the appropriate noiselet .mex files.  This step is only required if you plan to use noiselets with the SpaRCS code. Furthermore, pre-installed .mex files for 32-bit and 64-bit Windows systems are included.  

2) ./algs/ is a folder that contains the main code for SpaRCS.  Future releases will also contain code for the CS-SVT and CS-APG algorithms described in the paper.

3) ./operators/ contains operator functions for identity, noiselets, and random sampling operators. These can serve as templates for designing new operators for use with SpaRCS.

4) ./PROPACK/ contains the Propack numerical library, which can be used for efficient SVD computation in SpaRCS.

5) ./utility/ contains a variety of helper routines for least squares reconstruction and computation of noiselets.