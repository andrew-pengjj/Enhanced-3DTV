--------------------------------------------------------
What this code do :

This code demonstrates compressed sensing of Hyperspectral  images in the presence of impulse nose.
It solves following optimization problem :
 min_X || Y-AX||_1 + lambda ||Dh*X||_1 + lamdba ||Dv*X||_1 + mu ||X||_*
 X: Hyperspectral image
 A: Sparse binary measurement matrix
 Y: Compressive measurements
 Dh, Dv: Horizontal and vertical finite difference operators
||X||_* : Nuclear norm of matrix X

---------------------------------------------------------------
Toolbox Depandency: 

Our code uses SPOT toolbox by Ewout van den Berg. It can be freely downloaded from here:
Link : https://github.com/mpf/spot

after downloading this code, run "spottests.m" to check that SPOT toolbox is working correctly.

--------------------------------
How to Run this code :

Just run the DemoHSI.m file. It takes around 15 minutes to show the output on 160x160x160  hyperspectral image. 

--------------------------------
File Description :

DemoHSI.m                  : Simply run this file to see how the code works.
funHSI.m                : It is the main function which solves above problem using split-Bregman technique.
HyperSpectralImage.mat  : This is the portion of Washington DC mall image downloaded from here:
                          link: https://engineering.purdue.edu/%7ebiehl/MultiSpec/hyperspectral.html 
--------------------------------------------------------------------
Contact Information:

This code is released just to promote reproducible research and is not very robust.
If you face difficulty in running this code then please feel free to contact us. 

Hemant Kumar Aggarwal( jnu.hemant@gmail.com )  
Snigdha Tariyal (snigdha1491@iiitd.ac.in)



