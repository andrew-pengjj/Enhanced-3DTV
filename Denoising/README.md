# Enhanced-3DTV

The motivation map of Enhanced-3DTV is as follows. 

![motivation of Enhanced 3DTV](https://github.com/andrew-pengjj/Enhanced-3DTV/blob/master/Img/3DTVandOurs_v2.pdf)

The matlab code of paper ''Enhanced 3DTV Regularization and Its Applications on Hyper-spectral Image Denoising and Compressed Sensing''
 
  #### denoising
  The structure of Denoising is:
  * compete code
    * ALM_RPCA
    * BM4D
    * LLRT
    * LRTDTV
    * LRTV
    * TDL
    * WNNM
    * WSNM_RPCA
  * Enhanced3DTV in the paper
    * EnhancedTV.m
    * TV_operator
    * quality assess
    * simulation_case1_demo.m
    * simulation_case2_demo.m
    * simulation_case3_demo.m
    * simulation_case4_demo.m
    * simulation_case5_demo.m
    * simulation_case6_demo.m
  * quality assess
  * Demo_simulation_case1.m
  * Demo_simulation_case3.m
  * RunAllMethod.m

API of all methods are list in "RunAllMethod.m"  
```bash
# run "Demo_simulation_case1.m" and "Demo_simulation_case3.m"to test all the code. For example,
$ Demo_simulation_case3.m
# run the code inside "Enhanced3DTV in the paper" to see the performances of Enhanced 3DTV in Denoise tasks. For example,
$ simulation_case3_demo.m
```

```bash
# tips
# There are two methods in EnhancedTV.m to solve the TV inverse problem, 
# one is PCG method and the other is FFT method.
# The performance of two method on simulation case 3 are:
# PCG  time: 141(s) mpsnr: 43.09 mssim: 0.9926
# FFT  time: 35(s)  mpsnr: 43.01 mssim: 0.9926
```
