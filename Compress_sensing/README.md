# Enhanced-3DTV for Compress sensing

The motivation map of Enhanced-3DTV is as follows. 

![motivation of Enhanced 3DTV](https://github.com/andrew-pengjj/Enhanced-3DTV/blob/master/Img/3DTVandOurs_v2.pdf)

The matlab code of paper ''Enhanced 3DTV Regularization and Its Applications on Hyper-spectral Image Denoising and Compressed Sensing''

### structure 
  #### Compress sensing
  The structure of Compress sensing is:
  * compared method 
    * JtenRe3DTV   
    * KSC
    * LR&TV
    * SLNTCS
    * SpaRCS
  * Enhanced3DTV in the paper
    * demo_EnhancedTV_CS.m
    * EnhancedTV_CS.m
    * fdWHtrans.mexw64
    * quality assess
  * quality assess
  * demo.m

```bash
# run "demo.m" in command line window of matlab to test all codes
$ demo.m
# run the code inside "Enhanced3DTV in the paper" to see the performances of Enhanced 3DTV in compress sensing tasks.
$ demo_EnhancedTV_CS.m
```

```bash
# tips:
# downsampling operator involves mixed programming of matlab and C++. 
Here we provide fdWHtrans.mexw64, which can only be run success in windows platform
if the reader want to run in linux and mac platform, 
please compile fdWHtrans.cpp into the corresponding edition. 
The version suffixes of linux and mac are .mexa64 and .mexmaci64 respectively.

The compile process has two step:
please type the following code to command line of matlab.
1): mex -setup C++
2): mex MyMEXCode.cpp
```