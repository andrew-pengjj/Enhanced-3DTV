Outline:

    Folder structure
    Usage
    Citation
    URL of third-party toolboxes and functions


Folder structure:

    tensor_dl\      : major part of the TensorDL method
        |-- dependency\     : dependent toolbox
        |       |-- fkmeans\        : k-means++ [1]
        |       |-- tensor_toolbox\ : toolbox for tensor operations [2]
        |       |-- poblano_toolbox\: dependency of tensor_toolbox [3]
        |-- utils\          : helper functions of TensorDL
        |-- TensorDL.m      : core function of TensorDL
    invansc\        : VST tool for pre- and post-processing MSI contaminated by Poisson noise [4]
    quality_assess\ : functions of quality assessment indices
        |-- ssim_index.m            : SSIM [5]
        |-- FeatureSIM.m            : FSIM [6]
        |-- ErrRelGlobAdimSyn.m     : ERGAS
        |-- SpectAngMapper.m        : SAM
        |-- MSIQA.m                 : interface for calculating PSNR and the four indices above
    tdl_demo\       : comparison and GUI demo
        |-- compete_methods : competing methods (down loaded or implemented based on reference papers)
        |       |-- ksvdbox\        : K-SVD toolbox [7]
        |       |-- ompbox\         : dependency of toolbox ksvdbox [8]
        |       |-- naonlm3d\       : 3D NLM toolbox [9]
        |       |-- BM3D\           : BM3D toolbox [10]
        |       |-- BM4D\           : BM4D toolbox [11]
        |       |-- NLM3D.m         : 3D NLM using toolbox naonlm3d
        |       |-- LRTA.m          : LRTA method
        |       |-- PARAFAC.m       : PARAFAC method
        |-- ui_utils\       : scripts used in GUI
        |-- Comparison.m	: function that applies the methods and calculates the QA indices
        |-- TDLdemo.m       : script that display the results
        |-- pompoms.mat     : an example MSI data (512 x 512 x 31)
        |-- results.mat     : results from Comparison.m (to be used by TDLdemo.m)


Usage:
    
    For MSI with Gaussian noise, you can simply use the function TensorDL as follows:
        [ clean_img, basis, idx, cores, info ] = TensorDL( noisy_img, params )
    Please type 'help TensorDL' to get more information.

    For MSI with Poisson noise or mixed Gaussian and Poisson noise, you can follow these steps:
      1. Apply VST via function GenAnscombe_forward with estimated Gaussian STD \sigma.
         After VST, the STD of noise \sigma_VST is 1.
      2. Re-arrange the MSI into [0, 1]. Remember to scale \sigma_VST as well. 
      3. Apply TensorDL as in the Gaussian noise case.
      4. Reverse the operation in step 2 to re-arrange the filtered MSI.
      5. Apply inverse VST by function GenAnscombe_inverse_exact_unbiased (or alternatives provided
         by toolbox invsansc) to get the final result.
    You may find example codes in file tdl_demo\Comparison.m

    Also, you can use the demo to see some comparison. You can:
      1. Type 'TDLdemo' to see the pre-computed results.
      2. Type 'Comparison' to run various methods (NOTE: this may take a while, since some of the selected
         methods are time-consuming).
      3. Use 'help Comparison' for more information.
      4. Change test MSI by simply modifying variable 'filename' in Comparison.m (NOTE: make sure your MSI
         meets the format requirements).
      5. Change noise level by modifying variables 'kappa' and 'sigma_ratio' in Comparison.m
      6. Select competing methods by turn on/off the enable-bits in Comparison.m


Citation:

    Yi Peng, Deyu Meng, Xongben Xu, Chenqiang Gao, Yi Yang and Biao Zhang. 
    Decomposable Nonlocal Tensor Dictionary Learning for Multispectral Image Denoising.
    CVPR 2014 (accepted).


URL of third-party toolboxes and functions:

    [ 1] fkmeans            http://www.mathworks.com/matlabcentral/fileexchange/31274-fast-k-means
    [ 2] tensor_toolbox     http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
    [ 3] poblano_toolbox    https://software.sandia.gov//trac/poblano/
    [ 4] invansc            http://www.cs.tut.fi/~foi/invansc/
    [ 5] ssim_index.m       https://ece.uwaterloo.ca/~z70wang/research/ssim/
    [ 6] FeatureSIM.m       http://www4.comp.polyu.edu.hk/~cslzhang/IQA/FSIM/FSIM.htm
    [ 7] ksvdbox            http://www.cs.technion.ac.il/~ronrubin/software.html
    [ 8] ompbox             http://www.cs.technion.ac.il/~ronrubin/software.html
    [ 9] naonlm3d           http://personales.upv.es/jmanjon/denoising/arnlm.html
    [10] BM3D               http://www.cs.tut.fi/~foi/GCF-BM3D/
    [11] BM4D               http://www.cs.tut.fi/~foi/GCF-BM3D/