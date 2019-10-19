/*************************************************************************
* Description:
//  function [tempx2] = RelatedGradient2(X2,G,lambda_x2,z2,beta2)
//  
//   [m,n]= size(X2);
//   tempx2 = zeros(m,n);
//   for nfrm = 1:n %% please speed up!!!
//       
//      tempx2(:,nfrm) = sum(G.*(lambda_x2{nfrm} - beta(2)*zx2{nfrm}),2) +...
//                                 beta(2)*sum(diag(sparse(X2(:,nfrm)))*G, 2);
//                       
//   end;
* 
* Usage: tempX2 = mexRelatedGradient2(X2,G,lambda_x2,z2,beta2); 
* Note here lambda_x2 = cell2mat(lambda_x2); z2 also the same.
* Author: Wenfei Cao, Yao Wang;
* E-mail: caowenf2006@163.com.
*************************************************************************/

#include<mex.h>
#include<matrix.h>

// void display(double * x,mwSize col,mwIndex *Ir,mwIndex *Jc)
// {
//     mexPrintf("********Print Sparse Double Matrix *************\n");
//     
//     int idx =0;
//     for(int j=0;j<col;j++)
//         for(int i=0;i<Jc[j+1]-Jc[j];i++)
//         {
//             mexPrintf("\t\t %d,  %4.4f\n",Ir[idx]+1,x[idx]);
//             idx +=1;
//         }
//     
//     mexPrintf("\t\t idx:= %d \n", idx);
// }
// 
// void display(mxLogical * x,mwSize col,mwIndex *Ir,mwIndex *Jc)
// {
//     mexPrintf("********Print Sparse Logical Matrix *************\n");
//     
//     int idx =0;
//     for(int j=0;j<col;j++)
//         for(int i=0;i<Jc[j+1]-Jc[j];i++)
//         {
//             mexPrintf("\t\t %d  %d\n",Ir[idx]+1,x[idx]);
//             idx +=1;
//         }
//     mexPrintf("\t\t idx:= %d \n", idx);
// } 
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* declare the variables */
    double *X2,*Lam2,*Z2,beta;
    mxLogical *G;
    mwSize m_x,n_x,m_g,n_g,m_lam,n_lam,n_z2,Nz_g;
    mwIndex *Ir_g,*Jc_g;
    
    mxArray *Temp;
    mwIndex *Ir_tp,*Jc_tp,*Jc_lam,*Jc_z2;
    double *Pr_tp;
    double *Out;
   
    /* check for inputs and outputs */
    if(nlhs!=1)
        mexErrMsgTxt("Only one output argument required.");
    if(nrhs!=5)
        mexErrMsgTxt("Five input arguments required.");
    
    if(!mxIsDouble(prhs[0]))
        mexErrMsgTxt("The first argument must be double.");
    if(mxGetNumberOfDimensions(prhs[0])!=2)
        mexErrMsgTxt("The first argument must be of two dimensions");
    
    if(!mxIsLogical(prhs[1]))
        mexErrMsgTxt("The second argument must be logical.");
    if(!mxIsSparse(prhs[1]))
        mexErrMsgTxt("The second argument must be sparse.");
    
    if(!mxIsDouble(prhs[2]))
        mexErrMsgTxt("The third argument must be double.");
    if(!mxIsSparse(prhs[2]))
        mexErrMsgTxt("The third argument must be sparse.");
    
    if(!mxIsDouble(prhs[3]))
        mexErrMsgTxt("The forth argument must be double.");
    if(!mxIsSparse(prhs[3]))
        mexErrMsgTxt("The forth argument must be sparse.");
    if(mxGetM(prhs[2])!=mxGetM(prhs[3]) || mxGetN(prhs[2])!= mxGetN(prhs[3]))
        mexErrMsgTxt("The size for the third and forth argument must be equal."); 
    if(mxGetN(prhs[2])!=mxGetN(prhs[0])*mxGetN(prhs[1]))
        mexErrMsgTxt("The Cell array must be spreaded in the horizontal direction.");
    
    /*link to the input arguments */
    X2 = (double*)mxGetPr(prhs[0]);
    m_x = mxGetM(prhs[0]);
    n_x = mxGetN(prhs[0]);
    
    G = mxGetLogicals(prhs[1]);
    m_g = mxGetM(prhs[1]);
    n_g = mxGetN(prhs[1]);
    Ir_g = mxGetIr(prhs[1]);
    Jc_g = mxGetJc(prhs[1]);
    Nz_g = mxGetNzmax(prhs[1]);
                                     
    Lam2 = (double*)mxGetPr(prhs[2]);
    m_lam= mxGetM(prhs[2]);
    n_lam= mxGetN(prhs[2]);
    Jc_lam = mxGetJc(prhs[2]);
    
    Z2 = (double*)mxGetPr(prhs[3]);
    n_z2 = mxGetN(prhs[3]);
    Jc_z2= mxGetJc(prhs[3]);
    beta = mxGetScalar(prhs[4]);
    
    /*prepare for the output arguments */
    plhs[0] = mxCreateDoubleMatrix(m_x,n_x,mxREAL);
    Out = (double*)mxGetPr(plhs[0]);
    Temp = mxCreateSparse(m_g,n_g,Nz_g,mxREAL);
    Ir_tp = mxGetIr(Temp);
    Jc_tp = mxGetJc(Temp);
    Pr_tp = mxGetPr(Temp);
    
    
    /* main subroutine */
    for(int i=0; i<n_x; i++)
    {
        // Compute Temp<- G.*(lam_i - beta*Z2_i) + beta*diag(X2_i)*G
        int idx =0;
        Jc_tp[0]=0;
        for(int col=0; col< n_g; col++)
        {
            for(int row=0; row<Jc_g[col+1] -Jc_g[col]; row++)
            {
                if( Jc_lam[n_lam]!=0 && Jc_z2[n_z2] !=0){
                    
                    Pr_tp[idx] = (Lam2[i*Nz_g+idx] -beta*Z2[i*Nz_g+idx]) + beta*X2[i*m_x+Ir_g[idx]];
                    Ir_tp[idx] = Ir_g[idx];
                    idx++;
                    
                }else if(Jc_lam[n_lam]!=0){
                    
                    Pr_tp[idx] = Lam2[i*Nz_g+idx] + beta*X2[i*m_x + Ir_g[idx]];
                    Ir_tp[idx] = Ir_g[idx];
                    idx++;
                    
                }else if(Jc_z2[n_z2]!=0){
                    
                    Pr_tp[idx] = -beta*Z2[i*Nz_g+idx] + beta*X2[i*m_x + Ir_g[idx]];
                    Ir_tp[idx] = Ir_g[idx];
                    idx++;
                    
                }else{
                    
                    Pr_tp[idx] = beta*X2[i*m_x + Ir_g[idx]];
                    Ir_tp[idx] = Ir_g[idx];
                    idx++;
                }
                    
            }//END FOR

          Jc_tp[col+1]=Jc_tp[col] + (Jc_g[col+1] - Jc_g[col]);
        }

        // Compute sum(temp,2)
        for(int j=0; j< m_x; j++)
        Out[i*m_x +j]=0;

        idx = 0;
        for(int col=0; col<n_g; col++)
         for(int row=0; row<Jc_tp[col+1] -Jc_tp[col]; row++)
         {
            Out[i*m_x + Ir_tp[idx] ] += Pr_tp[idx];
            idx++;
         }
    
    }//END FOR
       
        
   /* clear up  */
   mxDestroyArray(Temp);
    
   return;        
}