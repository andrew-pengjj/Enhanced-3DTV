/*************************************************************************
 * Description: 
 * Update lambda_x2:
// function  [lambda_x2,norm_diff2] = mexUpdate_lambda_x2(X2,lambda_x2,zx2,...
//                                                G,beta2,gamma,isCont)
//    nfrm =size(X2,2);
//    norm_diff2 = 0;
// 
//    for k1 =1 : nfrm
// 
//        diagX2 = diag(sparse(X2(:,k1)));
// 
//        lambda_x2{k1} = lambda_x2{k1} - gamma*beta2*(zx2{k1} -diagX2*G);
// 
//        if isCont
//            norm_diff2 = norm_diff2 + (norm(zx2{k1} - diagX2*G,'fro'))^2;
//        end
//    end
// end
 *
 * Usage: [lambda_x2,norm_diff2] = mexUpdate_lambda_x2(X2,lambda_x2,zx2,...
 *                                                   G,beta2,gamma,isCont); 
 * Author: Wenfei Cao and Yao Wang.
 * E-mail:  caowenf2006@163.com
 ************************************************************************/

#include"mex.h"
#include"matrix.h"
#include<stdlib.h>
#include<math.h>
// #include "omp.h"

void display(double * x,mwSize col,mwIndex *Ir,mwIndex *Jc)
{
    mexPrintf("********Print Sparse Double Matrix *************\n");
    
    int idx =0;
    for(int j=0;j<col;j++)
        for(int i=0;i<Jc[j+1]-Jc[j];i++)
        {
            mexPrintf("\t\t %d,  %4.4f\n",Ir[idx]+1,x[idx]);
            idx +=1;
        }
    
    mexPrintf("\t\t idx:= %d \n", idx);
}

void display(mxLogical * x,mwSize col,mwIndex *Ir,mwIndex *Jc)
{
    mexPrintf("********Print Sparse Logical Matrix *************\n");
    
    int idx =0;
    for(int j=0;j<col;j++)
        for(int i=0;i<Jc[j+1]-Jc[j];i++)
        {
            mexPrintf("\t\t %d  %d\n",Ir[idx]+1,x[idx]);
            idx +=1;
        }
    mexPrintf("\t\t idx:= %d \n", idx);
} 
    
    
void callFunction(mxArray *outCell1,double *sum_diff, 
                  double *X,mwSize m,mwSize n,const mxArray *inCell1,const mxArray *inCell2,const mxArray *inPr3,double beta,double gamma,double isCont)
{  // callFunction(plhs[0],sum_diff,x2,m,n,prhs[1],prhs[2],prhs[3],beta2,gamma,isCont);
    
    mxArray *inPr1,*inPr2;
    double *out_lambda,*in_lambda,*z;
    mxLogical *g;
    double aux_diff=0.0;
    mwSize row1,col1,row2,col2,row3,col3;
    mwSize m_g,n_lambda,n_z,n_g,Nzmax_g;
    mwIndex *Ir_g, *Jc_g, *Jc_lambda, *Jc_z;
    mwIndex *out_Ir,*out_Jc;
    mwIndex idx;
    
    for(mwIndex index=0; index<n; index++){
        
//         mexPrintf("loop : ----------------%d ------------------\n",index+1);
        
        // Extract the mxarray pointer from cell array
        inPr1 = mxGetCell(inCell1,index);
        inPr2 = mxGetCell(inCell2,index);
     
        row1 = mxGetM(inPr1);
        col1 = mxGetN(inPr1);
        row2 = mxGetM(inPr2);
        col2 = mxGetN(inPr2);
        row3 = mxGetM(inPr3);
        col3 = mxGetN(inPr3);
        if(row1!=row3 || col1 !=col3 || row2 !=row3 || col2!=col3)
            mexErrMsgTxt("The size in arguments from the second to fifth must be equal!");
        
        if(!mxIsSparse(inPr1) ||!mxIsSparse(inPr2))
            mexErrMsgTxt("Each array contains one sparse matrix.");

        
        // link for out_lambda,in_lambda,z,g
        in_lambda = (double*)mxGetPr(inPr1);
        n_lambda = mxGetN(inPr1);
        Jc_lambda = mxGetJc(inPr1);
        
        z = (double*)mxGetPr(inPr2);
        n_z = mxGetN(inPr2);
        Jc_z = mxGetJc(inPr2);
        
        g = mxGetLogicals(inPr3);
        m_g = mxGetM(inPr3);
        n_g = mxGetN(inPr3);
        Nzmax_g = mxGetNzmax(inPr3);
        Ir_g = mxGetIr(inPr3);
        Jc_g = mxGetJc(inPr3);
//         display(g,n_g,Ir_g,Jc_g);
//         display(z,n_g,Ir_g,Jc_g);
//         display(in_lambda,n_g,Ir_g,Jc_g);
        
        // create sparse matrix for out_lambda
        mxArray *temp = mxCreateSparse(m_g,n_g,Nzmax_g,mxREAL);
        out_lambda = (double*)mxGetPr(temp);
        out_Ir = mxGetIr(temp);
        out_Jc = mxGetJc(temp);
      
        
        //lambda_out<-lambda_in - gamma*beta(z-Xn*g)
        out_Jc[0]=0;
        idx =0;
//         omp_set_num_threads(7);
//         #pragma omp parallel for reduction(+: aux_diff)
        for(int j=0; j<n_g; j++)
        {
            for(int i=0; i<Jc_g[j+1] - Jc_g[j]; i++) 
            {
                if(Jc_lambda[n_lambda]!=0 && Jc_z[n_z] !=0)
                 {
                    out_lambda[idx] = in_lambda[idx] - gamma*beta*(z[idx] - X[m*index+Ir_g[idx]]);
                    if(isCont)
                       aux_diff = aux_diff + (z[idx] - X[m*index+Ir_g[idx]])*(z[idx] - X[m*index+Ir_g[idx]]);
                    
                 }else if(Jc_lambda[n_lambda]!=0){
                     
                    out_lambda[idx] = in_lambda[idx] + gamma*beta*X[m*index+Ir_g[idx]];
                    if(isCont)
                       aux_diff = aux_diff + (X[m*index+Ir_g[idx]])*(X[m*index+Ir_g[idx]]);
                    
                 }else if(Jc_z[n_z] !=0){
                     
                    out_lambda[idx] = -gamma*beta*(z[idx] - X[m*index+Ir_g[idx]]);
                    if(isCont)
                       aux_diff = aux_diff + (z[idx] - X[m*index+Ir_g[idx]])*(z[idx] - X[m*index+Ir_g[idx]]);
                    
                 }else{
                    out_lambda[idx] = gamma*beta*(X[m*index+Ir_g[idx]]);
                    if(isCont)
                       aux_diff = aux_diff + (X[m*index+Ir_g[idx]])*(X[m*index+Ir_g[idx]]);
                 }// END IF 
                 
                out_Ir[idx] = Ir_g[idx];  

                idx += 1;
           }//END IF
            
       out_Jc[j+1] = Jc_g[j+1];
      }//END FOR
        
   // copy temp to outCell1
   mxSetCell(outCell1,index,temp);
       
   } // END FOR   
    
    
sum_diff[0] = aux_diff;
  
}//END FUN.



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{   
    /* decare the variables */
    double *x2, beta2, gamma, isCont, *sum_diff; 
    mwSize m, n, ndim;
    const mwSize *dims;    
    
    
    /* Check for input and output arguments */
    if(nrhs!=7)
        mexErrMsgTxt("Please check the correct number of input argument.");
    
    if(nlhs!=2)
        mexErrMsgTxt("Please check the correct number of output argument.");
    
    if(!(mxIsDouble(prhs[0])))
        mexErrMsgTxt("The first input argument must be of typer double.");
  
    if(mxGetNumberOfDimensions(prhs[0])!=2)
        mexErrMsgTxt("The first input argument must be two dimensional.");
    
    if(!mxIsCell(prhs[1]) ||!mxIsCell(prhs[2]))
        mexErrMsgTxt("The second and third arguments should be cell array.");
    
    if(mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]))
        mexErrMsgTxt("The size in the second and third argument must be equal.");
    
    if(!mxIsSparse(prhs[3]))
        mexErrMsgTxt("The forth argument should be sparse.");
    
    if(!mxIsLogical(prhs[3]))
        mexErrMsgTxt("The forth argument should be sparse.");
  
    
   /* Get the size and Pointer to the input arguments */
    //x2
    x2 =(double*)mxGetPr(prhs[0]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    //lambda_x2, zx2
    ndim = mxGetNumberOfDimensions(prhs[1]);
    dims = mxGetDimensions(prhs[1]);
    
    //beta2,gamma,isCont
    beta2 =mxGetScalar(prhs[4]);
    gamma =mxGetScalar(prhs[5]);
    isCont =mxGetScalar(prhs[6]);    

    
    /* Create space for output arguments */
    plhs[0] = mxCreateCellArray(ndim, dims);
    plhs[1] = mxCreateDoubleScalar(0.0);
    sum_diff = (double*)mxGetPr(plhs[1]);
       
    /* Main-subroutine  */
    callFunction(plhs[0],sum_diff,x2,m,n,prhs[1],prhs[2],prhs[3],beta2,gamma,isCont);
    
  return;
}//END FUNCTION
  
