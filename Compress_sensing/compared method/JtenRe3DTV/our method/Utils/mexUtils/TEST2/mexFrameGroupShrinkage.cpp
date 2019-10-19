/*************************************************************************
 * Group Shrinkage :
 *          Z_i = Thres(X_iG+(1/beta)*lambda_i, (1/beta))
 *
 *
 * Usage: Z = mexGroupShrinkage(X,G,lambda,beta);
 * Author: CaoWenfei Yao Wang @Xi'an Jiaotong Univerisity
 * E-mail: caowenf2006@163.com  
 ************************************************************************/

#include<mex.h>
#include<matrix.h>
#include<cstdlib>
#include<math.h>
#include <iostream>

 /* 
void display(double * x,mwSize col,mwIndex *Ir,mwIndex *Jc,mwSize Nz)
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
    mexPrintf("\t\t num:=%d \n", Jc[col]);
    
   idx =0;
   for(int k=0; k<Nz; k++)
   {
       mexPrintf("\t\t %d, %4.4f\n", Ir[idx]+1,x[idx]);
       idx ++;
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
} */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Declare the variables */
    double *X, *lambda, beta;
    mwSize m, n, m_g, n_g, m_lam, n_lam, Nz_g, Nz_lam;
    mxLogical *G;
    mwIndex *Ir_g, *Jc_g, *Ir_lam, *Jc_lam;
    mxArray *tempZ;
    mwIndex *Ir_z,*Jc_z,*Ir_out,*Jc_out;
    double *Pr_z,*Pr_out;
    
    /*Check the input and output arguments */
    if(nlhs!=1)
        mexErrMsgTxt("The number of output arguments must set to one.");
    if(nrhs!=4)
        mexErrMsgTxt("The number of input arguments must be set to four.");
    
    if(!mxIsDouble(prhs[0]))
        mexErrMsgTxt("The first argument must be double.");
    if(mxGetNumberOfDimensions(prhs[0])!=2)
        mexErrMsgTxt("The first input argument must be two dimensional.");
    
    if(!mxIsSparse(prhs[1]))
       mexErrMsgTxt("The second input argument must be sparse.");
    if(!mxIsLogical(prhs[1]))
       mexErrMsgTxt("The second argument should be logical.");
    
    if(!mxIsDouble(prhs[2]))
       mexErrMsgTxt("The third argument must be double.");
    if(mxGetNumberOfDimensions(prhs[2])!=2)
        mexErrMsgTxt("The third input argument must be two dimensional.");
    if(!mxIsSparse(prhs[2]))
        mexErrMsgTxt("The third input argument must be sparse.");
    if(mxGetN(prhs[2])!=mxGetN(prhs[0])*mxGetN(prhs[1]))
        mexErrMsgTxt("The Cell array must be spreaded in the horizontal direction.");
    
    /*Link to the input */
    X = (double*)mxGetPr(prhs[0]);
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    
    G = mxGetLogicals(prhs[1]);
    m_g = mxGetM(prhs[1]);
    n_g = mxGetN(prhs[1]);
    Ir_g = mxGetIr(prhs[1]);
    Jc_g = mxGetJc(prhs[1]);
    Nz_g = mxGetNzmax(prhs[1]);
//     display(G,n_g,Ir_g,Jc_g);

    
    lambda = (double*)mxGetPr(prhs[2]);
    m_lam = mxGetM(prhs[2]);
    n_lam = mxGetN(prhs[2]);
    Nz_lam = mxGetNzmax(prhs[2]);
    Jc_lam = mxGetJc(prhs[2]);
    Ir_lam = mxGetIr(prhs[2]);
//     display(lambda,n_lam,Ir_lam,Jc_lam,Nz_lam);
//     mexPrintf("m_lam:%d, n_lam:%d \n",m_lam,n_lam);
    beta = mxGetScalar(prhs[3]);

    /*temp argument for output */
    tempZ = mxCreateSparse(m_lam,n_lam,Nz_g*n,mxREAL); // Nz_g*n Not Nz_lam
    Pr_z = (double*)mxGetPr(tempZ);
    Ir_z = mxGetIr(tempZ);
    Jc_z = mxGetJc(tempZ);
//     display(Pr_z,n_lam,Ir_z,Jc_z,Nz_g*n);
    
    /* main subroutine */
    int idxZ = 0;
    int idxlam1 =0;
    int idxlam2 =0;
    int idxg1 =0;
    int idxg2 =0;
    int ncol =0;
    Jc_z[0]=0;
    for(int i =0; i< n; i++)// for i'th cell
    {
        for(int j=0; j<n_g; j++)// for the j'th col in the i'th cell
        {
//             mexPrintf("******************************%d'th col**********************************\n",ncol+1);
//             mexPrintf("numofLam: %d \n",Jc_lam[n_lam]);
            
            //compute the norm for the (i*n_g+j)'s th col
            double norm = 0.0;
            idxlam2 = idxlam1;
            idxg2 = idxg1;
            for(int k =0; k<Jc_g[j+1]-Jc_g[j]; k++)
            {
                if( Jc_lam[n_lam]!= 0 )
                {   
                    norm += ( X[i*m+Ir_g[idxg1]] + (1/beta)*lambda[idxlam1] ) * ( X[i*m + Ir_g[idxg1]] + (1/beta)*lambda[idxlam1] );
//                     mexPrintf("Cell No:%d, Col No:%d, Row No:%d, x=:%4.4f, lambda:=%4.4f,sum:=%4.4f\n",i+1,j+1,Ir_g[idxg1]+1,X[i*m+Ir_g[idxg1]],lambda[idxlam1],(X[i*m + Ir_g[idxg1]] + (1/beta)*lambda[idxlam1]) );
                    idxg1 ++;
                    idxlam1 ++;
                }
                else{
                    norm += ( X[i*m + Ir_g[idxg1] ] * X[i*m + Ir_g[idxg1] ]);
                    idxg1 ++;
                }
                    
            }
            norm = sqrt(norm);
//             mexPrintf("---norm:%4.4f \n", norm);  
            
            //shrink the (i*n_g +j)'s th col
            double thres = 1/beta;
            if(norm > thres)
            {
                //copy col
//                 Jc_z[ncol+1] = Jc_z[ncol] + ( Jc_lam[ncol+1] - Jc_lam[ncol]);
                Jc_z[ncol+1] = Jc_z[ncol] + ( Jc_g[j+1] - Jc_g[j] ); 
                ncol++;
                
                //shrinkage
                double scale = (norm -thres)/norm;
                for(int k=0; k<Jc_g[j+1] - Jc_g[j]; k++)
                {
                    if( Jc_lam[n_lam]!= 0 )
                    {
                        //element-wise shrikage
                        Pr_z[idxZ] = ( X[i*m+Ir_g[idxg2]] + (1/beta)*lambda[idxlam2] ) * scale;
//                         mexPrintf("row:%d,  %4.4f \n",Ir_g[idxg2]+1, Pr_z[idxZ]);
                        //copy row
                        Ir_z[idxZ] = Ir_g[idxg2];
                        //update count
                        idxZ++;
                        idxg2++;
                        idxlam2++;
                    }else{
                        Pr_z[idxZ] = ( X[i*m + Ir_g[idxg2]] ) * scale;
//                         mexPrintf("row:%d,  %4.4f \n",Ir_g[idxg2]+1, Pr_z[idxZ]);
                        Ir_z[idxZ] = Ir_g[idxg2];
                        idxZ++;
                        idxg2++;
                    }
                }//END FOR
                
            }else{
                
                Jc_z[ncol+1] = Jc_z[ncol];
//                 mexPrintf("goooooooooooooooooooooooooooooooooooooogle\n");
                ncol++;
                
            } // END IF 
            
        }//END FOR

    // reset idxg1, idxg2
    idxg1 = 0;
    idxg2 = 0;
            
    }//END FOR
    
//     mexPrintf("ncol:%d, idxZ:%d, numofcol:%d\n",ncol,idxZ,Jc_z[ncol]);
//     display(Pr_z,n_lam,Ir_z,Jc_z,idxZ);
  
    /* prepare for output */
    plhs[0] = mxCreateSparse(m_lam,n_lam,idxZ,mxREAL);
    Pr_out = mxGetPr(plhs[0]);
    Ir_out = mxGetIr(plhs[0]);
    Jc_out = mxGetJc(plhs[0]);
    
    memcpy(Ir_out, Ir_z, idxZ * sizeof(mwIndex) );
    memcpy(Jc_out, Jc_z, (n_lam +1)* sizeof(mwIndex) );
    memcpy(Pr_out, Pr_z, idxZ * sizeof(double) );
//     display(Pr_out,n_lam,Ir_out,Jc_out,idxZ);
 
    /* clean up */
    mxDestroyArray(tempZ);
    
    return;
    
} 
 