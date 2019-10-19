/* this script aims to perform kronecker unfolding of one matrix */
#include<stdio.h>
#include<mex.h>
#include<omp.h>

//C_flat  = kronUnFlod(C,szgA,szgB);
void mexFunction(int nlhs, mxArray *plhs[], \
                 int nrhs, const mxArray *prhs[])
{ 
    double *inPr1, *inPr2,*inPr3, *outPr;
	size_t I1,I2,J1,J2,Row,Col;
	size_t i1,i2,j1,j2,id1,id2,i,j,i_f,j_f;
	
	Row   = mxGetM(prhs[0]); //C
	Col   = mxGetN(prhs[0]);
	inPr1 = (double*)mxGetPr(prhs[0]);
	
    inPr2 = (double*)mxGetPr(prhs[1]); //A
	I1   = (size_t)inPr2[0]; 
	I2   = (size_t)inPr2[1];
	
    inPr3 = (double*)mxGetPr(prhs[2]);//B
	J1   = (size_t)inPr3[0];
	J2   = (size_t)inPr3[1];
    //mexPrintf("I1:%d, I2:%d, J1:%d, J2:%d \n", I1,I2,J1,J2);
    
	if(Row != I1*J1 || Col != I2*J2) //  A_I1xI2 B_J1xJ2, C_row x col = A \kron B
	   mexErrMsgTxt("The size of C is not equal to A kron B!");
		
	plhs[0] = mxCreateDoubleMatrix(I1*I2, J1*J2, mxREAL); // C_flat
	outPr   = (double*)mxGetPr(plhs[0]);
	
 #pragma omp parallel for shared(inPr1,outPr) private(i1,i2,j1,j2,i_f,j_f,i,j,id1,id2)
     for(i1=0;i1<I1;i1++)
	  for(i2=0;i2<I2;i2++)
	    for(j1=0;j1<J1;j1++)
		  for(j2=0;j2<J2;j2++)
		     {
			   i   =  i1*J1 + j1;
			   j   =  i2*J2 + j2;
	           id1 =  j*Row + i;
			   
			   i_f =  i2*I1 + i1;
			   j_f =  j2*J1 + j1;
			   id2 =  j_f*(I1*I2) + i_f;
			   
			   outPr[id2] = inPr1[id1];
			 }
			   
 return;
}
	   
	
	
          