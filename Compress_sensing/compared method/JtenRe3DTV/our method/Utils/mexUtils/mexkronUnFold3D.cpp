#include<stdio.h>
#include<mex.h>
#include<omp.h>

//X_kronflat     = mexkronUnFold3D(X,idx_A,idx_B)
void mexFunction(int nlhs, mxArray *plhs[], \
                 int nrhs, const mxArray *prhs[])
{ 
    double *X, *iA, *iB, *X_f;
	size_t I1,I2,I3,J1,J2,J3,Row,Col,Dep;
	size_t i1,i2,i3,j1,j2,j3,i,j,k,i_f,j_f,id1,id2;
	
	// input
	X   = (double*)mxGetPr(prhs[0]);
	size_t numofdim    = (size_t)mxGetNumberOfDimensions(prhs[0]);
	if(numofdim !=3)
	   mexErrMsgTxt("This script aims to perform kron unfolding of 3-order tensor!");	  
	const mwSize *dimX = mxGetDimensions(prhs[0]);
	Row = dimX[0]; Col = dimX[1]; Dep = dimX[2];
	
	iA  = (double*)mxGetPr(prhs[1]);
	I1  = (size_t)iA[0];
	I2  = (size_t)iA[1];
	I3  = (size_t)iA[2];
	
	iB  = (double*)mxGetPr(prhs[2]);
	J1  = (size_t)iB[0];
	J2  = (size_t)iB[1];
	J3  = (size_t)iB[2];
	
	if(Row !=I1*J1 || Col != I2*J2 || Dep != I3*J3)
	   mexErrMsgTxt("The size is not consistent!");
	   
	//output
	plhs[0] = mxCreateDoubleMatrix(I1*I2*I3,J1*J2*J3,mxREAL);
	X_f     = (double*)mxGetPr(plhs[0]);
	
	//core block
 #pragma omp parallel for shared(X_f,X) private(i1,i2,i3,j1,j2,j3,i_f,j_f,i,j,k,id1,id2)
    for(i1=0;i1<I1;i1++)
	   for(i2=0;i2<I2;i2++)
	     for(i3=0;i3<I3;i3++)
		   for(j1=0;j1<J1;j1++)
		    for(j2=0;j2<J2;j2++)
			  for(j3=0;j3<J3;j3++)
			     {
				    i   = i1*J1 + j1;
					j   = i2*J2 + j2;
					k   = j3*J3 + j3;
					id1 = k*Row*Col + j*Row + i;
					
					i_f = i3*I2*I1 + i2*I1 + i1;
					j_f = j3*J2*J1 + j2*J1 + j1;
					id2 = j_f*(I1*I2*I3)   + i_f;
					
					X_f[id2] = X[id1];
			    }
 return;
}
	  





