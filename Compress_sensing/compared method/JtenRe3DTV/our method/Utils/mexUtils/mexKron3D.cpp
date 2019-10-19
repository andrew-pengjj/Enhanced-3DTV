#include<stdio.h>
#include<mex.h>
#include<omp.h>



void kron3D(double *C,double *A,double*B,size_t I1,size_t I2,size_t I3,\
           size_t J1,size_t J2,size_t J3,size_t K1,size_t K2,size_t K3){

     size_t i1,i2,i3,j1,j2,j3,k1,k2,k3,id1,id2,id3;
	 
 #pragma omp parallel for shared(C,A,B) private(i1,i2,i3,j1,j2,j3,k1,k2,k3,id1,id2,id3)
	 for(i1 = 0; i1<I1; i1++)
	   for(i2=0; i2<I2; i2++)
	      for(i3=0; i3<I3; i3++)
		      for(j1=0;j1<J1;j1++)
			    for(j2=0;j2<J2;j2++)
				   for(j3=0;j3<J3;j3++)
				      {
					    k1 = i1*J1 + j1;
                        k2 = i2*J2 + j2;
                        k3 = i3*J3 + j3;
                        id3 = k3*K1*K2 + k2*K1 + k1;
						
						id1 = i3*I2*I1 + i2*I1 + i1;
						id2 = j3*J2*J1 + j2*J1 + j1;
						
                        C[id3] = A[id1]*B[id2];		
                      }						

}


//C = mexkron3D(A,B) C = A \kron B
void mexFunction(int nlhs, mxArray *plhs[], \
                 int nrhs, const mxArray *prhs[]){
	
    double *A, *B, *C;
    size_t I1, I2, I3, J1, J2, J3, K1,K2,K3;

	A = (double*) mxGetPr(prhs[0]);
	const mwSize *dim1 = mxGetDimensions(prhs[0]);
	B = (double*) mxGetPr(prhs[1]);
	const mwSize *dim2 = mxGetDimensions(prhs[1]);

    size_t numofdim1 = (size_t)mxGetNumberOfDimensions(prhs[0]);
	size_t numofdim2 = (size_t)mxGetNumberOfDimensions(prhs[1]);
	
	if(numofdim1 == 3){
	    I1 = (size_t)dim1[0];
		I2 = (size_t)dim1[1];
		I3 = (size_t)dim1[2];
    }else{
	    mexPrintf("The dimension of A in A kron B should be...%d\n",(int) numofdim1);
		return;
    }
	
	if(numofdim2 ==3){      //A
	    J1 = (size_t)dim2[0];
		J2 = (size_t)dim2[1];
		J3 = (size_t)dim2[2];
    }else if(numofdim2 == 2){//B
	    J1 = (size_t)dim2[0];
		J2 = (size_t)dim2[1];
		J3 = 1;
    }else{
	    mexErrMsgTxt(" The dimension of B in A kron B should be 2 or 3.");
    }
	
	K1 = I1*J1; K2 = I2*J2; K3 = I3*J3;
	int *dim = (int*)malloc(3*sizeof(int));
	dim[0] = K1; dim[1] = K2; dim[2] = K3;
    plhs[0] = mxCreateDoubleMatrix(K1,K2*K3,mxREAL);
    mxSetDimensions(plhs[0], reinterpret_cast<const int *>(dim), 3);
    C       = (double*) mxGetPr(plhs[0]);	
	
	kron3D(C,A,B,I1,I2,I3,J1,J2,J3,K1,K2,K3);
	free(dim);
}






				 