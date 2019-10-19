#include<mex.h>
#include<omp.h>
#include<stdio.h>


void TV1D_denoise(double* output, double* input, const int width, const double lambda) {
	if (width>0) {				/*to avoid invalid memory access to input[0]*/
		int k=0, k0=0;			/*k: current sample location, k0: beginning of current segment*/
		double umin=lambda, umax=-lambda;	/*u is the dual variable*/
		double vmin=input[0]-lambda, vmax=input[0]+lambda;	/*bounds for the segment's value*/
		int kplus=0, kminus=0; 	/*last positions where umax=-lambda, umin=lambda, respectively*/
		const double twolambda=2.0*lambda;	/*auxiliary variable*/
		const double minlambda=-lambda;		/*auxiliary variable*/
		for (;;) {				/*simple loop, the exit test is inside*/
			while (k==width-1) {	/*we use the right boundary condition*/
				if (umin<0.0) {			/*vmin is too high -> negative jump necessary*/
					do output[k0++]=vmin; while (k0<=kminus);
					umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
				} else if (umax>0.0) {	/*vmax is too low -> positive jump necessary*/
					do output[k0++]=vmax; while (k0<=kplus);
					umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
				} else {
					vmin+=umin/(k-k0+1); 
					do output[k0++]=vmin; while(k0<=k); 
					return;
				}
			}

			if(k>width-1)
				break;

			if ((umin+=input[k+1]-vmin)<minlambda) {		/*negative jump necessary*/
				do output[k0++]=vmin; while (k0<=kminus);
				vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
				umin=lambda; umax=minlambda;
			} else if ((umax+=input[k+1]-vmax)>lambda) {	/*positive jump necessary*/
				do output[k0++]=vmax; while (k0<=kplus);
				vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
				umin=lambda; umax=minlambda;
			} else { 	/*no jump necessary, we continue*/
				k++;
				if (umin>=lambda) {		/*update of vmin*/
					vmin+=(umin-lambda)/((kminus=k)-k0+1);
					umin=lambda;
				} 
				if (umax<=minlambda) {	/*update of vmax*/
					vmax+=(umax+lambda)/((kplus=k)-k0+1);
					umax=minlambda;
				} 	
			}
		}
	}
}

/* 0.5*||Z_i-X_i||_F^2 + mu*||X_i||_TV, i=1,2,3. */
void anIsotrTV(double *output,double *input, int m, int n, int d, const int mode, double shrinkCoef){
   
    int i,j,k;
    int idx;
    double *tempX, *tempY;
    
    // mode-unfold
    tempX = (double*)malloc(m*n*d*sizeof(double));
    tempY = (double*)malloc(m*n*d*sizeof(double));
    if(mode == 2)
    {
     #pragma omp parallel for shared(input,tempX) private(i,j,k,idx)
        for(k =0; k<d; k++)
            for(j=0; j<n; j++)
                for(i=0; i<m; i++)
                {
                    idx = (m*n)*k + i*n + j;
                    tempX[idx] = input[(m*n)*k + m*j +i];  //m x n x d
                }
    }else if(mode == 3)
     {   
      #pragma omp parallel for shared(input,tempX) private(i,j,k,idx)
         for(k = 0; k < d; k++)
             for(j = 0; j<n; j++)
                 for(i =0; i<m; i++)
                 {
                     idx = (m*d)*j + i*d + k;
                     tempX[idx] = input[(m*n)*k + m*j +i];
                 }
     }
     
     //TV denoise for mode-unfold
     if(mode == 1)
      {
		   #pragma omp parallel for shared(input,output) private(i)
		   for(i=0;i<d*n;i++)
		       TV1D_denoise(&output[m*i], &input[m*i], (int)m, shrinkCoef);
               
      }else if(mode==2){
     
           #pragma omp parallel for shared(tempX,tempY) private(j)
           for(j=0; j<m*d; j++)
               TV1D_denoise(&tempY[n*j], &tempX[n*j], (int)n, shrinkCoef);
     
          //output <= tempY;
          #pragma omp parallel for shared(tempY,output) private(i,j,k)
          for(k=0; k<d; k++)
              for(j=0; j<n; j++)
                  for(i=0; i<m; i++)
                      output[(m*n)*k + m*j +i] = tempY[(m*n)*k + i*n + j];        
     
      }else if(mode == 3){
     
          #pragma omp parallel for shared(tempX,tempY) private(k) 
           for(k=0; k<m*n; k++)
               TV1D_denoise(&tempY[k*d], &tempX[k*d], (int)d, shrinkCoef);
     
          //output <= tempY;
          #pragma omp parallel for shared(tempY,output) private(i,j,k)
          for(k=0; k<d; k++)
              for(j=0; j<n; j++)
                  for(i=0; i<m; i++)
                      output[(m*n)*k + m*j +i] = tempY[(m*d)*j + i*d + k];  
     
      }//END IF
     
     free(tempX);
     free(tempY);
}
 
/*  output = anisotrTV(input,mu,mode) */
void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    /*set up input arguments */
    double* input    =       mxGetPr(prhs[0]);
    const mwSize *dim =      mxGetDimensions(prhs[0]);
    double mu =              mxGetScalar(prhs[1]); 
    int mode =            (int)mxGetScalar(prhs[2]);  
    
    size_t numofdim =  (size_t) mxGetNumberOfDimensions(prhs[0]);   
    size_t m, n, d;
    if(numofdim == 3){
       m = (size_t) dim[0];
       n = (size_t) dim[1];
       d = (size_t) dim[2];
//        printf("The dimensions are %d x %d x %d\n", (int) m, (int) n, (int) d);
    }else{
        printf("The dimension should be 3...%d\n", (int) numofdim);
        return;
    }
    
    double *output;
    plhs[0] = mxCreateDoubleMatrix(m, n*d, mxREAL);
    mxSetDimensions(plhs[0], dim, 3);
    output = mxGetPr(plhs[0]);
    
    anIsotrTV(output,input,m,n,d, mode,mu);
}