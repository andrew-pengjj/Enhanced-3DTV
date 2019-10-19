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
void anIsotrTV(double *output,double *input, int m, int n, int d1, int d2, const int mode, double shrinkCoef){
   
    int k1,k2,k3,k4,i;
    int idx1,idx2;
    double *tempX, *tempY;
    
    // mode-unfold
    tempX = (double*)malloc(m*n*d1*d2*sizeof(double));
    tempY = (double*)malloc(m*n*d1*d2*sizeof(double));
    if(mode == 2)
    {
     #pragma omp parallel for shared(input,tempX) private(k1,k2,k3,k4,idx1,idx2)
	 for(k4 = 0; k4<d2; k4++)
        for(k3 =0; k3<d1; k3++)
            for(k2=0; k2<n; k2++)
                for(k1=0; k1<m; k1++)
                {
					idx1 = (n*m*d1)*k4 + (n*m)*k3 + n*k1 + k2;
					idx2 = (m*n*d1)*k4 + (m*n)*k3 + m*k2 + k1;
                    tempX[idx1] = input[idx2];  //m x n x d
                }
    }else if(mode == 3)
     {   
      #pragma omp parallel for shared(input,tempX) private(k1,k2,k3,k4,idx1,idx2)
	 for(k4 = 0; k4<d2; k4++)
        for(k3 =0; k3<d1; k3++)
            for(k2=0; k2<n; k2++)
                for(k1=0; k1<m; k1++)
                 {
       
					 idx1 = (m*n*d1)*k4 + (m*d1)*k2 + m*k1 + k3;
					 idx2 = (m*n*d1)*k4 + (m*n)*k3 + m*k2 + k1;
                     tempX[idx1] = input[idx2];
                 }
     }
	 else if(mode == 4)
	 {
	 #pragma omp parallel for shared(input,tempX) private(k1,k2,k3,k4,idx1,idx2)
	  for(k4 = 0; k4<d2; k4++)
        for(k3 =0; k3<d1; k3++)
            for(k2=0; k2<n; k2++)
                for(k1=0; k1<m; k1++)
                 {
       
					 idx1 = (m*n*d2)*k3 + (m*d2)*k2 + m*k1 + k4;
					 idx2 =  (m*n*d1)*k4 + (m*n)*k3 + m*k2 + k1;
                     tempX[idx1] = input[idx2];
                 }
	
	 }
     
     //TV denoise for mode-unfold
     if(mode == 1)
      {
		   #pragma omp parallel for shared(input,output) private(i)
		   for(i=0;i<d1*d2*n;i++)
		       TV1D_denoise(&output[m*i], &input[m*i], (int)m, shrinkCoef);
               
      }else if(mode==2){
     
           #pragma omp parallel for shared(tempX,tempY) private(i)
           for(i=0; i<m*d1*d2; i++)
               TV1D_denoise(&tempY[n*i], &tempX[n*i], (int)n, shrinkCoef);
     
          //output <= tempY;
          #pragma omp parallel for shared(tempY,output) private(k1,k2,k3,k4)
          for(k4=0; k4<d2; k4++)
		    for(k3=0; k3<d1; k3++)
               for(k2=0; k2<n; k2++)
                  for(k1=0; k1<m; k1++)
                       output[(m*n*d1)*k4 + (m*n)*k3+m*k2 +k1] = tempY[(m*n*d1)*k4 + (m*n)*k3 + n*k1 + k2];        
     
      }else if(mode == 3){
     
          #pragma omp parallel for shared(tempX,tempY) private(i)
           for(i=0; i<m*n*d2; i++)
               TV1D_denoise(&tempY[i*d1], &tempX[i*d1], (int)d1, shrinkCoef);
     
          //output <= tempY;
          #pragma omp parallel for shared(tempY,output) private(k1,k2,k3,k4)
          for(k4=0; k4<d2; k4++)
		    for(k3=0; k3<d1; k3++)
               for(k2=0; k2<n; k2++)
                  for(k1=0; k1<m; k1++)
                      output[(m*n*d1)*k4 + (m*n)*k3+m*k2 +k1]  = tempY[(d1*m*n)*k4 + (d1*m)*k2 + d1*k1 + k3];  
     
      }else if (mode == 4){
	       
		 #pragma omp parallel for shared(tempX,tempY) private(i)
           for(i=0; i<m*n*d1; i++)
               TV1D_denoise(&tempY[i*d2], &tempX[i*d2], (int)d2, shrinkCoef);
     
          //output <= tempY;
          #pragma omp parallel for shared(tempY,output) private(k1,k2,k3,k4)
          for(k4=0; k4<d2; k4++)
		    for(k3=0; k3<d1; k3++)
               for(k2=0; k2<n; k2++)
                  for(k1=0; k1<m; k1++)
                      output[(m*n*d1)*k4 + (m*n)*k3+m*k2 +k1]  = tempY[(d2*m*n)*k3 + (d2*n)*k2 + d2*k1 + k4];   
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
    size_t m, n, d1, d2;
    if(numofdim == 4){
       m = (size_t)dim[0];
       n = (size_t)dim[1];
       d1= (size_t)dim[2];
	   d2= (size_t)dim[3];
//        printf("The dimensions are %d x %d x %d\n", (int) m, (int) n, (int) d);
    }else{
        printf("The dimension should be ...%d\n", (int) numofdim);
        return;
    }
    
    double *output;
    plhs[0] = mxCreateDoubleMatrix(m, n*d1*d2, mxREAL);
    mxSetDimensions(plhs[0], dim, 4);
    output = mxGetPr(plhs[0]);
    
    anIsotrTV(output,input,m,n,d1,d2,mode,mu);
}