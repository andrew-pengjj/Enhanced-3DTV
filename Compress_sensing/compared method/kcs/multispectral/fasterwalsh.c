/**
 * C implementation of Gylson Thomas' fastwalsh. It is actually the inverse walsh transform.
 * Use 
 *     fasterwalsh(x, bitrevorder(1:length(x))) ./ length(x) 
 * for the forward transform.
 *
 * Use 
 *     fasterwalsh(x, bitrevorder(1:length(x))) 
 * for the inverse transform.
 */

#include <stdio.h>
#include <string.h>
#include "mex.h"


void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    int N;
    double *data, *idx, *x;

    if (nlhs != 1 || nrhs != 2 || !mxIsDouble (prhs[0]) || !mxIsDouble (prhs[1]))
        mexErrMsgTxt ("Usage: y = fastwalsh(x, idx)\n"
                      "  where x is a vector and \n"
                      "  idx is a permutation, for ordered transform use bitrevorder(1:length(x))");

    if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS || mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("Vector data class should be double.");

    N = mxGetNumberOfElements(prhs[0]);
    if (mxGetNumberOfElements(prhs[1]) != N)
        mexErrMsgTxt("The two vectors should have the same size!");

    data = mxGetPr(prhs[0]);
    idx = mxGetPr(prhs[1]);
            
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    x = mxGetPr(plhs[0]);

    {
        int i;
        for (i = 0; i < N; i++)
            x[i] = data[((int) idx[i]) - 1];
    }
    {
        int L, L1, k1, k2, k3, i1, i2, i3;
        L = 1;
        while ((1 << L) < N) L++;
        k1 = N, k2 = 1, k3 = N/2;
        for (i1 = 1; i1 <= L; i1++)
        {
            L1 = 1;
            for (i2 = 1; i2 <= k2; i2++)
            {
                for (i3 = 1; i3 <= k3; i3++)
                {
                    int i = i3 + L1 - 2;
                    int j = i + k3;
                    double temp1 = x[i],  temp2 = x[j];
                    if (i2 % 2 == 0)
                    {
                        x[i] = temp1 - temp2;
                        x[j] = temp1 + temp2;
                    } else
                    {
                        x[i] = temp1 + temp2;
                        x[j] = temp1 - temp2;
                    }
                }
                L1 = L1 + k1;
            }
            k1 = k1/2;
            k2 = k2*2;
            k3 = k3/2;
        }
    }
}
