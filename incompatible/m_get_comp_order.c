#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 /* #include "perplex.h" */
#include "math.h"
#include "mex.h"

#define p_cname_len 6

/* int ini_phaseq(char *inputfile); */

EXTERN_C void mexFunction(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[])
{
    char *order, *spcloc;
    char compname[p_cname_len];
    int i, j, ret, n, bytes_to_copy;
    int *charLengthArr;
    mxChar *dataPtr;

    if (nrhs > 0) {
        mexErrMsgIdAndTxt("PrplxWrap:m_get_comp_order:nrhs", "Too many arguments");
        return;
    }
    
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("PrplxWrap:m_get_comp_order:nlhs", "Expect one output argument");
        return;
    }
    
    ret = get_comp_order(&order);
    n = ret / p_cname_len;
    
    mexPrintf("%d \n %s\n", n, order);
    
    charLengthArr = (int *)malloc(sizeof(int)*n);
    for (i = 0; i < n; i++) {
        charLengthArr[i] = p_cname_len-1; /* one char is reserved for NUL */
    }
    
    plhs[0] = mxCreateCharArray(n, (const int *)charLengthArr);
    
    dataPtr = (mxChar *)mxGetData(plhs[0]);

    for (i = 0; i < n; i++) {
        memcpy(dataPtr + i * p_cname_len, order + i * p_cname_len, 1);
    }
    /*bytes_to_copy = n * p_cname_; /* n * mxGetElementSize(plhs[0]); 
    memcpy(dataPtr, order, bytes_to_copy);*/

    free(order);
}



