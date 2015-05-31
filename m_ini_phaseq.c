#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "mex.h"
#include "perplex.h"

/* int ini_phaseq(char *inputfile); */

static int ini_done = 0;
static int ncomp;
static double *wtphases, *cphases, *sysprop;
static char *namephases;
void checkIniDone();

EXTERN_C void mexFunction(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[])
{
    char *inputfile;
    double action, P, T;
    double *num_of_components, *comp, *ptr;
    int nphases, dbgprint, ret, i, j, k;

    if (nrhs < 1) {
        mexErrMsgIdAndTxt("PrplxWrap:m_ini_phaseq:nrhs", "Need at least one argument: action");
        return;
    }

    action = mxGetScalar(prhs[0]);

    if (action == 0) {
        /* initialization */
        if (nrhs != 2) {
            mexErrMsgIdAndTxt("PrplxWrap:m_ini_phaseq:nrhs", "Need one argument: 0 inputfile");
            return;
        }
        inputfile = mxArrayToString(prhs[1]);
        printf("%s\n", inputfile);
        ini_phaseq(inputfile);
        nlhs = 1;
        ncomp = number_of_components();
        plhs[0] = mxCreateDoubleScalar((double)ncomp);
        ini_done = 1;

        /* alloc the memory for later use (multiple times) when calling phaseq() */
        wtphases = (double *)mxMalloc(sizeof(double) * p_size_phases);
        cphases = (double *)mxMalloc(sizeof(double) * p_size_phases * p_size_components);
        sysprop = (double *)mxMalloc(sizeof(double) * p_size_sysprops);
        namephases = (char *)mxMalloc(sizeof(char) * p_size_phases * p_pname_len);
        mexMakeMemoryPersistent((void *)wtphases);
        mexMakeMemoryPersistent((void *)cphases);
        mexMakeMemoryPersistent((void *)sysprop);
        mexMakeMemoryPersistent((void *)namephases);

    } else if (action == 1) {
        checkIniDone();
        /* print comps */
        print_comp_order();
        nlhs = 0;
    } else if (action == 2) {
        checkIniDone();
        /* do minimization */
        if (nrhs == 4) {
           dbgprint = 0;
        } else if (nrhs == 5) {
           dbgprint = mxGetScalar(prhs[4]);
        } else {
            mexErrMsgIdAndTxt("PrplxWrap:m_ini_phaseq:nrhs", "Need four to five arguments: 2 P T composition [dbg]");
        }
        /* int phaseq(double P, double T, int ncomp, double *comp, int *nphases,
                    double *wtphases, double *cphases, double *sysprop, char *namephases, int dbgprint) */
        int output_array_sizes[2] = {1, ncomp};

        P = mxGetScalar(prhs[1]);
        T = mxGetScalar(prhs[2]);
        comp = mxGetPr(prhs[3]);

        ret = phaseq(P, T, ncomp, comp, &nphases, wtphases, cphases, sysprop, namephases, dbgprint);

        nlhs = 2;
        
        plhs[0] = mxCreateDoubleScalar((double)nphases);
        
        plhs[1] = mxCreateDoubleMatrix(1, nphases, mxREAL);
        ptr = mxGetPr(plhs[1]);
        for (i = 0; i < nphases; i++) {
            ptr[i] = wtphases[i];
        }
    }

    return;
}

void checkIniDone() {
    if (ini_done) return;
    mexErrMsgIdAndTxt("PrplxWrap:m_ini_phaseq:initialization", "Do init first");
    return;  /* should never get here */
}


