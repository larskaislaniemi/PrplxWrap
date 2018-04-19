#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include <string.h>
#include "mex.h"
#include "perplex.h"

/* int ini_phaseq(char *inputfile); */

static int ini_done = 0;
static int ncomp;
static double *wtphases, *cphases, *sysprop;
static char *namephases, *namephases2;
void checkIniDone();

EXTERN_C void mexFunction(int nlhs, mxArray *plhs[],
                          int nrhs, const mxArray *prhs[])
{
    char *inputfile;
    double action, P, T;
    double *num_of_components, *comp, *pind1, *pind2;
    int nphases, dbgprint, ret, iphase, icomp;

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
        ini_done=1;

        /* alloc the memory for later use (multiple times) when calling phaseq() */
        wtphases = (double *)mxMalloc(sizeof(double) * p_size_phases);
        pind1 = (double *)mxMalloc(sizeof(double) * p_size_phases);
        cphases = (double *)mxMalloc(sizeof(double) * p_size_phases * p_size_components);
        pind2 = (double *)mxMalloc(sizeof(double) * p_size_phases * p_size_components);
        sysprop = (double *)mxMalloc(sizeof(double) * p_size_sysprops);
        namephases = (char *)mxMalloc(sizeof(char) * p_size_phases * p_pname_len);
        mexMakeMemoryPersistent((void *)wtphases);
        mexMakeMemoryPersistent((void *)pind1);
        mexMakeMemoryPersistent((void *)cphases);
        mexMakeMemoryPersistent((void *)pind2);
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
        P = mxGetScalar(prhs[1]);
        T = mxGetScalar(prhs[2]);
        comp = mxGetPr(prhs[3]);

        ret = phaseq(P, T, ncomp, comp, &nphases, wtphases, cphases, sysprop, namephases, dbgprint);
        nlhs = 4;
        plhs[0] = mxCreateDoubleScalar((double)nphases);
        plhs[1] = mxCreateDoubleMatrix(nphases,1,mxREAL);
           pind1 = mxGetPr(plhs[1]);
           for (iphase=0;iphase<nphases;iphase++)
           {  pind1[iphase] = wtphases[iphase]; 
           }
        plhs[2] = mxCreateDoubleMatrix(nphases*ncomp,1,mxREAL);
           pind2 = mxGetPr(plhs[2]);
           for (iphase=0;iphase<nphases;iphase++)
           {  for (icomp=0;icomp<ncomp;icomp++)
              {  pind2[iphase*ncomp+icomp] = cphases[iphase*ncomp+icomp]; 
              }
           }
        plhs[3] = mxCreateString(namephases);
	plhs[4]= mxCreateDoubleScalar((double)ncomp);
    }
    return;
}

void checkIniDone() {
    if (ini_done) return;
    mexErrMsgIdAndTxt("PrplxWrap:m_ini_phaseq:initialization", "Do init first");
    return;  /* should never get here */
}
