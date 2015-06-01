#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "perplex.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

SEXP R_ini_phaseq(SEXP filename){
  const char *instr;
  char str[p_max_filename_len];
  int ret;
  SEXP result;

  instr = CHAR(STRING_ELT(filename,0));
  strncpy(str, instr, p_max_filename_len);

  ret = ini_phaseq(str);

  result = PROTECT(allocVector(INTSXP, 1));
  memcpy(INTEGER(result), &ret, sizeof(int) * 1);
  UNPROTECT(1);
  return result; 
}

SEXP R_print_comp_order() {
  SEXP result;
  int ret = 1;
  print_comp_order();
  result = PROTECT(allocVector(INTSXP, 1));
  memcpy(INTEGER(result), &ret, sizeof(int) * 1);
  UNPROTECT(1);
  return result;
}

SEXP R_get_comp_order() {
  char *order, *spcloc;
  char compname[p_cname_len];
  int i, ret, n;
  SEXP result;

  ret = get_comp_order(&order);
  n = ret / p_cname_len;

  result = PROTECT(allocVector(STRSXP, n));

  for (i = 0; i < n; i++) {
    strncpy(compname, order + sizeof(char) * i * p_cname_len, p_cname_len);
    spcloc = strstr(compname, " ");
    if (spcloc != NULL) {
      *spcloc = (char)0;
    }
    SET_STRING_ELT(result, i, mkChar(compname));
  }

  UNPROTECT(1);

  return result;
}
    
SEXP R_phaseq(SEXP P, SEXP T, SEXP comp) {
    int retval, ncomp, ncomp_proj;
    int dbgprint = 1;
    double rP, rT, *rcomp;
    
    int nphases;
    double *wtphases;
    double *cphases;
    double *sysprop;
    double *phsysprop;
    char *namephases;
   
    int i, ic, iprop;
    
    SEXP result, ret_wtphases, ret_retval, ret_namephases, ret_cphases, ret_sysprop, ret_phsysprop;
    
    ncomp = length(comp);
    ncomp_proj = length(R_get_comp_order());
    
    if (ncomp != ncomp_proj) {
        // wrong num of components in input
        fprintf(stderr, "!! Warning: number of components\n");
    }
    
    rP = REAL(P)[0];
    rT = REAL(T)[0];
    
    rcomp = REAL(comp); 
    
    wtphases = (double *)malloc(sizeof(double) * p_size_phases);
    cphases = (double *)malloc(sizeof(double) * p_size_phases * p_size_components);
    sysprop = (double *)malloc(sizeof(double) * p_size_sysprops);
    phsysprop = (double *)malloc(sizeof(double) * p_size_sysprops * p_size_phases);
    namephases = (char *)malloc(sizeof(char) * p_size_phases * p_pname_len);
    
    if (wtphases == NULL ||
        cphases == NULL ||
        sysprop == NULL ||
        phsysprop == NULL ||
        namephases == NULL) {
        fprintf(stderr, "err!\n");
    }
    
    //fprintf(stderr, "Calling phaseq() ... ");
    retval = phaseq(rP, rT, ncomp, rcomp, &nphases, wtphases, cphases, sysprop, phsysprop, namephases, dbgprint);
    //fprintf(stderr, "Done.\n");

    if (retval != 0) {
        /* phaseq error */
        free(wtphases);
        free(cphases);
        free(sysprop);
        free(namephases);
        result = PROTECT(allocVector(VECSXP, 1));
        ret_retval = PROTECT(allocVector(REALSXP, 1));
        REAL(ret_retval)[0] = (double)retval;
        SET_VECTOR_ELT(result, 0, ret_retval);
        UNPROTECT(2);
        return(result);
    }
    
    //fprintf(stderr, "Formatting results ...\n");
    
    // convert space paddings to zeros
    spc2null(namephases, p_size_phases*p_pname_len);

    // allocate a list of vectors for: 1) ret val, 2) wt% of phases, 3) names of phases,
    // 4) comp. of phases, 5) system properties, 6) system properties for each phase
    result = PROTECT(allocVector(VECSXP, 6));
    
    ret_retval = PROTECT(allocVector(REALSXP, 1));
    ret_wtphases = PROTECT(allocVector(REALSXP, nphases));
	ret_namephases = PROTECT(allocVector(STRSXP, nphases));
	ret_cphases = PROTECT(allocVector(VECSXP, nphases));
    ret_sysprop = PROTECT(allocVector(REALSXP, p_size_sysprops));
    ret_phsysprop = PROTECT(allocVector(VECSXP, nphases));
    
    REAL(ret_retval)[0] = (double)retval;
    
    for (i = 0; i < nphases; i++) {
        REAL(ret_wtphases)[i] = wtphases[i];
		SET_STRING_ELT(ret_namephases, i, mkChar((const char *)(namephases + i*p_pname_len)));
		SET_VECTOR_ELT(ret_cphases, i, allocVector(REALSXP, ncomp));
        SET_VECTOR_ELT(ret_phsysprop, i, allocVector(REALSXP, p_size_sysprops));
        for (iprop = 0; iprop < p_size_sysprops; iprop++) {
            REAL(VECTOR_ELT(ret_phsysprop, i))[iprop] = phsysprop[i*p_size_sysprops + iprop];
        }
		for (ic = 0; ic < ncomp; ic++) {
			REAL(VECTOR_ELT(ret_cphases, i))[ic] = cphases[i*ncomp + ic];
		}
    }
    
    for (iprop = 0; iprop < p_size_sysprops; iprop++) {
        REAL(ret_sysprop)[iprop] = sysprop[iprop];
    }
    
    SET_VECTOR_ELT(result, 0, ret_retval);
    SET_VECTOR_ELT(result, 1, ret_wtphases);
    SET_VECTOR_ELT(result, 2, ret_namephases);
    SET_VECTOR_ELT(result, 3, ret_cphases);
    SET_VECTOR_ELT(result, 4, ret_sysprop);
    SET_VECTOR_ELT(result, 5, ret_phsysprop);
    
    free(wtphases);
    free(cphases);
    free(sysprop);
    free(namephases);
    free(phsysprop);
    
    UNPROTECT(7);
    
    //fprintf(stderr, "Returning.\n");
    return result;
}
