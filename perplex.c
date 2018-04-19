#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "perplex.h"

/*
 * C wrapper for Perple_X v6.6.8
 * - Perple_X by Jamie Connolly (see http://www.perplex.ethz.ch/)
 * - This wrapper by Lars Kaislaniemi (lars.kaislaniemi@iki.fi)
 */

/* Some general notes on C-Fortran integration:
 * (Environment: GNU C, GNU Fortran, 64bit Linux)
 *
 *  = Fortran 'logical' is C 'int'
 *  = Fortran arr(i,j) is C arr[j-1][i-1]
 *    (Fortran starts arrays from 1 by default,
 *    C from 0 always). Thus, Fortran expression
 *    'arr1(integer_pointerarray_to_element(i))' is
 *    in C 'arr1[integer_pointerarray_to_element[i-1]-1]'.
 *  = Any common block definitions passed between Fortran
 *    subroutines should be defined as (extern) structs in C
 *    (see perplex.h for example) AND variables must be 
 *    in the same order inside the struct and the common block
 *    AND of same byte size. (CRUCIAL!)
 *  = Fortran does not end strings in \0, instead they are
 *    padded with spaces. So you need to know the string length
 *    in C in order to print it properly (otherwise you always
 *    print all of the remaining chars in the array).
 *  = Fortran subroutines and common blocks are referred in C
 *    by their names added with an underscore (_).
 *  = All Fortran subroutines are void functions in C,
 *    and all parameters must be passed as pointers.
 *  = With GNU environment, at least, you need to use option
 *    '-lgfortran' to link C object files with Fortran object files
 *
 * Some notes on usage of this wrapper:
 *  = For now, the input file needs to be built with BUILD
 *    before calling the wrapper. Pretty straightforward, and probably
 *    the best option is to define "fractionation calculation" option
 *    (num 5) and define a non-existent input file for P-T path.
 *    NB! Bulk composition in input file does not matter, BUT
 *    the values need to be greater than zero (i.e. bulk composition
 *    sum needs to be positive).
 *  = If vertex has been run for the same input file this wrapper
 *    is using, PerpleX asks (interactively, as would MEEMUM) whether 
 *    the calculated data should be used. Currently, all that can
 *    be done about this is to remove the output data of vertex.
 */


/* Changes:
 * - 2013-06-11: "Warm start" option added (see ini_phaseq()).
 *               Added print_comp_order()
 */

/* whether to set cst111_.istart to 1 before calling lpopt */
int lpopt_warmstart;

int ini_phaseq(char *inputfile) {
	/*
	 * ini_phaseq(): call Perple_X initialization subroutines.
	 * This needs to be called every time the input file 
	 * (by BUILD) has been changed.
	 * 
	 * input:
	 *  - inputfile	input file name (created by Perple_X BUILD)
	 *
	 * return value:
	 *  -  0: OK
	 *  - !0: error
	 */

	int first, output, err;

	/* Following is pretty much a copy&paste from
	 * Perple_X program MEEMUM. Instead of calling 
	 * iniprp subroutine (that asks the user interactively
	 * to give filenames etc.) we set the corresponding 
	 * parameters here and don't call iniprp at all.
	 */

	/*vrsion_();*/

	cst4_.iam = 2;	/* we fake us as being Perple_X program MEEMUM.
	                 * TODO: create new value for this wrapper in
	                 * Perple_X (to suppress extra output)
	                 */

	cxt20_.rxn = 0;
	first = 1;
	output = 1;
	err = 0;
	memcpy(cst228_.prject, inputfile, strlen(inputfile));
	input1_(&first, &output, &err);  // NB. 'err' is a new option 
	                                 // introduced in PerpleX 6.6.8.
	opts_.iopt[5] = 0; // no auto-refine
	input2_(&first);
	setau1_(&output);
	input9_(&first, &output);
	initlp_();
	
    /* This flag defines whether we try to do a "warm start" with
     * PerpleX. Unfortunately, it doesn't seem to work, or have
     * any positive effect in terms of computation time. (2013-07-04)
     * This flag will be set to 1 after the first call to phaseq().
     * If you want to disable "warm start" completely, switch 
     * the flag USE_LPOPT_WARM_START in makefile.
     */
    lpopt_warmstart = 0;
    
	
	/* TODO: check whether initialization succeeded or not,
	 * and adjust return value accordingly */

	return 0;
}
	

void print_comp_order() {
	/* This function prints the order of components. This is not 
	 * always the same as the order in input file. The composition
	 * ('comp' parameter for phaseq()) should be given using
	 * this order.
	 */
    int i;
    fprintf(stdout, "Expecting compositions in the following order:\n");
    for (i = 1; i <= cst6_.icomp; i++) {
        fprintf(stderr, "\t%.5s\n", csta4_.cname[i-1]);
    }
    return;
}

int get_comp_order(char **order) {
    char *tmporder;
    int n = number_of_components();
    int i;
    *order = (char *)malloc(n * p_cname_len * sizeof(char));
    if (*order == NULL) {
        return -1;
    }
    tmporder = *order;
    for (i = 1; i <= n; i++) {
        memcpy((void *)(tmporder + (i-1)*p_cname_len*sizeof(char)), 
               (void *)(csta4_.cname[i-1]),
               sizeof(char) * (p_cname_len-1));
        ((char *)(tmporder + i*p_cname_len - 1))[0] = '\0';
    }
    return n * p_cname_len;
}

int number_of_components() {
    return cst6_.icomp;
}

int number_of_sysprops() {
	return p_i8;
}

int phaseq(double P, double T, int ncomp, double *comp, int *nphases, 
	double *wtphases, double *cphases, double *sysprop, char *namephases, int dbgprint) {
	/*
	 * phaseq(): Calculate phase equilibria using Perple_X subroutines
	 * 
	 * input: 
	 *  - P     	  pressure in bars
	 *  - T     	  temperature in Kelvins
	 *  - ncomp 	  number of components
	 *  - comp  	  array of size ncomp: component amounts
	 *  - dbgprint	  debug printing on/off (1/0) ?
	 * 
	 * output:
	 *  - nphases	  number of stable phases
	 *  - wtphases	  array of size nphases: phase amounts 
	 *  - cphases	  array of size nphases*ncomp: composition of 
	 *                phases
	 *  - sysprop     an array of system properties (see below)
	 *  - namephases  names of the stable phases returned
	 * 
	 * return value:  0: OK
	 *               >0: Perple_X failed minimization
	 *               <0: Other error
	 * 
	 * The system properties array has length if p_i8 (see defs in perplex.h)
	 * and has following properties in it:
	 *  0 V J/bar*
	 *  1 H J*
	 *  2 3 Gruneisen_T
	 *  3 Ks bar
	 *  4 Mu bar
	 *  5 V0 km/s
	 *  6 Vp km/s
	 *  7 Vs km/s
	 *  8 ?
	 *  9 Rho kg/m3
	 *  10 ?
	 *  11 Cp J/K*
	 *  12 alpha 1/K
	 *  13 beta 1/bar
	 *  14 S J/K*
	 *  15 ?
	 *  16 N g*
	 * (*) = molar property
	 * 
	 */

	
	cst4_.iam = 2;	/* we fake us as being Perple_X program MEEMUM.
	                 * TODO: create new value for this wrapper in
	                 * Perple_X (to suppress extra output)
	                 */
	int i, j, l;
	int result;
	double sum;
	int same_comp = 1;
	
	/* some return values from Perple_X subroutines */
	int itri[4] = { 0, 0, 0, 0 };
	int jtri[4] = { 0, 0, 0, 0 };
	int ijpt = 0;
	double wt[3] = { 0, 0, 0 };
	int nodata = 0;

	int numofphases;
	
	if (*csta2_.vname[cst24_.jv[0]-1] == 'T' &&
        *csta2_.vname[cst24_.jv[1]-1] == 'P') {
		cst5_.v[cst24_.jv[0]-1] = T;
		cst5_.v[cst24_.jv[1]-1] = P;
	} else if (*csta2_.vname[cst24_.jv[0]-1] == 'P' &&
        *csta2_.vname[cst24_.jv[1]-1] == 'T') {
		cst5_.v[cst24_.jv[0]-1] = P;
		cst5_.v[cst24_.jv[1]-1] = T;
	} else {
		fprintf(stderr, "!! Independent variables in configuration file are not P & T, please fix.\n");
		return -1;
	}
	
	if (ncomp != cst6_.icp) {
		fprintf(stderr, "!! Number of components (ncomp=%d) disagrees with number or components in Perple_X input file (=%d).\n",
			ncomp, cst6_.icp);
		return -1;
	} else if (ncomp > p_k5) {
		fprintf(stderr, "!! Number of components (ncomp=%d) exceeds maximum number of components (=%d).\n",
			ncomp, p_k5);
		return -1;
	}
	

	for (i = 0; i < ncomp; i++) {
        cst300_.cblk[i] = comp[i];
	}
	for ( ; i < p_k5; i++) {
        cst300_.cblk[i] = 0.0;
	}
	

	/* Convert wt% in cblk to molar% and normalize to one 
	 * in cst313_.b[]
	 */
	for (i = 0; i < cst6_.icp; i++) {
		cst300_.cblk[i] = cst300_.cblk[i] / cst45_.atwt[i];
	}	
	sum = 0.0;
	for (i = 0; i < cst6_.icp; i++) {
		sum += cst300_.cblk[i];
	}
	for (i = 0; i < cst6_.icp; i++) {
        double val;
        val = cst300_.cblk[i] / sum;
        if (cst313_.b[i] != val) {
    		cst313_.b[i] = val;
            same_comp = 0;
        }
	}
	for ( ; i < p_k5; i++) {
		cst313_.b[i] = 0.0;
	}
	
	if (dbgprint) {
		fprintf(stderr, "ipot is %d\n", cst24_.ipot);
		fprintf(stderr, "jbulk is %d\n", cst300_.jbulk);
	

		/*fprintf(stderr, "Pointer for 0 is %d\n", cst24_.jv[1]);*/
		fprintf(stderr, "Pot are %s\n", csta2_.vname[cst24_.jv[0]-1]);
		/*fprintf(stderr, "Molar mass of %s is %f\n", csta4_.cname[0], cst45_.atwt[0]);*/
		fprintf(stderr, "Potential value of 1 is %f\n", cst5_.v[cst24_.jv[1]-1]);
		fprintf(stderr, "Molar composition (sum=%g: ", sum);
		for (i = 0; i < cst6_.icp; i++) {
			fprintf(stderr, " %f", cst313_.b[i]);
		}
		fprintf(stderr, "\n");
	}
	
	/* Perple_X phase eq optimization routine */
#ifdef USE_LPOPT_WARM_START
    if (lpopt_warmstart && same_comp) {
        result = -1;
        cst111_.istart = 2;
        fprintf(stderr, "1 lpopt0\n");
        lpopt0_(&result);
    } else {
        result = 0;
        cst111_.istart = 0;
        fprintf(stderr, "2 lpopt0\n");
        lpopt0_(&result);
    }
#else
    result = 0;
    cst111_.istart = 0;
    fprintf(stderr, "3 lpopt0\n");
    lpopt0_(&result);
#endif

	if (result > 0) {
		if (dbgprint) {
			fprintf(stderr, "!! Minimization failed (%d)\n", result);
		}
		return result;
	} else {
        if (dbgprint) {
            fprintf(stderr, "!! Minimization good (%d)\n", result);
        }
    }

	/* Perple_X subroutine to calculate derivative properties.
	 * This might be unnecessary if only phase composition data
	 * is needed.
	 */
	getloc_(itri, jtri, &ijpt, wt, &nodata);
	
	numofphases = cxt15_.ntot;
	*nphases = numofphases;


	/* Commented out: */
	/* This version assumes memory for the pointer arguments has been reserved by 
	 * the caller and that they are big enough to hold all the numbers... */
	/*freearr((void**)wtphases);
	freearr((void**)cphases);
	freearr((void**)namephases);
	freearr((void**)sysprop);
	*wtphases = (double *)calloc(numofphases, sizeof(double));
	*cphases = (double *)calloc(numofphases * ncomp, sizeof(double));
	*namephases = (char *)calloc(numofphases * p_pname_len, sizeof(char));
	*sysprop = (double *)calloc(p_i8, sizeof(double));*/

	for (i = 0; i < numofphases; i++) {
		strncpy(namephases + i*p_pname_len, (const char *)(&cxt21a_.pname[i]), p_pname_len-1);
		(namephases + (i+1)*p_pname_len - 1)[0] = ' '; /* JvH 20160806 changed this from '\0' */

		(wtphases)[i] = 100 * cxt22_.props[i][16] * cxt22_.props[i][15] / cxt22_.psys[16];

		for (j = 0; j < ncomp; j++) {
			(cphases)[i*ncomp + j] = cst324_.pcomp[i][j];
		}
	}
	for (i = 0; i < p_i8; i++) {
		(sysprop)[i] = cxt22_.psys[i];
	}

	if (dbgprint) {	
        fprintf(stderr, "\n\nSystem:\n\n");
        for (l = 0; l < 17; l++) fprintf(stderr, " %f ", cxt22_.psys[l]);

		for (l = 0; l < cxt15_.ntot; l++) {
			for (i = 0; i < cst6_.icomp; i++) {
				fprintf(stderr, "%f\t", cst324_.pcomp[l][i]);
			} 
			fprintf(stderr, "\n");
		}
		for (l = 0; l < cxt15_.ntot; l++) {
			for (i = 0; i < cst6_.icomp; i++) {
				fprintf(stderr, "%f\t", (cphases)[l*ncomp + i]);
			} 
			fprintf(stderr, "\n");
		}
	}
		


	if (dbgprint) {
		for (i = 0; i < cst24_.ipot; i++) {
			fprintf(stderr, "%.8s: %f\n", csta2_.vname[cst24_.jv[i]-1], cst5_.v[cst24_.jv[i]-1]);
		}

		fprintf(stderr, "Stable phases (%d), molar (0) / weight (1) proportions: %d:\n", numofphases, opts_.iopt[1]);
		for (i = 0; i < cst6_.icomp; i++) {
			fprintf(stderr, "%.5s\n", csta4_.cname[i]);
		}
		fprintf(stderr, "   Phase         \twt%%");

		for (i = 0; i < cst6_.icomp; i++) {
			fprintf(stderr, "\t%.5s", csta4_.cname[i]);
		}
		fprintf(stderr, "\n");
		for (i = 0; i < numofphases; i++) {
			fprintf(stderr, "%d: %.14s\t%f", i, cxt21a_.pname[i], 100*cxt22_.props[i][16]*cxt22_.props[i][15]/cxt22_.psys[16]);
			for (j = 0; j < cst6_.icomp; j++) {
				fprintf(stderr, "\t%f", cst324_.pcomp[i][j]);
			}
			fprintf(stderr, "\n");
		}
		
		fprintf(stderr, "\n\n!! DUMP\n\n");

		for (i = 0; i < p_k5; i++) {
			for (j = 0; j < p_i8; j++) {
				fprintf(stderr, "%f  ", cxt22_.props[i][j]);
			}
			fprintf(stderr, "\n");
		}
		
		for (i = 1; i <= cst6_.icomp; i++) {
			fprintf(stderr, "%.5s \t %f \n", csta4_.cname[i-1], cxt81_.fbulk[i-1]);
		}
		
		int val = 6;
		
		/* This is the original Perple_X print output subroutine. */
		calpr0_(&val);
	}

    /* Next time we can do warm start */
    lpopt_warmstart = 1;

	return 0;
}

void freearr(void **p) {
	/* "Safe free()". Frees the memory allocated for a pointer, but
	 * survives even if called with a NULL pointer (i.e. already free'd)
	 */
	if (p != NULL && *p != NULL) {
		free(*p);
		*p = NULL;
	}
}

