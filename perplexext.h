/*
 * C wrapper for Perple_X v6.6.8
 * - Perple_X by Jamie Connolly (see http://www.perplex.ethz.ch/)
 * - This wrapper by Lars Kaislaniemi (lars.kaislaniemi@iki.fi)
 */

/* Changes:
 * - 2013-06-11: Changes to parameters and common block definitions to 
 *   accommodate changes in PerpleX 6.6.8.
 * - 2012-03-23: Changed p_k21 from 500000 to 1000000, according 
 *   to changes in Perple_X on 2012-03-07.
 */

/* Phase name length. Not defined in PerpleX, but used as 
 * numeric constant 14. We add one to include the 
 * trailing \0 in C 
 */
#define p_pname_len 15


/* Functions of the wrapper */
int ini_phaseq(char *);
int phaseq(double, double, int, double*, int*, double**, double**, double**, char**, int);
void freearr(void **p);
void print_comp_order();

/* Variables of the wrapper */
static int lpopt_warmstart;   /* set cst111_.istart to 1 before calling lpopt */


