/*
 * C wrapper for Perple_X v6.7.2
 * - Perple_X by Jamie Connolly (see http://www.perplex.ethz.ch/)
 * - This wrapper by Lars Kaislaniemi (lars.kaislaniemi@iki.fi)
 */

/* Changes:
 * - 2015-10-19: Changes to make compatible with Perple_X v6.7.2
 * - 2013-06-11: Changes to parameters and common block definitions to 
 *   accommodate changes in PerpleX 6.6.8.
 * - 2012-03-23: Changed p_k21 from 500000 to 1000000, according 
 *   to changes in Perple_X on 2012-03-07.
 */

/* PerpleX constant parameters (see explanation below or in Perple_X 
 * sources perplex_parameters.h. These MUST HAVE the same values as
 * in Perple_X code (perplex_parameters.h), otherwise resulting
 * in array size mismatch and overall havoc.
 */
#define p_n10 7
#define p_n11 8
#define p_n1 11
#define p_n2 12
#define p_n3 13
#define p_n4 14
#define p_n5 15
#define p_n6 16
#define p_n7 17
#define p_n8 18
#define p_n9 19
#define p_n0 30	
#define p_mst 2
#define p_mdim 8
#define p_msp (p_mdim+1)
#define p_ms1 (p_msp-1)
#define p_h5 5
#define p_h6 500
#define p_h8 200
#define p_h9 30
#define p_i6 2
#define p_i7 20
#define p_i8 28
#define p_i9 200
#define p_i10 40
#define p_i11 80
#define p_j3 4
#define p_j4 5
#define p_j5 8
#define p_j6 8
#define p_j9 160000
#define p_k0 25
#define p_k1 1000000
#define p_k2 100000
#define p_k3 500
#define p_k4 32
#define p_k5 12
#define p_k7 (p_k5+1)
#define p_k8 (p_k5+2)
#define p_k9 30
#define p_k10 300
#define p_k12 15
#define p_k13 (p_mdim*p_k1)
#define p_k14 18
#define p_k15 6
#define p_k16 50
#define p_k17 7
#define p_k18 (29*p_k1)
#define p_k19 (2*p_k5+14)
#define p_k21 1700000
#define p_k20 (p_mdim*p_k21)
#define p_k22 (p_mdim*p_k19)
#define p_k23 5
#define p_k24 1
#define p_l2 5
#define p_l3 (p_l2+2)
#define p_l5 1000
#define p_l6 1000
#define p_l7 2048
#define p_l8 10
#define p_m0 8
#define p_m1 36
#define p_m2 8
#define p_m3 3
#define p_m4 15
#define p_m6 5
#define p_m7 15
#define p_m8 9
#define p_m9 10
#define p_m10 5
#define p_m11 5
#define p_m12 4
#define p_m13 8
#define p_m14 2
#define p_m15 9
#define p_m16 5
#define p_m17 5
#define p_m18 6
#define p_nsp 16
#define p_nx 500
#define p_ny 500
#define p_kd2 (p_k8*31)

/* Maximum length of the filename passed to ini_phaseq / in
 * cst228.prject. Including trailing \0.
 */
#define p_max_filename_len 100

/* Phase name length. Not defined in PerpleX, but used as 
 * numeric constant 14. We add one to include the 
 * trailing \0 in C 
 */
#define p_pname_len 15

/* Component name length. Not defined in Perplex, but used as
 * numeric constant 5. We add one to include the 
 * trailing \0 in C
 */
#define p_cname_len 6

#define p_size_sysprops p_i8   // num of sysprops
#define p_size_phases p_k0     // (max) num of phases ("database components")
#define p_size_components p_k5 // (max) num of components


/* PerpleX (Fortran) subroutines that read initial setup
 * and initialize arrays for calculations
 */
void input1_(int *first, int *output, int *err);
void input2_(int *first);
void setau1_(int *output);
void input9_(int *first, int *output);
void initlp_();
// void vrsion_(); // not used
// void fopen1_(); // not used
// void iniprp_(); // not used


/* PerpleX subroutines for phase equilibria calculations, 
 * derivative property calculations, and printing out
 * the results (calpr0, can be used for debugging here)
 */
void lpopt0_(int *result);
void getloc_(int *itri, int *jtri, int *ijpt, double *wt, int *nodata);
void calpr0_(int *lu);
// void rebulk_(int *isstatic); // not used
// void reopt_(int *result);    // not used
// void yclos1_(double *clamda, double *x, int *is, int *jphct, int *quit); // not used


/* Functions of the wrapper */
int ini_phaseq(char *);
int ini_phaseq_lt(char *, int, char *);
int phaseq(double, double, int, double*, int*, double*, double*, double*, double*, char*, int);
void freearr(void **p);
void print_comp_order();
int get_comp_order(char **);
int number_of_components();
void spc2null(char *arr, int len);
int match_comp_order(char *, int *, int *);


/* Lookup table definitions */
#define p_lt_bulk_comp_accuracy 5
int p_use_lookup_table;
char *p_lt_filename;
char *p_lt_meltphasename;
int p_lt_dataread;
int p_lt_ncomp;
double *p_lt_comps;
double p_lt_pmin, p_lt_pmax, p_lt_tmin, p_lt_tmax;
int p_lt_nppts, p_lt_ntpts;
double *p_lt_pdata, *p_lt_tdata, *p_lt_wtdata, *p_lt_compdata;


/* Variables of the wrapper */
/* int lpopt_warmstart; */   /* set cst111_.istart to 1 before calling lpopt */




/* The common block definitions of Perple_X. Always keep the order
 * of variables inside one struct the same as in Perple_X. 
 * The total size of struct should match the common block size in Fortran.
 * This may cause troubles between platforms (e.g. integer size 32 / 64 bit?)
 */
 
/* Modified 2013-06-11: Here, in principle, we only need to define those structs
 * that are read/written in the C code. Rest of them are now commented out.
 */

extern struct {
	int iam;
} cst4_;

extern struct {
	double goodc[3];
	double badc[3];
} cst20_;

extern struct {
	int gflu, aflu, fluid[p_k5], shear, lflu, volume, rxn;
} cxt20_;

extern struct {
	double nopt[p_i10];
	int iopt[p_i10];
	int lopt[p_i10];
} opts_;
      
extern struct {
	char prject[p_max_filename_len];
	char tfname[p_max_filename_len];
} cst228_;

extern struct {
      int icomp,istct,iphct,icp;
} cst6_;

extern struct {
	char xname[p_k5][8];
	char vname[p_l2][8];
} csta2_;

extern struct {
      double cblk[p_k5];
      int jbulk;
} cst300_;

extern struct {
      double atwt[p_k0];
} cst45_;

extern struct {
      double a[p_k1][p_k5], b[p_k5], c[p_k1];
} cst313_;

extern struct {
	int ipot;
	int jv[p_l2];
	int iv[p_l2];  /* Note: In meemum.f this is defined like here,
	                * in some other *.f files of Perple_X, this is defined
	                * with iv1,iv2,iv3,... etc. Hope this is right...
	                */
} cst24_;

extern struct {
      double v[p_l2];
      double tr,pr,r,ps;
} cst5_;

extern struct {
	int jphct, istart;
} cst111_;

extern struct {
	double cp3[p_k5][p_k0];
	double amt[p_k5];
	int kkp[p_k5], np, ncpd, ntot;
} cxt15_;

extern struct {
	char pname[p_k5][p_pname_len-1];  /* p_pname_len has space for NULL */
} cxt21a_;

extern struct {
	double props[p_k5][p_i8], psys[p_i8], psys1[p_i8], pgeo[p_i8], pgeo1[p_i8];
} cxt22_; 

extern struct {
	double pcomp[p_k5][p_k0];
} cst324_;

extern struct {
	char cname[p_k5][p_cname_len-1];  /* p_cname_len has space for NULL */
} csta4_;

extern struct {
	double gtot, fbulk[p_k0], gtot1, fbulk1[p_k0];
} cxt81_;

extern struct {
	int iwt;
} cst209_;

extern struct {
      double cptot[p_k5];
      double ctotal;
      int jdv[p_k19];
      int npt;
      int fulrnk;
} cst78_;

extern struct {
      int io3,io4,io9;
} cst41_;

/* These should not be needed here:
   ... and they date from version 6.6.6(??) so check before use... 

extern struct {
	double var[p_l3], dvr[p_l3], vmn[p_l3], vmx[p_l3];	
	int jvar;
} cxt18_; 

extern struct {
	char vnm[p_l3][8];
} cxt18a_; 

      
extern struct {
	double wmach[9];
} ax02za_;
      
extern struct {
	int ldt, ldq;
} be04nb_; 
      
extern struct {
	double vlaar[p_m4][p_m3];
	int jsmod;
} cst221_; 


extern struct {
	int jmsol[p_mst][p_m4], kdsol[p_m4];
} cst142_; 
      
extern struct {
	double wg[p_m3][p_m1];
	double xmn[p_msp][p_mst];
	double xmx[p_msp][p_mst];
	double xnc[p_msp][p_mst];
    double reach;                // added for 6.6.8
	int iend[p_m4];
	int isub[2][p_m2][p_m1];
	int imd[p_mst][p_msp];
	int insp[p_m4];
	int ist[p_mst];
	int isp[p_mst];
	int isite;
	int iterm;
	int iord;
	int istot;
	int jstot;
	int kstot;
} cst108_;

extern struct {
	char mname[p_m4][8];
} cst18a_;
      
extern struct {
	double xy[p_k1][p_mdim], y[p_k1][p_mst][p_ms1];
	int ntot, npairs;
} cst159_;
      
extern struct {
	int refine;
} cxt26_;
      
extern struct {
	char cfname[100];
} cst227_;
     
extern struct {
	int hkp[p_k21], mkp[p_k19];
} cst72_;
      
extern struct {
	double g2[p_k21], cp2[p_k21][p_k5];
    int jphct;          // order changed for 6.6.8 (typo in prev version?)
} cxt12_;
     
extern struct {
	double g[p_k1];
} cst2_;

extern struct {
	double thermo[p_k10][p_k4], uf[2], us[p_h5];
} cst1_;

extern struct {
	double cp[p_k1][p_k5];
} cst12_;

      
extern struct {
	double emod[p_k10][p_k15];
	int smod[p_h9], pmod[p_k10];
	int iemod[p_k10],kmod;
} cst319_;

      
extern struct {
	int imaf[p_i6];
	int idaf[p_i6];
} cst33_;

extern struct {
	char title[4][162];
} csta8_;

extern struct {
	double vmax[p_l2];
	double vmin[p_l2];
	double dv[p_l2];
} cst9_;

extern struct {
	char fname[p_h9][10];
} csta7_;



extern struct {
	int icp2;
} cst81_;

extern struct {
	char zname[5];
} cst209a_;

extern struct {
	char tcname[p_k0][5];
	char xcmpnt[p_k0][5];	
} csta9_;

extern struct {
	double buf[5];
} cst112_;


extern struct {
	int ivfl;
} cst102_;

extern struct {
	int jfct, jmct, jprct;
} cst307_;

extern struct {
	double c0, c1, c2, c3, c4, c5;
	int iind, idep;
} cst316_;

extern struct {
      double dblk[p_k5][3];
      double cx[2];
      int icont;
} cst314_;

extern struct {
      double dlnfo2,elag,gz,gy,gx;
      int ibuf,hu,hv,hw,hx;
} cst100_;

extern struct {
      double ctrans[p_k5][p_k0];
      int ictr[p_k5];
      int itrans;
} cst207_;

extern struct {
      int iff[2];
      int idss[p_h5];
      int ifug,ifyn,isyn;
} cst10_;

extern struct {
      int isoct;
} cst79_;


extern struct {
	double x3[p_msp][p_mst][p_k21];
} cxt16_; 
      
extern struct {
	int ifp[p_k1];
} cxt32_;

extern struct {
	int istg[p_h9], ispg[p_mst][p_h9], imlt[p_mst][p_h9], imdg[p_h9][p_mst][p_ms1];
} cxt6i_; 

extern struct {
	int idasls[p_k3][p_k5], iavar[p_k3][3], iasct, ias;
} cst75_;
      
extern struct {
	int ipa[p_k2], ibulk;
} cst74_; 

extern struct {
	int igrd[p_l7][p_l7];
} cst311_;

extern struct {
	double xcoor[p_k18];
	int icoor[p_k1];
} cxt10_; 

extern struct {
	double bg[p_k2][p_k5];
} cxt19_; 

extern struct {
	double mus[p_k2][p_k8];
} cst48_;

extern struct {
	double hsb[4][p_i8];
	int hs2p[6];
} cst84_;
      
      

extern struct {
      int ifct,idfl;
} cst208_;

extern struct {
      int ixct,iexyn,ifact;
} cst37_;

extern struct {
	char exname[p_h8][8];
	char afname[2][8];
} cst36_;

extern struct {
	int ids[p_h6][p_h5];
	int isct[p_h5];
	int icp1,isat,io2;
} cst40_;

extern struct {
	double sxs[p_k13], exces[p_m3][p_k1];
	int ixp[p_k1];
} cst304_;

extern struct {
	int jend[p_k12][p_h9];
} cxt23_;

extern struct {
	double y[p_m4], x[p_m4], pa[p_m4], p0a[p_m4], z[p_msp][p_mst], w[p_m1];
} cxt7_;

extern struct {
	int lorder[p_h9], lexces[p_h9], llaar[p_h9], lrecip[p_h9];
} cxt27_;

extern struct {
      int lstot[p_h9], mstot[p_h9], nstot[p_h9], ndep[p_h9], nord[p_h9];
} cxt25_;

extern struct {
	int ksmod[p_h9], ksite[p_h9], kmsol[p_mst][p_m4][p_h9], knsp[p_h9][p_m4];
} cxt0_;

extern struct {
      int isec,icopt,ifull,imsg,io3p;
} cst103_;


extern struct {
	int ivarrx[p_k2];
	int ivarip[p_k2];
	int isudo,ivar;
} cst62_;

extern struct {
      int jlow,jlev,loopx,loopy,jinc;
} cst312_;

extern struct {
      int oned;
} cst82_;


extern struct {
	int pindex, tindex;
	int usv;
} cst54_;

extern struct {
      int hcp,idv[p_k7];
} cst52_;

extern struct {
      char eoscmp[2][8];
} cst98_;
      
            
extern struct {
	int ncol;
	int nrow;
	int fileio;
} cst226_;
      





extern struct {
      double vnumu[p_k10][p_i6];
} cst44_;

extern struct {
	double zcoor[p_k20];
	int jcoor[p_k21], jkp[p_k21], jcoct;
} cxt13_; 

extern struct {
	double ctot[p_k1];
} cst3_;

extern struct {
	double ycoor[p_k22];
	int lcoor[p_k19], lkp[p_k19];
} cxt14_;
      
extern struct {
	char names[p_k1][8];
} cst8_;

extern struct {
      char name[8];
} csta6_;

extern struct {
      int ic[p_k0];
} cst42_;

extern struct {
      double comp[p_k0],tot;
      int icout[p_k0],ikind,icmpn,ieos; // changed for 6.6.8
} cst43_;

extern struct {
      int cl[p_k0];
      char cmpnt[p_k0][5], dname[80];
} csta5_;

extern struct {
      double tm[p_m6][p_m7],td[p_m8];
      int ilam,jlam,idiso,lamin,idsin; // changed for 6.6.8
} cst202_; 

extern struct {
      int ipoint,imyn;
} cst60_;

extern struct {
      double  mcomp[p_k0][p_k16];
      int nmak;
      int mksat[p_k16];
      char mknam[p_k17][p_k16][8];
} cst333_;

extern struct {
      double mkcoef[p_k17][p_k16], mdqf[p_k17][p_k16];
      int mkind[p_k17][p_k16], mknum[p_k16];
} cst334_;

extern struct {
      int make[p_k10];
} cst335_;

extern struct {
      int eos[p_k10];
} cst303_;

extern struct {
      int ikp[p_k1];
} cst61_; 

extern struct {
	int iam[p_k1], jam[p_k1], tloop, ploop;
} cst55_;

extern struct {
	double mu[p_k8];
} cst330_;

extern struct {
	int jtest, jpot;
} debug_;

extern struct {
	double dcp[p_h8][p_k5], soltol;
} cst57_;

*/
