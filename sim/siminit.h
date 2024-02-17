/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   initializing tables of species, sites, potentials  molecule descriptors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

void initNo(void);
void makemolgol(int nmol);
void initpot(void);
void initmolecules(int corr);
void initgammas(void);
double cutcorr(double LJcutoff,int corr); 
/* legacy: param LJcutoff not needed in COOK, should be removed */
double setL(vector L,double rho);
double initcutoff(double cutoff,vector L);

int initfix(char *fn);
#ifdef ANCHOR
void initanchor(char *fn,int init,double drmax);
#endif /*# ANCHOR */
void fixsites(ToIntPtr X);
int isfixed(int i);
void checkfixed(void);
void printmasses(void);

#ifdef WIDOM
double Widomcutcorr(int spreal,int spvirt,double LJcutoff);
#endif /*# WIDOM */

#if PARALLEL
void initparallel(void);
#endif /*# PARALLEL */

#ifdef ECC
double rescalecharges(int ecc,double epsf);
#endif /*# ECC */
