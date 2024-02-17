#include "prec.h"

typedef REAL sqfun(REAL *x);

#define MINIMIZE_SD -1
#define MINIMIZE_AMOEBA 0
#define MINIMIZE_CG 1
#define MINIMIZE_NR 2
#define MINIMIZE_MC 3

/* generic method: */
int Minimize(int method, /* -1: steepest descent
			     0: amoeba
			     1: conjugate gradient
			     2: Newton-Raphson
			     3: Monte Carlo */
	      int N,     /* # of variables */       
	      REAL *xx,  /* vector of variables */
	      sqfun ff,  /* function to minimize */
	      int maxit, /* (max) number of steps; negative = quiet */
	      REAL eps,  /* accuracy (irrelevant for MC) */
	      REAL D,    /* step to calc. num. deriv. (SD,NR,CG) */
	      REAL par); /* `nonsphericity' valey param. D>=1 (MC)
                            par<=0: repeats w. extrapolation */
/* 0 returned on success */

void minimize(int N, REAL *xx, sqfun ff, REAL D,int sd,int cg);
void MCminimize(int N, REAL *xx, sqfun ff, int mc, REAL acc,REAL rej,REAL trace);
int amoeba (int N, REAL *xx, sqfun ff, REAL eps, int maxit);

extern struct keys_s MinimizeKey[];
