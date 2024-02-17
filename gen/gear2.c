/* C2c gear2.C; cppnest gear2.c

  #####################################################
  USE VERSION WITH EXTENSION .C, then run C2c + cppnest
  #####################################################


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        one step t --> t+h of the Gear Predictor-Corrector method           %
%                 for second order differential equations                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 d2a/dt2 = b(t,a,da/dt) (#ifdef VELOCITY, used by MACSIMUS)

   or

 d2a/dt2 = b(t,a) (#ifndef VELOCITY, not tested recently))

   is solved, where a,b are vectors of neq real elements and the right
   hand side, b(t,a[,da/dt]), is computed by procedure rhs(b,a[,da/dt])

 Globals required:
   real t=time, h=time step
   struct gear_s gear (see simgear.h)
 typedef required:
   ToIntPtr=pointer to anything containing a contiguous vector of neq reals
 Macro or function giving pointer to the first real of the vector:
   D(a)
 If VELOCITY is #defined, rhs is assumed to depend on da/dt

 The contents of  a  and  da/dt  should not be changed except small
 corrections of the order of the method (e.g., to correct constraints).

 typedef real must be float or double
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
see gear.c for equation da/dt = b(t,a)
*/

/* this is duplicated in simopt.h - why? */
#define MAXGEARORDER 9

#ifdef POLAR
void Gear2pol(int neq,int nscf,int M,int K,ToIntPtr *a) /********** Gear2pol */
#else /*# POLAR */
void Gear(int neq,ToIntPtr *a) /*************************************** Gear */
#endif /*#!POLAR */
{
  int i,j,k,bin;
  real hh2=h*h/2,s,d,ss,*r,*hdr,*pred,*rhspred;
  static real *gc=NULL;
#ifdef VELOCITY
	real *preddr;
#endif /*# VELOCITY */

#ifdef POLAR
  int history=M+K; /* history length (M=order, K=additional lemgth) */
  static int knownhistory; /* at start less may be known */
#endif /*# POLAR */

  /* standard Gear's version gear.init=1 */
  static real GearCoef[MAXGEARORDER-2][MAXGEARORDER] = {
    /*3*/ {0.0,          1.0,         1.0},
    /*4*/ {1.0/6,        5.0/6,       1.0, 1.0/3},
    /*5*/ {19.0/120,     0.75,        1.0, 0.5,     1.0/12},
    /*6*/ {3.0/20,       251.0/360,   1.0, 11.0/18, 1.0/6,    1.0/60},
    /*7*/ {863.0/6048,   95.0/144,    1.0, 25./36,  35./144,  1./24,   1.0/360},
    /*8*/ {275.0/2016,   19087./30240,1.0,137./180, 5./16,    17./240, 1.0/120, 1./2520},
    /*9*/ {33953./259200,5257./8640,  1.0,49./60,   203.0/540,49.0/480,7.0/432, 1./720,1./20160} };

  /* Gear's velocity version gear.init=2 */
  static real GearCoefv[MAXGEARORDER-2][MAXGEARORDER] = {
    /*3*/ {0.0,          1.0,         1.0},
    /*4*/ {1.0/6,        5.0/6,       1.0, 1.0/3},
    /*5*/ {19.0/90,      0.75,        1.0, 0.5,      1.0/12},
    /*6*/ {3.0/16,       251.0/360,   1.0, 11.0/18,  1.0/6,    1.0/60},
    /*7*/ {863.0/5040,   95.0/144,    1.0, 25./36,   35./144,  1./24,    1./360},
    /*8*/ {275.0/1728,   19087./30240,1.0, 137./180, 5./16,    17./240,  1./120,  1./2520},
    /*9*/ {33953./226800,5257./8640,  1.0, 49./60,   203.0/540,49.0/480, 7.0/432, 1./720, 1./20160} };

  /* Janek's version with better energy conservation gear.init=3 */
  /* see DOI: 10.1080/00268976.2019.1674937 */
  static real GearCoefJ[MAXGEARORDER-2][MAXGEARORDER] = {
    /*3*/ {0.0,        1.0,         1.0},
    /*4*/ {0.0,        2.0/3,       1.0, 1.0/3},
    /*5*/ {1.0/12,     0.75,        1.0, 0.5,      1.0/12},
    /*6*/ {1.0/30,     23.0/36,     1.0, 11.0/18,  1.0/6,    1.0/60},
    /*7*/ {7.0/27,     95.0/144,    1.0, 25./36,   35./144,  1./24,    1./360},
    /*8*/ {43.0/840,   217.0/360,   1.0, 137./180, 5./16,    17./240,  1./120,  1./2520},
    /*9*/ {859./8640,  5257./8640,  1.0, 49./60,   203.0/540,49.0/480, 7.0/432, 1./720, 1./20160} };

#ifdef POLAR
#  define MAXK 3 /* maximum predictor length (K=k) in addition to order */
  static real A[3][MAXK+1][7]={ /* A[M-2][K] */
    /* A[0]: M=2 (ASPC); a general code is implemented with Verlet */
    { {2,   -1}, /* K=k=0 */
      {2.5, -2,     0.5},
      {2.8, -2.8,   1.2,    -0.2},
      {3,   -24./7, 27./14, -4./7, 1./14} },
    /* A[1]: M=3 (3rd order) */
    { {3,    -3,      1},
      {2.4,  -1.2,    -0.8,  0.6},
      {1.94, -0.46/3, -1.18, 0.06, 1./3},
      {0 /* n.a. */} },
    /* A[2]: M=4 (4th order, not useful) */
    { {4,     -6,     4,     -1},
      {10./3, -10./3, 5./3,  -2./3, 0},
      {2.71,  -1.34,  -1.74, 1.16,  0.71, -0.5},
      {2.15,  -0.05,  -1.80, -0.30, 0.95, 0.55, -0.5 } } };

  static real *scf[MAXGEARORDER+MAXK+1],*scfpred;
  static int nscf0=-1;
#endif /*# POLAR */

  if (gear.order<3 || gear.order>MAXGEARORDER) ERROR(("bad gear.order (option -m)"))
  if (gear.lastorder && gear.order!=gear.lastorder) ERROR(("change of gear.order = option -m not supported"))
  gear.lastorder=gear.order;

  r=D(a[0]);
  hdr=D(a[2]);
  pred=D(a[gear.order]);
#ifdef VELOCITY
	preddr=D(a[gear.order+2]);
#endif /*# VELOCITY */
  rhspred=D(a[gear.order+1]);

  if (gear.init) {
    /* performed once */
    char *info="?";

    switch (gear.init) {
      case 1: gc=GearCoefv[gear.order-3]; info="version without velocities"; break;
      case 2: gc=GearCoef [gear.order-3]; info="version with velocities"; break;
      case 3: gc=GearCoefJ[gear.order-3]; info="modification by J. Janek (lower order, better energy conservation)"; break;
      default: ERROR(("gear.init=%d is out of range",gear.init)); }

    prt("Gear m=%d %s\nCorrector coefficients: ",gear.order,info);
    loop (i,0,gear.order) {
      if (gear.C[i]>-8e9) gc[i]=gear.C[i]; }
      loop (j,0,gear.order) prt_("%14.10f",gc[j]);
      _n _n
      gear.init=0; }

#ifndef POLAR
  if (gear.order==4) {
    /* optimized code (speedup 1% for Ar fluid) */

    /* predictor for r and dr/dt */
    loop (i,0,neq) {
      pred[i] = (D(a[3]))[i] + (D(a[2]))[i] + (D(a[1]))[i] + (D(a[0]))[i];
#  ifdef VELOCITY
	preddr[i] = (3*(D(a[3]))[i] + 2*(D(a[2]))[i] + (D(a[1]))[i])/h;
#  endif /*# VELOCITY */
    }

    /* r.h.s. */
    t += h;

#  ifdef VELOCITY
	rhs(a[5],a[4],a[6]);
#  else /*# VELOCITY */
	rhs(a[5],a[4]);
#  endif /*#!VELOCITY */

    /* corrector and the rest of predictor */
    loop (i,0,neq) {
      ss = hh2*rhspred[i];
      d = ss-hdr[i]-3*(D(a[3]))[i];
      r[i] = pred[i] + d*gc[0];
#  ifdef VELOCITY
	(D(a[1]))[i] = h*preddr[i] + d*gc[1];
#  else /*# VELOCITY */
	(D(a[1]))[i] = d*gc[1] + 3*(D(a[3]))[i] + 2*(D(a[2]))[i] + (D(a[1]))[i];
#  endif /*#!VELOCITY */
      (D(a[3]))[i] = d*gc[3] + (D(a[3]))[i];
      hdr[i] = ss; }
  } else {
#endif /*# POLAR */

  /* general case - slow code */

  /* Gear's predictor */
  loop (i,0,neq) {
    s = 0;
    loop (j,0,gear.order) s += (D(a[j]))[i];
    pred[i] = s;
#ifdef VELOCITY
	s = 0; loop (j,1,gear.order) s += j*(D(a[j]))[i]; preddr[i] = s/h;
#endif /*# VELOCITY */
  }

#ifdef POLAR
  /* predictor - self-consistent field */

  if (M>4 || M<1 || K>4 || K<0 || A[M-2][K]==0)
    ERROR(("polar predictor order=%d, predictor length=%d: no such method",M,K))

  if (nscf0>=0 && nscf!=nscf0)
    ERROR(("# of polar sites changed from %d to %d",nscf0,nscf))

  loop (k,0,gear.order) scf[k]=(D(a[k]))+neq;
  scfpred=pred+neq; /* Drude after Gear's predictor at cfg[gear.order] */

  if (knownhistory==0) {
    loopto (k,gear.order,history) allocarrayzero(scf[k],nscf); /* active if scf history > gear.order */
    knownhistory=gear.order; /* C-style */
    if (knownhistory>history) knownhistory=history;
    prt_("Final POLAR predictor m=%d k=%d:\n ",M,K);
    loop (k,0,history) prt_(" %14.10f",A[M-2][K][k]); _n
    nscf0=nscf; }

  loop (i,0,nscf) {
    scfpred[i]=0;
    loop (k,0,knownhistory)
      scfpred[i] += scf[k][i]*A[M-2][knownhistory-M][k]; }
  if (knownhistory<history) knownhistory++;
#endif /*# POLAR */

    /* r.h.s. */
    t += h;

    /* NB: Drude positions after selffield (in rhs) is returned in a[gear.order] */
#ifdef VELOCITY
	rhs(a[gear.order+1],a[gear.order],a[gear.order+2]);
#else /*# VELOCITY */
	rhs(a[gear.order+1],a[gear.order]);
#endif /*#!VELOCITY */

    /* corrector and the rest of predictor */
    loop (i,0,neq) {
      ss = hh2*rhspred[i];
      d = ss-hdr[i];
      loop (j,3,gear.order) d -= j*(j-1)/2*(D(a[j]))[i];
      r[i] = pred[i] + d*gc[0];
#ifdef VELOCITY
	(D(a[1]))[i] = h*preddr[i] + d*gc[1];
#else /*# VELOCITY */
	s = d*gc[1]; loop (j,1,gear.order) s += j*(D(a[j]))[i]; (D(a[1]))[i]=s;
#endif /*#!VELOCITY */
      loop (j,3,gear.order) {
        s = d*gc[j]; bin = 1;
        loop (k,j,gear.order) {
          s += (D(a[k]))[i]*bin;
          bin = bin*(k+1)/(k+1-j); }
        (D(a[j]))[i] = s; }
      hdr[i] = ss; }
#ifndef POLAR
  } /* else optimized gear.order=4 */
#endif /*# POLAR */

#ifdef POLAR
  /* shift scf history */
  for (j=history-1; j>0; j--)
    memcpy(scf[j],scf[j-1],nscf*sizeof(scf[0][0]));
  memcpy(scf[0],(D(a[gear.order]))+neq,nscf*sizeof(scf[0][0]));
#endif /*# POLAR */
}
