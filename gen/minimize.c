/* different minimization methods */

/* 1/2001 update: ret. value of Newton-Raphson and Cartesian metric to
   calc. the difference
   3/2013 verbose switch added
   1/2019 small changes
*/

#define DISTURB 1  /* for amoeba : ndef: orig. algorithm
                                   1:improved,
                                   2:improved and verbose */

#include "ground.h"
#include "rndgen.h"
#include "minimize.h"

extern volatile int sig; /* to be set if ^C is pressed */

#define GOLD 0.6180339887498948

/*** this is mainly for SD and CG ***/

/* step to calculate numerical derivatives for SD and CG */
/* see also NRminimize */

static REAL h=1e-4; /* step for SD+CG, also search radius for MCminimize */

static int sN,verbose=1;
static sqfun *sff;
static REAL difh;

static REAL Upot(REAL *r,REAL *f)
/* forces from potential */
{
  int i;
  REAL plus,minus,old;

  loop (i,0,sN) {
    old=r[i]; r[i]+=difh;
    plus=sff(r);
    r[i]-=2*difh;
    minus=sff(r);
    f[i]=(minus-plus)/(2*difh);
    r[i]=old; }

  return sff(r);
}

/***************************************
steepest descent and conjugate gradients
rewritten from project blend (see blendmin.c)
inefficient here because the derivatives are calculated numerically
which is expensive
sd=# of SD steps
cg=# of CG steps
***************************************/

void minimize(int N, REAL *xx, sqfun ff, REAL D,int sd,int cg) /***** minimize */
{
  int size=N*sizeof(REAL),it,i;
  REAL *r,*f,*r0,*f0,*G,*H;
  REAL U,U0;

  difh=D;
  sN=N; sff=ff;

  ralloc(r0,size);
  ralloc(f0,size);

  /*** initial configuration from the molecule ***/
  copy(r0,xx,size);
  U0=Upot(r0,f0);

  put(U0);

  ralloc(r,size);
  ralloc(f,size);
  ralloc(G,size); /* conjugate grad only */
  ralloc(H,size); /* conjugate grad only */

  /***********/

  h*=3;
  U0=Upot(r0,f0); /* for sure...*/

#if 0
  /* re-consider step to prevent large moves when overlap etc. */
  U=0; loop (i,0,N) Max(U,fabs(f0[i]))
  if (U*h>0.5) h=0.5/U; /* max move in the 1st step = 0.5 A */
  if (U*h<0.001) h=0.001/U; /* min move in the 1st step = 0.001 A */
#endif

  if (sd) {
    /*** steepest descent method ***/
    if (verbose) prt("! %i steps of steepest descent:",sd);


    for (it=1; !sig; it++) {

      do {
        h*=0.55; /* magic value: works fine in a narrow valley */

        if (h<1e-12 || it>=sd) goto sddone;

        loop (i,0,N) r[i]=r0[i]+h*f0[i];
        U=Upot(r,f);

        } while (U>U0);

      copy(r0,r,size);
      copy(f0,f,size);
      U0=U;

/*      prt("%i %g %g",it,U,h);*/

      h*=2.5; /* magic: h should not shrink to zero */ }

    }

 sddone:

  if (cg) {
    /* this is conjugate gradient, rewritten from Numerical Recipes */
    REAL dgg,gg,h0;
    /* brrrr! U0=fp r0=p  f0=xi */
    int j,maxj=24;

    if (verbose) prt("! %i steps of conjugate gradients:",cg);
    /* from init, known r0 f0 U0 */
    h*=3;
    copy(G,f0,size);
    copy(H,f0,size);

    for (it=1; !sig && it<=cg; it++) {

      /***
      now we are looking for h>0 so that r=r0+h*f0 is minimum
      in most cases close to minimum we need only two evaluations of
      forces in one step
      ***/
      gg=0;
      loop (i,0,N) gg+=Sqr(f0[i]);
      if (gg==0) break;

      h0=h;
      for (j=0; ;j++) {
        loop (i,0,N) r[i]=r0[i]+h*f0[i];
        U=Upot(r,f);
        if (U<U0) break;
        h*=0.5;
        if (j>=maxj) {
          /* unsuccessful - reset CG (should occur rarely!) */
          U0=Upot(r0,f0);
          if (maxj==20) goto cgdone; /* unsuccessful again! */
          h=h0; maxj=20;
          copy(G,f0,size);
          copy(H,f0,size);
          goto OK; }
        }
      maxj=5;

      dgg=0;
      loop (i,0,N) dgg+=f[i]*f0[i];
      dgg/=gg;
      if (dgg<0.8) {
        /* try one step of the secant method - efficient close to minimum */
        h*=1/(1-dgg);
        loop (i,0,N) r0[i]+=h*f0[i];
        U0=Upot(r0,f0);
        if (U0<U) goto OK; }

      /***
      not better or too extrapolated - we must be happy with r,f,U
      this occurs once a while
      ***/

      h*=4; maxj=8; /* might be caused by too low h */
      U0=U;
      copy(f0,f,size);
      copy(r0,r,size);

    OK:
/*      prt("%i %g %g",it,U,h);*/

      gg=dgg=0;
      loop (j,0,N) {
        gg+=Sqr(G[j]);
        dgg+=(f0[j]-G[j])*f0[j]; }
      dgg/=gg;

      copy(G,f0,size);
      loop (j,0,N) H[j]=G[j]+dgg*H[j];
      copy(f0,H,size); }

    }

 cgdone:

  U0=Upot(r0,f0);
  copy(xx,r0,size);
  put(U0);

  release(r0);
}

/**************************************
AMOEBA: function ff(xx[N]) is minimized
Improved C++ version by Z.Wagner
**************************************/

int amoeba(int N, REAL *xx, sqfun ff, REAL eps,int maxit) /***** amoeba */
{
  int N1=N+1, jbest,jworst,j,i, iter = 0, ngood;
  int logscale=0; /* 0: xx[] both positive and negative; 1: of one sign */
#ifdef DISTURB
  int disturb=0;
  REAL qdisturb=-1.0/N;
#endif
  unsigned vecsize = N * sizeof(REAL);
  REAL *xguess,*xc,*xworst,*xbest,*y,**x;
  REAL ybest, yworst, yguess, yc, rms, q,qq, rmsold=pow(fabs(eps),GOLD);

  if (verbose) prt("amoeba %i steps eps=%g",maxit,Val(eps));

  if (eps<0) eps=0,logscale++;

  alloc(xguess,vecsize);
  alloc(xc,vecsize);
  alloc(xworst,vecsize);
  alloc(xbest,vecsize);
  alloc(y,vecsize+sizeof(REAL));

  alloc(x,N1*sizeof(REAL*));
  loop (j,0,N1) alloc(x[j],vecsize);

  /* initial simplex */

  loop (j,0,N1) {
    copy(x[j], xx, vecsize);
    if (j<N) {
      if (logscale) x[j][j] *= exp(rmsold);
      else  x[j][j] += rmsold; }
    y[j] = ff(x[j]); }

  do { /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

    /* Determination of min & max & centroid calculation */

    copy(xc, x[0], vecsize);
    jbest = jworst = 0;
    ybest = yworst = y[0];
    loop (j,1,N1) {
      if (y[j] < ybest) ybest = y[jbest=j];
      if (y[j] > yworst) yworst = y[jworst=j];
      loop (i,0,N) xc[i] += x[j][i]; }
    copy(xworst, x[jworst], vecsize); /* max */
    copy(xbest, x[jbest], vecsize); /* min */
    loop (i,0,N) xc[i] = (xc[i]-xworst[i])/N; /* centroid=center excl. worst */

    /* Reflection of the worst through xc */

    rms=0;
    loop (i,0,N) {
      xguess[i] = xc[i] + (q = xc[i] - xworst[i]);
      rms+=Sqr(q); }
    rms=sqrt(rms); /* to test error */
    yguess = ff(xguess);

    /* ngood = # of points from the simplex that are worse than the new guess */

    ngood=0;
    loop (j,0,N1) ngood+=yguess<=y[j];

    if (ngood>N1/3+1) { /* the guess is fairly good */

      if (ngood==N1) {
        /* the guess is the best: try expand once again
           (it does not pay to repeat it) */
        loop (i,0,N) xc[i] += (xc[i] - xworst[i])*1.6;
        if ((yc = ff(xc)) < ybest) {
          copy(x[jworst], xc, vecsize);  y[jworst] = yc; }
        else
          yguess = yworst; /* just to satisfy the following if */ }

      if (yguess >= ybest) {
        /* replace the worst point by the (improved) guess */
        copy(x[jworst], xguess, vecsize);  y[jworst] = yguess; } }

    else { /* the guess is too bad */

      if (yguess <= yworst) {
        /* the guess is still better than the worse case: replace it */
        copy(xworst, xguess, vecsize);  yworst = yguess; }

      /* contraction of the worst point towards the centroid */
      loop (i,0,N) xguess[i] = 0.5*(xworst[i] + xc[i]);

      if ((yguess = ff(xguess)) < yworst) {
        /* successfully contracted */
        copy(x[jworst], xguess, vecsize);  y[jworst] = yguess; }

      else {
        /* failed: shrink the whole simplex towards the best point */
        /* this occurs very rarely */
        loop (j,0,N1) {
          q=GOLD; qq=1-q;
          loop (i,0,N) x[j][i] = q*x[j][i] + qq*xbest[i];
          y[j] = ff(x[j]); } }

#ifdef DISTURB
      /* try to increase the volume of the simplex */
      /* of limited value in certain special cases */
      loop (i,0,N) xc[i]=0;
      loop (j,0,N1) loop (i,0,N) xc[i] += x[j][i];
      disturb=(disturb+1)%N1;
      q=qdisturb/N1;
      qq=1-q;
      q/=N1;
      loop (i,0,N) xc[i]=q*xc[i]+qq*x[disturb][i];
      yc=ff(xc);
      if (yc<y[disturb]) {
        copy(x[disturb],xc,vecsize); y[disturb]=yc; qdisturb*=2;
#if DISTURB==2
        if (verbose) prt("+");
#endif
      }
#if DISTURB==2
      else if (verbose) prt("-");
#endif
      qdisturb*=0.7;
#endif
    }

    iter++; rmsold=(N*rmsold+rms)/N1;
  }
  while ( sig==0 && rmsold>eps && iter<maxit );
  /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

  copy(xx, x[ybest < yworst ? jbest : jworst], vecsize);

  loop (j,0,N1) free(x[j]);
  free(x); free(y); free(xbest);
  free(xworst); free(xc); free(xguess);

  return -(iter>=maxit)-2*sig;
}


/***********************
improved MC minimization
************************
for complicated functions of many variables which do not contain long valleys,
esp. to obtain initial estimate for Newton-Raphson

parameters:
  N=vector length
  xx=vector to minimize
  ff=function to minimize
  mc=# of steps
  acc=multiplicative factor of the search radius on success, must be acc>1
  rej=multiplicative factor of the search radius on failure,
      must be 0<rej<1 and acc*rej>1
  trace=weight to sum previous points, must be 0<trace<1
      trace=0: no extrapolation in narrow valley
      trace->1: long trace, efficient in very narrow valley
global parameter:
  h=typical step size --> search radius

if acc<=0 then default acc=1.3
if rej<=0 then default rej=0.98
if trace<=0 then default trace=0.9

Method:
Points are shot at random with approximate Gaussian distribution of radius
(standard deviation) h around `shifted point'
If the new point is better, it is accepted and the radius is increaded (acc)
If the new point is worse, it is rejected and the radius is decreaded (rej)
even step: `shifted point' = point
odd  step: `shifted point' = extrapolated point in minus direction of
                             a weighted centre of previous points
*/

void MCminimize(int N, REAL *xx, sqfun ff, int mc, REAL acc,REAL rej,REAL trace)
{
  int size=N*sizeof(REAL),i;
  REAL *x,*c;
  REAL U,U0;

  if (acc<=0) acc=1.3;
  if (rej<=0) rej=0.98;
  if (trace<=0) trace=0.9;

  if (h<1e-8) h=1e-7;
  h*=Sqr(N); /* adjusting search radius size to SD/CG methods */
  if (verbose) prt("MC minimize %i steps accept=%.4f reject=%.4f trace=%.4f h=%g",
      mc,Val(acc),Val(rej),Val(trace),Val(h));
  if (rej>=1 || acc*rej<1 || trace>=1)
    ERROR(("bad values of MC minimization parms"))
  alloc(x,size);
  alloc(c,size);
  copy(x,xx,size);
  copy(c,xx,size);
  sN=N; sff=ff;

  U0=sff(xx);

  while (mc--) {

    /* x:=xx+Gauss(h) */
    loop (i,0,N) x[i]=xx[i]+(rnd()-rnd()+rnd()-rnd())*h;

    if (mc&1) {
      /* extrapolation */
      U=rnd();
      loop (i,0,N) x[i]+=U*(xx[i]-c[i]); }

    U=sff(x);
    if (U<U0) {
      /* accept */
      copy(xx,x,size); U0=U; h*=acc;
      loop (i,0,N) c[i]=c[i]*trace+x[i]*(1-trace); }
    else
      /* reject */
      h*=rej;

    if (h<1e-14 || sig) break;
  }

  free(x);
  free(c);
  if (verbose)
    prt("MC minimize finished: fmin=%g h=%g -> %g",Val(U0),Val(h),Val(h)/Sqr(N));
  h/=Sqr(N);
}

/**********************************
this is Newton-Raphson minimization, rewritten from old Pascal code 1989
Method: Newton-Raphson, derivatives calculated numerically
If fails, minimum is looked for in the direction of the obtained NR point
If not successful, simply the best from all calculated points is taken
as the next iteration
eps=accuracy
D=step for numerical derivatives
nr=max # of steps
**********************************/

#include "gjlineq.h"

static REAL *b,*g,*xmin,*y,*sx;
static REAL fx,fxmin;
static int size;

#define FandMIN { fx=sff(y); if (fx<fxmin) { fxmin=fx; copy(xmin,y,size); } }

static REAL fdir(REAL q)
{
  /* =fun(x+q*g), i.e. q from the point x in the direction g */
  int i;

  loop (i,0,sN) y[i]=sx[i]+q*g[i];
  FandMIN

  return fx;
}

REAL NRminimize(int N,REAL *x,sqfun ff,REAL D,REAL eps,int nr)
/* looks for a vector x which minimizes ff(x).
  The minimum is returned as a function value */
{
  REAL **a;
  REAL *d;
  int isg,i,j,lq;
  REAL f0,f1,f2,s,r,q;

  size=sizeof(REAL)*N;

  ralloc(b,size);
  ralloc(g,size);
  ralloc(xmin,size);
  ralloc(y,size);
  ralloc(d,size);
  loop (i,0,N) d[i]=D; /* hoo! */
  ralloc(a,N*sizeof(REAL*));
  loop (i,0,N) ralloc(a[i],size);


  sN=N; sff=ff; sx=x;
  fxmin=1e35;
  do {
    copy(y,x,size);
    loop (i,0,N) loop (j,i+1,N) {

      /* cross 2nd derivatives */
      y[i]+=d[i];     y[j]+=d[j];     FandMIN s=fx;
                      y[j]=x[j]-d[j]; FandMIN s-=fx;
      y[i]=x[i]-d[i];                 FandMIN s+=fx;
                      y[j]=x[j]+d[j]; FandMIN
      a[i][j]=(s-fx)/d[i]/d[j]/4; a[j][i]=a[i][j];
      y[i]=x[i]; y[j]=x[j]; }

    FandMIN f0=fx; isg=N; s=f0*2;

    loop (i,0,N) {
      /* 1st and 2nd derivatives */
      y[i]-=d[i]; FandMIN f1=fx;
      y[i]=x[i]+d[i]; FandMIN f2=fx; y[i]=x[i];
      a[i][i]=(f2-s+f1)/Sqr(d[i]);
      b[i]=(f1-f2)/2/d[i];   /* minus 1st deriv */
      if (a[i][i]>0) isg--; }

#undef FandMIN

    r=1;
    lq=!GJlineq(N,a,b,g);

    if (lq && (isg==0 /*q-form is +*/ )) {
      q=0; f1=fdir(r); }
    else {
      if (!lq) {
        if (verbose) prt("lineq KO");
        copy(g,d,size); };
      q=-r; f1=f0; f0=fdir(q); }

    if (f0<f1) {
      f2=f0; f0=f1; f1=f2; q+=r; r=-r; }
    /* now f0>f1 */
  again: f2=fdir(r+r+q);
    if (f2<f1) { f1=f2; r+=r; f2=fdir(r+r+q); goto again; };
    if ((s=f2-f1-f1+f0)!=0) if (fabs(s=(f0-f2)/s)<5) f0=fdir(q+r+r/2*s);
    r=0;
    loop (i,0,N) r+=Sqr(x[i]-xmin[i]);
    r=sqrt(r);
    copy(x,xmin,size);
    if (verbose) prt("signature(0=convex)=%i   fmin=%g  dx in 1 iter=%g",
        isg,Val(fxmin),Val(r));
  } while (r>eps && !sig && --nr);

  release(b);
  /*.....return fxmin;*/
  return r+isg;
}


/*** generic method: ***/

static int oneMinimize(int method, /* -1: steepest descent
			     0: amoeba
			     1: conjugate gradient
			     2: Newton-Raphson
			     3: Monte Carlo */
	      int N,      /* # of variables */
	      REAL *xx,   /* vector of variables */
	      sqfun ff,   /* function to minimize */
	      int maxit,  /* (max) number of steps */
	      REAL eps,   /* accuracy (irrelevant for MC) */
	      REAL D,     /* step to calc. num. deriv. (SD,NR,CG) */
	      REAL par)   /* `nonsphericity' valey param. D>=1 (MC) */
{
  REAL x;

  switch (method) {
    case 0: return     amoeba(N,xx,ff,eps,maxit);
    case 1:          minimize(N,xx,ff,D,0,maxit); break;
    case -1:         minimize(N,xx,ff,D,maxit,0); break;
    case 2: return NRminimize(N,xx,ff,D,eps,maxit)>eps;
    case 3: if (par==0) par=10;
      x=pow(par,-0.3333333333);
      MCminimize(N,xx,ff,maxit*2*N,1+0.7*x,1-0.3/par,1-x*x);
      break;
    default: ERROR(("unknown minimization method %d",method))  }

  return 0;
}

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
	      REAL par)  /* `nonsphericity' valey param. D>=1 (MC)
                            par<=0: repeats w. extrapolation */
/* 0 returned on success */
{
  REAL x,e;
  REAL *x0,ff0,ffx;
  int i,rc;
  static REAL rnd=0;
  int iter=0,exiter=0;

  if (maxit==0) maxit=1000;
  if (maxit<0) verbose=0,maxit=-maxit;

  rc=oneMinimize(method,N,xx,ff,maxit,eps,D,fabs(par));

  if (sig) return 1;

  if (par>0) return rc;

  alloc(x0,sizeof(REAL)*N);
  copy(x0,xx,sizeof(REAL)*N);

  ff0=ff(x0);

  for (;;) {
    rc=oneMinimize(method,N,xx,ff,maxit,eps,D,fabs(par));

    ffx=ff(xx);
    if (ffx==ff0) {
      rc=2; goto ex; }

    if (ffx>ff0) {
      x=ffx; ffx=ff0; ff0=x;
      loop (i,0,N) {
        x=xx[i]; xx[i]=x0[i]; x0[i]=x; }
      if (++iter>6) goto ex; }
    else
      iter=0;

    /* PATCH: dirty end */
    if (exiter++>100) {
      rc=8; goto ex; }

    if (sig) goto ex;

    rnd+=GOLD; if (rnd>1) rnd-=1;
    e=(1+rnd)*ffx/(ff0-ffx); /* ideal factor is 2 */
    rnd+=GOLD; if (rnd>1) rnd-=1;
    /*.....if ((rnd+=GOLD)>1) rnd-=1; DOES NOT WORK WITH C++ */
    Min(e,5+rnd)
    Max(e,-2)

    if (verbose) prt("##### minimization extrapolation value=%f",Val(e));

    loop (i,0,N) {
      x=xx[i];
      xx[i]+=e*(xx[i]-x0[i]);
      x0[i]=x; }

    ff0=ffx; }

 ex: free(x0);
  return rc;
}

struct keys_s MinimizeKey[]={
  {"SD",-1},{"steepest descent",-1},
  {"amoeba",0},
  {"CG",1},{"conjugate grad",1},
  {"NR",2},{"Newton-Raphson",2},
  {"MC",3},{"Monte Carlo",3},
  {NULL,0}
};
