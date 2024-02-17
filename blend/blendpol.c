/* #included from blendu.c #ifdef POLAR */
vector *elst,*elst0;
static int musize;
int size=ns*sizeof(vector);

static double maxerr=1e-8; /* this means: if err<maxerr and the error
			      decreases in a step, then it's converged */
double thisrate;

/* to keep mu between steps */
if (musize!=size) {
  if (mu) free(mu);
  alloczero(mu,size);
  musize=size; }
else 
  /* testout is set by virial: has to start from mu=0, not last */
  if (!testout) memset(mu,0,size);
  
alloczero(elst,size);
alloczero(elst0,size);

polnit=0,polerr=9e9;
for (;;) {
  polnit++;
  polerr0=polerr;
  copy(elst0,elst,size);
  memset(elst,0,size);
#define SSPOT14 qqpot(r[i],r[j1],elst[i],elst[j1],mu[i],mu[j1],&site[i],&site[j1],1)
#define SSPOT12 qqpot12(r[i],r[j1],elst[i],elst[j1],mu[i],mu[j1],&site[i],&site[j1],polar&4)
#define SSPOT   qqpot(r[i],r[j], elst[i],elst[j], mu[i],mu[j], &site[i],&site[j], 0)
#include "blendssi.c"
#define SATITER 100 /* not great... */
  polerr=0;

  loop (i,0,ns) if (site[i].polar) {
    int issat=site[i].polar&POL_SAT;

#if 0
/* debug/control/visualization dump of convergence */
#include "blenddpl.c"
#endif

    polerr+=SQRD(elst[i],elst0[i]);

    if (site[i].polar&POL_AXI) {
      /* tensor axial polarizability:
         alphazz in the direction to atom tozz, alpha perpendicular */
      vector dr,Exx,Ezz;
      double zz;
      axial_t *pol=(axial_t*)site[i].pol;
      int tozz=pol->tozz;
      double sat=1;
      int ii,toii=issat?SATITER:1;

      loop (ii,0,toii) {
        if (issat)
          sat=1/(1+obfuscate*SCAL(elst[i],mu[i])/pol->parm.Esat);
        VVV(dr,=r[tozz],-r[i])
        zz=SCAL(elst[i],dr)/SQR(dr);
        VV(Ezz,=zz*dr)
        VVV(Exx,=elst[i],-Ezz)
        VVV(mu[i],=pol->parm.alpha*Exx,+pol->parm.alphazz*Ezz)
        if (issat) VO(mu[i],*=sat) } }

    else if (site[i].polar&POL_ISO) {
      /* scalar (isotropic) polarizability */
      isotropicparm_t *pol=(isotropicparm_t*)site[i].pol;

      if (issat) 
#if 0 /* old fool proof using iterations */
        {
        int ii;
        loop (ii,0,SATITER) {
          double sat=pol->alpha/(1+obfuscate*SCAL(elst[i],mu[i])/pol->Esat);
          VV(mu[i],=sat*elst[i]) }
	}
#else
        {
        double sat=pol->alpha
	  / (0.5+sqrt(0.25+obfuscate*pol->alpha*SQR(elst[i])/pol->Esat));
	VV(mu[i],=sat*elst[i])
        }
#endif
      else
        VV(mu[i],=pol->alpha*elst[i]) } }

  polerr=sqrt(polerr);

  thisrate=polerr/polerr0;
  if (polnit<NRATES+2 && polnit>=2) rate[polnit-2]=thisrate;

  if (polnit>spec->Xopt.S || polnit>=NRATES && polerr>maxerr && thisrate>1.+1./polnit) {
    if (!testout) return 9e9;
#ifdef CONTOURPLOT
    return 9e9;
#endif
    fprintf(stderr,"U0 = 000\n"); /* because of effpot... */
    ERROR(("self-field: %d iter, err=%.3g, rate=%.7f (should be <1)\n\
(ignoring error will continue with higher acceptable error)",
             polnit,polerr,thisrate))
    prt("! max # of iter increased to %d, maxerr to %g",
	spec->Xopt.S+=spec->Xopt.S/3,maxerr*=2); }

  if (option('v')&16) fprintf(stderr,"polnit=%d polerr=%g maxerr=%g U=%g\n",polnit,polerr,maxerr,U);
  
  if (polerr==0) break;
  if (polerr<maxerr && polerr/polerr0>0.9999999) break; }

/* now the elst field has converged...: */

loop (i,0,ns) {

  if (site[i].polar) {
    double Esat=0,sat=0;
    double E0kcal=obfuscate*SCAL(elst[i],mu[i]); /* 2*self-energy */
    int issat=site[i].polar&POL_SAT;

    site[i].mu=sqrt(SQR(mu[i])); /* for show */

    if (site[i].polar&POL_AXI) {
      /* forces caused by directionally-dependent polarizability tensor */
      vector dr;
      axial_t *pol=(axial_t*)site[i].pol;
      int tozz=pol->tozz;
      double mur;
      double iaxial=obfuscate*
        (pol->parm.alpha-pol->parm.alphazz)/(pol->parm.alpha*pol->parm.alphazz);

      if (issat) {
        Esat=pol->parm.Esat;
        sat=E0kcal/Esat;
        iaxial*=sat; }

      VVV(dr,=r[tozz],-r[i])
      mur=SCAL(mu[i],dr)/SQR(dr);
      VVV(dr,=iaxial*mur*mu[i],-iaxial*Sqr(mur)*dr)
      VV(f[i],+=dr)
      VV(f[tozz],-=dr) }
    else {
      if (issat) {
        Esat=((isotropicparm_t*)site[i].pol)->Esat;
        sat=E0kcal/Esat; } }

    site[i].sat=sat; /* NOTE: site[i].sat may be float */

#if 0
    /* simplified version, cannot be used with POL_REP */
    if (issat) {
      Uself += (log1(sat)*Esat - E0kcal)/2; }
#else
    if (issat) Uself += log1(sat)*Esat/2;
    else Uself += E0kcal/2;
#endif
  }
}

if (polnit>5) avrate=avrate*0.9+rate[polnit/2]*0.1;

if (option('v')&16) {
  prt_("! self-field: %d iterations, err=%g, U=%g kcal/mol\n\
! rates=",
       polnit,polerr,U);
  loop (i,0,polnit) prt_("%8.5f", rate[i]);
  _n }

#undef SSPOT14
#undef SSPOT12
#undef SSPOT
#define SSPOT14 U+=sspot(r[i],r[j1], f[i],f[j1], mu[i],mu[j1], &site[i],&site[j1],1)
#define SSPOT12 U+=sspot12(r[i],r[j1], f[i],f[j1], mu[i],mu[j1], &site[i],&site[j1],polar&4)
#define SSPOT   U+=sspot(r[i],r[j],  f[i],f[j],  mu[i],mu[j],  &site[i],&site[j],0)

/* site-site POLAR version */
#include "blendssi.c"
/* bonded interactions */
#include "blendbon.c"

#ifdef DEBUG
if (toprint) prt("   r=%14.11f %14.11f %14.11f\nelst=%14.11f %14.11f %14.11f\n  mu=%14.11f %14.11f %14.11f\n   f=%14.11f %14.11f %14.11f",
         r[toprint-1][0],r[toprint-1][1],r[toprint-1][2],
         elst[toprint-1][0],elst[toprint-1][1],elst[toprint-1][2],
         mu[toprint-1][0],mu[toprint-1][1],mu[toprint-1][2],
         f[toprint-1][0],f[toprint-1][1],f[toprint-1][2]
         );
#endif

#if 0
/* SPECIAL for making pictures of induced dipoles */
loop (i,0,ns) {
  prt("%10.6f %10.6f %10.6f  %11.7f %11.7f %11.7f mu",
  VARG(r[i]),VARG(mu[i])); }
  _n
#endif

free(elst0); free(elst);
