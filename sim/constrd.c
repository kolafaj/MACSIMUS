/*
  Lagrangian constraint and SHAKE
  global #defines:
    SHAKE: 4 versions (since V2.6d should be #defined)
    POLAR: selffield calculations
    SHEAR: shear viscosity support
    VERLET: with SHAKE: kin.energy algorithm
  local #defines:
    PRECOND DEBUG
  removed (see old+misc/scf+lambda.c):
    LAMBDA
*/

/*
#define DEBUG
  debug energy terms, pressure, scaling */

#include "ground.h"
#include "sds.h"
#include "statics.h"
#include "simglob.h"
#include "simcg.h"
#include "norm.h"
#include "constrd.h"
#include "forces.h" /* the same for LINKCELL */
#include "ewald.h" /* POLAR: epspol */
#include "cputime.h"
#include "units.h"
#include "simils.h"
#include "asksig.h"
#include "simdef.h" /* because of sitedef_t */
#define REAL double
#include "gjlineq.h"
#undef REAL
#ifdef LINKCELL
#  include "linklist.h"
#endif /*# LINKCELL */

#ifdef POLAR
#  include "intermac.h"
#  if POLAR&2
#    include "log1.c"
#  endif /*# POLAR&2 */
#  ifndef FREEBC
#    include "elst.h"
#  endif /*# FREEBC */
#endif /*# POLAR */

#ifdef SHEAR
#  include "units.h"
#endif /*# SHEAR */

#ifdef ANCHOR
#    define FORCEUNIT (massunit*lengthunit/Sqr(timeunit))
#    define TORQUEUNIT (massunit*Sqr(lengthunit/timeunit))

static void printanchor(double unit,double *var,int mask) /***** printanchor */
/*
   unit = one of 1,FORCEUNIT,TORQUEUNIT
   var = vector or double* to record
   mask = vector components to record (ASCII file anchor.f: also zero components recorded)
*/
{
  int k;

  if (anchor.f) fprintf(anchor.f," %12g %12g %12g ",unit*var[0],unit*var[1],unit*var[2]);

  loop (k,0,3) if ((1<<k)&mask) {
    if ((-option('k')) & 3)
      if (anchor.i>=anchor.col) ERROR(("printanchor: number of recorded variables exceeds the precalculated value"))
    if ((-option('k')) & 1) StaAdd(anchor.rec[anchor.i].info,unit*var[k]);
    if ((-option('k')) & 2) anchor.rec[anchor.i].sum+=unit*var[k];
    anchor.i++; }
}

static void normalizer(vector cm) /****************************** normalizer */
/* local version of normalize for 1 vector, see also normalize in norm.c */
{
  return;
#  ifndef FREEBC
  if (cm[0]<-box.Lh[0]) cm[0]+=box.L[0]; else if (cm[0]>box.Lh[0]) cm[0]-=box.L[0];
  if (cm[1]<-box.Lh[1]) cm[1]+=box.L[1]; else if (cm[1]>box.Lh[1]) cm[1]-=box.L[1];
#    ifndef SLIT
  if (cm[2]<-box.Lh[2]) cm[2]+=box.L[2]; else if (cm[2]>box.Lh[2]) cm[2]-=box.L[2];
#    endif /*# SLIT */
#  endif /*# FREEBC */
}
#endif /*# ANCHOR */

#if PRESSURETENSOR&PT_VIR
static double PTvirc[PT_DIM];  /* of constraint forces */
#endif /*# PRESSURETENSOR&PT_VIR */

static double blogT(double E,double T) /****************************** blogT */
/* 1/2*ln(E/kT), but always kept in [-1,1] */
{
  const double exp2=7.38905609893065022723;

  if (E>T*exp2) return 1;
  if (E*exp2<T) return -1;

  return log(E/T)/2;
}

#ifdef ECC
#  include "ecc.c"
#else /*# ECC */
#  define ecc(X) /* void */
#endif /*#!ECC */

static void measurePconstraints(void) /***************** measurePconstraints */
/*
   - finish pressure calculations and constraint statistics
   - if measure=1: once per cycle, with NPT every cycle
   - returns En.P, En.trPt, En.Pcfg (according to `virial')
   - En.PdV calculated elsewhere
*/
{
  int i;

  //  StaSet(0,lag.ierr,2,lag.in); ?
  StaSet(0,lag.err,2,lag.n);
  if (option('v')&4) {
    StaAdd("Pvirc [Pa]",En.virc*(Punit/DIM)/box.V);
    StaAdd("Pvir pair excl. el. part [Pa]",En.vir*(Punit/DIM)/box.V); }

#ifdef POLAR
  StaAdd("Eel(excl. self-term)",En.el);
  StaAdd("Eel(incl. self-term)",En.el+En.self);
  En.vir+=En.virc+En.el-En.self*2; /* P=En.vir/3V (in main.c) */
  /* ... patched below for GAUSSIANCHARGES */

  /* cf. also VVV(En.Pvir,-=rpol[i],*fpol[i]) below */
  StaAdd("Eself",En.self);
  /* 9/2001: -En.self*2 added instead of the virial of forces on aux
     charges that incorrectly did not include Ewald corrections.
     WARNING: not tested for axial and saturated polarizability !!! */
#else /*# POLAR */
  ecc(0); /* ECC only */
  StaAdd("Eel",En.el);
#  ifdef ECC
  /* NB: NOT En.ECC_virnel=En.vir+En.virc if combined with ECC for point dipoles
     "full" version of Pecc -> correct value for short enough dipoles only */
  En.ECC_virnel=En.vir;
#  endif   /*# ECC */
  En.vir+=En.virc+En.el; /* P=En.vir/3V (in main.c) */
  /* ... patched below for GAUSSIANCHARGES */
#endif /*#!POLAR */

#if PRESSURETENSOR&PT_KIN
  {
    static int pass=0;
    double x=SUM(En.Pkin)/2-En.kin/2; /* still doubled */

    if (!pass && fabs(x)>1000) {
      prt("WARNING: Ekin=%g K  tr(Pkin)/2=%g K  dif=%g K\n\
          (Pkin uncorrected here; more warnings suppressed)",
          En.kin/2,SUM(En.Pkin)/2,x);
      pass=1; }
  }
  loop (i,0,PT_DIM) En.Pkin[i]*=No.Pkinq/box.V;
#  if PRESSURETENSOR&PT_MOL
  {
    static int pass=0;
    double x=SUM(En.PKin)/2-En.kin_tr/2; /* still doubled */

    if (!pass && fabs(x)>1000) {
      prt("WARNING: Ekininter=%g K  tr(PKin)=%g K  dif=%g K\n\
          (Pkin uncorrected here; more warnings suppressed)",
          En.kin_tr/2,SUM(En.PKin)/2,x);
      pass=1; }
  }
  loop (i,0,PT_DIM) En.PKin[i]*=No.Pkinq/box.V;
#  endif /*# PRESSURETENSOR&PT_MOL */
#endif /*# PRESSURETENSOR&PT_KIN */

#if PRESSURETENSOR&PT_VIR
  loop (i,0,PT_DIM) {
    En.Pvir[i]+=PTvirc[i];
    En.Pvir[i]/=box.V;
#  if PRESSURETENSOR&PT_KIN
    En.Ptens[i]=En.Pvir[i]+En.Pkin[i];
#  endif /*# PRESSURETENSOR&PT_KIN */
  }

  /* NB: there was patch En.vir=SUM(En.Pvir)/(Punit/box.V); here */
  /* to add (again) comparison of pressures ? */
  //    StaSet(0,2,2,lag.n);
  //    StaAdd("Pvirerr [Pa]",x);

#  ifdef SLAB
  //    StaSet(0,lag.ierr,2,lag.in); ?
    En.Ptv=(En.Pvir[0]+En.Pvir[1]-2*En.Pvir[2])/DIM;
    StaSet(0,2,2,lag.n);
    StaAdd("Pvir xy-zz [Pa]",En.Ptv*Punit);
#    if PRESSURETENSOR&PT_KIN
    StaSet(0,lag.err,2,lag.n);
    En.Pt=((En.Pvir[0]+En.Pvir[1]-2*En.Pvir[2])
          +(En.Pkin[0]+En.Pkin[1]-2*En.Pkin[2]))/DIM;
    StaAdd("Pt [Pa]",En.Pt*Punit);
#    endif /*# PRESSURETENSOR&PT_KIN */
#  endif /*# SLAB */
    StaSet(0,2,2,lag.n);

    /* pressure based on the trace of the pressure tensor */
    En.trPt=SUM(En.Ptens)/DIM+En.corr/Sqr(box.V)
#ifdef ECC
      + En.ECC_Pcorr
#endif
      ;

#endif /*# PRESSURETENSOR&PT_VIR */

#if PRESSURETENSOR&PT_ANY
    loop (i,0,PT_DIM) {
      char *xx="xx\0yy\0zz\0yz\0zx\0xy"+i*3;

#  if PRESSURETENSOR&PT_VIR
      StaAdd(string("Pvir.%s [Pa]",xx),En.Pvir[i]*Punit);
#  endif /*# PRESSURETENSOR&PT_VIR */
#  if PRESSURETENSOR&PT_KIN
      StaAdd(string("Pkin.%s [Pa]",xx),En.Pkin[i]*Punit);
#  endif /*# PRESSURETENSOR&PT_KIN */
#  if PRESSURETENSOR&PT_MOL
      StaAdd(string("PKin.%s [Pa]",xx),En.PKin[i]*Punit);
#  endif /*# PRESSURETENSOR&PT_MOL */
#  if PRESSURETENSOR&PT_MOM
      if (i<3) {
        StaAdd(string("PKin2.%s [Pa^2]",xx),En.PKin2[i]*Sqr(Punit*No.Pkinq/box.V));
        StaAdd(string("PKin3.%s [Pa^3]",xx),En.PKin3[i]*Cub(Punit*No.Pkinq/box.V)); }
#  endif /*# PRESSURETENSOR&PT_MOM */
#  if (PRESSURETENSOR&PT_ANY) == PT_ANY
      StaAdd(string("Ptens.%s [Pa]",xx),En.Ptens[i]*Punit);
#  endif /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */
    }
#endif /*# PRESSURETENSOR&PT_ANY */

  if (No.f_in) En.T_in=En.kin_in/No.f_in;
  else En.T_in=0;
  if (No.f_tr>0) En.T_tr=En.kin_tr/No.f_tr;
  else En.T_tr=0;

/* presure based on elst_virial = -elst_energy */
  En.P=(En.kin*No.Pkinq+En.vir)/DIM/box.V + En.corr/Sqr(box.V)
#ifdef ECC
      + En.ECC_Pcorr
#endif
    ;
  
  /* @91  (NB: No.Pkinq=1 assumed for NPT/MTK) */
  switch (virial) {
    case 1: En.Pcfg=En.P; break;
    case 2: En.Pcfg=En.PdV; break;
#if PRESSURETENSOR&PT_VIR                                                     
    case 3: En.Pcfg=En.trPt; break;
#endif
    default: ERROR(("pressure type virial=%d not supported",virial))
  }
} /* measurePconstraints */

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             constraint dynamics                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formulae:
% ^^^^^^^^^
% M[a,b] = SUM{i} grad[i]c[a].grad[i]c[b] / m[i]
%
% G[a] = SUM{i} p[i].grad[i]c[a] / m[i]         (Hamilton, removed in V 2.6d)
%
% G[a] = SUM{i} F[i].grad[i]c[a] / m[i]
%      + SUM{i,j} V[i].grad[i]grad[j]c[a].V[j]  (Lagrange formalism)
%
% where constraint c[a] = (r[ia]-r[ja])**2/2 - BondLength**2/2
% ( ia = si[a].pair[0], ja = si[a].pair[1], BondLength**2 = si[a].bondq )
%
% F[i] = -grad[i]U is force, V[i] is velocity, . is dot (scalar) product
%
% Eq. for Lagrange multipliers: (note: g(Lagrange) = dg(Hamilton)/dt )
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% SUM{b} M[a,b] g[b] = G[a]
%
% Eqs. of motion:
% ^^^^^^^^^^^^^^^
% m[i] dr[i]/dt = p[i] - SUM{a} g[a]grad[i]c[a]
% dp[i]/dt      = F[i] + SUM{a}SUM{j}g[a]grad[i]grad[j].dr[i]/dt (Hamilton)
%
% m[i] ddr[i]/dt/dt = F[i] - SUM{a} g[a]grad[i]c[a]              (Lagrange)
%
% Using matrix M to correct constraints: (see Lcorrect)
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Dr[i] = SUM{b} e[b] 1/m[i]*grad[i]c[b]
% where  e  solves
% SUM{b} M[a,b] e[b] = -c[a]
%
% Velocity constraints in the Lagrangian formalism:
% Dv[i] = SUM{b} ve[b] 1/m[i]*grad[i]c[b]
% where  ve  solves
% SUM{b} M[a,b] ve[b] = -SUM{i}v[i]*grad[i]c[a]
%
% Thermostats (Lagrange only):
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Friction:
%   ddr[i]/dt/dt -= v[i] ln(Tkin/T)/2/tau.T
% Nose:
%   ddr[i]/dt/dt -= v[i] dlogs/dt
%   ddlogs/dt/dt = (Tkin/T-1)/tau.T^2
% Since ddr[i]/dt/dt is a `friction' force, the corresponding term appears
% also in G[]
%
% #define PRECOND :
%   uses preconditioned matrix M for conjugate gradients:
%     SUM{b} M*[a,b] g*[b] = G*[a]
%   where M*=diag.M.diag, G*=diag.G and the solution is g=diag.g*
%   ( . is matrix multiplication)
%   and diag=diag(M)^(-1/2)
%   It appears that about 2 iterations are spared but the control operations
%     cost about the same.
%   Not implemented for constraints corrections (option -c1)
%   (sometimes referred to as the Cholesky preconditioning)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*** Lagrange equations of motion *********************** constraintdynamics */
void constraintdynamics(ToIntPtr B, ToIntPtr A, ToIntPtr V)
{
  int n,a,b,nc,ns,sp,i,j,k;
  static int irhs=0;
  vector *v,*f,*r,*ra=NULL,*rp=A->rp;
#ifdef PRECOND
  double *diag;
#endif /*# PRECOND */
  molecule_t *mn;
  vector auxv;
  /* warning: even #ifdef FLOAT, the constraint calculations are in double */
  double *g=NULL,*aux,*G=NULL;
  double **gg=NULL;
  double z, Vlogs=0, Vlogs_mol=0;
  M_t *M,*Msp;
  siteinfo_t *si;
  vector v_mol;
  double m_mol,mi;

  En.virc=0;

  if (No.pred) allocarray(gg,No.pred);

  if (No.maxc) {
    rallocarray(ra,No.maxc); /* = r[i]-r[j] for constraint a=(i,j) */
    rallocarray(G,No.maxc);  /* rhs of constraint equation */
#ifdef PRECOND
    rallocarray(diag,No.maxc);
#endif /*# PRECOND */
    if (No.pred==0) rallocarrayzero(g,No.maxc); }

  /* ! all kin.energies (En/kin, En.kin_tr, En.kin_in) are doubled here !*/
  switch (thermostat) {
    case T_NOSE: /* Nose */
      Vlogs=V->logs; break;
    case T_BERENDSEN: /* friction thermostat: must know En.kin=twice kin.en. in advance! */
      if (En.kin==0) Vlogs=0;
      else Vlogs=blogT(En.kin/No.f,T)/tau.T;
      break;
    case T_BUSSI: /* canonical sampling through velocity rescaling */
      if (En.kin==0) Vlogs=0;
      // EXPERIMENTAL!!!
      else ERROR(("not implemented"))
      break;
    case T_LANGEVIN: case T_LANGEVIN_CM:
      Vlogs=0.5/tau.T;
      break;
    case T_TR: /* friction thermostat - intermolecular only */
      Vlogs_mol=blogT(En.kin_tr/No.f_tr,T)/tau.T;
      break;
    case T_IN: /* friction thermostat - intramolecular only */
      Vlogs=blogT(En.kin_in/No.f_in,T)/tau.T;
      break;
    case T_FRICTIONS: { /* decoupled friction thermostats, T_tr/T_in=T_tr_in */
      double T_tr=(No.f_in+No.f_tr)*T/(T_tr_in*No.f_in+No.f_tr);
      double T_in=(No.f_in+No.f_tr)*T/(No.f_in+No.f_tr/T_tr_in);

      Vlogs=blogT(En.kin_in/No.f_in,T_in)/tau.T;
      Vlogs_mol=blogT(En.kin_tr/No.f_tr,T_tr)/tau.T; } }

  Vlogs_mol-=Vlogs;

  /*** constraints: computing M and P ***/
  loop (n,FROM,No.N) { /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
    mn=molec+n;
    sp=mn->sp;
    ns=mn->ns;
    si=spec[sp]->si;
    r=rof(mn,rp);
    v=rof(mn,V->rp);
    f=rof(mn,B->rp);

    if ( (nc=mn->nc) ) {

      if (No.pred) {
        loop (k,0,No.pred) gg[k]=gs[k]+mn->ig;
        g=gg[0]; }

      M=Msp=Ms[option('c')&1 ? n : sp];

      loop (a,0,nc) {
        i=si[a].pair[0]; j=si[a].pair[1];
        VVV(ra[a],=r[i],-r[j])

        M->Mab = M->imass*SQR(ra[a]); /* a==b */
        // <V3.1l:   M->Mab = (si[i].imass+si[j].imass)*SQR(ra[a]); /* a==b */
        while ( (b=(++M)->b) < a )
          M->Mab = M->imass*SCAL(ra[a],ra[b]);

        VVV(auxv, =si[i].imass*f[i], -si[j].imass*f[j])
        G[a] = SCAL(auxv,ra[a]);
        VVV(auxv, =v[i], -v[j])
        G[a] += SQR(auxv);
        if (tau.T) G[a] -= 2*Vlogs*SCAL(auxv,ra[a]);

        /*** predictor of Lagrange multipliers g ***/
        if (No.pred>0) {
          z=0;
          loop (i,0,No.pred) z+=binom[i]*gg[i][a];
          gg[No.pred-1][a]=z; }
      } /* a */

      if (No.pred) g=gg[No.pred-1];
      else memset(g,0,nc*sizeof(g[0]));

#ifdef PRECOND
      /* find the diagonal */
      if (option('c')&1) Error("-c&1 not allowed with PRECOND");
      /* ...should pass diag to Lcorrect (norm) and use it in the same way */
      M=Msp;
      loop (a,0,nc) {
        diag[a]=1./sqrt(M->Mab); /* a==b */
        while ( (b=(++M)->b) < a ); }

      /* calculating diag.M.diag and diag.P */
      M=Msp;
      loop (a,0,nc) {
        M->Mab=1; /* a==b */
        G[a]*=diag[a];
        g[a]/=diag[a]; /* to fix the predictor of g */
        while ( (b=(++M)->b) < a ) M->Mab*=diag[a]*diag[b]; }
#endif /*# PRECOND */

      /*** solving the eqs. for Lagrange multipliers ***/
      constrit[sp].nit[LAGR_MULT] += omega==0
        ? conjgrad  (nc,g,Mvec,Msp,G, irhs<No.pred && eps>=1?1e-12:eps)
        : directiter(nc,g,Msp,G,omega,irhs<No.pred && eps>=1?1e-12:eps);

#ifdef PRECOND
      loop (a,0,nc) g[a]*=diag[a];
#endif /*# PRECOND */

    } /* if (nc) */

    /*** equations of motion ***/
    loop (a,0,nc) {
      i=si[a].pair[0]; j=si[a].pair[1];

      /*** forces equiv to constraint motions, virial ***/
      VV(auxv, =g[a]*ra[a])
      VV(f[i], -=auxv)
      VV(f[j], +=auxv)
      En.virc -= g[a]*si[a].bondq;
#if PRESSURETENSOR&PT_VIR
      En.Pvir[0] -= g[a]*ra[a][0]*ra[a][0];
      En.Pvir[1] -= g[a]*ra[a][1]*ra[a][1];
      En.Pvir[2] -= g[a]*ra[a][2]*ra[a][2];
#  if PRESSURETENSOR&PT_OFF
      En.Pvir[3] -= g[a]*ra[a][1]*ra[a][2];
      En.Pvir[4] -= g[a]*ra[a][0]*ra[a][2];
      En.Pvir[5] -= g[a]*ra[a][0]*ra[a][1];
#  endif /*# PRESSURETENSOR&PT_OFF */
#endif /*# PRESSURETENSOR&PT_VIR */
    } /*a*/

    /*** thermostats & acceleration from forces (kin. energy now in rhs) ***/

    if (thermostat==0 || tau.T==0)
      /*** no thermostat ***/
      loop (i,0,ns)
        VO(f[i], *=si[i].imass)

    else if (thermostat<=T_NOSE || thermostat==T_BUSSI)
      /*** Nose, NPT, or simple Berendsen (friction) thermostat ***/
      loop (i,0,ns)
        VVV(f[i], =si[i].imass*f[i], -Vlogs*v[i])

    else if (thermostat==T_LANGEVIN)
      /*** atom-based Langevin thermostat ***/
      /* WARNING: OK with long tau.T only */
      loop (i,0,ns) {
        double x=sqrt(T*si[i].imass/(tau.T*h)); /* TO BE OPTIMIZED */

        VVVO(f[i], =si[i].imass*f[i], -Vlogs*v[i], +x*rndgauss()) }

    else if (thermostat==T_LANGEVIN_CM) {
      /*** molecule-based Langevin thermostat ***/
      /* WARNING: OK with long tau.T only */
      double x=sqrt(T/(spec[sp]->mass*tau.T*h)); /* TO BE OPTIMIZED */
      vector rf;

      VO(rf,=x*rndgauss())
      loop (i,0,ns)
        VVVV(f[i], =si[i].imass*f[i], -Vlogs*v[i], +rf) }

    else if (thermostat>=T_TR) {
      /*** decoupled inter/intra thermostats ***/
      VO(v_mol,=0)
      m_mol=0;
      loop (i,0,ns) {
        mi=si[i].mass;
        m_mol+=mi;
        VV(v_mol,+=mi*v[i]) }

      VO(v_mol,/=m_mol) /* velocity of CM of a molecule */

      loop (i,0,ns)
        VVVV(f[i], =si[i].imass*f[i], -Vlogs*v[i], -Vlogs_mol*v_mol) }

  } /* n >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

  if (No.pred) {
    /*** dirty trick: pointer shuffling ***/
    k=No.pred-1; aux=gs[k];
    while (k) { gs[k]=gs[k-1]; k--; }
    gs[0]=aux;
    free(gg); }

  if (No.maxc) release(ra);

  /* final energies and statistics: repeated in Shake */
  box.V=PROD(box.L);
  if (measure==1) measurePconstraints();

  /* NB: No.Pkinq=1 assumed for NPT/MTK */
  /* this En.Pcfg not needed anyway since Gear+NPT/MTK not implemented (yet);
     cf. thermo.c */
  En.Pcfg = ((En.kin*No.Pkinq+En.vir)/DIM+En.corr/box.V)/box.V
#ifdef ECC
        +En.ECC_Pcorr
#endif /*# ECC */
        ;

  En.T=En.kin/No.f; /* NB: En.kin still doubled here ! */
  En.kin/=2;
#ifdef POLAR
  En.pot += En.self;
#endif /*# POLAR */
  En.pot += En.el; /* do not move to forces(), but is in scforces() */
  En.U = En.Unc = En.kin+En.pot;
#ifndef FREEBC
  En.U += En.corr/box.V;
#endif /*# FREEBC */
  En.tot = En.U; /* cutoff corrections added in V2.7 because of the barostat */

  /*.....if (measure) fprintf(stderr,"%12.2f %12.2f %12.2f %12.2f %12.2f %12.2f\n", En.tot,En.kin,En.pot,En.el,En.intra,En.self);*/

  if (thermostat==T_NOSE) {
    En.ext = No.f*T*(A->logs+Sqr(tau.T)/2*Sqr(Vlogs));
    En.tot += En.ext;
    En.logs=A->logs;
    En.dlogs=Vlogs;
    B->logs = (En.T/T-1)/Sqr(tau.T); }

  irhs++;

  if (FROM) /* to keep 1st FROM molecules on place */
    loop (n,0,FROM) {
      mn=molec+n;
      ns=mn->ns;
      memset(rof(mn,V->rp),0,ns*sizeof(vector));
      memset(rof(mn,B->rp),0,ns*sizeof(vector)); }

  CPUtime("constr");

#ifdef POLAR
  if (option('p')/10%10==0) {
    if (init>=0 && init<3)
      ERROR(("Gear + Car-Parrinello NOT SUPPORTED"))
    else
      WARNING(("Gear + Car-Parrinello: BRUTE FORCE INITIALIZATION"))

    loop (n,0/*FROM ?*/,No.N) {
      mn=molec+n;
      sp=mn->sp;
      ns=mn->ns;
      si=spec[sp]->si;
      r=polarrof(mn,rp);
      v=polarrof(mn,V->rp);
      f=polarrof(mn,B->rp);
      loop (i,0,ns) VVVO(r[i],=v[i],=f[i],=0); } }
#endif /*# POLAR */

} /* constraintdynamics */


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   SHAKE                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef SHAKE
#  error SHAKE is not #defined, should be one of {-1,-2,1,2}
#endif /*# SHAKE */

void Shake(double eps, double prob) /********************************* Shake */
/*
  Verlet integration (equivalent to leap-frog)
  with SHAKE constraint dynamics for bonds and similar external constraints
  eps = accuracy (if eps<1), see parameter epsc
      = # of Shake iteration sweeps over all sites (if eps>=1, |SHAKE|=1 only)
  prob = probability of changing velocity for Andersen and Maxwell-Boltzmann
         thermostats, see verlet.c
  omega (see omegac) = relaxation parameter (omega=1: direct iteration)
  data:
    rof(molecule*,cfg[0]->rp) ->r = r(t)
    rof(molecule*,cfg[1]->rp) ->r1= r(t)-r(t-h) = h*v(t-h/2)
    rof(molecule*,cfg[2]->rp) ->p = forces; used to keep r(t+h) during shaking

  POLAR:
    polarrof(molecule*,cfg[0]->rp) = dr(t) (rel.aux.site of ind. dipole)
    polarrof(molecule*,cfg[1]->rp) = dr(t-h)
    polarrof(molecule*,cfg[2]->rp) = used in forces calculation

  Compile-time versions:
    SHAKE>0: simplified formula, equivalent upto h^2,
             usu more efficient for complex systems
    SHAKE<0: `exact' formula using the scalar product,
             perhaps better for diatomics
    |SHAKE|=1: update in sweeps (generally better)
    |SHAKE|=2: update for unprecise bonds only
               (better for e.g. mixture of complex molecule + diatomics)

  History:
    2012 MTK barostat
    2010 Nose-Hoover, several VERLET versions
    2008 pressure tensor, automatic omega optimization
    2003 ANCHOR constraints, |SHAKE|=1 only
    2001 Maxwell/Andersen thermostat
    1995 revisited
    1992 first version
*/
{
  vector rij,pij;
  int al,n,i,j,sp,nc,ns=-1,it,maxit0=-16*log(eps*0.1);
  molecule_t *mn;
  vector *r,*p,*rp=cfg[0]->rp,*rp1=cfg[1]->rp,*r1,*rp2=cfg[2]->rp,*v;
  siteinfo_t *si;
  double rxp,bondq,z,ze,mi,mj,sq;
  vector v_mol;
  vector v_mol1,v_mol2;
  double m_mol;
  double hh=h*h,Thh=T*hh,hhh=hh/2,hscale;
  double omega; /* = omegac/2 */
  int CMbasedthermostat=thermostat==T_ANDERSEN_CM || thermostat==T_MAXWELL_CM;

  /* Nose-Hoover; arrays are allocated forever
     indices: [0] = logs (= denoted also xi in the algorithms)
              [1,2,3] = lambda (vector)
              [RPOFFSET..] = rp (RPOFFSET=4)
  */
  static double **Vhist; /* Vhist[0]=X(t)-X(t-h) (pointer passed); X=logs,lambda,r
                            Vhist[1]=X(t-h)-X(t-2h) ... */
  static int k; /* predictor extra length, actual value, growing up to No.pred */
  static double *b; /* coefficients */
  static double *Vpred; /* Vpred[0]=h*d X(t)/dt, predicted */
  static double Vfactor;
  static vector *vpred; /* vpred=pointer to predicted h*v[0] */
  static vector scale={1,1,1}; /* bond length scaling in T_NPT and tau.rho<1 */

#if SHAKE==1 || SHAKE==-1
  int ifit=eps>=1;
  int maxit = (int)eps;
  double maxsq;
#else /*# SHAKE==1 || SHAKE==-1 */
  int *moved=NULL,*moving,done;
#endif /*#!SHAKE==1 || SHAKE==-1 */

#ifdef ANCHOR
  vector cf={0,0,0};
  int notinplace=0;

  anchor.i=0; /* this is called first */
#endif /*# ANCHOR */

  t+=h;
  En.virc=0;

  if (thermostat>=T_NPT) {
    /* virial pressure needed for the T_NPT algorithm at every step */
    measure=1;
    /* NB: box.L[] is the box throughout cook; cfg[0]->lambda[] is here
       calculated, integrated, and finally converted back to box.L[] */
    // @1
    loop (j,0,DIM) cfg[0]->lambda[j]=log(box.L[j]); }

  VV(box.Lh,=0.5*box.L)
  box.V=PROD(box.L); // @2

  box.cq=Sqr(box.cutoff);
  En.kin=En.kin_tr=0;

  if (option('c')&2) {
    En.r1=constrainterror(cfg[0],cfg[1]);
    En.v1=vconstrainterror/h; }
  depend_r(cfg[0],0);

#ifdef POLAR
  /* pred1.c is old version with -@%=A[0] -_%=A[1]: WARNING - options changed */
#  include "pred2.c"
  scf.nit=0;
  if (option('p')/10%10==0)
    mechpolar(cfg[2],cfg[0]);
  else do { /* iterate self-field */
#endif /*# POLAR */
      //    fprintf(stderr,"%d====MD==== %g %g\n",scf.nit,cfg[0]->rp[0][0],cfg[2]->rp[0][0]);
    zeroEn();
    forces(cfg[2],cfg[0]); // @3
#ifdef POLAR
  } while (!selffield(cfg[2],cfg[0],scf.eps,scf.omega,1));
  //  fprintf(stderr,"%d====MD==== %g %g\n",scf.nit,cfg[0]->rp[0][0],cfg[2]->rp[0][0]);

  StaAdd("polar no of iter",En.polnit=scf.nit);
  StaAdd("polar Drude maxdr",En.polmaxdr=scf.maxdr);
#endif /*# POLAR */

#if PRESSURETENSOR&PT_VIR
  if (measure) memset(PTvirc,0,sizeof(PTvirc));
#endif /*# PRESSURETENSOR&PT_VIR */

  depend_f(cfg[0],cfg[2]);

#ifdef SHEAR
  Shear(cfg[2],cfg[0],cfg[1],h);
#endif /*# SHEAR */

  /*** independent pressure measurement via <dU/dV>: pass 1 */
  if (measure) measureP(1);

#ifdef POLAR
  /* SCF test vs. full iteration */
  if (Eext.isB) ERROR(("Magnetic field (el.B) not implemented in POLAR"))
  if (measure) testSCF();
#endif /*# POLAR */

  /* Time-Reversible Velocity Predictor + scaling predictor for T_NPT */
#include "trvpscale.c"

  /* <<<<< loop over all molecules w. Verlet+SHAKE starts here <<<<< */

  loop (n,FROM,No.N) {
    mn=molec+n;
    sp=mn->sp;
    r=rof(mn,rp);
    p=rof(mn,rp2);
    v=rof(mn,vpred);
    r1=rof(mn,rp1);
    si=spec[sp]->si;
    nc=mn->nc;
    omega=constrit[sp].omega; /* HALF omegac, the relaxation parameter */

    if (Eext.isB)
      loop (i,0,mn->ns) {
        vector F;
        
        VECT(F,v[i],Eext.Bh)
        /* BUG: virial ignored */
        VV(p[i],+=si[i].charge*F) }
    
#if SHAKE>0 /* simplified and in general more efficient algorithm */
#  define CALCULATE_RXP rxp=bondq;
#else /* more complicated algorithm using angle between p & r */ /*# SHAKE>0 */
#  define CALCULATE_RXP rxp=SCAL(rij,pij) \
         if (rxp<bondq*1e-6) ERROR(( \
           "constraint failure\n" \
           "it=%d %d[%d-%d]:(%.3f %.3f %.3f)->(%.3f %.3f %.3f)", \
           it,n,i,j,pij[0],pij[1],pij[2],rij[0],rij[1],rij[2]));
#endif /*#!SHAKE>0 */

#if SHAKE==1 || SHAKE==-1 /* constraints updated in sweeps */
    ns=mn->ns;

#  ifdef ANCHOR
   /* actions independent on anchoring:
      - measure force on center-of-mass and torque
      - set zeros to cf of groups */
#    include "anchorm.c"
#  endif /*# ANCHOR */

/* Verlet step (for 1 molecule) + possible Andersen/Maxwell thermostat */
#  include "verlet.c"

    if (!ifit) maxit=maxit0*(nc+4);

    it=0;
    do { it++;
      maxsq=0;

#  ifdef ANCHOR
      /* anchor site: */
#    include "anchors.c"
      /* anchor CM of group(s): */
#    include "anchorg.c"
#  endif /*# ANCHOR */

      loop (al,0,nc) {
        i=si[al].pair[0]; j=si[al].pair[1];
        sq=bondq=si[al].bondq;
        VVV(pij,=p[i],-p[j])
        VV(pij,*=scale) // @7 INEFFICIENT - should optimize
        sq -= SQR(pij);
        Max(maxsq,fabs(sq/bondq))

        VVV(rij,=r[i],-r[j])
        CALCULATE_RXP
        mi=si[i].imass; mj=si[j].imass;
        sq *= omega/(mi+mj);
#  include "virc.c"
        sq /= rxp;
        VO(rij,*=sq)
        VV(p[i],+=mi*rij)
        VV(p[j],-=mj*rij) } /*al*/

    } while ((ifit || maxsq>eps) && it<maxit);

#  ifdef ANCHOR
    /* constrain CM, axes, and print site/cm constraint force */
#    include "anchora.c"
#  endif /*# ANCHOR */

   if (!ifit && it==maxit) {
     int ii;

     ERROR(("Shake: species=%d molecule=%d: too many iterations (%d)\n\
*** maxsq=%g eps=%g omega/2=%g\n\
*** (selecting (i)gnore in interactive run will dump the molecule)",
            mn->sp,n,it,
            maxsq,eps,omega))
     prt("%d\n",mn->ns);
     loop (ii,0,ns) prt("%4s %g %g %g",sitedef[spec[mn->sp]->si[ii].st].name,VARG(r[ii])); }

#else /*# SHAKE==1 || SHAKE==-1 */

/* SHAKE=+-2: info about moved sites is kept, eps<1 only */
#  ifdef ANCHOR
#    error ANCHOR not supported with this version of SHAKE
#  endif /*# ANCHOR */

    if (ns<mn->ns) {
      ns=mn->ns;
      if (moved!=NULL) { free(moved); free(moving); }
      allocarray(moved,ns);
      allocarray(moving,ns); }

/* Verlet step (for 1 molecule) + possible Andersen/Maxwell thermostat */
#  include "verlet.c"

    loop (i,0,ns) { moving[i]=0; moved[i]=1; }

    done=it=0;
    do {
      done=1;

      if (it>maxit0*(nc+4)) {
        int ii;

        ERROR(("Shake: species=%d molecule=%d: too many iterations (%d)\n\
*** maxsq=%g eps=%g omega/2=%g\n                                        \
*** (selecting (i)gnore in interactive run will dump the molecule)",
               mn->sp,n,it,
               maxsq,eps,omega))
        prt("%d\n",mn->ns);
        loop (ii,0,ns) prt("%4s %g %g %g",sitedef[spec[mn->sp]->si[ii].st].name,VARG(r[ii])); }

      loop (al,0,nc) { it++;
        i=si[al].pair[0]; j=si[al].pair[1];

        if (moved[i] || moved[j]) {

          sq=bondq=si[al].bondq;
          VVV(pij,=p[i],-p[j])
          VV(pij,*=scale) // @7 INEFFICIENT - should optimize
          sq -= SQR(pij);

          if (fabs(sq) > eps*bondq) {

            VVV(rij,=r[i],-r[j])
            CALCULATE_RXP
            mi=si[i].imass; mj=si[j].imass;
            sq *= omega/(mi+mj);
#  include "virc.c"
            sq /= rxp;
            VO(rij,*=sq)
            VV(p[i],+=mi*rij)
            VV(p[j],-=mj*rij)

            moving[i]=moving[j]=1; done=0; }
          } /* i and j moved */
        } /*a*/
      loop (i,0,ns) { moved[i]=moving[i]; moving[i]=0; }
    } while (!done);

#endif /*#!SHAKE==1 || SHAKE==-1 */
    //    prt("%d %g %g %g",n,VARG(rof(molec+1,rp)[0]));

    if (nc) constrit[sp].nit[COR_BONDS] += it;

#ifdef ANCHOR
  } /* Shake n-loop over all molecules ends here (ANCHOR)  */

/* constrain CM of not-anchored molecules - is global, must be out of n-loop */
#  include "anchorc.c"

  /* another n-loop needed to finish */
  loop (n,FROM,No.N) {
    mn=molec+n;
    ns=mn->ns;
    sp=mn->sp;
    r=rof(mn,rp);
    p=rof(mn,rp2);
    v=rof(mn,vpred);
    r1=rof(mn,rp1);
    si=spec[sp]->si;
#endif /*# ANCHOR */

/*
  Kinetic energy, pressure tensor (if measure),
  and final cfg[0], cfg[1] update (except rescaling for T_NPT and tau.rho<0)
  Several versions using different formulas for velocity:
  recommended: VERLET=3 (leap-frog)
  ususally less accurate: VERLET=1 (velocity Verlet)
  more options: 2 (harmonic), 5 (best equipartition (?)), 6 (aver.~Beeman)
  development/not recommended: VERLET=0,4,9,30
*/
#define V1 r1[i] /* [r(t)-r(t-h)]/h */
#if VERLET==0
#  include "shakev0.c"
#elif VERLET==1
#  include "shakev1.c"
#elif VERLET==2
#  include "shakev2.c"
#elif VERLET==3
#  include "shakev3.c"
#elif VERLET==4
#  include "shakev4.c"
#elif VERLET==5
#  include "shakev5.c"
#elif VERLET==6
#  include "shakev6.c"
#elif VERLET==30
#  include "shakev30.c"
#else /*#!VERLET==0!VERLET==1!VERLET==2!VERLET==3!VERLET==4!VERLET==5!VERLET==6!VERLET==30 */
#  error "wrong value of compile-time switch VERLET (one of 0,1,2,3,4,5,9,30 expected)"
#endif /*#!VERLET==0!VERLET==1!VERLET==2!VERLET==3!VERLET==4!VERLET==5!VERLET==6!VERLET==30 */
#undef V1
  } /* ANCHOR: end of 2nd n-loop, not ANCHOR: of the Shake n-loop */

#     include "anchork0.c"

  //prt("%d %g %g %g",n,VARG(rof(molec+1,rp)[0]));
  cfg[0]->dep=0;

#ifdef POLAR
#  include "shakecp.c"
#endif /*# POLAR */

#if SHAKE==2 || SHAKE==-2
  if (moved!=NULL) { free(moved); free(moving); }
#endif /*# SHAKE==2 || SHAKE==-2 */

  /* En.kin etc = TWICE the respective kinetic energies: */
#if VERLET==3
  hscale = 0.5/hh;
#elif VERLET==1
  hscale = 0.25/hh;
#else /*#!VERLET==3!VERLET==1 */
  hscale = 1./hh;
#endif /*#!VERLET==3!VERLET==1 */

  if (h==0) hscale=0;

  En.kin *= hscale;
  En.kin_tr *= hscale;

#if VERLET==2
#  define HSCALE (i<DIM?hscale:hscale/2)
#else /*# VERLET==2 */
#  define HSCALE hscale
#endif /*#!VERLET==2 */

#if PRESSURETENSOR&PT_KIN
  loop (i,0,PT_DIM) En.Pkin[i]*=HSCALE;

#  if VERLET==30
  hscale*=4;  /* ? */
#  endif /*# VERLET==30 */

#  if PRESSURETENSOR&PT_MOL
  loop (i,0,PT_DIM) En.PKin[i]*=HSCALE;
#    if PRESSURETENSOR&PT_MOM
  loop (i,0,DIM) {
    En.PKin2[i]*=Sqr(hscale);
    En.PKin3[i]*=Cub(hscale); }
#    endif /*# PRESSURETENSOR&PT_MOM */
#  endif /*# PRESSURETENSOR&PT_MOL */
#endif /*# PRESSURETENSOR&PT_KIN */

  En.kin_in=En.kin-En.kin_tr;

#ifdef SHEAR
  /*
    systematic flow kinetic energy removed
    NOTE: in some cases, the rest of the kinetic energy is calculated later
  */
  En.kin-=En.kinshear;
  En.kin_tr-=En.kinshear;
#endif /*# SHEAR */

  En.T=En.kin/No.f; /* NB: En.kin still doubled here ! */

  if (measure) measureP(2);

  /* SHAKE bond length correction -> virial of constraint forces */
  En.virc/=hh;
#if PRESSURETENSOR&PT_VIR
  loop (n,0,PT_DIM) PTvirc[n]/=hh;
#endif /*# PRESSURETENSOR&PT_VIR */

  /* note that measure=1 for MTTK NPT every step */
  if (measure==1) measurePconstraints();

  /* (finishing) thermostats and tau.rho<0 */
#include "thermo.c"

  En.kin /= 2;

#ifdef POLAR
  En.pot += En.self;
#endif /*# POLAR */
  En.pot += En.el; /* do not move to forces() */

  En.U=En.Unc=En.kin+En.pot;
#ifndef FREEBC
  En.U+=En.corr/box.V;
#endif /*# FREEBC */
  En.tot=En.U+En.ext; /* cutoff corrections added in V2.7 because of the barostat */

  //  if (measure) prt("%.9g %.9g %.9g %.9g NOSE",En.kin,En.pot,No.f*T*cfg[0]->logs,No.f*T*Sqr(tau.T*Vpred[0]/h)/2);

  if (measure) if (!finite(En.tot))
    ERROR(("numeric problem: En.kin=%g En.pot=%g",En.kin,En.pot))

#ifdef ANCHOR
  if (anchor.f) {
    if (notinplace) {
      fprintf(anchor.f," not in place or too large force/torque");
      fprintf(stderr,"ANCHOR: not in place or too large force/torque, %dx at t=%g\n",notinplace,t); }
    fprintf(anchor.f," %g\n",t); }
#endif /*# ANCHOR */

  CPUtime("Verlet+Shake");
} /* Shake */


#ifdef POLAR
#  include "scf.c"
#endif /*# POLAR */

#ifdef SHEAR
void Shear(ToIntPtr B, ToIntPtr A, ToIntPtr VH,double H) /************ Shear */
/*
  (1) adds forces f_i[x] = f_i[y] = m_i*Kf*cos(2 PI z_i/L[z])
      where Kf=Cf/sqrt(2)  (because forces are in x and y directions)
  (2) calculates sumfc = SUM_i f_i[x] cos(2 PI z_i/L[z])
  (3) calculates sumpc = SUM_i m_i (v_i[x]+v_i[y]) cos(2 PI z_i/L[z])
      and Kv = sqrt(2)*Cv = sumpc*Kf/sumfc
  (4) records  num = sumfc*No.mass*L[2]/(2*PI*PI*L[0]*L[1])
      and      den = sumpc
      Then, <viscosity> = <num>/<den>
  (5) also Cv=sumpc/sumfc*Cf/2 is recorded in the convergence profile
      (amplitude of the velocity profile)
  (6) Twice the systematic kinetic energy contribution from the flow is
      En.kinshear = Kv^2*sumfc/2/Kf

  NOTES:
  * VH is a product of velocities and H (normally H=1, for Shake H=h)
  * For Shake, the velocities are by h/2 delayed (does not matter here)
  * Shear should be included before any thermostat (the systematic term
      En.Ekinshear has to be subtracted form the total [and
      intermolecular] energy)
*/
{
  molecule_t *mn;
  siteinfo_t *si;
  int ns,sp,n,i;
  vector *r,*f,*v;
  double
    cosz,m,fx,
    rho,
    phi,Kf,
    aux,
    sumfc,sumpc,sumfx;
  static int pass;

  if (shear==0) return;
  if (shear<0) ERROR(("negative shear not supported"))

  phi=2*PI/box.L[2];
  if (fabs(box.L[1]/box.L[0]-1)>1e-6)
    ERROR(("Shear implementation limitation: L[0]!=L[1]"))

  /* eta*dT/dt/rho^2 in program units: */
  rho=No.mass*rhounit/(box.L[0]*box.L[1]*box.L[2]);
  aux=shear/Sqr(rho)*(massunit/Pow5(lengthunit)*Sqr(timeunit));

  /* Kf=Cf/sqrt(2) because in x and y directions */
  Kf=PI*sqrt(2*aux*No.f/(Cub(box.L[2])*box.L[0]*box.L[1]));

  if (pass==0) {
    double Cf=Kf*sqrt(2.);
    double viscunit=Punit*timeunit;

    pass++;
    put3(shear,Kf,Cf);
    put2(viscunit,No.mass);
    put2(box.L[2],box.L[0]*box.L[1]*box.L[2]);
    prt("expected heating rate in prog.units [K/ps, etc.]: dE/dt = %g/eta",
        sqr(Cf*No.mass*box.L[2]/4/PI)/(box.L[0]*box.L[1]*box.L[2])); }

  /*
     NOTE to the estimate of Cf [Kf]: we use the `kinetic' estimate of
     heat capacity = No.f/2 (No.f would be probably better)
     because the same estimate is used in friction thermostats;
     (the time-scaling errors cancel)
  */

  sumfc=sumpc=sumfx=0;

  /* shear force + measurements */
  loop (n,FROM,No.N) {
    mn=molec+n;
    ns=mn->ns;
    sp=mn->sp;
    si=spec[sp]->si;
    f=rof(mn,B->rp);
    r=rof(mn,A->rp);
    v=rof(mn,VH->rp);
    loop (i,0,ns) {
      m=si[i].mass;
      cosz=cos(phi*r[i][2]);
      fx=Kf*(m*cosz);
      sumfx+=fx;
      f[i][0]+=fx; f[i][1]+=fx;
      sumfc+=fx*cosz;
      sumpc+=(v[i][0]+v[i][1])*(m*cosz); } }

  /* correction to total acceleration in x,y directions
     (exact momentum conservation) */
  fx=sumfx/No.mass;
  loop (n,FROM,No.N) {
    mn=molec+n;
    ns=mn->ns;
    sp=mn->sp;
    si=spec[sp]->si;
    f=rof(mn,B->rp);
    loop (i,0,ns) {
      m=si[i].mass;
      f[i][0]-=fx*m; f[i][1]-=fx*m; } }

  sumpc/=H;

  if (measure) {
    StaSet(0,2,2,lag.n);
    StaAdd("shear:num",No.mass*box.L[2]/(2*PI*PI*box.L[0]*box.L[1])*sumfc);
    StaAdd("shear:den",sumpc); }
  En.Cv=sumpc/sumfc*Kf*0.7071067811865476;
  En.kinshear=Sqr(En.Cv)*sumfc/Kf; /* twice! */

  /* kinetic T contribution equivalent to global cos-like velocity profile */
  StaSet(0,2,2,0);
  StaAdd("Tkinshear",En.kinshear/No.f);
}
#endif /*# SHEAR */

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   POLAR ERRORS + VIRTUAL VOLUME CHANGE                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

struct constrd_s constrd;
#ifndef FREEBC
static ToIntPtr locA,locB;
#endif /*# FREEBC */

void measureP(int pass) /***************************************** measureP */
/*
  Virtual volume/area method.

  A selected combination of diagonal components of the pressure tensor are
  calculated using the virtual change of volume or shape of a rectangular
  simulation box.  These include: total pressure, pressures in x,y,z (see
  variable `rescale').

  If densprof.slab&2 (=constrd.mode&RESCALE_SLAB), then x,y,z are rescaled so
  that volume does not change (virtual area change by Gloor GJ, Jackson G, Blas
  FJ, de Miguel E: "Test-area simulation method for the direct determination of
  the interfacial tension of systems with continuous or discontinuous
  potentials", JCP, 123, 134703 (2005))

  Limitations:
  - Only one rescaling in one call (i.e., cannot calculate all components of
    the pressure tensor at once).
  - Cannot calculate off-diagonal elements.
  POLAR:
    polar precision epspx must be small because of numerical derivative
  NOTES:
  - constrd.dV is in units of V (relative change)
  - uses two evaluations of forces (sparing one would mean to include this
    code to the place where the forces are normally calculated or exporting
    them which would be mess)
  - must work in 2 passes because Ekin need not be known while pass 1:
    measureP(1) : when positions are known
    measureP(2) : when Ekin is known
*/
{
#ifndef FREEBC
  static int warn;
  static double Pkin;
  double Ekin,Ep,Em,vir=0;
  vector L0,QM0={0,0,0}; /* initialized to suppress compiler warning */
  En_t En0;
#  ifdef ECC
  double U2p,U2m,U3p,U3m;
#  endif /*# ECC */
#  ifdef DEBUG
  En_t Enp;
#  endif /*# DEBUG */

  if (!constrd.dV) return;

  box.V=PROD(box.L);
  box.cq=Sqr(box.cutoff);

  /* ECC only: remember default box */
  ecc(-1);

  if (pass==1) {
    /* numerical derivative dU/dV */

    if (No.c && !(rescale&RESCALE_CM) && !warn++)
      WARNING(("%d constrained bonds and atom-based rescaling:\n\
*** pressure by virtual volume change will be wrong\n\
*** - for small molecules use CM-based rescaling: rescale=8|(more options)\n\
*** - for large molecules use -u9999 (vibrating bonds) and check the thermostat",No.c))

    En0=En;
    VV(L0,=box.L)
    if (Q) VV(QM0,=Q->M)

    if (!locA) sdsalloc(locA,cfg[0]->size);
    if (!locB) sdsalloc(locB,cfg[0]->size);

    sdscopy(locA,cfg[0]);

    /* enlarged box/area */
    box.V=rescalecfg(locA,constrd.mode,exp(constrd.dV/No.ncoord),NULL);
    if (!(constrd.mode&RESCALE_CM)) depend_r(locA,1); /* meaningfull for Rowlinson only */
#  ifdef POLAR
    scforces(locB,locA);
    StaAdd("selffieldx iter",scf.nit);
    StaAdd("selffieldx maxdr",scf.maxdr);
    Ep=En.pot; /* +En.el+En.self is in scforces() */
#  else /*# POLAR */
    zeroEn();
    forces(locB,locA);
#    ifdef ECC
    ecc(1);
    U2p=En.ECC_U2;
    U3p=En.ECC_U3;
#    endif /*# ECC */
    Ep=En.pot+En.el;
#  endif /*#!POLAR */
#  ifdef DEBUG
    Enp=En;
#  endif /*# DEBUG */

    /* shrunk box/area */
    box.V=rescalecfg(locA,constrd.mode,exp(-2*constrd.dV/No.ncoord),NULL);
    if (!(constrd.mode&RESCALE_CM)) depend_r(locA,1); /* meaningfull for Rowlinson only */
#  ifdef POLAR
    scforces(locB,locA);
    StaAdd("measureP:scf.nit",scf.nit);
    StaAdd("measureP:scf.maxdr",scf.maxdr);
    Em=En.pot; /* +En.el+En.self is in scforces() */
#  else /*# POLAR */
    zeroEn();
    forces(locB,locA);
#    ifdef ECC
    ecc(1);
    U2m=En.ECC_U2;
    U3m=En.ECC_U3;
#    endif /*# ECC */
    Em=En.pot+En.el;
#  endif /*#!POLAR */

#  include "constrddebug1.c"

    /* restore previous state obtained while integration */
    En=En0;
    if (Q) VV(Q->M,=QM0)
    VV(box.L,=L0)
    VV(box.Lh,=0.5*box.L)
    box.V=PROD(box.L);

    vir=(Em-Ep)/(2*constrd.dV);
#  ifdef ECC
    StaAdd("ECC P2fulldV [Pa]",(U2m-U2p)/(2*box.V*constrd.dV)*Punit);
    StaAdd("ECC P3fulldV [Pa]",(U3m-U3p)/(2*box.V*constrd.dV)*Punit);
#  endif /*# ECC */
  } /* pass=1 */

  else if (pass==2) {
    /* kinetic contribution (NB: En.kin, En.kin_tr doubled here) */
    if (constrd.mode & RESCALE_CM) Ekin=En.kin_tr*No.PkinqCM;
    else Ekin=En.kin*No.Pkinq;

    if (constrd.mode & RESCALE_SLAB) Pkin=0;
    else Pkin=Ekin/(DIM*box.V);

#  if PRESSURETENSOR & PT_KIN && PRESSURETENSOR & PT_MOL
    /*
       center-of-mass rescaling: use pressure tensor component, not total Ekin
       (is the same on average anyway...)
    */
    /* UNFINISHED - probably not important */
#  endif /*# PRESSURETENSOR & PT_KIN && PRESSURETENSOR & PT_MOL */

    En.PdVnc=Pkin+vir/box.V;
    En.PdV=Pkin+(vir/box.V+En.corr/Sqr(box.V));
    // StaSet(0,lag.ierr,2,lag.in); ?
    StaSet(0,lag.err,2,lag.n);
    StaAdd(constrd.PdVname, En.PdV*Punit);
    StaSet(0,lag.err,2,2);
    StaAdd("PdVkin [Pa]",Pkin*Punit);
    StaAdd("dVvir",vir);
    StaAdd("PdVvir [Pa]",vir/box.V*Punit);

#  include "constrddebug2.c"
  }
  else
    ERROR(("measureP: bad pass=%d",pass))
#endif /*# FREEBC */
}

#ifdef POLAR
void testSCF(void) /************************************************ testSCF */
/*
  For testing quality of SCF integration: full iterations
  FREEBC also supported
  See also measureP() where some test are also performed (except FREEBC)
  WARNING: 1st evalation of forces unnecessarily calculated twice
*/
{
  En_t En0;
  FILE *frun=NULL,*fex=NULL,*ferr=NULL;
  int n,i;

  if (!scf.test) return;

  En0=En;

  box.V=PROD(box.L);
  box.cq=Sqr(box.cutoff);

  if (!locA) sdsalloc(locA,cfg[0]->size);
  if (!locB) sdsalloc(locB,cfg[0]->size);

  sdscopy(locA,cfg[0]);

  if (!(constrd.mode&RESCALE_CM)) depend_r(locA,1); /* meaningfull for Rowlinson only */

  scforces(locB,locA);
  depend_f(locA,locB); /* warning: param order ! */
  StaAdd("scf.nit",scf.nit);
  StaAdd("scf.maxdr",scf.maxdr);

  En=En0;
  En.Polmaxerr=En.Polstderr=0;
  En.Fstderr=En.Fmaxerr=0;

  if (option('v')&32) {
    char *fmode=init<2 || (init_append&2) ? "at" : "wt";

    frun=fopen(Fn("run.pol"),fmode);
    fex=fopen(Fn("ex.pol"),fmode);
    ferr=fopen(Fn("err.pol"),fmode);
    fprintf(frun,"\
#  t=%.4f  Running induced dipole moments in p.u. (%.11g D)\n\
#  dip_x              dip_y              dip_z               |dip|\n",t,Debye);
    fprintf(fex,"\
#  t=%.4f  Iterated induced dipole moments in p.u. (%.11g D)\n\
#  dip_x              dip_y              dip_z               |dip|\n",t,Debye);
    fprintf(ferr,"\
#  t=%.4f  Errors in induced dipole moments in p.u. (%.11g D)\n \
#  dip_x              dip_y              dip_z               |dip|\n",t,Debye); }

  loop (n,FROM,No.N) {
    molecule_t *mn=molec+n;
    int ns=mn->ns;
    siteinfo_t *si=spec[mn->sp]->si;
    vector *rpolRun=polarrof(mn,cfg[0]->rp);
    vector *rpolA=polarrof(mn,locA->rp);
    vector *fRun=rof(mn,cfg[2]->rp);
    vector *fB=rof(mn,locB->rp);
    vector dd;
    double f;

    loop (i,0,ns) if (si[i].chargepol) {

      VVV(dd,=rpolRun[i],-rpolA[i])
      f=SQR(dd)*si[i].qqpol;
      if (frun) fprintf(frun,"%18.12f %18.12f %18.12f  %18.12f\n",
                        rpolRun[i][0]*si[i].chargepol,
                        rpolRun[i][1]*si[i].chargepol,
                        rpolRun[i][2]*si[i].chargepol,
                        sqrt(SQR(rpolRun[i])*si[i].qqpol));
      if (fex) fprintf(fex,"%18.12f %18.12f %18.12f  %18.12f\n",
                       rpolA[i][0]*si[i].chargepol,
                       rpolA[i][1]*si[i].chargepol,
                       rpolA[i][2]*si[i].chargepol,
                       sqrt(SQR(rpolA[i])*si[i].qqpol));
      if (ferr) fprintf(ferr,"%18.12f %18.12f %18.12f  %18.12f\n",
                        dd[0]*si[i].chargepol,
                        dd[1]*si[i].chargepol,
                        dd[2]*si[i].chargepol,
                        sqrt(f));
      En.Polstderr+=f; Max(En.Polmaxerr,f)

      f=SQRD(fRun[i],fB[i]);
      En.Fstderr+=f; Max(En.Fmaxerr,f) } }

  if (frun) { fprintf(frun,"\n"); fclose(frun); }
  if (fex) { fprintf(fex,"\n"); fclose(fex); }
  if (ferr) { fprintf(ferr,"\n"); fclose(ferr); }

  // StaSet(0,lag.ierr,2,lag.in); ?
  StaSet(0,lag.err,2,lag.n);

  En.Polstderr=sqrt(En.Polstderr/No.pol);
  En.Polmaxerr=sqrt(En.Polmaxerr);
  StaAdd("Polar stderr",En.Polstderr);
  StaAdd("Polar maxerr",En.Polmaxerr);

  En.Fstderr=sqrt(En.Fstderr/No.pol);
  En.Fmaxerr=sqrt(En.Fmaxerr);
  StaAdd("Polar force stderr [N]",En.Fstderr*(massunit*lengthunit/Sqr(timeunit)));
  StaAdd("Polar force maxerr [N]",En.Fmaxerr*(massunit*lengthunit/Sqr(timeunit)));
} /* testSCF() */

void measureepspol(void) /************************************ measureepspol */
/*
   to determine eps^inf (optical permittivity)
   Ewald needed
*/
{
#  ifndef FREEBC
  //  static ToIntPtr locA=NULL,locB; // not local now
  vector M0,dM;
  int i,nit0;
  double E,aux;
  static int phase; /* x,-x,y,-y,z,-z */
  static int stride=0; /* not every cycle */

  /* save previous state */
  struct Eext_s Eext0;
  En_t En0;
  vector QM0;

  if (scf.Estride && stride++%scf.Estride) return;

  En0=En;
  Eext0=Eext;
  VV(QM0,=Q->M)

  box.cq=Sqr(box.cutoff);

  if (!locA) sdsalloc(locA,cfg[0]->size);
  if (!locB) sdsalloc(locB,cfg[0]->size);

  sdscopy(locA,cfg[0]);

  scforces(locB,locA);
  nit0=scf.nit;

  VVV(M0,=dM,=Q->M)

  E=Eext.Eepspol*(1-2*(phase&1)); /* in prog. units, w. sign */

  Eext.E[phase/2]+=E;
  Eext.isE=1;

  scforces(locB,locA);

  VVV(dM,=Q->M,-M0)
  aux=4*PI/PROD(box.L)/E;
  if (option('v')&4) prt("%g %g %g  %c E=%g p.u. it=%d,%d HIGHFREQ",
      dM[0]*aux,dM[1]*aux,dM[2]*aux,
      'x'+phase/2,
      E,nit0,scf.nit);

  StaSet(0,2,2,0);
  StaAdd("hf:scf.nit(0)",nit0);
  StaAdd("hf:scf.nit(E)",scf.nit);

  StaSet(0,2,2,lag.in);
  En.chi_hf=dM[phase/2]*aux;
  StaAdd("high-frequency chi",En.chi_hf);
  aux=4*PI/(3*T*PROD(box.L));
  StaAdd("total chi",En.chi_hf+aux*SQR(M0));

  phase=(phase+1)%6;
  Eext=Eext0;
  /* to fix total chi (original Eext) */
  if (Eext.isE && scf.Estride>1) {
    aux=sqrt(aux);
    putv(Eext.E)
    loop (i,0,3) if (Eext.E[i]) StaAdd(string("xchi_%c",'x'+i),aux*M0[i]); }
  En=En0;
  VV(Q->M,=QM0)
#  endif /*# FREEBC */
}
#endif /*# POLAR */

#include "jacobi.c"
#include "jacobins.c"

typedef struct sort_s {
  double l;
  int i;
} sort_t;

static int sortnmf(const void *a,const void *b)
{
  sort_t *A=(sort_t *)a;
  sort_t *B=(sort_t *)b;
  if (A->l<B->l) return -1;
  else return A->l>B->l;
}

/*                dr  eps  ampl mem/GiB modes frames zero key*/
#ifdef FREEBC
struct nm_s nm = {0, 1e-9, 0.3, 13.,    0,    7,     6,   0};
#else /*# FREEBC */
struct nm_s nm = {0, 1e-9, 0.3, 13.,    0,    7,     3,   0};
#endif /*#!FREEBC */

#include "nm.c"
#include "nmc.c"
