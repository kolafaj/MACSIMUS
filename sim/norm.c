#include "ground.h"
#include "simglob.h"
#include "simdef.h" /* nspec needed because of constrit stuff... */
#include "simcg.h"
#include "norm.h"

/*
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %             errors of constraints : measuring and correcting            %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/* errors of constraints and SHAKE omegac optimization */
constrit_t *constrit; /* [nspec] */

void initconstrit(double omegac) /**************************** initconstrit */
{
  if (!constrit) {
    /* we do not store optimized omegac, but they are remembered
       within one job for any init */
    int sp;

    allocarrayzero(constrit,nspec);
    loop (sp,0,nspec) {
      constrit[sp].omegac=fabs(omegac);
      constrit[sp].omega=fabs(omegac)/2;
      constrit[sp].d=0.03125;
      constrit[sp].weight=0.25/4*sqrt(min(16,spec[sp]->N*spec[sp]->nc)); }
    prt("constraints statistics%s initialized",omegac<0?" and omegac optimization":""); }
}

void doconstrit(double omegac) /********************************* doconstrit */
{
  int sp,i;
  static double often=1;
  static int nsteps=0;

  loop (sp,0,nspec) {

    if (omegac<0) {
      /* omegac optimization */
      if (option('v')&4)
        prt("constrit: sp=%d omegac=%f nit=%d",
                       sp,constrit[sp].omegac,constrit[sp].nit[COR_BONDS]);
      if (constrit[sp].omegaclast) {
        /* omegac was trial */
        if (constrit[sp].nit[COR_BONDS] && constrit[sp].nit[COR_LAST]) {
          /* two omegac's have been measured */
          if (constrit[sp].nit[COR_BONDS]<constrit[sp].nit[COR_LAST]) {
            /* move closer to the trial=better */
            constrit[sp].omegac=constrit[sp].omegaclast*(1-constrit[sp].weight)+constrit[sp].omegac*constrit[sp].weight;
            constrit[sp].omega=constrit[sp].omegac/2; }
          else
            /* keep the old=better */
            constrit[sp].omegac=constrit[sp].omegaclast; }
        constrit[sp].omegaclast=0; }
      else if (rnd()<often) {
        /* schedule trial omegac */
        constrit[sp].omegaclast=constrit[sp].omegac; /* store last */
        constrit[sp].nit[COR_LAST]=constrit[sp].nit[COR_BONDS];
        constrit[sp].omegac+=rnd()*constrit[sp].d; /* trial omegac */
        constrit[sp].omega=constrit[sp].omegac/2;
        constrit[sp].d=-constrit[sp].d; }
      else
        constrit[sp].omegaclast=0; }

    loop (i,0,COR_LAST) {
/* do not report here (every step), but after every cycle ***
      if (constrit[sp].nit[i]) StaAdd(string("constrit%d spec%d",i,sp),
        (double)constrit[sp].nit[i]/spec[sp]->N);
*/
      constrit[sp].nit[i]=0; } }

  often=(1-1./200)*often+0.1/200; /* from 1 -> 0.1 within ~200 steps*/
  if (nsteps++==200) loop (sp,0,nspec) {
    constrit[sp].weight/=4;
    constrit[sp].d/=2; }
}

double vconstrainterror=-1;

double constrainterror(ToIntPtr A,ToIntPtr V) /* ========== constrainterror */
/***
  Standard deviation of relative constraint error defined by
    sqrt ( SUM{a} (r[a]^2/bond^2-1)^2 / No.c )
  where
    r[a] = r[site0 of constraint a] - r[site1 of constraint a]
  is returned.
  In addition, the following measure of the error of velocity is calculated:
    vconstrainterror = sqrt ( SUM{a} (v[a].r[a]/bond)^2 / SUM{a} v[a]^2 )
  if not measured then -1 is returned
***/
{
  int n,nc,a,i,j;
  molecule_t *m;
  siteinfo_t *si;
  vector *r,*v;
  vector ra,va;
  double sq,sumsq,sumsc,sumvv;
  double z;
  static double errlim=0.5;

  if (measure==0 || No.c==0) return -1;

  sumsc=sumsq=0;
  sumvv=1e-33;

  loop (n,FROM,No.N) {
    m=molec+n;
    si=spec[m->sp]->si;
    r=rof(m,A->rp);
    v=rof(m,V->rp);
    nc=m->nc;

    loop (a,0,nc) {
      i=si[a].pair[0]; j=si[a].pair[1];
      VVV(ra,=r[i],-r[j])
      sq=SQR(ra)/si[a].bondq-1;
      if (fabs(sq)>errlim) {
        ERROR(("bad constraint %d [%d-%d] in molecule %d (%.3f instead of %.3f)",
               a,i,j,n,
               sqrt(SQR(ra)), sqrt(si[a].bondq) ))
#if 0
        errlim=fabs(sq);
        prt("error limit increased to %f",errlim);
#endif       /*# 0 */
        }
      sumsq += sq*sq;
      VVV(va,=v[i],-v[j])
      sumvv += SQR(va);
      sq=SCAL(ra,va); sumsc += sq*sq/si[a].bondq;
      } }

  z=sqrt(sumsq/No.c);
  vconstrainterror=sqrt(sumsc/sumvv);
  if (z>1) Error("constraint failure");

  return z;
} /* constrainterror */

int Scorrect( /* ================================================== Scorrect */
  molecule_t *m,
  vector *r,    /* = rof(m,rp), NOT rp for sites */
  vector *rpvel,/* = rp, NOT rof(m,rp) for velocities */
  double eps)
/***
  Positions of sites and velocities of molecule  m  are corrected to satisfy
  the constraints by the direct iteration method similar to SHAKE (RATTLE).
  The center of mass and total momentum do not change.

  If  eps<1  then iterations stop once relative errors of all constraints
  are less than eps.
  If eps>=1  then  (int)eps  iterations are performed irrespective of errors.
  If eps<=0  then nothing is done (see init==4 and initcfg)
  Number of iterations (sweeps over all constraints) is returned.
  NEW: Errors of velocity constraints are now also used to control the
  number of iterations.
***/
{
  vector *v;
  vector ra,dr;
  int nc=m->nc,i,j,a,it=0;
  int ifit=eps>=1;
  int maxit = ifit ? (int)eps : nc*100;
  siteinfo_t *si=spec[m->sp]->si;
  double sq,maxsq,sc,maxsc,mi,mj,epsv,omegav,omega;

  if (eps<=0 || !nc) return 0;

  omegav = constrit[m->sp].omegac;
  eps *= 2; /* because  (ra^2/bond^2-1) ~ 2*|ra/bond-1| */
  omega = constrit[m->sp].omega; /* omegac/2 - the same reason */

  v=rof(m,rpvel);
  /* calculate epsv = error of vel. constr. incl. typical dimension */
  sq=0;
  loop (a,0,nc) {
    i=si[a].pair[0]; j=si[a].pair[1];
    VVV(ra,=r[i],-r[j])
    VVV(dr,=v[i],-v[j])
    sq+=SQR(ra)*SQR(dr); }
  sq=sqrt(sq/nc);
  epsv=sq*eps;

  do { it++;
    maxsq=maxsc=0;
    loop (a,0,nc) {
      i=si[a].pair[0]; j=si[a].pair[1];
      VVV(ra,=r[i],-r[j])
#ifdef TWODIM
      sc=ra[0]*(v[i][0]-v[j][0])+ra[1]*(v[i][1]-v[j][1]);
#else /*# TWODIM */
      sc=ra[0]*(v[i][0]-v[j][0])+ra[1]*(v[i][1]-v[j][1])+ra[2]*(v[i][2]-v[j][2]);
#endif /*#!TWODIM */
      if (fabs(sc)>maxsc) maxsc=fabs(sc);
      sq = SQR(ra)/si[a].bondq-1;
      /* added 6/2001 to catch diverging constraints when very bad input */
      if (sq>100) ERROR(("cannot correct constraints"))
      if (fabs(sq)>maxsq) maxsq=fabs(sq);
      mi=si[i].imass; mj=si[j].imass;
      sq *= omega/(mi+mj); /* weighting to conserve the center of mass */
      VV(dr,=sq*ra)
      VV(r[i],-=mi*dr)
      VV(r[j],+=mj*dr)
      sc *= omegav/(si[a].bondq*(mi+mj));
      VO(ra,*=sc)
      VV(v[i],-=mi*ra)
      VV(v[j],+=mj*ra) } /*a*/
  } while ((ifit || maxsq>eps || maxsc>epsv) && it<maxit);

  if (!ifit && it==maxit) Error("Scorrect: too many iterations");

  return it;
} /* Scorrect */

int shake_r(ToIntPtr A) /* ========================================= shake_r */
/***
  the same as Scorrect without velocities, for the whole configuration
***/
{
  vector ra,dr,*r;
  molecule_t *mn;
  int n,nc,i,j,a,it=0;
  int sumit=0,ifit=eps>=1,maxit;
  siteinfo_t *si;
  double sq,maxsq,mi,mj;
  double eps=1e-14; /* TO CHANGE */
  double omega=0.5; /* conservative fixed value - no dynamic optimization*/

  if (eps<=0) return 0;

  eps *= 2; /* because  (ra^2/bond^2-1) ~ 2*|ra/bond-1| */

  loop (n,FROM,No.N) {
    mn=molec+n;

    if ( (nc=mn->nc) ) {

      maxit = ifit ? (int)eps : nc*100;
      si=spec[mn->sp]->si;
      r=rof(mn,A->rp);

      if (nc>No.maxc) Error("nc");

      it=0;

      do { it++;
        maxsq=0;
        loop (a,0,nc) {
          i=si[a].pair[0]; j=si[a].pair[1];
          VVV(ra,=r[i],-r[j])
          sq = SQR(ra)/si[a].bondq-1;
          if (sq>100) ERROR(("cannot correct constraints"))
          if (fabs(sq)>maxsq) maxsq=fabs(sq);
          mi=si[i].imass; mj=si[j].imass;
          sq *= omega/(mi+mj); /* weighting to conserve the center of mass */
          VV(dr,=sq*ra)
          VV(r[i],-=mi*dr)
          VV(r[j],+=mj*dr) } /* a */
      } while ((ifit || maxsq>eps) && it<maxit);

      sumit+=it;
      if (!ifit && it==maxit) Error("shake_r: too many iterations");
    }
  }
  return sumit;
} /* shake_r */


void Lcorrect(ToIntPtr B, ToIntPtr V) /* ========================== Lcorrect */
/***
  Positions of sites of all molecules and their velocities are corrected
  to satisfy the constraints by solving linearized equations for the
  corrections.
  Knowledge of all matrices  M  computed in the previous call to rhs
  (using the predictor) is required.
  See also comments on the constraint dynamics later!
***/
{
  int n,a,nc,i,j;
  vector *v,*r,*ra=NULL;
  molecule_t *mn;
  vector auxv;
  double *dc=NULL,*e=NULL,*dv=NULL,*ve=NULL;
  double z,sc,sq;
  M_t *Msp;
  siteinfo_t *si;

  if (No.maxc) {
    rallocarray(ra,No.maxc);  /* = r[i]-r[j] for constraint a=(i,j) */
    rallocarray(dc,No.maxc);  /* errors of constraints */
    rallocarray(e,No.maxc);   /* multipliers to correct constraints */
    rallocarray(dv,No.maxc);  /* errors of velocity constraints */
    rallocarray(ve,No.maxc);  /* multipl. to correct vel. constraints */
  }

  loop (n,FROM,No.N) {
    mn=molec+n;

    if ( (nc=mn->nc) ) {
      si=spec[mn->sp]->si;
      r=rof(mn,B->rp);
      v=rof(mn,V->rp);

      Msp=Ms[n];

      if (nc>No.maxc) Error("nc");

      loop (a,0,nc) {
        i=si[a].pair[0]; j=si[a].pair[1];
        ra[a][0]=z=r[i][0]-r[j][0]; sq=si[a].bondq-z*z; sc= -z*(v[i][0]-v[j][0]);
        ra[a][1]=z=r[i][1]-r[j][1]; sq-=z*z;            sc-= z*(v[i][1]-v[j][1]);
#ifndef TWODIM
        ra[a][2]=z=r[i][2]-r[j][2]; sq-=z*z;            sc-= z*(v[i][2]-v[j][2]);
#endif /*# TWODIM */
        dv[a]=sc;
        dc[a]=sq/2; } /* a */

      memset(e,0,nc*sizeof(e[0]));
      constrit[mn->sp].nit[COR_BONDS] += omega==0
        ? conjgrad(nc,e,Mvec,Msp,dc,epsc)
        : directiter(nc,e,Msp,dc,omegac,epsc);

      memset(ve,0,nc*sizeof(ve[0]));
      constrit[mn->sp].nit[COR_VELOC] += omega==0
        ? conjgrad(nc,ve,Mvec,Msp,dv,epsc)
        : directiter(nc,ve,Msp,dv,omegac,epsc);

      loop (a,0,nc) {
        i=si[a].pair[0]; j=si[a].pair[1];

        /* length constraints */
        VV(auxv, =e[a]*ra[a])
        VV(r[i], +=si[i].imass*auxv)
        VV(r[j], -=si[j].imass*auxv)

        /* velocity constraints */
        VV(auxv, =ve[a]*ra[a])
        VV(v[i], +=si[i].imass*auxv)
        VV(v[j], -=si[j].imass*auxv) } /*a*/

    } /* if (nc) */

  } /* n */

  if (No.maxc) release(ra);
} /* Lcorrect */


/*
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %                 normalization of positions and momenta                  %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

int normalize(int m) /******************************************** normalize */
/***
 * Positions of molecule(s) are normalized so that all coordinates are
   positive and there is at least one site in the periodic box
   [0,Lx)*[0,Ly)*[0,Ly).
 * It is assumed that no more than one L must be added or subtracted.
 * No nearest-image convention needs to be applied to calculate distances
   within one molecule.
 * A molecule must not extend over [0,3/2*L] ([0,5/2*L] if WORM is #defined)
 * m=-2: Normalizes all molecules, verbose.
   m=-1: Normalizes all molecules, quiet.
   m>=0: Normalizes molecule m.
 * Void #ifdef FREEBC.
***/
{
#ifdef FREEBC
  return 0;
#else /*# FREEBC */
  int mshifted=0,pshifted=0;
  int i,n,k,ns,pn,mn;
  vector rmin,rmax,RMIN,RMAX;
  vector *r,*rp=cfg[0]->rp;
  double dL;
  vector rlimit;
  siteinfo_t *si;

#  ifdef SLIT
  int dim=DIM-1; /* max dim of periodic b.c. */
#  else /*# SLIT */
  int dim=DIM; /* max dim of periodic b.c. */
#  endif /*#!SLIT */

  int from=0,to=No.N;

#  ifdef LINKCELL
  vector LR0;

  VV(LR0,=box.r0+box.L)
#    define R0 box.r0
#  else /*# LINKCELL */
#    define R0 0.0
#    define LR0 box.L
#  endif /*#!LINKCELL */

#  ifdef WORM
  VV(rlimit,=2.5*box.L)
#  else /*# WORM */
  VV(rlimit,=1.5*box.L)
#  endif /*#!WORM */

  if (m>=0) from=m,to=m+1;
  if (m<-1) { VO(RMIN,=9e9) VO(RMAX,=-9e9) }

  //  if (!a[0]->dep) {
  //    static int pass;
  //    depend_r(a[0]);
  //    if (pass<30) WARNING(("internal: normalization without dependants, pls report this incident!"))
  //    pass++; }

  loop (n,from,to) {
    r=rof(molec+n,rp);
    ns=molec[n].ns;
    VO(rmin,=9e9)
    VO(rmax,=-9e9)
    si=spec[molec[n].sp]->si;
    /* WARNING: only for masses sites (dependant may be unknown)
       NB: L-dependants require box.rmin set */
    loop (i,0,ns) if (si[i].mass) loop (k,0,dim) {
      if (r[i][k]<rmin[k]) rmin[k]=r[i][k];
      if (r[i][k]>rmax[k]) rmax[k]=r[i][k]; }

    if (m<-1) loop (k,0,dim) {
      if (rmin[k]<RMIN[k]) RMIN[k]=rmin[k];
      if (rmax[k]>RMAX[k]) RMAX[k]=rmax[k]; }
    loop (k,0,dim) {
      dL=0;
      if (rmin[k]<R0) {
        dL=box.L[k]; rmax[k]+=box.L[k];
        mn=n;
        mshifted++; }
      else if (rmin[k]>=LR0[k]) {
        dL= -box.L[k]; rmax[k]-=box.L[k];
        pn=n;
        pshifted++; }
      if (rmax[k]>rlimit[k])
        ERROR(("t=%g, molecule[%d]: %c=%g out of box or too long (check #define WORM)",
             t, n,k+'x',rmax[k]))
      if (dL!=0) loop (i,0,ns) r[i][k]+=dL; } } /* n */

  if (m<-1) {
    underline("configuration normalization");
    prt("minimum and maximum coordinates:\n\
     MIN = %9.6f %9.6f %9.6f\n\
     MAX = %9.6f %9.6f %9.6f\n\
   MAX-L = %9.6f %9.6f %9.6f",
                VARG(RMIN),
                VARG(RMAX),
                RMAX[0]-box.L[0], RMAX[1]-box.L[1],RMAX[2]-box.L[2]);
    if (mshifted) prt("NOTE: %d molecules*coordinates < %g   (last n=%d)",mshifted,R0,mn);
    if (pshifted) prt("NOTE: %d molecules*coordinates > L+%g (last n=%d)",pshifted,R0,pn);
    if (mshifted+pshifted) prt("INFO: unnormalized molecules may occur if:\n\
- a configuration was loaded with a different value of box.r0 (if LINKCELL)\n\
- plb-file (with a centered configuration) has been loaded\n\
- LINKCELL and direct nearest image versions were mixed\n\
- L-dependants are present\n\
NOTE: the configuration has been normalized");
    _n }

  return pshifted+mshifted;
#endif /*#!FREEBC */
} /* normalize */

void unsplitpbc(int xyz) /*************************************** unsplitpbc */
/***
   Normalize input coordinates such that a molecule is not split over periodic
   b.c. if it is on input.
   To be called before normalize().
***/
{
#ifdef FREEBC
  return;
#else /*# FREEBC */
  vector *r,*rp=cfg[0]->rp;
  siteinfo_t *si;
  int n,i,k,ns;

#  ifdef WORM
  WARNING(("WORM (long molecules) and unsplitpbc: unsplit=%d",xyz))
#  endif /*#!WORM */

  loop (n,0,No.N) {
    r=rof(molec+n,rp);
    ns=molec[n].ns;
    si=spec[molec[n].sp]->si;
    /* WARNING: only for masses sites (dependant may be unknown)
       NB: L-dependants require box.rmin set */
    loop (i,1,ns) if (si[i].mass) loop (k,0,3) if (xyz&(1<<k)) {
      if (r[i][k]>r[i-1][k]+box.Lh[k]) r[i][k]-=box.L[k];
      if (r[i][k]<r[i-1][k]-box.Lh[k]) r[i][k]+=box.L[k]; }
  }
#endif /*#!FREEBC */
} /* unsplitpbc() */

char *prtxyz(int xyz) /********************************************** prtxyz */
/* print 1=x,2=y,4=z */
{
  static char s[8];

  s[0]=0;
  if (xyz&1) strcat(s,"x");
  if (xyz&2) strcat(s,"y");
  if (xyz&4) strcat(s,"z");
  if (s[0]==0) strcpy(s,"none");

  return s;
}

static void prtdrift(int drift) /********************************** prtdrift */
{
#ifdef FREEBC
  if (drift&7) prt("* drift of center-of-mass from center corrected for %s",prtxyz(drift&7));
#else /*# FREEBC */
  if (drift&7) prt("* drift of center-of-mass from box center corrected for %s",prtxyz(drift&7));
#endif /*#!FREEBC */
  if (drift&7*8) prt("* velocity drift corrected for %s",prtxyz((drift&7*8)/8));
  if (drift&7*64) prt("* angular momentum drift corrected for %s",prtxyz((drift&7*64)/64));
  if ((drift&511)==0) prt("* no conserved quantity to correct");
}

void setdrift(void) /********************************************** setdrift */
{
  int k;
#ifdef SLIT
  int locdrift=DRIFT_X|DRIFT_Y|DRIFT_VX|DRIFT_VY|DRIFT_AZ;
#else /*# SLIT */
  int locdrift=511;
#endif /*#!SLIT */

#ifndef FREEBC
  locdrift&=56; /* vdrift only (no CoM nor angular momentum)  */
#endif /*# FREEBC */

  if (FROM || fixsites0 || (thermostat>=T_ANDERSEN && thermostat<=T_LANGEVIN_CM)) locdrift=0;

  loop (k,0,3) if (center.cmK[k]!=0 || center.K[k]!=0)
    locdrift &= ~(DRIFT_VX<<k);
#ifdef SLAB
  loop (k,0,NCENTER) if (slab.Kz[k]) locdrift &= ~DRIFT_VZ;
  if (wall.is) locdrift &= ~DRIFT_VZ;
#endif /*# SLAB */

  underline("integrals of motion");

  prt("\
In dependence on boundary conditions, external forces, and thermostat, some\n\
quantities are conserved. To avoid cumulation of errors, they are set to zero.\n\
You may specify these quantities by variable \"drift\", or rely on the\n\
automatic setup (default). The number of conserved quantities is also\n\
subtracted from the number of degrees of freedom used for calculating the\n\
kinetic temperature and pressure; the automatically determined value can be\n\
changed by variable \"conserved\".\n\
\n\
The automatic setup fails for elliptical harmonic forces (center.K),\n\
fixed points, etc.\n\
\n\
Conservation of energy (Hamiltonian) is treated in a different module\n\
and does not count here.\n\
");

  if (drift&DRIFT_AUTO) {
    drift = (drift & 0x7ffffe00) | locdrift;
    prt("Automatically set drift = %d = %s:",drift,int2sumbin(drift));
    prtdrift(drift);
    /* NB: Eext.isE and Eext.isB not set yet */
#ifdef FREEBC
    if (SQR(el.E)+SQR(el.B)>0) ERROR(("\
Implementation limitation: I cannot set the drift if electric or magnetic\n\
*** field is present, do it manually!"))
#else
    if (SQR(el.E)+SQR(el.B)>0) WARNING(("\
Implementation limitation: Electric or magnetic field is present,\n\
*** double check the value of drift!"))
#endif
    }
  else {
    prt("Requested drift=%d:",drift);
    prtdrift(drift);
    if ((drift&511)!=(locdrift&511)) {
      prt("WARNING: this differs from the suggested drift=%d:",locdrift);
      prtdrift(locdrift);  } }

  if (drift&63) {
    prt_("* corrections will be performed at start (init/load) ");
    switch (drift&DRIFT_WHEN) {
      case DRIFT_CYCLE: prt("and every cycle."); break;
      case DRIFT_STEP: prt("and every step."); break;
      case DRIFT_SAVE: prt("and every save."); break;
      case DRIFT_START: prt("only."); break;
      default: ERROR(("drift&DRIFT_WHEN")) } }
}

void CoM(vector CM,ToIntPtr A) /**************************************** CoM */
/*
   center of mass of the whole configuration w.r.t. (0,0,0)
   box.center should be subtracted if needed (not in simmeasx.c)
*/
{
  int n,i;

  VO(CM,=0)

  loop (n,0,No.N) {
    molecule_t *mn=molec+n;
    siteinfo_t *si=spec[mn->sp]->si;
    vector *r;

    r=rof(mn,A->rp);

    loop (i,0,mn->ns) VV(CM,+=si[i].mass*r[i]) }

  VO(CM,/=No.mass)
  // NOT VV(CM,-=box.center): to be applied elsewhere
}

static char *remstr(int mask) /************************************** remstr */
{
  static char ret[4];
  int d=drift/mask,i;

  loop (i,0,3) {
    if (d&(1<<i)) ret[i]='R';
    else ret[i]='-'; }

  return ret;
}

void removedrifts(int pr) /************************************ removedrifts */
/*
  remove drifts in:
  - center of mass w.r.t box.center[]
  - total linear momentum
  - angular momentum w.r.t box.center[]
  if the respective flags in `drift' are set, the drift is removed
  see simglob.h for the flag values
  pr = print (verbose) flag
*/
{
  int n,i,ns,isCM=0;
  molecule_t *mn;
  vector *r,*v;
  vector CM;
  siteinfo_t *si;
  real mi;

  /* center of mass */
  if (drift&(DRIFT_X|DRIFT_Y|DRIFT_Z)) {

    CoM(CM,cfg[0]);
    isCM++;
    VV(CM,-=box.center)
    if (pr)
      prt("Remove CM drifts [%s] = [ %.9g %.9g %.9g ]",
          remstr(DRIFT_X), VARG(CM));

    if (!(drift&DRIFT_X)) CM[0]=0;
    if (!(drift&DRIFT_Y)) CM[1]=0;
    if (!(drift&DRIFT_Z)) CM[2]=0;
    En.CMdrift=sqrt(SQR(CM)); /* for CP */

    loop (n,0,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      ns=mn->ns;
      r=rof(mn,cfg[0]->rp);

      loop (i,0,ns) VV(r[i],-=CM) }

    VV(CM,=box.center) // new value, to be used in angular momenum drift

    normalize(-1); /* inefficient: normalizes always all coordinates */ }

  /* linear momentum (velocity drift) */
  if (drift&(DRIFT_VX|DRIFT_VY|DRIFT_VZ)) {
    vector Vdrift;

    VO(Vdrift,=0)
    loop (n,FROM,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      ns=mn->ns;
      v=rof(mn,cfg[1]->rp);
      loop (i,0,ns) VV(Vdrift,+=si[i].mass*v[i]) }

    VO(Vdrift,/=No.free_mass)
    if (pr)
      prt("Remove v drifts [%s] = [ %.9g %.9g %.9g ]",
          remstr(DRIFT_VX), Vdrift[0]/h,Vdrift[1]/h,Vdrift[2]/h);

    if (!(drift&DRIFT_VX)) Vdrift[0]=0;
    if (!(drift&DRIFT_VY)) Vdrift[1]=0;
    if (!(drift&DRIFT_VZ)) Vdrift[2]=0;
    En.vdrift=sqrt(SQR(Vdrift))/h;

    loop (n,FROM,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      ns=mn->ns;
      v=rof(mn,cfg[1]->rp);
      loop (i,0,ns) VV(v[i],-=Vdrift) } }

  /* angular momentum drift */
  if (drift&(DRIFT_AX|DRIFT_AY|DRIFT_AZ)) {
    vector I[3];
    vector angularM,rv,omega;
    real det,rr;
    int j,k;

    /*
       the angular momentum and the inertia matrix I
       with respect to box.center
    */
    VO(angularM,=0)
    memset(I,0,sizeof(I));

    if (!isCM) {
      CoM(CM,cfg[0]);
      VV(CM,-=box.center) }
    loop (n,0,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      ns=mn->ns;
      r=rof(mn,cfg[0]->rp);
      v=rof(mn,cfg[1]->rp);

      loop (i,0,ns) {
        mi=si[i].mass;
        VVV(omega,=r[i],-CM) /* omega=temporary vector here = r wrt center */
        VECT(rv,omega,v[i])
        VV(angularM,+=mi*rv)
        rr=SQR(omega);
        loop (j,0,DIM) loop (k,j,DIM)
          /* matrix is symetric: j<=k */
          I[j][k]+=mi*((j==k)*rr-omega[j]*omega[k]); } } /* n */

    VO(angularM,/=h)
    if (pr)
      prt("Remove angularM drifts [%s] = [ %.9g %.9g %.9g ]",
          remstr(DRIFT_X), VARG(angularM));

    if (!(drift&DRIFT_AX)) angularM[0]=0;
    if (!(drift&DRIFT_AY)) angularM[1]=0;
    if (!(drift&DRIFT_AZ)) angularM[2]=0;
    En.AMdrift=sqrt(SQR(angularM));

    // see AM2omega.mw for double check
    det =
       I[0][0]*I[1][1]*I[2][2] + 2*I[0][1]*I[1][2]*I[0][2]
     - I[1][1]*Sqr(I[0][2]) - I[2][2]*Sqr(I[0][1]) - I[0][0]*Sqr(I[1][2]);
    omega[0] =
      ( angularM[0]*(I[1][1]*I[2][2]-Sqr(I[1][2]))
      + angularM[1]*(I[1][2]*I[0][2]-I[0][1]*I[2][2])
      + angularM[2]*(I[0][1]*I[1][2]-I[1][1]*I[0][2]) ) / det;
    omega[1] =
      ( angularM[1]*(I[2][2]*I[0][0]-Sqr(I[0][2]))
      + angularM[2]*(I[0][2]*I[0][1]-I[1][2]*I[0][0])
      + angularM[0]*(I[1][2]*I[0][2]-I[2][2]*I[0][1]) ) / det;
    omega[2] =
      ( angularM[2]*(I[0][0]*I[1][1]-Sqr(I[0][1]))
      + angularM[0]*(I[0][1]*I[1][2]-I[0][2]*I[1][1])
      + angularM[1]*(I[0][2]*I[0][1]-I[0][0]*I[1][2]) ) / det;

    if (pr)
      prt("Remove omega drifts [%s] = [ %.9g %.9g %.9g ]",
          remstr(DRIFT_AX), VARG(omega));

    if (!(drift&DRIFT_AX)) omega[0]=0;
    if (!(drift&DRIFT_AY)) omega[1]=0;
    if (!(drift&DRIFT_AZ)) omega[2]=0;
    En.omegadrift=sqrt(SQR(omega));

    if (!finite(En.omegadrift))
      ERROR(("numeric problem: omega drift (consider changing drift)"))

    loop (n,0,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      ns=mn->ns;
      r=rof(mn,cfg[0]->rp);
      v=rof(mn,cfg[1]->rp);

      loop (i,0,ns) {
        VECT(rv,omega,r[i])
        VV(v[i],-=h*rv) } } }
}

#define HOMOGENEITY 0.1

void homogeneity(int corr) /************************************ homogeneity */
{
  vector s={0,0,0},c={0,0,0};
  double a,sumq=0,suma=0;
  vector K;
  int n,i,k;

  VV(K,=2*PI/box.L)

  loop (n,0,No.N) {
    molecule_t *mn=molec+n;
    siteinfo_t *si=spec[mn->sp]->si;
    vector *r=rof(mn,cfg[0]->rp);

    loop (i,0,mn->ns) {
      loop (k,0,3) {
        c[k]+=si[i].mass*cos(r[i][k]*K[k]);
        s[k]+=si[i].mass*sin(r[i][k]*K[k]); }
      sumq+=si[i].mass; } }

  prt_("\nhomogeneity index = ");
  loop (k,0,3) {
    a=sqrt((Sqr(s[k])+Sqr(c[k])))/sumq;
    suma+=a*a;
    prt_("%c%.3g","[  "[k],a); }
  prt("]  total=%.3g",sqrt(suma));

  if (suma>HOMOGENEITY && corr&3) WARNING(("The system is inhomogeneous (index>%g) and corr=%d is set.\n\
*** Doublecheck whether the homogeneous cutoff corrections are appropriate!\n\
*** In the SLAB version, consider slab.K and slab.range.",HOMOGENEITY,corr))

}

#if defined(SLAB) && SLAB & 2
void cleavescaling(double *Pscale,double tau,int noint,double factor,double maxscale)
                                                        /***** cleavescaling */
/*
   calculates scalings of both z-compartments (Pscale[0,1])
   from Pscale[2] in the z-direction
*/
{
  double s=exp((tau<0?noint:-1)/tau*factor*2);
  double dz=cleave.z[1]-cleave.z[0];

  if (s<1/maxscale || s>maxscale) {
    WARNING(("cleavescaling = %g (tau=%g, factor=%g) out of range [%g,%g]",
             s,tau,factor,1/maxscale,maxscale))
    if (s<1) s=1/maxscale;
    else s=maxscale;
    prt("*** scaling changed to %g",s); }

  Pscale[0]=Pscale[2]/(1+(s-1)*dz);
  Pscale[1]=s*Pscale[0];
}

double cleaverescale(ToIntPtr A,int mode,double q0,double q1)
/*
   version for cleaving with two cleaving potentials
   |  <-q0->  |  <-q1->  |  <-q0-> |
   0          z0         z1        Lz
*/
{
  int n,i;
  molecule_t *mn;
  siteinfo_t *si;
  vector *r;
  double cmz;
  double z0=cleave.z[0]*box.L[2];
  double z1=cleave.z[1]*box.L[2];

  //  prt("%11.9f %11.9f %11.9f cleaverescale", q0,q1, q0+(cleave.z[0]-cleave.z[1])*(q0-q1)/box.L[2]);

  if (mode & RESCALE_CM) /* center-of-mass based scaling */
    loop (n,0,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      r=rof(mn,A->rp);

      cmz=0;
      loop (i,0,mn->ns) cmz+=si[i].mass*r[i][2];

      cmz/=spec[mn->sp]->mass;

      if (cmz<z0) cmz*=q0-1;
      else if (cmz<z1) cmz=(q1-1)*cmz+z0*(q0-q1);
      else cmz=(q0-1)*cmz+(z0-z1)*(q0-q1);

      loop (i,0,mn->ns) r[i][2]+=cmz; }

  else /* site-based rescaling */
    loop (n,0,No.N) {
      mn=molec+n;
      si=spec[mn->sp]->si;
      r=rof(mn,A->rp);

      loop (i,0,mn->ns) {
        if (r[i][2]<z0) r[i][2]*=q0;
        else if (r[i][2]<z1) r[i][2]=q1*r[i][2]+z0*(q0-q1);
        else r[i][2]=q0*r[i][2]+(z0-z1)*(q0-q1); } }

  /* L scaling */
  box.L[2]=box.L[2]*q0+(z0-z1)*(q0-q1);

  cleave.z[0]*=q0;
  cleave.z[1]=(z0*(q0-q1)+q1*z1)/box.L[2];
  //  prt("%g %g %g cleave2",cleave.z[0],cleave.z[1]*box.L[2],box.L[2]);
  VV(box.Lh,=0.5*box.L)

#  ifdef TWODIM
  return box.L[0]*box.L[1];
#  else /*# TWODIM */
  return box.L[0]*box.L[1]*box.L[2];
#  endif /*#!TWODIM */
} /* cleaverescale */
#endif /*# defined(SLAB) && SLAB & 2 */

double rescalecfg(ToIntPtr A,int mode, double q1,double *q)
                                              /****************** rescalecfg */
/*
   Rescale the whole configuration according to bits set in mode.
   q==NULL: use factor q1 for all coordinates the bits set in mode.
   q!=NULL: use the respective q[i] for all coordinates the bits set in mode
            (i.e., for mode=2 and q[]={2,2,2}, only y is rescaled)
   The rescaling mode is a sum of binary digits:
   ========================================================================
   bit 2^bit function
   ------------------------------------------------------------------------
    0    1   rescale x coordinates by qx
    1    2   rescale y coordinates by qy  (all: RESCALE_XYZ=7)
    2    4   rescale z coordinates by qz
    3    8   rescaling is center-of-mass based (molecules not distorted)
             if not set, rescales all atom positions independently
             (RESCALE_CM)
    4   16   (in rescale only: separate x,y,z scaling; RESCALE_PTENS)
    5   32   (in rescale only: x=y,z scaling; RESCALE_XisY)
    6   64   uniform (xscale=yscale=zscale) (RESCALE_XisYisZ)
    7  128   rescale also logs and lambda (RESCALE_H)
    8  256   rescale also corresponding L (RESCALE_L)
             (e.g., mode=17 = rescale x of all sites and also L[0])
    9* 512   REMOVED
   10* 1024  SLAB: for surface-tension by virtual area change (densprof.slab&2)
             scaling factors qy:=qx, qz:=1/qx^2
             ERROR if RESCALE_L not set (RESCALE_SLAB)
    ========================================================================
    *lower bits are ignored

   returns the volume
   dependants are not recalculated
*/
{
  int n,i;
  molecule_t *mn;
  siteinfo_t *si;
  vector *r,cm;
  real M;
  double qx,qy,qz;

#if 0
  if (q)
    fprintf(stderr,"DEBUG t=%g rescaling cfg q=%.12g %.12g %.12g-times\n", t,VARG(q));
  else
    fprintf(stderr,"DEBUG t=%g rescaling cfg q1=%.12g-times\n", t,q1);
#endif

#if defined(SLAB) && SLAB & 2
  if (mode & RESCALE_CLEAVE) return cleaverescale(A,mode,qx,qy);
#endif /*# defined(SLAB) && SLAB & 2 */

  if (q) qx=q[0],qy=q[1],qz=q[2];
  else qx=qy=qz=q1;

#ifdef SLAB
  if (mode & RESCALE_SLAB) {
    if (!(mode & RESCALE_L)) ERROR(("RESCALE_SLAB without RESCALE_L\n\
*** (both cfg+L should be rescaled in the virtual area method)"))
    qy=qx;
    qz=1./Sqr(qx); }
#endif /*# SLAB */

  if (qx!=1 || qy!=1 || qz!=1) {

    if (mode & RESCALE_CM) /* center-of-mass based scaling */
      loop (n,0,No.N) {
        mn=molec+n;
        si=spec[mn->sp]->si;
        r=rof(mn,A->rp);
        VO(cm,=0)
        M=spec[mn->sp]->mass;
        loop (i,0,mn->ns) VV(cm,+=si[i].mass*r[i])
        cm[0]*=(qx-1)/M, cm[1]*=(qy-1)/M, cm[2]*=(qz-1)/M;
        if ((mode & RESCALE_XYZ)==RESCALE_XYZ)
          loop (i,0,mn->ns) VV(r[i],+=cm)
        else {
          if (mode & RESCALE_X) loop (i,0,mn->ns) r[i][0]+=cm[0];
          if (mode & RESCALE_Y) loop (i,0,mn->ns) r[i][1]+=cm[1];
          if (mode & RESCALE_Z) loop (i,0,mn->ns) r[i][2]+=cm[2]; } }

    else { /* site-based rescaling */
      loop (n,0,No.N) {
        mn=molec+n;
        si=spec[mn->sp]->si;
        r=rof(mn,A->rp);

        if ((mode & RESCALE_XYZ)==RESCALE_XYZ)
          loop (i,0,mn->ns) r[i][0]*=qx,r[i][1]*=qy,r[i][2]*=qz;
        else {
          if (mode & RESCALE_X) loop (i,0,mn->ns) r[i][0] *= qx;
          if (mode & RESCALE_Y) loop (i,0,mn->ns) r[i][1] *= qy;
          if (mode & RESCALE_Z) loop (i,0,mn->ns) r[i][2] *= qz; } } }

    /* L scaling */
    if (mode & RESCALE_L) {
      if (mode & RESCALE_X) box.L[0] *= qx;
      if (mode & RESCALE_Y) box.L[1] *= qy;
      if (mode & RESCALE_Z) box.L[2] *= qz; }

    /* termostat/barostat data scaling (normally qx=qy=qz) */
    if (mode & RESCALE_H) {
      A->logs*=qx;
      A->lambda[0]*=qx;
      A->lambda[1]*=qy;
      A->lambda[2]*=qz; }
  }

  VV(box.Lh,=0.5*box.L)

  return PROD(box.L);
} /* rescalecfg */

int rhorescale(int rescale,double tau_rho,vector finalL,double halftcyc)
/* rescaling to reach given rho */                /************** rhorescale */
/* 1 returned if the target box has been reached */
{
  vector rscale;
  // old: rescale not taken into account
  //  int ii,mode=RESCALE_L | (rescale&RESCALE_CM);
  // new: rescale works only as requested by 'rescale'
  int ii,mode=RESCALE_L | rescale;

  loop (ii,0,DIM) if (finalL[ii]) {
    double d=log(finalL[ii]/box.L[ii]);
    double x=halftcyc/tau_rho;
    //double x=(1+(d>0))*halftcyc/tau_rho;

    if (fabs(d)>5e-15 && (1<<ii & rescale)) {
      mode|=1<<ii;
      rscale[ii]=exp((fabs(d)<x ? d : d/fabs(d)*x)); }
    else rscale[ii]=1; }

  rescalecfg(cfg[0],mode,0,rscale);

  if (option('v')&4) putv(box.L)

  loop (ii,0,DIM) if (finalL[ii]) {
    double d=log(finalL[ii]/box.L[ii]);
    if (fabs(d)>5e-15 && (1<<ii & rescale)) return 0;  }

  return 1;
}

double scaling(double tau,int noint,double factor,double maxscale)
                                                              /***** scaling */
{
  double s=exp((tau<0?noint:-1)/tau*factor);

  if (s<1/maxscale || s>maxscale) {
    WARNING(("scaling = %g (tau=%g, factor=%g) out of range [%g,%g]",
             s,tau,factor,1/maxscale,maxscale))
    if (s<1) s=1/maxscale;
    else s=maxscale;
    prt("*** scaling changed to %g",s); }

  return s;
}

void distancecheck(void) /********************************** distancecheck */
/*
  debugging tool
  finds min and max distance of sites in different b.c.
*/
{
  int i,j,k;
  vector *r=cfg[0]->rp,dr;
  double rrmin[2],rrmax[2],rr;
  int imin[2]={0,0},imax[2]={0,0},jmin[2]={0,0},jmax[2]={0,0}; // init to suppress warning

  if (!cfg[0]) return;

#ifndef FREEBC
  if (box.L[0]<=0 || box.L[1]<=0 || box.L[2]<=0) ERROR(("zero or negative box"))
  loop (k,0,3)
    if (fabs(box.Lh[k]/box.L[k]-0.5)>1e-6)
      ERROR(("L[%d]=%g Lh[%d]=%g",k,box.L[k],k,box.Lh[k]))
#endif /*# FREEBC */

  loop (k,0,2) { rrmin[k]=3e33; rrmax[k]=0; }

#define check(K) \
  if (rr<rrmin[K]) { rrmin[K]=rr; imin[K]=i; jmin[K]=j; } \
  if (rr>rrmax[K]) { rrmax[K]=rr; imax[K]=i; jmax[K]=j; }

#if 0 /* very special hook */
loop (i,0,No.s) {
  rr=SQR(r[i]);
  if (rr>250000) {
    VO(cfg[1]->rp[i],=0)
    WARNING(("%d r=%f velocity set to zero",i,sqrt(rr))) }
  }
#endif /*# 0 very special hook */

  loop (i,0,No.s)
    loop (j,0,i) {
      VVV(dr,=r[i],-r[j])
      rr=SQR(dr);
      check(0)
#ifndef FREEBC
        loop (k,0,3) {
          while (dr[k]>box.Lh[k]) dr[k]-=box.L[k];
          while (dr[k]<-box.Lh[k]) dr[k]+=box.L[k]; }
      rr=SQR(dr);
      check(1)
#endif /*# FREEBC */
    }

#ifdef FREEBC
  k=0;
#else /*# FREEBC */
  loop (k,0,2)
#endif /*#!FREEBC */

    prt("%s min: |%d-%d|=%g  max: |%d-%d|=%g",k?"n.i.":"abs.",
        imin[k],jmin[k],sqrt(rrmin[k]),
        imax[k],jmax[k],sqrt(rrmax[k]));

#ifdef GOLD
  rr=9e9;
  loop (i,0,No.s) Min(rr,r[i][2])
    prt("GOLD min z=%g",rr);
#endif /*# GOLD */
}

void zeroEn(void) /************************************************* zeroEn */
{
#if PRESSURETENSOR&PT_VIR
  memset(En.Pvir,0,sizeof(En.Pvir));
#endif /*# PRESSURETENSOR&PT_VIR */
#if PRESSURETENSOR&PT_KIN
  memset(En.Pkin,0,sizeof(En.Pkin));
#endif /*# PRESSURETENSOR&PT_KIN */
#if PRESSURETENSOR&PT_MOL
  memset(En.PKin,0,sizeof(En.PKin));
#endif /*# PRESSURETENSOR&PT_MOL */
#if PRESSURETENSOR&PT_MOM
  memset(En.PKin2,0,sizeof(En.PKin2));
  memset(En.PKin3,0,sizeof(En.PKin3));
#endif /*# PRESSURETENSOR&PT_MOM */
#ifdef LOG
  En.pot0=En.el0=En.LJ0=En.bonded0=
  En.potX=En.elX=En.LJX=En.intra=En.LJ=
#endif /*# LOG */
  En.pot=En.vir=En.el=En.bonded=0;
#if SLAB
  En.Pwall[0]=En.Pwall[1]=0;
#endif /*# SLAB */
}

/*** dependants ***/

void depend_r(ToIntPtr A,int always) /***************************** depend_r */
/*
  calculate positions of dependent sites
  always=1 : always
  always=0 : only if not A->dep (this is marked during simulation)
*/
{
  molecule_t *mn;
  int n,sp;
  depend_t *d;
  vector *r;
  vector xx,yy,zz;
  int i,k;
  double rr;
  struct depitem_s *dep;
  real *depr;

  if (always) A->dep=0;

  if (A->dep) return;

  loop (n,FROM,No.N) {
    mn=molec+n;
    sp=mn->sp;
    r=rof(mn,A->rp);

    looplist (d,spec[sp]->dependants) {
      depr=r[d->indx];

      if (d->type<=DEP_M) {
        /* "Middle" (linear) dependant */
        VO(depr,=0)
        for (dep=d->dep; dep<d->dep+d->n; dep++)
          VV(depr,+=dep->w*r[dep->i]) }

      else if (d->type==DEP_R) {
        /* "Rowlinson" dependant: L-D0 is perpendicular to plane (D0,D1,D2)
            (as in Rowlinson flexible water)
            bonds/angles in D0-D1-D2 may be flexible
            L            d->indx
            |            |
            D0---D1      dep[0].i---dep[1].i
              \           \
               D2          dep[2].i                   */
        vector h1,h2,l;
        double norm;

        dep=d->dep;
        VVV(h1,=r[dep[1].i],-r[dep[0].i])
        VVV(h2,=r[dep[2].i],-r[dep[0].i])
        VECT(l,h1,h2)
        norm=d->wz/sqrt(SQR(l));
        VVV(depr,=r[dep[0].i],+norm*l) }

      else {
        /* "Lone" (out-of-plane) dependant, (D0,D1,D2) must be rigid */
        dep=d->dep;
        VVO(xx,=yy,=0)
        VO(depr,=0)

        loop (k,0,3) {
          i=dep[k].i;
          VV(depr,+=dep[k].w*r[i])
          VV(xx,+=d->x[k]*r[i])
          VV(yy,+=d->y[k]*r[i]) }

        VECT(zz,xx,yy)
        rr=d->wz/sqrt(SQR(zz));

        VV(depr,+=rr*zz) } } }

  A->dep=1; /* dependants calculated */
}

#if PRESSURETENSOR&PT_VIR
#  if PRESSURETENSOR&PT_OFF
#    if PRESSURETENSOR&PT_OFF_TR
/* tranposed to the standard version (debug only: all versions are equivalent) */
#      define CALCULATEPRESSURETENSOR(Q,G,F) \
    En.Pvir[0]+=Q*G[0]*F[0]; \
    En.Pvir[1]+=Q*G[1]*F[1]; \
    En.Pvir[2]+=Q*G[2]*F[2]; \
    En.Pvir[3]+=Q*G[2]*F[1]; \
    En.Pvir[4]+=Q*G[0]*F[2]; \
    En.Pvir[5]+=Q*G[1]*F[0];
#    elif PRESSURETENSOR&PT_OFF_AV
/* average of both version (debug only: all versions are equivalent) */
#      define CALCULATEPRESSURETENSOR(Q,F,G) \
    En.Pvir[0]+=Q*G[0]*F[0]; \
    En.Pvir[1]+=Q*G[1]*F[1]; \
    En.Pvir[2]+=Q*G[2]*F[2]; \
    En.Pvir[3]+=Q/2*(G[2]*F[1]+G[1]*F[2]); \
    En.Pvir[4]+=Q/2*(G[0]*F[2]+G[2]*F[0]); \
    En.Pvir[5]+=Q/2*(G[1]*F[0]+G[0]*F[1]);
#    else  /*#!PRESSURETENSOR&PT_OFF_TR!PRESSURETENSOR&PT_OFF_AV */
/* standard version (equivalent to the above versions) */
#      define CALCULATEPRESSURETENSOR(Q,F,G) \
    En.Pvir[0]+=Q*G[0]*F[0]; \
    En.Pvir[1]+=Q*G[1]*F[1]; \
    En.Pvir[2]+=Q*G[2]*F[2]; \
    En.Pvir[3]+=Q*G[1]*F[2]; \
    En.Pvir[4]+=Q*G[2]*F[0]; \
    En.Pvir[5]+=Q*G[0]*F[1];
// prt("%g %g %g  %g %g %g PTyz",G[1],F[2],G[1]*F[2], F[1],G[2],F[1]*G[2]);
#    endif  /*#!PRESSURETENSOR&PT_OFF_TR!PRESSURETENSOR&PT_OFF_AV */
#  else /*# PRESSURETENSOR&PT_OFF */
#    define CALCULATEPRESSURETENSOR(Q,G,F) \
    En.Pvir[0]+=Q*G[0]*F[0]; \
    En.Pvir[1]+=Q*G[1]*F[1]; \
    En.Pvir[2]+=Q*G[2]*F[2];
#  endif /*#!PRESSURETENSOR&PT_OFF */
#else /*# PRESSURETENSOR&PT_VIR */
#  define CALCULATEPRESSURETENSOR(Q,G,F) /* empty */
#endif /*#!PRESSURETENSOR&PT_VIR */

void depend_f(ToIntPtr A,ToIntPtr B) /***************************** depend_f */
/* distribute forces from dependent to independent sites */
/* Lagrange only */
{
  molecule_t *mn;
  depend_t *d;
  int n,sp;
  vector *r,*f,dr;
  real *depr,*depf; /* dependant position, force to distribute */
  int k;
  double rr;
  struct depitem_s *dep;

  loop (n,FROM,No.N) {
    mn=molec+n;
    sp=mn->sp;
    r=rof(mn,A->rp);
    f=rof(mn,B->rp);

    looplist (d,spec[sp]->dependants) {

      depf=f[d->indx];

      if (d->type<=DEP_M) {
        /* "Middle" (linear) dependant (old and new versions) */
        for (dep=d->dep; dep<d->dep+d->n; dep++)
          VV(f[dep->i],+=dep->w*depf) }

      else if (d->type==DEP_R) {
        /* "Rowlinson" dependant: L-D0 is perpendicular to plane (D0,D1,D2)
            (as in Rowlinson flexible water)
            bonds/angles in D0-D1-D2 may be flexible
            L            d->indx                L
            |            |                      l
            D0---D1      dep[0].i---dep[1].i    O---h1--H1
              \           \                      h2
               D2          dep[2].i               H2      */
        vector h1,h2,l,M,Mp;
        double h1h1,h2h2,h1h2,det,a1,a2;

        dep=d->dep;
        VVV(h1,=r[dep[1].i],-r[dep[0].i])
        VVV(h2,=r[dep[2].i],-r[dep[0].i])
        VVV(l,=r[d->indx],-r[dep[0].i])
        VECT(M,l,depf)
        VECT(Mp,l,M)
        h1h1=SQR(h1);
        h2h2=SQR(h2);
        h1h2=SCAL(h1,h2);
        det=(h1h1*h2h2-Sqr(h1h2))*SQR(l);
        a1=(SCAL(Mp,h1)*h2h2-SCAL(Mp,h2)*h1h2)/det;
        a2=(SCAL(Mp,h2)*h1h1-SCAL(Mp,h1)*h1h2)/det;
        VVV(f[dep[0].i],+=depf,-(a1+a2)*l)
        VV(f[dep[1].i],+=a1*l)
        VV(f[dep[2].i],+=a2*l)
        /* triangle components */
        CALCULATEPRESSURETENSOR(a1,l,h1)
        CALCULATEPRESSURETENSOR(a2,l,h2)

        CALCULATEPRESSURETENSOR(-1,l,depf) /* correction for moving f from L to O */
        En.vir-=SCAL(l,depf); }
      else {
        /* "Lone" (out-of-plane) dependant, (D0,D1,D2) must be rigid */
        vector xx,yy,zz; /* local orthonormal molecule coordinate system */
        double fx,fy,fz; /* f_dependant=(fx,fy,fz) in the local system */
        double alpha;

        dep=d->dep;
        VVO(xx,=yy,=0)

        loop (k,0,3) {
          VV(xx,+=d->x[k]*r[dep[k].i])
          VV(yy,+=d->y[k]*r[dep[k].i]) }
        VECT(zz,xx,yy)

        /* another 1/sqrt(SQR(zz)) because will be multiplied by zz */
        /* this is a small correction anyway */
        rr=SQR(zz);
        fx=SCAL(xx,depf)/sqrt(SQR(xx)*rr);
        fy=SCAL(yy,depf)/sqrt(SQR(yy)*rr);
        fz=SCAL(zz,depf)/rr;

        depr=r[d->indx];
        rr=sqrt(rr);
        loop (k,0,3) {
          alpha=fx*d->tx[k]+fy*d->ty[k];

          VVV(f[dep[k].i],+=dep[k].w*depf,+alpha*zz)
          VVV(dr,=r[dep[k].i],-depr)
          /* sum=Tr(Pt)=0 => we omit the virial contribution */
          //          En.vir += alpha*SCAL(dr,zz);
          /* DEBUG!: check OFF order */
          CALCULATEPRESSURETENSOR(alpha,zz,dr) }
        En.vir-=d->wz/rr*fz; // nonzero contribution; divide by rr? (~1)
        CALCULATEPRESSURETENSOR(-d->wz/rr,depf,zz) }

      /* not needed: */
      /*  VO(depf,=0) */
    }
  }
}

int iscube(void) /*************************************************** iscube */
/* returns 1 if the box is a cube (with machine precision) */
{
  int i;

  loopto (i,1,2)
    if (fabs(box.L[0]-box.L[i])/(box.L[0]+box.L[i])>No.eps) return 0;

  return 1;
}

int centergroups(ToIntPtr A,int valence) /********************* centergroups */
/*
  centers all atoms of given valence to its neighbors
  to be used (with valence=option('h') = 4) to remove wrong pyramidal
  configurations of e.g. methane
*/
{
  int n,i,j,a;
  molecule_t *mn;
  siteinfo_t *si;
  vector *r,c;
  double rr;
  int nch=0,nnbr;

  if (valence<2) return 0;

  loop (n,0,No.N) {
    mn=molec+n;
    si=spec[mn->sp]->si;
    r=rof(mn,A->rp);

    loop (i,0,mn->ns) {
      nnbr=0;

      loop (a,0,mn->nc)
        nnbr+=(i==si[a].pair[0] || i==si[a].pair[1]);

      if (nnbr==valence) {
        /*
          if CoM of the neigbors if too far from the central atom,
          move to central atom to the CoM
          (good for wrong methane etc.)
        */
        VO(c,=0)
        rr=0;
        loop (a,0,mn->nc) {
          j=-1;
          if (i==si[a].pair[0]) j=si[a].pair[1];
          if (i==si[a].pair[1]) j=si[a].pair[0];
          if (j>=0) {
            VV(c,+=r[j])
            rr+=SQRD(r[i],r[j]); } }
        VO(c,/=valence)
        if (SQRD(c,r[i])>rr/valence/(valence-1)*0.03) {
          VV(r[i],=c)
          nch++; }

        /*
          if the angles X-center-Y are too small, try to fix this
        */
        loop (a,0,mn->nc) {
          j=-1;
          if (i==si[a].pair[0]) j=si[a].pair[1];
          if (i==si[a].pair[1]) j=si[a].pair[0];
          if (j>=0) {
            int b,k;

            loop (b,0,mn->nc) {
              k=-1;
              if (i==si[b].pair[0]) k=si[b].pair[1];
              if (i==si[b].pair[1]) k=si[b].pair[0];
              if (k>=0 && j!=k) {
                vector drj,drk;
                double cosa;

                VVV(drj,=r[j],-r[i])
                VVV(drk,=r[k],-r[i])
                cosa=SCAL(drj,drk)/sqrt(SQR(drj)*SQR(drk));
                // put(cosa)
                if (cosa>1-1./valence) {
                  nch++;
                  VVV(r[j],=r[i],-drj) } } } } }
        } }
    }

  return nch;
} /* centergroups */

void addshift(int key,int nshift,vector shift) /******************* addshift */
/* key: 0=positions
        1=velocities
        2=accelerations (not used)
   nshift: number of affected molecules
           >0 The shift (push) applies to the first nshift molecules (max No.N)
           <0 The shift (push) applies to the first No.N+nshift molecules,
              the remaining -nshift molecules are shifted (pushed)
              in the opposite direction
*/
{
  int n,i;
  static char *q[3]={"positions","velocities","accelerations"};
  molecule_t *mn;
  vector *r;

  if (shift[0]==0 && shift[1]==0 && shift[2]==0) return;
  if (nshift==0) return; /* to be used for something else? */

  VO(shift,*=powi(h,key))

  loop (n,0,No.N) {
    mn=molec+n;

    r=rof(mn,cfg[key]->rp);

    if (nshift>0) {
      if (n<nshift) loop (i,0,mn->ns) VV(r[i],+=shift) }
    else {
      if (n<nshift+No.N) loop (i,0,mn->ns) VV(r[i],+=shift)
      else loop (i,0,mn->ns) VV(r[i],-=shift) } }


  if (nshift>0) {
    n=min(nshift,No.N);
    prt("%s of %d molecules shifted by [%g %g %g]",q[key],nshift,VARG(shift)); }
  else {
    n=min(nshift+No.N,No.N);
    prt("%s of the first %d molecules shifted by [%g %g %g]\n\
positions of the remaining molecules shifted by -[%g %g %g]\n\
NOTE: shift:=0 (not shifted next time unless you specify the shift again)",
        q[key],nshift,VARG(shift),VARG(shift)); }
  VO(shift,=0)
  if (key==0) normalize(-1);
}

void bounddr(double drmax) /* ************************************** bounddr */
/*
  if the displacement of an atom in one integration step exceeds drmax,
    it is lowered and a warning is printed
  useful for initial equilibrating at high temperature
  turned off for drmax=0, or drmax>0 and init=0,1,2
  negative: always turned on
  should be called BEFORE removedrifts
*/
{
  int n,i,ns;
  molecule_t *mn;
  vector *p;
  double vvmax,q,maxq=9e99;
  int nvv=0;
  static int first=1;

  if (!drmax) return;

  if (drmax>0) {
    if (init<=2 || !thermostat) return; }

  vvmax=Sqr(drmax);
  /* NOTE: cfg[1] contains the difference in 1 step, not the speed */

  loop (n,FROM,No.N) {
    mn=molec+n;
    ns=mn->ns;
    p=rof(mn,cfg[1]->rp);
    loop (i,0,ns) {
      double vv=SQR(p[i]);

      if (vv>vvmax) {
        nvv++;
        q=sqrt(vvmax/vv);
        Min(maxq,q)
        VO(p[i],*=q) } }
    }
  if (nvv) {
    prt("t=%.3f: dr in 1 step reduced for %d site(s), min %g x",t,nvv,maxq);
    if (first) WARNING(("A displacement in 1 MD step was reduced because it exceeded limit |drmax|\n\
*** DO NOT IGNORE THIS MESSAGE, THE EQUATIONS OF MOTION ARE SEVERELY AFFECTED:\n\
*** - kinetic energy is decreased\n\
*** - energy (Hamiltonian) conservation is affected\n\
*** THIS CONDITION IS ACCEPTABLE ONLY DURING UNEQUILIBRATED START AND\n\
*** NEVER SHOULD APPEAR IN A PRODUCTIVE RUN.\n\
*** Remedy:\n\
*** - decrease the timestep h\n\
*** - if you are sure, increase variable drmax or set drmax=0 to turn off\n\
***   the step reduction algorithm\n\
*** Current value:\n\
*** drmax=%g (drmax>0 means drmax:=0 after a sweep ended by ; in the data)",drmax))
    first=0; }

} /* bounddr */


void sortmolecules(int sort) /******************************** sortmolecules */
/***
  sorts molecular positions (not velocities etc.) according to the
  position (|sort|=1=x, 2=y, 3=z) of the center-of-mass;
  sort>0: ascending, sort<0:descending
***/
{
  int i,n,ns,ig;
  molecule_t *mn1,*mn2;
  siteinfo_t *si;
  vector *r1,*r2;
  double sg=sort>0?1:-1;
  double d=-1e-16,q=1+1./No.N,DCMz;
  int KO,pass=-1;

  sort=abs(sort);
  if (sort==0 || sort>3) return;
  sort--; /* now 0=x,1=y,2=z */

  do {
    KO=0;
    pass++;
    d*=q;
    loop (n,0,No.N-1) {

      mn1=molec+n;
      mn2=molec+n+1;
      if (mn1->sp==mn2->sp) {
        ns=mn1->ns;
        si=spec[mn1->sp]->si;
        r1=rof(mn1,cfg[0]->rp);
        r2=rof(mn2,cfg[0]->rp);

        DCMz=0;
        loop (i,0,ns) DCMz+=si[i].mass*(r2[i][sort]-r1[i][sort]);
        if (sg*DCMz<d) {
          loop (ig,0,option('m')) {
            vector *R1=rof(mn1,cfg[ig]->rp);
            vector *R2=rof(mn2,cfg[ig]->rp);
            vector X;

            loop (i,0,ns) {
              VV(X,=R1[i]) VV(R1[i],=R2[i]) VV(R2[i],=X) } }
          KO++; } } }
  } while (KO);

  if (pass) {
    prt("%d molecules of %d species sorted according to %c of CM in %d passes\n",
         No.N,           nspec,                         sort+'x',   pass);
    if (nspec>1) prt("NOTE: each species is sorted separately"); }
  else
    prt("configuration has already been sorted");
}
