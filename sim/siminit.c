/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              initializing tables of species, sites, potentials
             molecule descriptors, configuration, and matrices M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "ground.h"
#include "sds.h"
#include "simglob.h"
#include "interpot.h"
#include "simdef.h"
#include "siminit.h"
#include "simils.h"
#include "simpot.h"
#include "units.h"
#include "simcg.h"
#include "cputime.h"

static double machprec(void) /************************************* machprec */
/*
  worst case one FPU operation relative resolution (error=1/2)
  see also soft/nsk/nsksys.c
*/
{
#define NN 999
  real x[NN];
  real err=0,d;
  int i;

  loop (i,0,NN) x[i]=(real)(i+2)/(2*i+3);
  loop (i,0,NN) {
    d=fabs(x[i]*(2*i+3)-(i+2))/(2*i+3);
    Max(err,d) }

  i=-log(d)/log(2);
  return powi(2,-i);
#undef NN
}

void initNo(void) /************************************************** initNo */
/*
  Calculates No.s, No.N, and No.mass
  To be called (in pass=1) after spec[sp]->N and spec[sp]->ns
  have been read from sysname.def
*/
{
  int sp,nd;
  int downN=FROM,jN,nq=0,nnq; // NB: FROM=0 unless #define LOG
  double molmass,molcharge,sumq=0,molq,err;
  double molalpha; /* POLAR or ECC */

#ifdef POLAR
  int npol;

  No.A=0;
#endif /*# POLAR */

  No.s=No.free_s=No.N=No.free_N=No.depend[0]=No.free_depend=0;
  No.mass=No.free_mass=0;
  No.charge=0;
  No.minmass=3e33;

  if (No.eps==0) {
    /* in 1st pass only */
    No.eps=machprec();
    if (el.epsq<=0) {
#ifdef POLAR
      el.epsq=50*pow(No.eps,0.7); // guess POLAR
#else /*# POLAR */
      el.epsq=pow(No.eps,0.75); // guess nonpolar
#endif /*#!POLAR */
    }
    prt("machine precision=%.3g, el.epsq=%.3g (for charge neutrality test)",No.eps,el.epsq);
    prt("\n\
The following initialization is based on the numbers of molecules\n\
defined optionally in the blend-file and typically in the def-file"); }
  else
    prt("\n\
The following re-initialization is based on the numbers of molecules possibly\n\
changed in the configuration initializer or read from the cfg-file\n\
(as allowed by load.N=%d).",load.N);

  loop (sp,0,nspec) {
    int ns=spec[sp]->ns;
    int N=spec[sp]->N,i;
    depend_t *d;

    molalpha=0;
#ifdef POLAR
    npol=0;
#endif /*# POLAR */
    molmass=molcharge=molq=0;

    loop (i,0,ns) {
      siteinfo_t *si=spec[sp]->si+i;
      sitedef_t *sd=sitedef+si->st;
#ifdef POLAR
      molalpha += sd->alphapol;
      if (sd->alphapol) npol++;
      molcharge += si->chargepol+si->charge;
      if (si->chargepol!=0) nq+=N; /* nq is incl. Drude charges */
#else /*# POLAR */
      molalpha += sd->LJ[0].alpha;
      molcharge += si->charge;
#endif /*#!POLAR */
      if (si->charge!=0) nq+=N;
      molq += Sqr(molcharge);
      if (si->mass) {
        Min(No.minmass,si->mass)
        si->imass=1/si->mass; /* may be changed for equalize<0 */
        No.massy_s+=N; }
      molmass += si->mass; }

    No.A += molalpha*N;
#ifdef POLAR
    No.pol += npol*N;
#endif /*# POLAR */
    No.mass += molmass*N;
    No.charge += molcharge*N;
    sumq += molq*N;
    No.s += ns*N;
    No.N += N;

    jN = N<downN ? 0 : downN>0 ? N-downN : N;
    No.free_N += jN;
#ifdef POLAR
    No.free_pol += npol*N;
    /* there should be also  No.free_A ... */
#endif /*# POLAR */
    No.free_s += ns*jN;
    No.free_mass += molmass*jN;

    /* dependants */
    nd=0;
    looplist (d,spec[sp]->dependants) nd++;
    if (nd && fixsites0) WARNING(("dependants + fixed sites"))
    No.depend[0] += nd*N;
    No.free_depend += nd*jN;

    downN -= N; } /* loop sp */

  prt("No.N=%d molecules  No.s=%d sites (of these, %d dependants)",
            No.N,              No.s,               No.depend[0]);
  prt("total mass of configuration M = %.4f g/mol = %g kg",
                                   No.mass*Munit,No.mass*massunit);
  prt("total charge Q = %g prog.u. = %g e",No.charge,No.charge/electron);
  prt("sum q^2 = %g prog.u. = %g e^2",sumq,sumq/Sqr(electron));
  prt("No.bonded=%d",No.bonded);

  put3(No.depend[1],No.depend[2],No.depend[3])

  err=No.charge/sqrt(sumq);
  if (fabs(err)>el.epsq) {
    prt("*** THE SYSTEM IS CHARGED ***\n\
    charge = %g e, charge/sqrt(sum q^2) = %g, detection limit el.epsq=%g",
        No.charge/electron,err,el.epsq);
    switch (el.bg) {
      case 0:
        ERROR(("The system is charged and el.bg=0: check el.bg and el.epsq."))
        break;
      case 1:
        WARNING(("The system is charged and el.bg=1:\n\
***   Ewald: \"background energy\" will be added\n\
***   Other: no action"
#ifdef SLAB
                 "\n*** The background is not included in charge profiles"
#endif /*# SLAB */
                 ))
        break;
      case 2: {
        double dq=-No.charge/nq;
#ifdef POLAR
        WARNING(("The system is charged and el.bg=2:\n\
*** A neutralizing charge %g will be added to all %d charged sites\n\
*** incl. Drude charges",nq,dq))
#else /*# POLAR         */
        WARNING(("The system is charged and el.bg=2, neutralizing charge %g e\n\
*** will be added to all %d charged sites",dq/electron,nq))
#endif /*#!POLAR         */
        if (nq<2)
          ERROR(("nq=%g is not enough charges to neutralize",nq))

        nnq=0;

        loop (sp,0,nspec) {
          int ns=spec[sp]->ns,i;

          loop (i,0,ns) {
#ifdef POLAR
            if (spec[sp]->si[i].chargepol) {
              spec[sp]->si[i].chargepol+=dq;
              nnq+=spec[sp]->N; }
#endif /*# POLAR */
            if (spec[sp]->si[i].charge) {
              spec[sp]->si[i].charge+=dq;
              nnq+=spec[sp]->N; } } }

        if (nq!=nnq)
          ERROR(("internal: # of charges changed %d -> %d",nq,nnq)) }
        break;
      default:
        ERROR(("el.bg=%d is invalid",el.bg))
    }
  }

  /*
     counting ions (after optional neutralizing system charge)
     the original No.charge is kept (the new should be zero)
     BUG: zero charge and nonzero Drude charge not treated correctly
  */
  nnq=0;
  No.ion=0;
  molq=0; /* used for total charge */
  loop (sp,0,nspec) {
    int ns=spec[sp]->ns;
    int N=spec[sp]->N,i;

    molcharge=0;

    loop (i,0,ns) {
      siteinfo_t *si=spec[sp]->si+i;
#ifdef POLAR
      molcharge += si->chargepol;
#endif /*# POLAR */
      molcharge += si->charge;
      if (si->charge) nnq+=N; }

    molq+=N*molcharge;

    if (fabs(molcharge)>=0.01*electron) No.ion+=N; }

  prt("number of ions No.ion=%d, system charge=%g e",No.ion,No.charge/electron);

  if (Sqr(molq)/sumq>Sqr(el.epsq) && el.bg==2)
    ERROR(("Has been neutralized, but still charged - check el.epsq"))

  if (FROM>0)
    prt("free %d molecules  %d sites (%d dependants) %.4f g/mol",
	No.N-FROM,No.free_s,No.free_depend,No.free_mass*Munit);
  if (No.N<=0) { ERROR(("no molecule - nothing to do")) exit(-1); }
  if (FROM>=No.N)
    ERROR(("No.first=%d > # of molecules N=%d",FROM,No.N))
#ifdef POLAR
  prt("total polarizability = %g AA^3",No.A);
#endif /*# POLAR */
#ifdef ECC
  prt("total polarizability= %g AA^3 (from combining-rule column)",No.A);
#endif /*# ECC */

  prt("Minimum mass = %g p.u., total mass = %g p.u., # of massy sites = %d",
      No.minmass, No.mass, No.massy_s);

  if (equalize.cfg) {
    /* Global equalization of masses of species equalize.sp
       or range 0..|equalize.sp|
       See simdef.c: readblend() for equalization in molecules */
    double M=0, msum=0;
    int nM=0;
    int spfrom=equalize.sp;
    int spto=equalize.sp;
    static int equalized; /* flag to make global equalization only once */

    if (equalized) ERROR(("IMPLEMENTATION LIMITATION\n\
*** equalize.cfg=%g is set and equalizecfg() is about to be called for the 2nd\n\
*** time. This happens if the number of molecules changes after reading the\n\
*** cfg-file or for init=\"crystal\". Use the correct number of molecules in\n\
*** thedef-file, or use equalize.mol instead.",equalize.cfg))

    if (equalize.sp<0) {
      /* big negative equalize.sp changed to -(nspec-1) -> whole cfg */
      spfrom=0;
      if (equalize.sp<=-nspec) equalize.sp=1-nspec;
      spto=-equalize.sp; }

    if (spfrom<0 || spto>=nspec)
      ERROR(("equalize.sp = %d for nspec = %d refers to species out of range",equalize.sp,nspec))

    if (equalize.cfg<0 || equalize.cfg>1)
      WARNING(("equalize.cfg = %g out of range [0,1]",equalize.cfg))

    /* calculate the sum and number of masses to equalize */
    loopto (sp,spfrom,spto) {
      int ns=spec[sp]->ns;
      int N=spec[sp]->N,i;

      loop (i,0,ns) {
        siteinfo_t *si=spec[sp]->si+i;

        if (si->mass) {
          M+=si->mass*N;
          nM+=N; } } }

    M/=nM;

    /* equalize */
    loopto (sp,spfrom,spto) {
      int ns=spec[sp]->ns;
      int N=spec[sp]->N,i;

      spec[sp]->mass=0;

      loop (i,0,ns) {
        siteinfo_t *si=spec[sp]->si+i;
        if (si->mass) {
          if (option('v')&2) prt_("%d.%d: %.4f", sp,i, si->mass*Munit);
          si->mass=equalize.cfg*M+(1-equalize.cfg)*si->mass;
          if (option('v')&2) prt(" --> %.4f g/mol", si->mass*Munit);
          spec[sp]->mass+=si->mass;
          si->imass=1/si->mass; }
        msum+=si->mass*N; } }
    put3(No.massy_s,No.mass,msum)

    equalized=1;
    prt("\
Masses of all atoms in species %d..%d (incl.) equalized with factor %g\n\
(equalize=0: original masses, equalize=1: all masses equal)\n\
Mean mass in the above range = %g g/mol\n\
The total mass of the configuration is unchanged\n\
WARNING: kinetic quantities are affected",spfrom,spto,equalize.cfg,M*Munit);

    /* minimum mass No.minmass and mass fool-proof check */
    msum=0;
    No.minmass=3e33;
    loop (sp,0,nspec) {
      int ns=spec[sp]->ns;
      int N=spec[sp]->N,i;

      loop (i,0,ns) {
        siteinfo_t *si=spec[sp]->si+i;
        if (si->mass) {
          Min(No.minmass,si->mass)
          msum+=si->mass*N; } } }
    put3(No.massy_s,No.mass,msum)

    if (fabs(No.mass-msum)/No.mass>No.eps*sqrt(No.s)*3)
      ERROR(("doublecheck of the total mass %.15g->%.15g after equalization failed",No.mass,msum))
  }
}

void printmasses(void) /**************************************** printmasses */
/*
  print all masses (in g/mol)
*/
{
  int n,i,ns;
  molecule_t *mn;
  siteinfo_t *si;

  prt("! masses of all sites\n\
! sp n: mass mass ... [g/mol]");

  loop (n,0,No.N) {
    mn=molec+n;
    si=spec[mn->sp]->si;
    ns=mn->ns;

    prt_("%d %d:",mn->sp,n);
    loop (i,0,ns) prt_(" %g",si[i].mass*Munit);
    _n }
}


#ifdef ECC
double rescalecharges(int ecc,double epsf) /***************** rescalecharges */
/* ECC only, EXPERIMENTAL */
{
  double scale;
  int sp;

  if (!ecc) return 1;

  underline("Electronic Continuum Correction");
  prt("WARNING: EXPERIMENTAL CODE");
  prt("el.ecc=%d (%s)",el.ecc,abs(el.ecc)==1?"ions":"dipolar molecules");
  prt("\
* Ewald required\n\
* combination of dipoles and ions is not supported\n\
* only NVT simulation is supported\n\
* charge scaling factor calculated from el.epsf is applied\n\
* Lorentz polarization energy for ions or dipoles is added to elst energy;\n\
  via the virial of force, this becomes a part of pressure\n\
* pressure corrections (difference frome naive formula) is printed separately");

  if (epsf)
    prt("ECC: el.epsf=%g specified",epsf);
  else {
    epsf=(4*PI/3)*No.A/box.V;
    epsf=(1+2*epsf)/(1-epsf);
    prt("ECC: el.epsf=%g calculated from total polarizability",epsf); }

  switch (abs(ecc)) {
    case 1: scale=1/sqrt(epsf); break;
    case 2: scale=(2+epsf)/sqrt(9*epsf); break;
    default: ERROR(("wrong el.ecc=%d",el.ecc)) }

  if (ecc<0)
    prt("ECC: scaling factor %g (already scaled or taken from ble-file)",scale);
  else {
    loop (sp,0,nspec) {
      int ns=spec[sp]->ns,i;

      loop (i,0,ns)
        if (spec[sp]->si[i].charge)
          spec[sp]->si[i].charge*=scale; }
    prt("ECC: all charges have been scaled by %g",scale); }
  _n

  return epsf;
}
#endif /*# ECC */

void makemolgol(int nmol) /************************************** makemolgol */
{
  int sp,err=0,nn;
  char cmd[2048];
  char *nm=simils.simname;

  if (nmol<=0) return;
  Min(nmol,No.N)

  strcpy(cmd,"molcfg ");
  loop (sp,0,nspec) {
    if (!strcmp(spec[sp]->name,simils.simname)) err++;
    nn=min(nmol,spec[sp]->N);
    if (nn) sprintf(strend(cmd),"-%d:%c%%03d %s ",
                    nn, sp+'a', spec[sp]->name);
    nmol-=nn; }

  if (err) {
    WARNING(("simulation name %s matches (one of) species name(s)\n\
*** CONFIG.mol and .gol will be generated instead of %s.mol",simils.simname,simils.simname))
    nm="CONFIG"; }

  strcat(cmd,nm);
  underline("show support");
  prt("Command to generate .mol and .gol files for show will be called:\n  %s",cmd);
  putenv("MOLCFG=1");
  sp=system(cmd);
  if (sp) WARNING(("%s.{mol,gol} not generated\n\
*** (missing .mol and .gol files for species?)",nm))
  else
    prt("Files %s.{mol,gol} have been generated.\n\
To show the trajectory, use command\n  show %s", nm,nm);
}

void initpot(void) /*********************************************** initpot */
{
  int i,j;

#ifdef LINKCELL
  ralloc(pot1,nspec*sizeof(pot1_t*));
#endif /*# LINKCELL */
  ralloc(pot,nspec*sizeof(pot2_t**));
  loop (i,0,nspec) {
#ifdef LINKCELL
    pot1[i]=intramol;
#endif /*# LINKCELL */
    ralloc(pot[i],nspec*sizeof(pot2_t*));
    /* old COOK: charmm2; */
    loop (j,0,nspec) pot[i][j]=Potential(i,j); }
}

void initmolecules(int corr) /******************************** initmolecules */
/*
  allocate configurations a[] with pointers to molecules
  allocate M for constraint dynamics
  initialize degrees of freedom
*/
{
  int sp,n=0,irp=0,i,Mlen=0,sub,tsub;

  if (RPOFFSET!=(real*)(cfg[0]->rp)-(real*)&(cfg[0]->logs))
    ERROR(("wrong RPOFFSET, real padding in structure ToInt, or other declaration problems"))

  No.c=No.maxc=No.maxs=No.free_c=0;

  rallocarrayzero(molec,No.N); /* 11/2004 because of ANCHOR */
  ralloc(Ms,(option('c')&1?No.N:nspec)*sizeof(M_t*));

  underline("initializing molecules");
  if (option('v')&2) header("species     N     ns     nc   Mlen");
  loop (sp,0,nspec) {
    specinfo_t *s=spec[sp];

    if (s->N) {
      if (!(option('c')&1)) Ms[sp]=defineM(s,&Mlen);

      loop (i,0,s->N) {
        molecule_t *m=molec+n;

        /* INEFFICIENT! should allocate (by molecules) in advance/bigger chunks */
        if (option('c')&1) Ms[n]=defineM(s,&Mlen);
        m->sp=sp;
        m->ns=s->ns;
        m->nc=s->nc;
        m->ir=sizeof(vector)*irp;
        m->ig=No.c;

        Max(No.maxs,m->ns)
        irp += m->ns;
        Max(No.maxc,m->nc)
        No.c += m->nc;
        if (n>=FROM) No.free_c += m->nc;
        n++; } }
    if (option('v')&2)
      prt("%6d %6d %6d %6d %6d", sp,s->N,s->ns,s->nc,Mlen); }
  if (option('v')&2) header("");

  if (irp!=No.s) ERROR(("internal: irp=%d != No.s=%d",irp,No.s))
  No.nreal=No.s*DIM;
#ifdef POLAR
  No.nreal*=2;
#endif /*# POLAR */

  No.eq=No.s*DIM+RPOFFSET; /* # of eqs of motion (incl. logs, lambda) */

/* DEGREES OF FREEDOM
   ^^^^^^^^^^^^^^^^^^
   No.f = DIM*No.s
         -No.c  constraints
         +1     additional degree of freedom for Nose
         -1     energy (Hamiltonian) conservation
         -DIM*(# of dependants) for all sites lin. dependent on other sites
         -(fixed degrees of freedom, see variables  conserved,drift)

   WARNING: not correct for a few fixed sites

   Intramolecular degrees of freedom (No.f_in) contain all moves relative to
   the center of mass (incl. rotations).

   Intermolecular degrees of freedom (No.f_tr) contain only translational
   moves of centers of mass and the above corrections apply BUT -1 for energy
   conservation if no Nose thermostat is used.

   The Maxwell-Boltzmann, Andersen, (Langevin) thermostats
   sample also translations and rotations.
*/

  sub=0;
  if ((drift&DRIFT_WHEN)!=DRIFT_START)
    for (i=DRIFT_VX; i<=DRIFT_AZ; i*=2) if (drift&i) sub++;
  prt("Calculated number of conserved degrees of freedom = %d\n\
(includes conserved momenta and angular momenta, as calculated from drift=%d;\n\
Hamiltonian conservation is treated separately).",sub,drift);
  if (conserved<0)
    prt("\
This value will be subtracted from the total number of degrees of freedom as\n\
well as intermolecular (based on center-of-mass) number of degrees of freedom.");
  else {
    if (sub!=conserved)
      prt("\
WARNING: specified conserved=%d does not match the calculated value,\n\
         I assume that you have calculated it correctly.",conserved);
    sub=conserved; }

  tsub=0; /* for Nose, Maxwell/Andersen */
  if (thermostat<=T_BERENDSEN || (thermostat>=T_TR && thermostat<=T_FRICTIONS)) {
    /* Berendsen-like or NVE */
    if (corr&64) {
      tsub=1; /* energy conservation */
      prt("\
INFO: corr&64 set (for %s):\n\
  One degree of freedom for energy conservation subracted: Tkin is affected.\n \
  This finite-size correction is usually less accurate - consult the manual.\n\
  in addition, the translational temperature (from centers-of-mass) Ttr and\n\
  the intramolecular and rotational temperature Tin for an equilibrated system,\n \
  will by ~O(1/N) differ from the total kinetic temperature Tkin",
          thermostat?"Berendsen thermostat":"NVE simulation"); }
    else
      prt("INFO: corr&64 is not set - usually a good choice"); }

  No.f0=DIM*(No.free_s-No.free_depend)-No.free_c;
#ifdef POLAR
  if (option('p')/10%10==0) {
    /* Car-Parrinello-like */
    No.f0+=DIM*No.free_pol;
    prt("Car-Parrinello-like: %d*%d mechanical degrees of freedom for dipoles added",
	No.free_pol,DIM); }
#endif /*# POLAR */
  No.f=No.f0-sub;
  /* No.f=degrees of freedom without energy conservation */

  put3(No.free_s,No.free_c,No.free_depend)
  if (No.f<=0) {
    WARNING(("%d degrees of freedom: replaced by 1, and E conservation not included",No.f))
    tsub=0;
    No.f=1; }

  /* inter- and intramolecular degrees of freedom (without energy conservation) */
  No.f0_tr=DIM*No.free_N;
  No.f_tr=No.f0_tr-sub;
  No.f_in=No.f-No.f_tr;

  No.f-=tsub; /* (rather heuristic or approximate) correction for energy conservation */

  No.conserved=sub+tsub;

  put3(No.s,No.c,No.f)
  put3(No.s*3,No.eq,No.nreal)
  put3(No.free_N,No.f_tr,No.f_in)
  put3(No.conserved,No.f0,No.f0_tr)

  No.fx=No.f;

  if (thermostat>=T_NPT) {
    No.NPT=1.+(thermostat==T_NPT)*3./No.f; /* kinetic pressure correction for MTK NPT */
    No.invf=(double)(thermostat==T_NPT)/No.f; /* Ekin multiplication factor */
    if (rescale&RESCALE_CM) {
      ERROR(("rescale=%d includes CM-based rescaling, which is not supported with NPT/MTK\n\
*** Hints:\n\
***   for isotropic barostat, use rescale=\"xyz\"\n\
***   to maintain Pzz=P, use rescale=\"Z\"",rescale))
      /* rescale-=RESCALE_CM; */
      }
    if (tau.P<.999e6) {
      if (rescale&RESCALE_PT) No.fx+=No.ncoord;
      else No.fx++; }
    else
      WARNING(("large tau.P>=1e6 forces Nose (no extra degrees of freedom added)\n\
*** but term PV is still included in total energy"))
    No.M_T=Sqr(tau.T)*T*No.f;
    No.M_Th=No.M_T/2;
    No.M_P=Sqr(tau.P)*T*(No.f+3)*No.ncoord;
    No.M_Ph=No.M_P/2;
    put3(No.ncoord,No.fx,No.NPT)
    put2(No.M_T,No.M_P) }

  if (tau.rho<0) {
    if (thermostat>=T_NPT) ERROR(("tau.rho=%g and NPT not allowed",tau.rho))
    if (thermostat) WARNING(("tau.rho=%g and thermostat unknown interaction",tau.rho))
    if (rescale&RESCALE_CM) {
      ERROR(("tau.rho=%g: rescale=%d includes CM-based rescaling, which is not supported\n\
*** with V(t) or L(t) control. Use rescale=\"xyz\", \"Z\", etc.",tau.rho,rescale))
      rescale-=RESCALE_CM; } }

/*** physical configuration allocation ***/
#ifdef POLAR
  i=sizeof(ToInt)-sizeof(vector)+No.s*2*sizeof(vector);
  polar_off=No.s*sizeof(vector);
#else /*# POLAR */
  i=sizeof(ToInt)-sizeof(vector)+No.s*sizeof(vector);
#endif /*#!POLAR */
  if (option('r')>option('m')) ERROR(("option -r%d exceeds option -m%d",option('r'),option('m')))
  /* probably too many of them - why? */
  loop (n,0,option('m')+3) sdsralloczero(cfg[n],i)
} /* initmolecules */


void initgammas(void) /****************************************** initgammas */
{
  int gssize=No.c*sizeof(gs[0][0])+(No.c==0);

  if (No.pred) {
    int i; int b= -1;

    ralloc(gs,No.pred*sizeof(gs[0]));
    ralloc(binom,No.pred*sizeof(binom[0]));
    loop (i,0,No.pred) {
      ralloczero(gs[i],gssize);
      binom[i] = b= -b*(No.pred-i)/(i+1); } }
} /* initgammas */


double cutcorr(double LJcutoff,int corr) /************************** cutcorr */
/***
    Standard bulk fluid cutoff correction of pair Lennard-Jones-like
    forces over all intermolecular pairs is calculated.
    Intramolecular pairs are not included!
    The result should be divided by V to obtain cutoff correction of Epot
    The result should be divided by V^2 to obtain cutoff correction of P
    In addition, # of pairs for rdf histograms are calculated
    (incl. intramolecular pairs more distant than 1--4)
    (elst forces not included)
***/
{
  int n,m,i,j,si,sj;
  double np;
  double c,fcorr=0;
  rdf_t *rdfij;

  if (LJcutoff==0) return 0;

  loop (n,0,nspec) {
    loopto (m,0,n) {
      if (m==n) {
        if (corr&32) np=(double)spec[n]->N*spec[n]->N/2.;
        else np=(double)spec[n]->N*(spec[n]->N-1)/2.; }
      else np=(double)spec[n]->N*spec[m]->N;
      c=0;
      loop (i,0,spec[n]->ns) loop (j,0,spec[m]->ns) {
	c += sstab [si=spec[n]->si[i].st] [sj=spec[m]->si[j].st] .corr;
	if ( (rdfij=rdf[si][sj]) ) rdfij->npair+=np;
        /* this is OK because rdf[si][sj] and rdf[sj][si]
           point to the same rdf_t */ }
      fcorr+=c*np; }

    /* now intramolecular pairs but 1-2,1-3, and 1-4 */
    /* cutoff corr not included because it is not clear how */

    loop (i,0,spec[n]->ns) {
      exception_t *exc=spec[n]->si[i].exc;
      int j0=0,j1;

      do {
	j1=exc->indx;

	if (exc->type==ONEFOUR)
	  if ( (rdfij=rdf[spec[n]->si[i].st][spec[n]->si[j1].st]) )
	    rdfij->npair += spec[n]->N;

	loop (j,j0,j1)
	  if ( (rdfij=rdf[spec[n]->si[i].st][spec[n]->si[j].st]) )
	    rdfij->npair += spec[n]->N;

	exc++;
      } while ( (j0=j1+1)<i );
    }
  } /* n */

#ifdef DEBUG
  {
    int n,m,sp;
    molecule_t *mn,*mm;
    double debugcorr=0;
    vector *auxf;

    measure=2;

    loop (n,0,No.N) {
      mn=molec+n;
      sp=mn->sp;

      allocarrayzero(auxf,mn->ns);

      loop (m,0,n) {
	mm=molec+m;
	debugcorr += (*pot[sp][mm->sp]) (auxf,auxf,mn,mm,cfg[0]->rp); }

      (*pot[sp][mn->sp]) (auxf,auxf,mn,mn,cfg[0]->rp); /* ugly, but it works */
      free(auxf); }

    put2(fcorr,debugcorr)
    if (fabs(fcorr/debugcorr-1)>1e-10) Error("");
  }
#endif /*# DEBUG */
  return fcorr;
} /* cutcorr */

#ifdef WIDOM
double Widomcutcorr(int spreal,int spvirt,double LJcutoff) /****************** Widomcutcorr */
/***
    Standard bulk fluid cutoff correction of pair Lennard-Jones-like
    forces over intermolecular pairs of virtually inserted molecule of
    species spvirt is calculated.
    If spreal>=0, the analogous term for species spreal is subtracted
    so that the correction for spreal->spvirt identity change is returned.
    Intramolecular pairs are not included!
    The result should be divided by V to obtain the cutoff correction
    of the chemical potential by the Widom method.
***/
{
  int m,i,j,si,sj,np;
  double c,corr=0;

  if (spvirt>=nspec) ERROR(("widom.sp undefined (>=nspec)"))
  if (spreal>=nspec) ERROR(("widom.spreal undefined (>=nspec)"))

  loop (m,0,nspec) {
    np=spec[m]->N;
    c=0;
    loop (i,0,spec[spvirt]->ns) loop (j,0,spec[m]->ns)
      c += sstab [si=spec[spvirt]->si[i].st] [sj=spec[m]->si[j].st] .corr;
    corr+=c*np; }

  if (spreal>=0) loop (m,0,nspec) {
    np=spec[m]->N;
    if (m==spreal) np--;
    c=0;
    loop (i,0,spec[spreal]->ns) loop (j,0,spec[m]->ns)
      c += sstab [si=spec[spvirt]->si[i].st] [sj=spec[m]->si[j].st] .corr;
    corr-=c*np; }

  if (spreal>=0)
    prt("\nWidom cutoff corection for species id-change [%d]->[%d] for a homogeneous box",
	spreal,spvirt);
  else
    prt("\nWidom cutoff corection for species [%d] for a homogeneous box",spvirt);

  prt("corr = %g/V = %g\n\
Widom correcting Boltzmann factor = %g\n\
(the latter two values will change if the box size changes)\n",
      corr,corr/(L[0]*L[1]*L[2]),
      exp(-corr/(L[0]*L[1]*L[2]*T)));

  return corr;
} /* Widomcutcorr */
#endif /*# WIDOM */

double setL(vector L,double rho) /************************************* setL */
/*
  rescale L, or calculate missing L[] from rho given
  returns final rho
*/
{
  int i,nzero=0;
  double x,p=1;

  loop (i,0,3) if (L[i]<0) ERROR(("negative L[%d]=%g",i,L[i]))

  if (rho*(L[0]*L[1]*L[2])!=0) {
    /* both rho[initrho] and L given: scaling */
    x=cbrt(No.mass*rhounit/(L[0]*L[1]*L[2])/rho);
    if (fabs(x-1)>5e-16) {
      prt("Current L=[%g %g %g] will be rescaled to reach %g kg/m3.",VARG(L),rho);
      VO(L,*=x) } }
  else {
    loop (i,0,3) if (L[i]==0) nzero++; else p*=L[i];

    if (nzero==0) {
      /* all L[] given ==> calculate rho */
      rho=No.mass*rhounit/(L[0]*L[1]*L[2]);
      prt("Density %g [kg/m3] was calculated from box sides",rho); }
    else {
      if (rho==0) {
#ifdef FREEBC
        WARNING(("The reference density was set to rho=1000\n\
*** It is used by the initializer only and does not affect calculations"))
#else /*# FREEBC */
        ERROR(("I cannot calculate the reference box because both density (rho) and\n\
*** %d box sides(s) (L[]) are undefined (zero).\n\
*** (This is no longer warning to force users to specify the density;\n\
*** nevertheless, ignoring this error will continue as if rho=1000.)",nzero))
#endif /*#!FREEBC */
        rho=1000; }
      /* calculate missing L[] from rho */
      x=pow(No.mass*rhounit/rho/p,1./nzero);
      loop (i,0,3) if (L[i]==0) L[i]=x;

      prt("%d missing box sides calculated from density rho=%g [kg/m3]",
           nzero,                                           rho); } }

  prt("Density doublecheck (calculated from new L[]): current rho=%g [kg/m3]",
    No.mass*rhounit/(L[0]*L[1]*L[2]));

  return rho;
}

double initcutoff(double cutoff,vector L) /********************** initcutoff */
/*
  L is the reference L
*/
{
  double minL=9e99;
  int i;

  prt("Truncation of electrostatics forces and maximum cutoff for LJ and RDF:");

#ifdef FREEBC
  if (cutoff<=0) cutoff=9e9;
  prt("cutoff = %g A",cutoff);
#else /*# FREEBC */

  loop (i,0,DIM) Min(minL,L[i])

  if (cutoff<=0) {
#  ifdef NIBC
    cutoff=9e9;
#  elif COULOMB<0
    /* Ewald: L/2 break-even for 1000 atoms (changed from 2000 in V2.7d)*/
    cutoff=cbrt(PROD(L)/8)*pow(No.s/1000.0,-1.0/6);
    prt("default cutoff for Ewald summation estimated:\n\
WARNING: may be far from optimum! (consider the test module, see el.test)");
#  else /*#!NIBC!COULOMB<0 */
    cutoff=12;
    prt("default cutoff for short-range electrostatics set:");
#  endif /*#!NIBC!COULOMB<0 */
#  ifndef PERSUM
    Min(cutoff,minL/2)
#  endif /*# PERSUM */
  }

  prt("cutoff = %.8f  cutoff/min(L)=%g",cutoff,cutoff/minL);

  loop (i,0,DIM)
    if (cutoff>L[i]/2)
      prt("WARNING: cutoff=%g > L%c/2=%g (by %.3g%%)",
          cutoff,'x'+i,L[i]/2,(cutoff/(L[i]/2)-1)*100.);
#endif /*#!FREEBC */

#ifdef PERSUM
  if (cutoff/minL<=0.5) WARNING(("cutoff/min(L)<=0.5: PERSUM is not necessary"))
  loop (i,0,DIM)
    No.nimg[i]=(int)((cutoff+2*No.molspan)/L[i]+1);
  prt("PERSUM setup: No.nimg = [ %d %d %d ]",No.nimg[0],No.nimg[1],No.nimg[2]);
#endif /*# PERSUM */

  return cutoff;
}

/**************************** fixed sites support ****************************/

int initfix(char *fn) /********************************************* initfix */
{
  static int from;
  FILE *fs=fopen(fn,"rt");
  char line[128],*s;
  fixsites_t *f;
  int i,to,isum=0;

  if (!fs) {
    ERROR(("no %s file -- no sites kept fixed, check option -k",lastFn))
    return 0; }

  if (FROM)
    ERROR(("combination of option -j and file %s not supported",lastFn))
      
  while (fgets(line,128,fs)) if (!strchr("!#",line[0])) {
    s=strtok(line," \t\n"); if (!s) continue;
    if (!isdigit(s[0])) ERROR(("%s: %s not a number",lastFn,s))
    from=atoi(s);
    to=from;
    s=strchr(s,'-');
    if (s) sscanf(s+1,"%d",&to);
    if (to<from || from<0)
      ERROR(("%s: %d-%d bad order or negative site",lastFn,from,to))

    alloczero(f,sizeof(fixsites_t));
    f->from=from;
    to++; /* last not incl. */
    f->to=to;
    
    while ( (s=strtok(NULL," \t\n")) && f->isr<3)
      f->r[f->isr++]=atof(s);
    if (f->isr%DIM)
      ERROR(("%s: format of optional site %d position",lastFn,from))

    loop (i,from,to) {
      if (isfixed(i)) prt("%s: %d multiply listed",lastFn,i);
      else isum++; }
    
    f->next=fixsites0; fixsites0=f; }

  fclose(fs);
  if (!isum) {
    WARNING(("%s: no site kept fixed",lastFn))
    from=-1;
    return 0; }
  prt("%s: %d sites kept fixed",lastFn,isum); 

  return from;
}

int isfixed(int i) /************************************************ isfixed */
{
  fixsites_t *f;

  for (f=fixsites0; f; f=f->next) if (i>=f->from && i<f->to) return 1;

  return 0;
}

void checkfixed(void) /****************************************** checkfixed */
{
  fixsites_t *f;

  for (f=fixsites0; f; f=f->next) if (f->to>No.s)
    ERROR(("fixed site %d out of range",f->to))
}

#ifdef ANCHOR
#  define ALINELEN 1024
void initanchor(char *fn,int init,double drmax) /**************** initanchor */
{
  FILE *fs;
  char line[ALINELEN];
  int i,n,k,iaxis[3],col;
  vector r;
  molecule_t *mn;
  int keep[2];
  int STAT[64]; /*to print statistic of different ANCHOR modes */
  static int pass;

  if (option('k')>=0 || option('k')<=-16) {
    WARNING(("This is the ANCHOR version of cook, but option -k%d implies\n\
*** %s",option('k'),option('k')>0?"harmonic springs instead of SHAKE-based constraints":option('k')?"no anchoring":"another algorithm and file format"))
    return; }

  if (option('m')!=2) ERROR(("Anchoring requires -m2 (Verlet+SHAKE), cf. option -k"))

  loop (i,0,64) STAT[i]=0;

  prt("reading %s",fn);

  loop (n,0,No.N) {
    molec[n].anchor=0;
    molec[n].xyz=7; /* probably not needed (now) */ }

  if ( !(fs=fopen(fn,"rt")) ) {
    ERROR(("no %s file requested by option -k",lastFn))
    return; }

  while (fgets(line,ALINELEN,fs)) if (isalpha(line[0])) {
    char *c=line;
    int anch=0,xyz=0;

    while (isalpha(*c)) c++;

    iaxis[0]=iaxis[1]=iaxis[2]=0;

    STAT[((unsigned char*)line)[0]&63]++;

    switch (line[0]) {
      case 'v': anch=ANCHOR_v; goto readn;
      case 'r': anch=ANCHOR_r; goto readn;
      case 'x': anch=ANCHOR_x; goto readn;
      case 'f': anch=ANCHOR_f; goto readn;
      case 'm':	anch=ANCHOR_m;
      readn:
        if (1!=sscanf(c,"%d",&n))
	  ERROR(("%s: missing molecule number in line:\n%s",lastFn,line))
	break;
      case 's':
        anch=ANCHOR_s; goto s;
      case 'a':
        anch=ANCHOR_a;
      s:
        if (5!=sscanf(c,"%d%d%lf%lf%lf",&n,iaxis,r,r+1,r+2))
          ERROR(("%s: not enough data in line:\n%s",lastFn,line))
        break;
      case 'p':
        anch=ANCHOR_p; /* bug fixed V3.0f */
        if (6!=sscanf(c,"%d%d%d%lf%lf%lf",&n,iaxis,iaxis+1,r,r+1,r+2))
          ERROR(("%s: not enough data in line:\n%s",lastFn,line))
        break;
      case 't':
        anch=ANCHOR_t;
        if (7!=sscanf(c,"%d%d%d%d%lf%lf%lf",&n,iaxis,iaxis+1,iaxis+2,r,r+1,r+2))
          ERROR(("%s: not enough data in line:\n%s",lastFn,line))
        break;
      case 'C':
        if (3!=sscanf(c,"%lf%lf%lf",anchor.r0,anchor.r0+1,anchor.r0+2))
          ERROR(("%s: not enough data in line:\n%s",lastFn,line))
        anchor.xyz=(anchor.r0[0]>=0) + (anchor.r0[1]>=0)*2 + (anchor.r0[2]>=0)*4;
        n=-0x7fffffff; /* no molecule */
        break;
      case 'c':
        anch=ANCHOR_c;
        goto c;
      case 'i':
        anch=ANCHOR_i;
      c:
        if (4!=sscanf(c,"%d%lf%lf%lf",&n,r,r+1,r+2))
          ERROR(("%s: not enough data in line:\n%s",lastFn,line))
        /* void for ANCHOR_i */
        xyz=(r[0]>=0) + (r[1]>=0)*2 + (r[2]>=0)*4;
        break;
      case 'g': {
        struct anchorgroup_s *group;
        char *tok=strtok(line," \t\n");
        int i,ns;

        tok=strtok(NULL," \t\n"); if (!tok) ERROR(("%s: group: not enough data",lastFn))
        n=atoi(tok);
        if (n<0 || n>=No.N)
          ERROR(("%s: group: %d is invalid molecule number",lastFn,n))
        mn=molec+n;
        tok=strtok(NULL," \t\n"); if (!tok) ERROR(("%s: group: not enough data",lastFn))
        ns=atoi(tok);
        if (ns<1 || ns>=mn->ns)
          ERROR(("%s: group: ns=%d out of range",lastFn,ns))

        anch=ANCHOR_g;

        ralloc(group,sizeof(*group)+(ns-1)*sizeof(group->site[0]));
	group->ns=ns;
        loop (i,0,ns) {
          tok=strtok(NULL," \t\n"); if (!tok) ERROR(("%s: group: not enough data",lastFn))
          group->site[i]=atoi(tok);
          if (group->site[i]<0 || group->site[i]>=mn->ns)
            ERROR(("%s: group: site %d out of range",lastFn,group->site[i])) }

        loop (i,0,DIM) {
          tok=strtok(NULL," \t\n"); if (!tok) ERROR(("%s: group: not enough data",lastFn))
          group->r0[i]=atof(tok); }

        group->next=mn->group;
        mn->group=group; }
	break;
      default:
        ERROR(("%s: wrong key in line:\n%s",lastFn,line)) }

    if (n!=-0x7fffffff) {
      if (n<0 || n>=No.N)
        ERROR(("%s: n=%d is invalid molecule number",lastFn,n))

      mn=molec+n;

      if ( (anch & ANCHOR_g) && (mn->anchor & (ANCHOR_pos|ANCHOR_axis))
       ||  (mn->anchor & ANCHOR_g) && (anch & (ANCHOR_pos|ANCHOR_axis)) )
        ERROR(("%s: constrain group and axis/cm/site combined for molecule %d",lastFn,n))
      if ( (anch & ANCHOR_pos) && (mn->anchor & ANCHOR_pos) ) {
        if ((-option('k'))&8)
          WARNING(("%s: constrain position doubly defined for molecule %d",lastFn,n))
        else
          ERROR(("%s: constrain position doubly defined for molecule %d (consider -k8)",lastFn,n)) }
      if ( (anch & ANCHOR_axis) && (mn->anchor & ANCHOR_axis) )
        ERROR(("%s: constrain direction doubly defined for molecule %d",lastFn,n))

      mn->anchor|=anch;
      if (anch==ANCHOR_c) mn->xyz=xyz;

      if ( (mn->anchor & ANCHOR_s) && (mn->anchor & ANCHOR_axis) )
        ERROR(("%s: combination of site and axis constraining for molecule %d",lastFn,n))

      if (anch & ANCHOR_pos) {
        if (ANCHOR_s) mn->iaxis[0]=iaxis[0];
        VV(mn->r0,=r) }
      else if (anch & ANCHOR_axis) {
        double x=sqrt(SQR(r));
        if (x==0)
          ERROR(("%s: zero vector for molecule %d",lastFn,n))
        copy(mn->iaxis,iaxis,sizeof(iaxis));
        VVO(mn->axis,=r,/x) } } }

  fclose(fs);

  /* setting bounds to moves caused by constraints - cheap */
  anchor.rr=Sqr(drmax)/40;
  anchor.sin=sqrt(anchor.rr);
  anchor.cos=1-anchor.rr;

  if (drmax==0 || (drmax>0 && (init<=2 || !thermostat))) {
    anchor.rr=9e99;
    anchor.sin=2; anchor.cos=-2; }

  n=0;
  loop (i,0,64) n+=STAT[i];
  if (n) {
    prt_("ANCHOR:");
    loop (i,0,64) if (STAT[i]) prt_(" %c=%d",i+'@',STAT[i]);
    _n }

  if ((-option('k'))&4) {
    anchor.f=fopen(Fn("anc"),init<2?"at":"wt");
    if (!anchor.f) ERROR(("cannot write to %s",lastFn)) }
  else
    anchor.f=NULL;

  if (init>=2) {
    /* active with ascii output (-k-4); for -k-3 see the block below */
    col=0;
    if (anchor.f) fprintf(anchor.f,"# UNITS: [force]=N, [torque]=N.m, [length]=AA, [time]=ps\n#");
    loop (n,0,No.N) {
      mn=molec+n;
      if (mn->anchor & ANCHOR_g) {
	struct anchorgroup_s *group;

	looplist (group,mn->group) {
	  if (anchor.f) fprintf(anchor.f,"--constr.f. mol.%3d group ns=%d(%2d..)---|",n,group->ns,group->site[0]);
	  col+=3; } }
      if (mn->anchor & ANCHOR_r) if (anchor.f) fprintf(anchor.f,"-------- position molecule%3d ---------|",n),col+=3;
      if (mn->anchor & ANCHOR_v) if (anchor.f) fprintf(anchor.f,"-------- velocity molecule%3d ---------|",n),col+=3;
      if (mn->anchor & ANCHOR_x) if (anchor.f) fprintf(anchor.f,"------ acceleration molecule%3d -------|",n),col+=3;
      if (mn->anchor & ANCHOR_f) if (anchor.f) fprintf(anchor.f,"---------- force molecule%3d ----------|",n),col+=3;
      if (mn->anchor & ANCHOR_m) if (anchor.f) fprintf(anchor.f,"-------- momentum molecule%3d ---------|",n),col+=3;

      /* error fixed (?) - needs doublecheck! */
      if (mn->anchor & ANCHOR_s) if (anchor.f) fprintf(anchor.f,"---- constr.force mol.%3d site %3d ---|",n,mn->iaxis[0]),col+=3;
      if (mn->anchor & ANCHOR_c) if (anchor.f) fprintf(anchor.f,"- constr.force mol=%3d center-of-mass -|",n),col+=3;
      if (mn->anchor & ANCHOR_axis) if (anchor.f) fprintf(anchor.f,"- constr.moment mol=%3d center-of-mass |",n),col+=3; }

    if (anchor.xyz) if (anchor.f) fprintf(anchor.f,"-   center-of-mass of the whole box    |"),col+=3;

    if (anchor.f) fprintf(anchor.f,"     t\n#");
    loop (n,1,col+2) if (anchor.f) fprintf(anchor.f,"     #%-7d%s",n,n%3?"":"|");
    if (anchor.f) fprintf(anchor.f,"\n"); }

#  define ADVANCEINFO sprintf(anchor.rec[col].info,"A%02d%c",col,"xyz"[k])

  if ((-option('k'))&3) {
    /* V3.6b: support for StaAdd and CP:
       as above AND taking into account only selected coordinates */
    if (!pass) loop (pass,0,2) {
      col=0;
      /* print info and allocate arrays */
      loop (n,0,No.N) {
        mn=molec+n;
        if (mn->anchor & ANCHOR_g) {
          struct anchorgroup_s *group;

          looplist (group,mn->group)
            loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s constrained force mol.%d group ns=%d(%d..)",
                  anchor.rec[col].info,n,group->ns,group->site[0]); }
            col++; } }

        if (mn->anchor & ANCHOR_r)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            prt("ANCHOR: %s position molecule %d",anchor.rec[col].info,n);
            col++; }

        if (mn->anchor & ANCHOR_v)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s velocity molecule %d",anchor.rec[col].info,n); }
            col++; }

        if (mn->anchor & ANCHOR_x)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s acceleration molecule%3d",anchor.rec[col].info,n); }
            col++; }

        if (mn->anchor & ANCHOR_f)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s force molecule %d",anchor.rec[col].info,n); }
            col++; }

        if (mn->anchor & ANCHOR_m)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s momentum molecule %d",anchor.rec[col].info,n); }
            col++; }

        /* error fixed (?) - needs doublecheck! */
        if (mn->anchor & ANCHOR_s)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s constrained force mol.%d site %d",anchor.rec[col].info,n,mn->iaxis[0]); }
            col++; }

        if (mn->anchor & ANCHOR_c)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s constrained force mol=%d wrt. center-of-mass",anchor.rec[col].info,n); }
            col++; }

        if (mn->anchor & ANCHOR_axis)
          loop (k,0,3) if ((1<<k)&mn->xyz) {
            if (pass) {
              ADVANCEINFO;
              prt("ANCHOR: %s constrained angular momentum mol=%d wrt. center-of-mass",anchor.rec[col].info,n); }
            col++; } } /* n */

      if (anchor.xyz)
        loop (k,0,3) if ((1<<k)&anchor.xyz) {
          if (pass) {
            ADVANCEINFO;
            prt("ANCHOR: %s constrained force center-of-mass of the whole box excl. above",anchor.rec[col].info); }
          col++; }

      if (!pass) {
        anchor.col=col;
        allocarrayzero(anchor.rec,col); } } }

  keep[0]=keep[1]=0;
  loop (n,0,No.N) {
    molecule_t *mn=molec+n;

    if (mn->anchor & ANCHOR_pos) keep[0]++;
    if (mn->anchor & ANCHOR_axis) keep[1]++; }

  
  prt("ANCHOR: %d points (sites, molecular centers of mass) anchored",keep[0]);
  prt("ANCHOR: %d axes anchored",keep[1]);
  if (anchor.xyz) prt("ANCHOR: CM of remaining molecules anchored");

  if (keep[0]+keep[1])
    WARNING(("There are anchored points or axes, but I am not smart enough to\n\
*** calculate the degrees of freedom and drift corrections.\n\
*** You must set variables `conserved' and `drift' by yourself!\n\
*** Ignore this message if you have already done so, otherwise\n\
*** DO NOT IGNORE THIS MESSAGE -- KINETIC TEMPERATURE MAY BE WRONG!"))
}
#endif /*# ANCHOR */

#if PARALLEL
void initparallel(void) /************************************** initparallel */
{
  static int pass;

#  if PARALLEL==1
#    define MAXTHREADS 24
  if (No.th<1 || No.th>MAXTHREADS)
    ERROR(("parallel: number of threads (%d) is out of permitted range [1,%d]",No.th,MAXTHREADS))
  prt("PARALLEL=%d (linked-cell list parallelized), %d threads",
      PARALLEL,No.th);

#  endif /*# PARALLEL==1 */

#  if PARALLEL==2
  if (No.th!=2) ERROR(("PARALLEL==2: the number of threads (%d) is not 2",No.th))
#  endif /*# PARALLEL==2 */

#  ifdef SERIAL
  WARNING(("SERIAL: this is debugging version with all %d threads running sequentially",No.th))

#  else /*# SERIAL */

/* Thread initialization. Taken from:
**  GNU Pth - The GNU Portable Threads
**  Copyright (c) 1999-2003 Ralf S. Engelschall <rse@engelschall.com>
*/
  if (pass==0) {
    if (pthread_attr_init(&No.thread_attr))
      ERROR(("parallel: pthread_attr_init"))
    if (pthread_attr_setdetachstate(&No.thread_attr, PTHREAD_CREATE_JOINABLE))
      ERROR(("parallel: pthread_attr_setdetachstate")) }
#  endif /*#!SERIAL */

  rallocarray(No.thread,No.th);
  fprintf(stderr,"initparallel: %d threads allocated\n",No.th);
  prt("initparallel: %d threads allocated",No.th);

  pass=1;

  if (option('t')) {
#  if PARALLEL==2
    partimes.nrspace=2;
#  else /*# PARALLEL==2 */
    partimes.nrspace=No.th;
#  endif /*#!PARALLEL==2 */
    allocarrayzero(partimes.rspace,partimes.nrspace);
    partimes.nkspace=No.th;
    allocarrayzero(partimes.kspace,partimes.nkspace); }
}

#endif /*# PARALLEL */
