#define VERSION "3.7b"

#if defined(LINKCELL) && defined(FREEBC)
#  error "LINKCELL not supported for FREEBC"
#endif /*# defined(LINKCELL) && defined(FREEBC) */

#include "ground.h"
#include "sds.h" // BJERRUM only
#include <time.h>

#include "options.h"
#include "optlist.c"

#include <signal.h>

#include "varfile.h"
#include "statics.h"
#include "linregr.h"

#include "simglob.h"

#include "asksig.h"

#include "units.h"
#include "simgear.h"
#include "simcg.h"

#include "interpot.h"
#include "simdef.h"
#include "forces.h"

#ifdef COULOMB
#  include "ewald.h"
#  include "elst.h"
#endif /*# COULOMB */

#include "simpot.h"
#include "siminit.h"
#include "simmeas.h"
#include "simmeasx.h"

#include "setss.h" /* widomrdf(), poteps */

#include "norm.h"
#ifdef LINKCELL
#  include "linklist.h"
#endif /*# LINKCELL */

#include "simils.h"
#include "rhs.h"
#include "constrd.h"
#include "cputime.h"

#ifdef MARK
// not tested recently
#  include "mark.h"
#endif /*# MARK */

#include "maindef.c"

static char *email=NULL,*host,*arg0;
static FILE *lock;
static char *lockname;
static void removelock(void) /* --------------------------------- removelock */
/*
  removes the lock file SIMNAME.loc
  registered by atexit() => to be called at exit()
*/
{
  if (lock) fclose(lock);
  if (lockname) remove(lockname);

  //  if (sig!=4 && sig!=2 && sig!=15)
  /* do not mail if stopped by .stp, kill-2, kill-15 */
  if (email) {
    if (system(string("echo \"%s %s %s on %s finished\" | mail -s \"%s_finished\" %s",
                      arg0,simils.sysname,simils.simname,host,simils.simname,email)))
    WARNING(("email command failed")) }
}

static int justnow(double dt,double half) /* ----------------------- justnow */
{
  dt=fabs(dt);
  if (dt==0) return 0;

  return (int)((t+half)/dt+1e4) > (int)((t-half)/dt+1e4);
}

static void prtsfill(char *s) /* ---------------------------------- prtsfill */
{
  char *c;
  int n=0;

  for (c=s; *c; c++) n+=*c=='\n';

  while (n++<22) _n

  prts(s);
}

static void removedriftssta(int pr) /* --------------------- removedriftssta */
{
  removedrifts(pr);
  StaSet(0,2,2,0);
  StaAdd("CM drift to remove",En.CMdrift);
  StaAdd("v drift to remove",En.vdrift);
  StaAdd("AM drift to remove",En.AMdrift);
  StaAdd("omega drift to remove",En.omegadrift);
}

int main(int narg,char **arg) /**************************************** main */
{
  int
    noint=10,     /* number of steps in 1 cycle = frequency of measurements */
    no=10000,     /* # of cycles in 1 main cycle */
    nomax=0,      /* max # of cycles from start of measurements (init>=1) */
    initno=0,     /* option -NUMBER (overrides no in the 1st cycle) */
    key=0,        /* "command" */
    sort=0,       /* sort by CM, coordinate: 1=x, etc. (deprecated) */
    checkfrac=1,  /* check fractional charges once, turned off by -q */
    pass=0,       /* pass of reading data */
    CPnbit=0,     /* resolution (1/2^CPnbit) for packed convergence profiles */
    rhoreached,   /* switch to check whether rho was reached */
    *group,       /* [nspec] groups of species */
    initvel=0,    /* # of molecules to initialize velocities */
    issta=0,      /* statistics (of "Tkin") exists */
    expectedbuffer=0, /* to check remaining disk space, see -w */
    plbopened=0,
    unsplit=0,    /* usplit periodically split molecules, in x=1,y=2,z=4 */
    isPkincorr,
    iarg,nsteps,i;

  volatile int
    icyc=0,iint=0;

  double
    stop=0,       /* max time from simulation start (init>=2) */
    Tstop=0,      /* stop if Tstop reached */
    Emax=1e4,     /* max energy of inserted molecule - see initcfg */
    drmax=0.2,    /* max displacement in one integration step (init>2 only)
                     drmax<0: applies for any init */
    bulkmodulus=0,/* for NPT */
    P=101325,     /* barostat pressure */
    maxscale=1.03,/* max scaling of box, applies with tau.P, tau.sig, (tau.rho) */
    pins=0.01,    /* limit probability of insertions - see initcfg */
    initrho=0,    /* initial density - see initcfg; 0=rho */
    dV=0,         /* volume-change and surface-change, cf. constrd.dV */
    DT,           /* cycle in s */
    oldt,newt,    /* to check changed t */
    lastspeed=0,  /* to predict next job time */
    LJcutoff=-3,  /* Lennard-Jones cutoff - see comments in simdef.c! */
    err;

  struct {
    double
      grid,       /* intervals/1A; negative:s-s rdf's merged; 0:none */
      cutoff;     /* range of rdf */
    int
      onefour;    /* 1 if 1-4 are to be included to RDF's, 0 otherwise */
  } rdf = {0,0,1};

  struct {
    double prt; /* frequency of printing t,T,E,... to out [in ps] (default=every) */
    double cfg; /* save configuration frequency [in ps] - see savecfg */
    double plb; /* frequency of saving playback file */
  } dt = { 1e-4,0,0 };

  struct {
    /* with -m0 or -m1: */
    int from;  /* first frame to read from plb or # in SIMNAME.# */
    int to;    /* last frame to read (incl.): must be set */
    int by;    /* stride of frames */
    int frame; /* private: current frame */
  } reread={1,0,1,0};

  int nshift=0x7ffffffe; /* shift mode, see addshift() */
  vector shift={0,0,0}; /* added to cfg. before run, then :=0 */
  vector vshift={0,0,0}; /* added to cfg. before run, then :=0 */

  struct removemol_s {
    int n; }
  removemol = {-1};

#ifdef WIDOM
  struct {
    int spreal; /* for scaled particle (id change): the real species */
    int sp; /* species to insert or be changed to virtually */
    int n; /* # of attempted insertions per cycle [SLIT,SLAB]; negative: excl. wall-mol pot. */
#  ifdef SLAB
    int mode; /* also symmetrized profile */
    double z0,z1,dz; /* define slabs */
#  endif /*# SLAB */
    double corr; /* cutoff correction without division by V */
  } widom={-1,0,0,2};
#endif /*# WIDOM */

  /* #ifdef XSECTION: xs used to be here, moved to simglob */
#ifdef RGYR
  struct {
    int end[2];   /* indices of `ends' */
    int cp;       /* -1: averaged Rgyr and end-to-end in CP, >=0 : 1 mol.val. */
  } rg={{0,-1},0};
#endif   /*# RGYR */

#ifdef DIHHIST
  struct {
    double res;  /* resolution, in deg */
    int grid;   /* # of intervals / 360 deg */
    int cp,dcp;
  } dih={0,0,-1,0};
#endif /*# DIHHIST */

#ifdef BJERRUM
  struct {
    int mode;   /* sum: 1=topology(?), 2=Gauss */
    double from;/* Gaussian width -- start */
    double to;  /* Gaussian width -- the last one, in cbrt(V) */
    double q;   /* Gaussian width quotient for loop */
    double eps; /* accuracy of minimization */
    ToIntPtr cfg; /* to average (bj.mode&4) */
  } bj={0, 4.0, 0.25, 1.09050773266525766, 1e-6};
#endif /*# BJERRUM */

  vector Pscale={1,1,1};
  double sigvdW=0; /* particle size to adjust: see tau.sig, tau.i, tau.j */
  double lastEtot,halftcyc;

  const char *mark; /* temporary; mark for release */
  double starttime,stoptime,stop_start; /* see cputime.h, #define CHEAPTIME */
  double rho=0; /* density [kg/m3] */
  double refrho;
  vector L={0,0,0}; /* box size (in input data) */
  double cutoff=0; /* cutoff (in input data) */
  /* int corr=3+16; made global in V3.6v */
  FILE *indef=NULL; /* to reread the def-file */

  /* keys for ID="mnemonics" on data and help moved to mainhlp.c */
#define INFO "*** " PROJECT " *** V" VERSION " *** (c) J.Kolafa 1991-now ***"
  //out=stderr;
  out=stdout;

#include "mainhlp.c"

#ifdef PERSUM
  option('x')=0; /* the default is no optimized waters because they are
                    wrong with PERSUM */
  corr|=32;
#endif  /*# PERSUM */


  /*** option analysis ***/
  /* the scaling factor and defaults first */
  loop (iarg,1,narg)
    if (arg[iarg][0]=='-' && arg[iarg][1]=='_') option('_')=atoi(arg[iarg]+2);
  option('j')=option('a')=option('q')=option('_');
  optionscaling=option('_');
  if (option('_')!=100) fprintf(stderr,"WARNING: denominator for options -j -q -^ -a is %d\n",option('_'));

  arg0=arg[0];
  loop (iarg,1,narg)
    if (arg[iarg][0]=='-')
      if (arg[iarg][1]>='0' && arg[iarg][1]<='9')
        initno=atoi(arg[iarg]+1);
      else {
        if (arg[iarg][1]=='@') email=arg[iarg]+2;
        getoption(arg[iarg])
        if (arg[iarg][1]=='q') checkfrac=0; }
#if CHECKHEAP==-2
    /* to debug memory problems: option +AllocRange */
    else if (arg[iarg][0]=='+') {
      AllocTrace++;
      AllocRange=atoi(arg[iarg]+1); }
#endif /*# CHECKHEAP==-2 */
    else { /* no -OPTION */
      char *c;

      for (c=arg[iarg]; *c; c++)
        if ((!isalnum(*c) && !strchr("%~_-=+,.@#/",*c)) || (unsigned)*c&128)
          prt("WARNING: character \'%c\' = %d dec in file name is not recommended",*c,*c);

      if (!simils.sysname)
        simils.sysname=arg[iarg];
      else if (!simils.simname)
        simils.simname=arg[iarg];
      else if (!simils.plbname)
        simils.plbname=simils.cfgname=arg[iarg];
      else
        ERROR(("max 3 file/name arguments allowed (SYSNAME, SIMNAME, PLBNAME), more found")) }

#if PARALLEL
  if (option('i')==-9) option('i')=1+(!option('s'));
  if (option('i')==0) WARNING(("option -i0 for PARALLEL"))
#  if PARALLEL==1 || PARALLEL==3
  if (sizeof(struct pll_global_s)%CACHELINE) DISASTER(("sizeof(pll_global_s)=%d is not an integer multiple of the cache line=%d\n\
*** adjust padding in `struct pll_global_s' in sim/interpot.h and recompile!",sizeof(struct pll_global_s),CACHELINE))
#  endif /*# PARALLEL==1 || PARALLEL==3 */
#else /*# PARALLEL */
  if (option('i')==-9) option('i')=2*!option('s');
#endif /*#!PARALLEL */

  /* install SIGINT (2) and SIGTERM (15) signal handler */
  if (option('i')>=0) {
    signal(SIGINT,asksig);
    /* SIGTERM not for MPI ...? */
    signal(SIGTERM,asksig); }

  if (!simils.sysname) ERROR(("no file/name argument (SYSNAME) found"))

  if ( (mark=getext(simils.sysname)) && !strcmp(mark,".ble") )
    simils.sysname=stripext(simils.sysname);

  if (simils.simname) {
    if ( (mark=getext(simils.simname)) && !strcmp(mark,".get") )
      simils.simname=stripext(simils.simname); }
  else {
    prt("WARNING: SIMNAME not specified, SIMNAME=SYSNAME will be used.");
    simils.simname=simils.sysname; }

  if (simils.plbname) {
    if ( (mark=getext(simils.plbname)) && !strcmp(mark,".plb"))
      simils.plbname=stripext(simils.plbname); }
  else {
    prt("WARNING: PLBNAME not specified, PLBNAME=SIMNAME will be used if needed.");
    simils.cfgname=simils.plbname=simils.simname; }

  if (option('m')<2 && !strcmp(simils.plbname,simils.simname))
    ERROR(("Third name/file parameter PLBNAME missing or is identical to SIMNAME\n\
*** which is not allowed in the playback mode (options -m0,-m1)"))

  lockname=strdup(Fn("loc"));
  lock=fopen(lockname,"r+");

  out=stdout;
  if (lock) {
    if (option('o')) ERROR(("Lock file \"%s\" exists\n\
(\"cook %s\" may be running or has not been stopped properly)",
      lastFn,simils.simname))
    else
      fprintf(stderr,"WARNING: Lock file exists but is ignored because of -o0\n"); }
  else
    lock=fopen(Fn("loc"),"w");

  atexit(removelock);

  if (option('s')==1) option('s')=256;
  initscroll(option('s')*1024);

  if (option('s')) {
    /* interactive mode (deprecated) */
    lastFn="STDOUT";
    out=stdout; }
  else {
    /* batch mode */
    backup("prt");
    out=fopen(Fn("prt"),"wt"); }

  if (!out) ERROR(("cannot write to %s",lastFn))

  prts(INFO);

  prts_("***");

  printdefines(out); /* out of the (deprecated) scroll system */

  { /* options, sex, environment */
    int Int=1  ;
    char *end=(char*)&Int;
    char *user=getenv("USERNAME");

    if (!user) user=getenv("USER");
    if (!user) user="UNKNOWN";

    host=getenv("HOSTNAME");
    if (!host) host=getenv("HOST");
    if (!host) host="UNKNOWN";

    prt("*** %s endian short=%d int=%d long=%d pointer=%d ***",
        *end==1?"little":*end==0?"BIG":"?",
        sizeof(short),sizeof(int),sizeof(long),sizeof(void*));

    loop (i,0,narg) {
      prt_("%s%c",arg[i]," \n"[i==narg-1]);
      fprintf(stderr,"%s%c",arg[i]," \n"[i==narg-1]); }
    prt("HOST=%s USER=%s PWD=%s",host,user,getenv("PWD"));
    fprintf(stderr,"HOST=%s USER=%s PWD=%s\n",host,user,getenv("PWD"));
    prt("Edited arg list:");

    if (simils.sysname) prt_("%s ",simils.sysname);
    if (simils.simname) prt_("%s",simils.simname);
    if (simils.simname || simils.sysname) _n

    loop (i,0,32) {
      if (i==17) _n
      prtoption(i,optionlist) }
    _n

    /* option('m'): 0=reread cfg, 1=reread plb, 2=Verlet+SHAKE, 3=Gear order */
    gear.order=option('m');
    if (gear.order==1) prt("playback mode (-m1): %s.plb will be read instead of simulating",simils.plbname);
    if (gear.order==0) prt("playback mode (-m0): %s.1,%s.2,.. will be read instead of simulating",simils.plbname,simils.plbname);
    if (gear.order>2 && option('c')&1) epsc=0.05;

    if (option('p')<0) ERROR(("negative option -p no longer supported"))
    No.pred=option('p')%10;
    if (No.pred==9) No.pred=gear.order;
  }

  /* remove the stp-file so that the simulation is not immediately stopped */
  remove(Fn("stp"));
  _n

#if PARALLEL
  if (option('i')==0) WARNING(("option -i0 for PARALLEL"))
  {
    char *env=getenv("NSLOTS");

    if (!env) ERROR(("PARALLEL: undefined number of threads\n\
*** (no environment variable NSLOTS)"))
    No.th=atoi(env);
  }
#endif /*# PARALLEL */

  if (option('n')&1) {
    if (option('n')&2 && gear.order!=2)
      WARNING(("SIMNAME.plb shifted by -h/2 but SIMNAME.vlb not, check options -n -m"))
    if (((option('n')&2)==0) && gear.order==2)
      WARNING(("SIMNAME.vlb shifted by -h/2 but SIMNAME.plb not, check options -n -m"))
    if (option('n')&4) prt("WARNING: centering SIMNAME.plb is an experimental feature"); }

  rndinit(7L/*low quality enough*/,option('z'));

  /* open sysname.def for pass=0 */
  if ( !(in=fopen(Fn("def"),"rt")) ) {
    if (strcmp(simils.simname,simils.sysname)) {
      /* simils.simname.def not found: trying sysname.def */
      char *swapname=simils.simname;

      simils.simname=simils.sysname; /* swap sysname/simils.simname temporarily */
      in=fopen(Fn("def"),"rt");
      simils.simname=swapname; /* normal sysname/simils.simname again */
      if (!in) ERROR(("no %s.def nor %s.def found",simils.sysname,simils.simname)) }
    else
      ERROR(("no %s.def",simils.sysname,simils.simname)) }

  _Id.ignoreid=1; /* NB: no check of data! */
  getdata
    /* reading def-file -- 1st pass */
    get(equalize.mol) get(equalize.cfg) get(equalize.sp)
    get(no) get(noint) get(h)
    get(E) get(rho) getvec(L,,DIM)
    get(tau.T) get(T)
#ifdef ECC
    get(el.ecc) get(el.epsf)
#endif /*# ECC */
  enddata
  box.rho=rho; /* shadowed */
  _Id.ignoreid=0; /* check of data again */

  rewind(in);
  indef=in;

  /* read ble-file (swap sysname/simname temporarily) */
  mark=simils.simname;
  simils.simname=simils.sysname;
  in=fopen(Fn("ble"),"rt");
  simils.simname=mark; /* normal sysname/simname again */
  if (!in) ERROR(("open system ble-file %s",simils.sysname))
  readblend();
  fclose(in);

  in=indef; /* for pass=0, rewound above */

  /* <<<<<<<<<<<<<< loop over "jobs" = ';' in input data <<<<<<<<<<<<<< */
  while (sig<2) {
    sig= -1;

    checkranges(1);

/* get data (from def- and get-files) */
#include "mainget.c"

    /* CAVEAT: a check missing for changing groups for init=0,1 */
    ngroups=0;
    loop (i,0,nspec) {
      prt("species %d is in group %d",i,group[i]);
      if (group[i]<0 || group[i]>=nspec) ERROR(("group[%d] out of range",i))
      spec[i]->group=group[i];
      Max(ngroups,group[i]) }
    ngroups++;
    prt("%d groups found",ngroups);

    lag.ierr=noint*lag.err; // by DT/noint, needed only once
    lag.in=lag.n+log(noint)/log(2)+0.5; // needed only once
    No.bulkmodulus=bulkmodulus/Punit;
    No.P=P/Punit;

    if (h==0) WARNING(("The timestep h=0, I hope you know what you are doing.\n\
*** Ekin etc. will be set to zero, some test suppressed."))

    // moved to be compatible with load.N
    init_append=0;

    if (tau.E!=0 && thermostat) {
      ERROR(("both tau.E and thermostat set: only one allowed"))
      thermostat=0;
      prt("thermostat=\"none\" set"); }

    if (tau.E==0 && !thermostat && gear.order==2 && epsc<1 && epsc>3e-8)
      WARNING(("adiabatic simulation and epsc=%g may not be sufficiently small",epsc))

#ifdef SLAB
    slab.outx=slab.out-slab.outc;
    wall.is=wall.g||wall.n;
#endif /*# SLAB */
    setdrift();

    if (T) {
      double maxh=0.03*sqrt(No.minmass/T);

      prt("\
minimum atom mass (after equalization, if any) = %g p.u. = %g g/mol\n\
suggested timestep h at T=%g K in (%.2g,%.2g) ps\n\
(based on thermal motion only, shorter h needed in case of large forces)",
        No.minmass, No.minmass*Munit, T, maxh*0.5,maxh);

      if (thermostat) {
        if (h>maxh*1.5)
          WARNING(("the timestep h=%g is probably too long\n\
*** bad energy conservation or crash may occur",h))
        else if (h>maxh) prt("\
WARNING: the timestep h=%g is rather long, the results may be inaccurate",h);
      } }

    if (rho==0 && PROD(L)==0)
      if (init>=3) ERROR(("undefined box: set box (L[0],L[1],L[2]) or density (rho)"))

    No.ncoord="\0\1\1\2\1\2\2\3"[rescale&RESCALE_XYZ];

    underline("rescaling of the simulation box");
    prt("applies if tau.P (barostat), tau.rho, or dV (virtual volume change)");
    prt("rescale = %d = 0x%x = %s",rescale,rescale,int2sumbin(rescale));
#if (PRESSURETENSOR&PT_ANY) != PT_ANY
    if (rescale&RESCALE_PT)
      ERROR(("\
Pressure tensor components are requested but this version of cook\n\
*** was not compiled so. Use at least PRESSURETENSOR=3."))
#endif /*# (PRESSURETENSOR&PT_ANY) != PT_ANY */
    if (rescale&RESCALE_XYZ) {
      prt("configuration rescaling in these %d coordinate(s): %s\n\
based on %s",
        No.ncoord,
        prtxyz(rescale&RESCALE_XYZ),
        rescale&RESCALE_PT?"pressure tensor components":"isotropic pressure"); }
    else {
      prt("no configuration rescaling allowed, dV=tau.P=tau.rho=0 set\n\
*** (check rescale if this is not what you want)");
      dV=tau.P=tau.rho=0; }

    prt("configuration rescaling will be based on %s",
      rescale&RESCALE_CM ? "molecular centers-of-mass" : "sites");
    if (dV<0) prt("WARNING: dV<0 no longer denotes molecular-based rescaling");

    prt("box.center=%g %g %g\n\
This is used for center of mass (CM) and angular momentum (AM) mesurements:\n\
lag.CM=%d  lag.LM=%d  lag.AM=%d\n\
NB: slab.sp (in SLAB version) sets box.center automatically",
        VARG(box.center),lag.CM,lag.LM,lag.AM);

    underline("pressure");

    if (virial>=1 && virial<=3)
      prt_("User-selected pressure formula:");
    else {
      prt_("This version default pressure formula:");
#if (PRESSURETENSOR&PT_ANY) == PT_ANY
      virial=3;
#else
      virial=1;
#endif
    }
    prt(" virial=%d = %s",virial,
        virial==1 ?   "Pevir = virial theorem (elst_virial = -elst_energy)"
        : virial==2 ? "PdV   = virtual volume (shape) change method"
        : virial==3 ? "Ptr   = trace of pressure tensor (from forces)" : "ERROR");
    prt("\
Will be reported as  Pcfg [Pa]  and shown in convergence profile\n\
(unless superseded by the virtual volume change method, see dV).\n\
Pressures available, cf. variable virial:\n\
  Pevir = uses the virial theorem (elst_virial = -elst_energy):\n\
          - good with Ewald point charges (cookew*)\n\
          - wrong for Gaussian charges (cookgc*)\n\
          - inaccurate for cutoff electrostatics (cookce*)\n\
  PdV = virtual volume change, good for debugging, see variables dV and rescale\n\
  Ptr = trace of pressure tensor (from forces):\n\
          - most accurate\n\
          - for Ewald point charges is a bit faster (if the full pressure tensor\n\
            calculation is turned off by #define PRESSURETENSOR=0 in simopt.h");
    if (virial==2) {
      if (dV==0) ERROR(("Virtual volume change method (virial=2) selected but dV=0"))
      if (thermostat>=T_NPT) ERROR(("thermostat=%d cannot use the virtual volume change method"))
      if (tau.P<0)
        WARNING(("Once-by-cycle Berendsen barostat (tau.P<0)\n\
*** has not been tested with the virtual volume change method."))
      if (tau.P>0)
        prt("WARNING: Every step Berendsen barostat will use pressure by virtual volume change"); }

    constrd.dV=dV;
    constrd.PdVname = rescale&RESCALE_CM ? "PdVmol [Pa]" : "PdVatom [Pa]";

#ifdef SLAB
#  if SLAB & 1
#    if (PRESSURETENSOR&PT_ANY) != PT_ANY
#      error "SLAB & 1 requires PRESSURETENSOR at least 3"
#    endif /*# (PRESSURETENSOR&PT_ANY) != PT_ANY */
    if (tau.L) {
      underline("MODIFICATION - extended slab control");
      prt("Target box sizes maintained with correlation time tau.L=%g:",tau.L);
      prt("Lx=%g + %g*T + %g*T^2 + %g*T^3",slab.Lx[0],slab.Lx[1],slab.Lx[2],slab.Lx[3]);
      prt("Ly=%g + %g*T + %g*T^2 + %g*T^3",slab.Ly[0],slab.Ly[1],slab.Ly[2],slab.Ly[3]);
      prt("Enthalpy (energy) E=%.12g K maintained with correlation time tau.E=%g",E,tau.E);
      if (!tau.E)
        WARNING(("tau.E=0 (no energy/enthalpy correction) is not recommended"))
      if ( (rescale & (RESCALE_XYZ|RESCALE_PT) ) != (RESCALE_XYZ|RESCALE_PT) )
        WARNING(("nonstandard rescale=%d (consider \"XYZCM\" or \"XYZ\")",rescale)) }
#  endif /*# SLAB & 1 */
#endif /*# SLAB */

    if (thermostat>=T_ANDERSEN && thermostat<=T_MAXWELL_CM) {
#if VERLET==1 || VERLET==2
      WARNING(("kinetic energy inacccurate with Andersen/Maxwell thermostat\n\
*** recompile with VERLET=3 (recommended) or 0 (less accurate)"))
#endif /*# VERLET==1 || VERLET==2 */
      if (gear.order>2)
        WARNING(("Gear integration inacccurate with Andersen/Maxwell thermostat")) }

    if (gear.order>2) {
      if (thermostat>=T_NPT) ERROR(("NPT via Nose/MTK not implemented with Gear"))
      if (tau.rho<0) ERROR(("Box size as a function of time not implemented with Gear.")) }

    if (tau.T && !thermostat)
      prt("WARNING no thermostat: nonzero tau.T=%g ignored\n\
(if you mean friction thermostat as in old cook, use thermostat=\"friction\")",tau.T);
    if (tau.T==0 && thermostat>0 && thermostat<T_NPT) {
      ERROR(("thermostat=%d was specified but tau.T=0. There is no default anymore!\n\
*** Hints for tau.T:\n\
***   0.1--0.3  Nose, NPT\n\
***   0.05--0.2 Berendsen/frictions start\n\
***   0.5--2    Berendsen/frictions productive runs\n\
***   0.2--2    Andersen,Maxwell\n\
***   >1        Langevin\n",thermostat))
      thermostat=0; }

    if (tau.P==0 && thermostat>=T_NPT) {
      ERROR(("thermostat=NPT (i.e., thermostat+barostat) specified but tau.P=0"))
      thermostat=0; }

    /* total charge check delayed to be interactive */
    if (checkfrac) {
      int sp, nfrac=0;
      double charge;

      checkfrac=0;
      loop (sp,0,nspec) {
        charge=spec[sp]->charge/electron;
        if (fabs(fmod(fabs(charge+0.5),1)-0.5) > 1e-5) {
          nfrac++;
          prt("charge =%10.6f e of species %d is fractional",charge,sp); } }
      if (nfrac)
        WARNING(("\
%d molecule(s) with fractional charges\n\
(check mol-file, check blend option -_PREC)",nfrac)) }

    /* reread.to is included */
    if (gear.order<2) {
      if (reread.to<1 ) ERROR(("reread.to=%d unset or bad with the playback mode (-m1 or -m0)",reread.to))
      no=(reread.to-reread.from)/reread.by+1;
      prt("playback mode: reread.from=%d reread.to=%d reread.by=%d => no=%d",reread.from,reread.to,reread.by,no); }

    if (initno) {
      prt("number of cycles no has been overridden by option -%d",initno);
      no=initno;
      initno=0; }

    underline("units");
    prt("\
Selected program units (p.u.):\n\
  length = 1 [AA] = 1e-10 [m]\n\
  time = 1 [ps] = 1e-12 [s]\n\
  mass = %.11g [kg] = %.11g [g/mol]\n\
  energy = 1 [K] (more precisely, (1 K)*k_B)\n\
         = %.11g [J] \"=\" %.11g [J/mol] = %.11g [cal/mol]\n\
           (where mol = 1 mol of simulation boxes)\n\
  pressure = %.11g [Pa]\n\
  charge = %.11g [C] (CGS-style energy = q^2/r)\n\
  dipole moment = %.11g [D] = %.11g [C m]\n\
                = %.11g e AA\n\
  el.field intensity = %.11g [V/m]\n\
A quantity printed without units is in p.u. or dimensionless,\n\
see MACSIMUS/sim/units.h for details.\n\
Units are in [] for compatibility with the MACSIMUS calculator \'evu\';\n\
also, [+UNIT] converts to p.u. and [-UNIT] backwards, see evu -> p(.)units.",
        massunit,Munit,
        energyunit,Eunit,Eunit/kcal*1e3,
        Punit,
        chargeunit,
        Debye,chargeunit*lengthunit,1./electron,
        massunit*lengthunit/Sqr(timeunit)/chargeunit);
    _n

    sig=0;

#ifdef DIHHIST
    if (dih.cp<-1) dih.cp=-1;
#endif /*# DIHHIST */

#ifdef COULOMB
#  if COULOMB<0
    if (el.kappa<=0 && gear.order>0) WARNING(("el.kappa=%g is suspicious for Ewald",el.kappa))
    if (el.corr) {
      if (el.corr==8) WARNING(("el.corr=8 requests partial Yeh-Berkowitz correction, but no coordinate given"))
      if (fabs(el.epsinf)<1e15) ERROR(("el.corr=%d: the Yeh-Berkowitz correction requires el.epsinf=infinity",el.corr))
      prt("el.corr=%d: %s Yeh-Berkowitz correction in these coordinates: %s",
          el.corr,
          el.corr&8?"partial":"full",
          prtxyz(el.corr));
      if (el.corr&8) prt("\
- only pressure tensor components and energy corrected\n\
- trajectory without correction => energy conservation violated"); }
# else
    if (Eext.isE && Eext.f) WARNING(("oscillating el. field: many features not available (Ewald needed)"))
#  endif /*#!COULOMB<=-1 */
# else
  if (Eext.isE && Eext.f) WARNING(("oscillating el. field: many features not available (Ewald needed)",Eext.f))
#endif /*# COULOMB */

#ifdef POLAR
    if ((option('p')/10)%10==0) if (tau.dip==0)
      ERROR(("Car-Parrinello-like method and tau.dip=0"))
    if (abs(option('^'))>999) option('^')=abs(option('^'))%1000;
    if (scf.omega==-9) { /* scf.omega not specified in data */
      if (option('^'))
        scf.omega=option('^')/100.;
      else
        /* no -^ nor scf.omega: defaults
           - for -m2 (ASPC) the default is postponed to pred2.c
           - for -m>2 (full-iteration) the default is scf.omega=0.9
        */
        if (gear.order!=2) {
          WARNING(("undefined scf.omega set to %g",scf.omega=0.9)) } }
    else {
      /* scf.omega specified in data */
      option('^')=100*scf.omega+0.5; /* not to be set to default in pred2.c */ }
    if (thermostat>=T_ANDERSEN && thermostat<=T_MAXWELL_CM) WARNING(("Maxwell-Boltzmann-based thermostat and POLAR:\n\
*** check the parameters (ASPC inaccurate, iterations inefficient)"))
#endif /*# POLAR */

    putkey(init,initkey)
    if (init>=3 && init<10) {
      putkey_(pins,pinskey) put2(initrho,Emax) }
    put2(el.Perr,el.epsinf)
    put2(rho,tau.rho)
    put2(E,tau.E)
    put_(T) put_(tau.T) putkey(thermostat,thermostatkey)
    if (thermostat>=T_TR && thermostat<=T_FRICTIONS) put(T_tr_in)
    put2(P,tau.P)
    if ((tau.P!=0 || tau.rho!=0) && tau.sig!=0)
      ERROR(("intermolecular parameter adjustement requires NVT"))
    if (rdf.grid) put2(rdf.grid,rdf.cutoff)

#ifdef DIHHIST
    if (dih.res!=0) {
      if (dih.grid)
        WARNING(("max one of dih.grid, dih.res may be specified (dih.grid will be ignored)"))
      dih.grid=(int)(360./dih.res+0.5);
      if (fabs(360./dih.res/dih.grid-1)>1e-6)
        WARNING(("dih.res rounded to %f",360./dih.grid));
      dih.res=0; }
#endif /*# DIHHIST */

/* rdf is parallelized with PARALLEL==1, dih is serial (intramolecular) */
#if PARALLEL==3
    if (rdf.grid) {
    prt("WARNING: rdf not parallelized, all measurements serial (check noint)");
    No.measureserial++; }
#  ifdef DIHHIST
    if (dih.grid) {
      WARNING(("WARNING: dih not parallelized, all measurements serial (check noint)"))
      No.measureserial++; }
#  endif /*# DIHHIST */
#endif /*# PARALLEL==3 */

    if (tau.P && tau.rho) {
      ERROR(("both tau.rho and tau.P set (nonzero): only one allowed"))
      tau.P=0;
      prt("tau.P=0 set"); }

    if (dV && tau.rho)
      WARNING(("Both dV and tau.rho set. Column 4 of the cp-file will be density."))

    /* put(omega) */
    if (!(option('c')&1)) put_(omegac)
    put2(epsc,eps)
    put3(no,noint,h)
    No.dt=noint*h;
    DT=No.dt*timeunit;
    prt("cycle (of noint MD steps): dt=%g [s]%s",DT,gear.order<2?" (void because in the re-read mode)":"");
    put3(dt.prt,dt.cfg,dt.plb)
    if (dt.cfg!=0 && gear.order==0) ERROR(("Trying to save configurations SIMNAME.# in the re-read mode (-m0)\n\
*** which reads the same configurations."))
    put2(lag.err,lag.n)
#if PRESSURETENSOR&PT_OFF
    put_(lag.Pt)
#endif /*# PRESSURETENSOR&PT_OFF */
#ifdef RGYR
    put_(lag.v)
#endif /*# RGYR */
    put(lag.J)
    StaSet(DT,lag.err,2,lag.n);
    put3(cutoff,rdf.cutoff,LJcutoff)
    if (LJcutoff==0) WARNING(("LJcutoff=0 (no Lennard-Jones interaction)"))

#ifdef RGYR
    put3(rg.cp,rg.end[0],rg.end[1])
#endif /*# RGYR */
#ifdef XSECTION
    prt("\
xs.grid=%g  xs.RvdW=%.8f  xs.Rscale=%.8f\n\
xs.mode=%d  xs.freq=%d    xs.sizelimit=%d",
xs.grid,    xs.RvdW,      xs.Rscale,
xs.mode,    xs.freq,   xs.sizelimit);
    if (xs.sizelimit>AllocSizeLim) {
ERROR(("xs.sizelimit=%d exceeds AllocSizeLim=%d\n\
*** (if you are sure, increase also AllocSizeLim)",xs.sizelimit,AllocSizeLim))
      AllocSizeLim=xs.sizelimit; }
#endif /*# XSECTION */
    put2(corr,cache)

    /*** repeatable initializations ***/
    ralloc(mark,1);

    underline("box size setting and determination of cutoff");
    refrho=setL(L,rho);
    VV(box.L,=L)

    prt("reference box = %.8f %.8f %.8f [AA]",VARG(L));
    prt("reference rho = %.4f [kg/m3]",refrho);
    prt("> The reference box and density are derived from input rho and L[].\n\
> They are used for setting cutoff and Ewald parameters el.alpha, el.kappa\n\
> (if these are not given explicitly).\n\
> The actual box may differ:\n\
> * in NPT simulation (see tau.P and rescale)\n\
> * if other box is loaded and change is not allowed (see tau.rho and rescale)");

    /* load configuration if exists - 1st pass */
    if ((init>=0 && init<3) || gear.order<2) {
      loadcfg(-1,&sigvdW); /* 1st pass */
      if (load.N==4) {
        int err=0;

        if (load.N&3) ERROR(("cannot combine load.N=4 with 1 or 2"));

        loop (i,0,nspec)
          if (spec[i]->N!=simils.spec[i].N) {
            prt("species %d: %d molecules loaded, %d specified",i,simils.spec[i].N,spec[i]->N);
            spec[i]->N=simils.spec[i].N;
            err++; }

        if (err) {
          WARNING(("\
The loaded numbers of molecules differ from the specified ones\n\
*** (in def-file), the latter will be ignored because of load.N=4.\n\
*** The re-initialization is not compatible with equalize.cfg."))
          initNo();
          if (No.N!=simils.N || No.s!=simils.Ns) ERROR(("internal")) } } }
    else if (init<0)
      readplayback(-init,0); /* L and NS only */
    else if (init>10)
      readasc(0,0); /* L only */
    else {
      /* init=3,4,5[,or more] */
      if (initrho && fabs(initrho/refrho-1)>5e-16) {
        prt("Initial density (initrho) differs from the reference density.");
        setL(box.L,initrho); }
      else
        prt("Initial density is given by rho/L specified"); }

    /* now: box.L=actual box loaded or to be initialized,
            L=reference or target box */

    VV(simils.changeL,=load.L) /* ugly control for historical reasons */
    replaceL(L);
    /* VO(load.L,=0) after initreplicate */
    rhoreached=0;

    if (tau.sat!=0 && No.ion!=0)
      WARNING(("Saturation of polarization (tau.sat=%g) and %d free ions:\n\
*** dielectric constant calculation turned on which does not make sense\n\
*** if the ions can freely move (as in liquid), may be OK in a crystal.",
               tau.sat,No.ion))

    err=sqrt(Sqr(log(box.L[0]/L[0]))+Sqr(log(box.L[1]/L[1]))+Sqr(log(box.L[2]/L[2])));
    if (err>1e-15) {
      prt("\
WARNING: The actual (loaded/initial) box size differs from the reference size\n\
         box.L=[%.12g %.12g %.12g]\n\
         V=%.12g  rho=%.12g [kg/m3]\n\
         (mean relative error = %g)",
          VARG(box.L),PROD(box.L),rhounit*No.mass/PROD(box.L),
          err);
      if (err>1e-7) {
        if (tau.P) {
          if (box.L[0]/L[0]-1<0 || box.L[1]/L[1]-1<0 || box.L[2]/L[2]-1<0)
            prt("WARNING: tau.P specified and the reference box is larger than actual\n\
         (check accuracy of Ewald summation)"); }
        else if (tau.rho==0 || (rescale&7)==0) {
          prt("WARNING: Box size change is not allowed because tau.rho==0 or (rescale&7)==0\n\
=> the reference box will be replaced by the actual box (and kept constant)");
          VV(L,=box.L) }
        else if (tau.rho>0) {
          rhoreached++;
          prt("Since tau.rho is specified, the box will be gradually rescaled until the\n\
reference box is reached."); } } }

    if (tau.P)
      prt("Since tau.P is specified, the box size will fluctuate; the reference box\n\
will be used for parameter setting.");

    VV(box.Lh,=0.5*box.L)

    box.cutoff=initcutoff(cutoff,L); /* set cutoff; needs No, L */

    if (rdf.cutoff>box.cutoff)
      WARNING(("rdf.cutoff=%g > cutoff=%g",rdf.cutoff,box.cutoff))
    if (rdf.cutoff==0) rdf.cutoff=fmin(box.cutoff,10);
    initrdf(rdf.grid,rdf.cutoff,Fn("s-s"));
#ifdef DIHHIST
    /* only histograms initialized - control structures already initialized */
    initdihhist(dih.grid);
#endif /*# DIHHIST */
#ifdef HARMONICS
    initharmonics(rdf.grid,rdf.cutoff,T,refrho);
#endif /*# HARMONICS */

    constrd.mode=rescale+RESCALE_L; /* force box rescale */

#ifdef SLAB
    if (slab.mode && !slab.grid)
      WARNING(("\
slab.mode specified, but slab.grid not given:\n\
*** surface tension and cutoff corrections will not be calculated\n"))
    if (slab.mode&2) constrd.mode|=RESCALE_SLAB;
    slab.qT=0;
    if (thermostat && slab.T>=0) {
      slab.qT=sqrt(slab.T/T);
      prt("separate thermostat T=%g for slab z in [%g,%g), otherwise T=%g",
          slab.T,slab.Tz0,slab.Tz1,T);
      if (slab.Tz0>=slab.Tz1) WARNING(("slab-based thermostat range is empty"))
      if (gear.order>2 || !strchr("\1\3\4\5\6",thermostat)) /* T_BERENDSEN, T_ANDERSEN .. T_MAXWELL_CM */
        WARNING(("slab thermostat not supported for -m%d and thermostat=%d",gear.order,thermostat))
        if (thermostat==T_BERENDSEN) WARNING(("The slab-based Berendsen thermostat may not work as you expect,\n\
*** please consult the manual.")) }

    initdpr(slab.grid,slab.max,Fn("dpr"));

#endif /*# SLAB */

#ifdef ECC
    el.epsf=rescalecharges(el.ecc,el.epsf);
    el.ecc=-abs(el.ecc);
#endif /*# ECC */

    /* moved here in V3.5i */
    if (option('k')>=0) initfix(Fn("fix"));
    checkfixed();

    initss(LJcutoff);
    if (poteps<0) {
      poteps=0;
      prt("poteps:=0 (potentials have been tested, further tests turned off)"); }

    initpot();
    initmolecules(corr); /* here cfg allocated */

#ifdef METAL
    eldensities(NULL); /* rho's are 0 for initcfg */
    eldens.rref=cfg[0]->rp; /* will be used by initcfg */
#endif /*# METAL */

#ifdef SLAB
    if (wall.is) initwall();
#endif /*# SLAB */

    initconstrit(omegac);

/* external forces initialization and info */
#include "mainextf.c"

    if (No.c && (thermostat==T_ANDERSEN || thermostat==T_MAXWELL))
      WARNING(("Constraint dynamics and atom-based Maxwell/Andersen thermostat:\n\
*** The results will be likely severely affected.\n\
*** Consider center-of-mass-based Andersen and Maxwell thermostats instead\n\
*** (thermostat=\"AndersenCM\" or thermostat=\"MaxwellCM\")."))

    En.corr=cutcorr(LJcutoff,corr); /* called because of counting rdf->npairs */
    if (corr&32) WARNING(("corr&32: the number of pairs of any atom is N^2/2 instead of N(N-1)/2\n\
*** this violates strict definition of RDF, but gives consistent cutoff\n\
*** corrections for replicated cell (e.g., with PERSUM)"))

    Hamaker();
                   
#include "mainkinc.c"

#ifdef COULOMB
    if (el.sf && (option('f') || gear.order>1)) ERROR(("\
IMPLEMENTATION LIMITATION: Structure factor (el.sf=%d) cannot be calculated\n\
*** during simulation or when forces (option -f) are requested.",el.sf))
#endif /*# COULOMB */

#ifdef SLAB
    if (corr&3 && slab.K) ERROR(("\
Homogeneous cutoff corrections (corr&3) and slab cutoff corrections\n\
*** cannot be used simultaneously."))
    if (slab.K && PROD(slab.n)>1) WARNING(("IMPLEMENTATION LIMITATION:\n\
*** Slab cutoff correction calculated incorrectly with box replication."))
#endif /*# SLAB */

#ifdef WIDOM
    if ( (widom.n || widom.spreal>=0) && thermostat==0) WARNING(("\
Widom insertion or scaling requires a thermostat because\n\
*** temperature T is used to calculate the insertion probability"))
    widom.corr=Widomcutcorr(widom.spreal,widom.sp,LJcutoff);
    if (!(corr&1)) {
      prt("Since !(corr&1), the above correction will not be included.\n");
      widom.corr=0; }
    else
      _n
    if (widom.spreal>=0 && widom.n)
      ERROR(("Widom insertion (widom.n>0) and Widom scaling (widom.spreal>=0) are exclusive"))
#endif /*# WIDOM */

#if defined(TWODIM) || defined(FREEBC)
    En.corr=0; /* not supported */
#endif /*# defined(TWODIM) || defined(FREEBC) */

    prt("Mass of the whole configuration in program units: ");
    put2(No.mass,No.free_mass)
    prt("mass = %.8g g/mol = %.8g kg", No.mass*Munit,No.free_mass*massunit);
    prt("rho = %.12g / V  (for rho in kg m^-3 and V in AA^3)", rhounit*No.mass);
    _n

    if (gear.order==1 && init>=0) init=-1;

    initgammas();

    if (init==0) {
      StaLoad(Fn("sta"));
      issta++;
      if (nomax>0) if (StaN("Tkin")>=nomax) {
        fprintf(stderr,"already nomax=%d cycles has finished - stop\n",nomax);
        prt("already nomax=%d cycles has finished - stop\n",nomax);
        exit(1); } }
    LRFree();

    if (init==0) loadrdf(rdf.grid);
#ifdef DIHHIST
    if (init==0) loaddih(dih.grid);
    if (dih.dcp) dihcp=fopen(Fn("dcp"),init<2?"at":"wt");
#endif /*# DIHHIST */
#ifdef HARMONICS
    if (init==0) loadharmonics(rdf.grid,rdf.cutoff);
#endif /*# HARMONICS */
#ifdef SLAB
    if (init==0) loaddpr(slab.grid);
#endif /*# SLAB */

    /* fix sites/axes by SHAKE constraint, measure forces/torques */
#ifdef ANCHOR
    initanchor(Fn("fix"),init,drmax);
#endif /*# ANCHOR */

    /*** load configuration 2nd pass or initialize ***/
    if (init>=3) {
      if (PROD(load.n)!=1 || load.N || SUM(load.L) || load.tr)
        WARNING(("load.n[], load.N, load.L[], load.tr: specified but ignored because init>=3"))
#ifdef COULOMB
      /* not clear why needed here ... see Ewaldtest() below */
      if (init<10) initrspace();
#endif /*# COULOMB */
      if (init>=10) {
        readasc(0,1);
        if (init>=11) readasc(1,1);
        if (init>=12) readasc(2,1);
        if (init>12) ERROR(("init=%d not implemented",init)) }
      else if (init==5)
        initcryst((int)pins,Emax);
      else {
        if (MC && option('y')) {
          openplayback(option('y'),0);
          plbopened=1; }
        initcfg(pins,Emax,
                option('y')?dt.plb==0?1:(int)(fabs(dt.plb)/h+0.5):0,
#ifdef SLAB
                slab.geom
#else /*# SLAB */
                option('d')<0?1-option('d'):0
#endif /*#!SLAB */
                ); } }
    else {
      int repl=initreplicate();

      if (init>=0) {
        loadcfg(-1,NULL); /* 2nd pass */
        if (repl>1 && simils.frommol) ERROR(("cannot replicate and add molecules at once"))
        if (init==2) t=0;
        if (stop>0 && t>stop-h/2) sig=3;
        if (simils.frommol) {
#ifdef COULOMB
          /* not clear why needed here ... see Ewaldtest() below */
          if (init<10) initrspace();
#endif /*# COULOMB */
#ifdef SLAB
          initcfg(0,Emax,0,slab.geom);
#else /*# SLAB */
          initcfg(0,Emax,0,option('d'));
#endif /*#!SLAB */
          simils.frommol=0; } }
      else {
        if (plbopened && simils.plbname)
          ERROR(("internal: attempt to open the same playback file for reading and writing"))
        readplayback(-init,1);
        prt("initialized from plb frame %d, continuing as with init=2",-init);
        init=2; }
      replicatecfg();
      VO(load.n,=1)
      VO(load.L,=0)
      load.N=load.tr=0; }

    if (option('y')>0 && (init<0 || init>1)) makemolgol(option('y'));

    if (oldt!=newt) {
      t=newt;
      prt("\nyou changed the running simulation time from t=%.3f to t=%.3f",oldt,newt);
      if (init==0 || init==1) WARNING(("changing t is strongly discouraged for init=\"cont\",\"append\"")) }

    if (gear.order>=2) {
      if (option('d')>=2) distancecheck();
      homogeneity(corr); }

    /* nullify some coordinates - only once */
    zerocfg(load.zero); load.zero=0;

    /* to suppress undefined/badly defined logs in case of thermostat changes */
    if (thermostat!=T_NOSE && thermostat<T_NPT) {
      loop (i,0,gear.order) if (cfg[i]->logs) {
        prt("NB: Nose extended variable^(%d):=0 because other thermostat",i);
        cfg[i]->logs=0; } }

    /* to suppress undefined/badly defined lambda in case of barostat change */
    if (thermostat<T_NPT) {
      loop (i,0,gear.order) if (cfg[i]->logs) {
        prt("NB: barostat extended variable:=0 (%d deriv)",i);
        VO(cfg[i]->lambda,=0) } }
#ifdef SLAB
      if (wall.n) {
        if (thermostat>=T_NPT)
          ERROR(("implementation limitation: cannot use the MTK barostat with wall(s)\n\
*** try Berendsen-like barostat with rescale=\"z\" or rescale=\"zCM\""))

        if (abs(wall.n)==3 && (constrd.mode&7)==4 && tau.P)
          prt("NPT in the slit pore: two walls, z-rescale, and tau.P is set:\n\
*** (Pwall[0]+Pwall[1])/2 will be used as the pressure to be maintained");
        else
          prt("WARNING: if a barostat were requested, the virial pressure would be used\n(which likely does not make sense)"); }

#  if SLAB & 2
    if (cleave.init&1) {
      cleave.init&=0x7ffffffe;
      cleave.z[0]=cleave.z0;
      cleave.z[1]=cleave.z1;
      prt("cleave.z0=%.9f cleave.z1=%.9f copied from data to cleave.z[0,1]\n\
cleave.init bit 0 unset (not to copy next time)"); }
    if (cleave.n==2) {
      if (fabs(cleave.z[0]-cleave.z[1])<1e-6)
        WARNING(("cleave.n==2 and cleave.z[0] == cleave.z[1] (with precision 1e-6)"))
      if (cleave.z[1]<cleave.z[0])
        ERROR(("cleave.n==2 and cleave.z[1]<cleave.z[0]")) }
#  endif /*# SLAB & 2 */
#endif /*# SLAB */

    if (init<0 || init>1) {
      En.tot=3e33; /* = will be ignored in the statistic of dEtot^2 */
      /* for Nose: restart En.tot */
      cfg[0]->logs=0; }

    /* fix sites by springs, -k0 hard, do not measure force */
    if (option('k')>=0) {
      makefixa();
      loadfixa();
      savefixa(); }

    if (initvel) Maxwell(0,initvel,1);
    initvel=0; /* just once */

    if (option('h')) {
      int valence=abs(option('h'));

      while (valence) {
        prt("WARNING: %d %d-groups or neighbors centered\n",
            centergroups(cfg[0],valence%10), valence%10);
        valence/=10; }
      if (option('h')>0) option('h')=0; }

    box.V=PROD(box.L);

    /* editing the configuration and velocities */
    if (key>=97) {
      int nn;
      char line[128];

      /* discard any text after ; to EOL */
      if (!fgets(line,128,in)) ERROR(("unexpected EOF while cfg edit expected"))

      loop (nn,0,key==97?2:1) {
        int n,nl=0,k,ned=0;

        if (key==98) n=1; else n=nn;
        prt_("reading/editing %s: ",n?"velocities":"positions");

        while (fgets(line,128,in)) {
          char *tok=strtok(line," \t");

          if (!tok) break;
          if (tok[0]==';') break;
          if (tok[0]=='.' && tok[1]<=' ') break;

          if (nl>=No.s) { ERROR(("number of atoms overflow nl=%d",nl)) nl=0; }
          loop (k,0,3) {
            double x=atof(tok);
            if (tok[0]!='=') {
              ned++;
              if (n) x/=h;
              cfg[n]->rp[nl][k]=x; }
            tok=strtok(NULL," \t"); }
          nl++; }
        prt("%d lines read, %d items replaced",nl,ned); } }

#ifdef LINKCELL
    lcsetup(L);
#endif /*# LINKCELL */

#ifdef SLIT
    if (wall.is && abs(wall.n)!=3)
      WARNING(("SLIT and less than 2 walls may cause error out of box"))
#endif /*# SLIT */

#if PARALLEL
    initparallel();
#  if PARALLEL==1 || PARALLEL==3
    parallelizerdf(1);
#  endif /*# PARALLEL==1 || PARALLEL==3 */
#endif /*# PARALLEL */
    addshift(0,nshift,shift);
    addshift(1,nshift,vshift);
    if (unsplit) {
      unsplitpbc(unsplit);
      unsplit=0;
      normalize(-2); }
    else
      if (init>=2 || init<0) normalize(-2);

    putv(L)
    putv(box.L)
    put2(box.cutoff,box.V)

    /* init=0,1: apppend playback; init=2,3: new playback */
    if (option('y')) {
      if (!plbopened) openplayback(option('y'),init<2&&init>=0);
      if (init>=2 && dt.plb>=0) writeplayback(); /* t=0: do'nt write for dt.plb<0*/ }

    /*    constrd.Pcorr=En.corr/Sqr(box.V)*Punit; removed in 2.7t */
    underline("cutoff corrections");
    prt("site-site potential cutoff=%g %s => cutoff correction for homogeneous fluid\n\
En.corr=%.9g [K AA3] (Ecorr=En.corr/V, Pcorr=En.corr/V^2)",
        fabs(LJcutoff),LJcutoff>0?"AA":"vdWsig",En.corr);
    if (tau.P)
      prt("NOTE: tau.P specified, the cutoff corrections in E and P depend on V\n\
      see variables Ecorr, Pcorr for their averaged values");
    else {
      prt("NVT: Ecorr=%g   Pcorr=%g  (in p.u.)",
          En.corr/box.V,En.corr/Sqr(box.V));
      prt("     Ecorr=%.9g [J/mol]  Pcorr=%.8g [Pa]",
          En.corr/box.V*Eunit,En.corr/Sqr(box.V)*Punit); }
    if ((corr&3) == 2 ) ERROR(("implementation limitation: corr=2 not supported, consult the manual!"))
    prt("Because of corr=%d, the above corrections are:\n\
  %sincluded in P, diagonal pressure tensor components, and E\n\
  %sincluded in P used by the barostat",corr,
        corr&1?"":"NOT ",
        corr&2?"":"NOT ");
    if (!(corr&1)) En.corr=0;

#ifdef SLAB
    if (slab.K) prt("Fourier transform slab cutoff correction for kz<slab.K=%d set",slab.K);
#endif /*# SLAB */

    halftcyc=h*noint/2;

    if (el.sf)
      prt("el.sf=%d: structure factor mode, elst forces not initialized");
    else {
#ifdef COULOMB
      underline("initialization of electrostatic forces");
      if (Ewaldtest(L)) goto TheEnd;
#  if COULOMB<0
      prt("Ewald final values: cutoff=%g el.alpha=%g el.kappa=%g Q->N=%d",
          box.cutoff,el.alpha,el.kappa,Q?Q->N:0);
#  else /*# COULOMB<0 */
      prt("no k-space (Ewald) electrostatics, final values:\n\
  cutoff=%g el.alpha=%g (el.kappa=%g ignored)",
          box.cutoff,el.alpha,el.kappa);
#  endif /*#!COULOMB<0 */
#endif /*# COULOMB */
    }

#ifdef MARK
    /* Not tested recently! */
    mark_do(no);
    /* Everything suppressed for !MARK */
#else /*# MARK */

    if (init%10>=3) {
      /* thus, forces are not calculated for init=10 */
      measure=1; /* debugging only! */
      /* necessary rhs if Verlet? */
#  ifdef POLAR
      { /* OOPS! scf.omega may be -9 here to be set later... */
        double o=scf.omega;
        scf.omega=0.9;
#    if POLAR&32
        initfqcharges();
#    endif /*# POLAR&32 */
        rhs(cfg[2],cfg[0],cfg[1]);
        scf.omega=o; }
#  else /*# POLAR */
      rhs(cfg[2],cfg[0],cfg[1]);
#  endif /*#!POLAR */
      rescalecfg(cfg[2],RESCALE_XYZ,h*h/2,NULL); /* for Gear only? */

#  ifdef POLAR
      /*.....  memset(polarrof(molec,cfg[2]->rp),0,No.s*sizeof(vector));*/
      loop (i,1,gear.order)
        copy(polarrof(molec,cfg[i]->rp),polarrof(molec,cfg[0]->rp),No.s*sizeof(vector));
#  endif /*# POLAR */
    }

#  if 0
#    include "foolproo.c"
#  endif /*# 0 */

    /* columns 4,5 of the convergence profile */
#  ifdef ECC
    if (tau.P==0 && tau.rho==0 && tau.sig==0)
      // prior V2.8c:      initCP(no,CPnbit,dV==0?(thermostat==T_NOSE?"Ep+k":No.first?"Ep0":"Ein"):constrd.PdVname,"P",corr);
      initCP(no,CPnbit,dV==0?(thermostat==T_NOSE?"Ep+k":FROM?"Ep0":"rho"):constrd.PdVname,"PECC",corr);
    else /* "P" changed to "Pcfg" in V3.6l, but in V3.7a changed to "Pref" */
      WARNING(("tau.P, tau.rho, tau.sig not supported with ECC except initialization"))
      initCP(no,CPnbit,tau.sig?"svdW":"rho",dV==0?"Pref":constrd.PdVname,corr);
#  else /*# ECC */
    if (tau.P==0 && tau.rho==0 && tau.sig==0)
      // prior V2.8c:      initCP(no,CPnbit,dV==0?(thermostat==T_NOSE?"Ep+k":No.first?"Ep0":"Ein"):constrd.PdVname,"P",corr);
      initCP(no,CPnbit,dV==0?(thermostat==T_NOSE?"Ep+k":FROM?"Ep0":"rho"):constrd.PdVname,"Pref",corr);
    else
      initCP(no,CPnbit,tau.sig?"svdW":"rho",dV==0?"Pref":constrd.PdVname,corr);
#  endif /*#!ECC */

    if (tau.sig==0) {
      if (sigvdW!=0) prt("\
WARNING: no sigvdW adjustment (tau.sig=0), but sigvdW=%g in cfg-file found.\n\
This value will be ignored and the force-field (ble-file) value will be used.\n\
sigvdW=0 was set for a subsequent write to the cfg-file, so that next time\n\
tau.sig is specified, sigvdW will be initialized from the ble-file",sigvdW);
      sigvdW=0; }
    else {
      nbfix_t *fix=fixij(tau.i,tau.j);

      sigvdWptr=&fix->onefour[0].sig;

      if (!fix) ERROR(("sites tau.i=%d tau.j=%d not found in the nbfixes table",tau.i,tau.j))
      prt("tau.i=%d tau.j=%d in the nbfixes table: sig=%g eps=%g",
          fix->indx[0],fix->indx[1],fix->onefour[0].sig,fix->onefour[0].eps);

      if (sigvdW==0) {
        sigvdW=*sigvdWptr;
        prt("sigvdW adjustment: sigvdW=0 (loaded or initialized) replaced by force field value %g",sigvdW); }
      else {
        prt("sigvdW adjustment: sigvdW=%g loaded from %s",sigvdW,Fn("cfg"));
        *sigvdWptr=sigvdW; } }

#  ifdef DIHHIST
    gauchetrans(-3);
#  endif /*# DIHHIST */

    starttime=inittime(option('t')?8:0);
    nsteps=0;

#  ifdef POLAR
    underline("POLAR: SCF iterations");
    prt_("Method:\n  ");
    switch (option('p')/10%10) {
      case 0: prt("Car-Parrinello-like without thermostatting to T=0 (do not use)"); break;
      case 1: prt("no prediction (last value)"); break;
      case 2: prt("Always Stable Predictor-Corrector (2nd order), recommended"); break;
      case 3: prt("Partially Stable Predictor-Corrector (3rd order), rarely useful"); break;
      case 4: prt("Partially Stable Predictor-Corrector (4th order), not recommended"); break;
      default: ERROR(("bad option -p-pKPC (value of P=%d)")) }

    prt("Additional predictor length: k=%d",(option('p')/100)%10);

    if (scf.eps>2e33) scf.eps=1;
    put3(scf.eps,scf.omega,scf.maxit)
    scf.eps0=scf.eps;
    if (scf.epsq>=1) ERROR(("scf.epsq=%g is invalid",scf.epsq))
    if (scf.eps>0) {
      if (scf.epsq>0) scf.eps *= 1-scf.epsq;
      if (scf.epsx>2e33) scf.epsx=scf.eps0*0.001; }
    else {
      scf.epsq=0;
      if (scf.epsx>2e33) scf.epsx=(int)(scf.eps*2-1); }
    put(scf.epsq)
    put2(scf.epsx,scf.omegax)
    if (tau.dip) prt("WARNING: tau.dip=%g: Car-Parrinello-like (extender Lagrangian) is deprecated",tau.dip);
#  endif /*# POLAR */

    {
      char *st=myctime(starttime);

      if (st[24]=='\n') st[24]=0;
      /* this string used by cutprt.c - do not change! */
      fprintf(stderr,"%d cycle%s started at %s, t=%.3f\n",no,"s"+(no==1),st,t);
      prt("\n*** %d cycle%s started at %s, t=%.3f",no,"s"+(no==1),st,t);
    }

    if (option('b')) {
      if (!option('s')) {
        prt("batch mode -b%d: check for %s (to stop) every %d %s",
            option('b'),Fn("stp"),abs(option('b')),option('b')<0?"s":"cycle(s)"); }
      else
        prt("interactive mode -b%d (beep level)",option('b')); }

    if ((drift&DRIFT_WHEN)==DRIFT_START) {
      if (drift&0777) WARNING(("\
No CoM/v/angM drift corrected anymore because drift&%d specified\n\
final drift = %d = %s",DRIFT_START,drift,int2sumbin(drift)))
      else prt("no drift correction");
      // drift&=DRIFT_DEPEND; DRIFT_DEPEND removed (always performed)
      drift=0; }

    initSF();
    if (MSD.mode&3) initdiff(gear.order==0?dt.cfg*reread.by:
                              gear.order==1?dt.plb*reread.by:
                              h*noint);
#  ifdef CLUSTERS
    readclusterdef(Fn("cli"));
#  endif /*# CLUSTERS */

    forDUMP=cfg[0]; /* debugging: see internp.C */

    if (lastspeed) {
      lastspeed=h*noint*no/lastspeed;
      prt("predicted duration of next job = %.5g h = %.5g days",
          lastspeed/3600., lastspeed/86400.); }
    No.t0=t;

    if (option('v')&4) printmasses();

    fflush(out);

    /* <<<<<<<<<<<<<<<<<<<<< main cycle = no sweeps <<<<<<<<<<<<<<<<<<<<< */
    for (icyc=0; icyc<no && sig==0; icyc++) {

      lastEtot=En.tot;

      if (option('d')==2) distancecheck();

      /* reading previously stored trajectory instead of one sweep */
      if (gear.order<2) {
        reread.frame=reread.from+icyc*reread.by;

        if (gear.order==0) {
          double dummy;
          if (dt.cfg<=0) ERROR(("Re-read mode (-m0) and dt.cfg<=0.\n\
*** Set dt.cfg to the stride between PLBNAME.# files."))
          /* read full cfgs instead of simulating */
          loadcfg(reread.frame,&dummy); /* 1st pass */
          loadcfg(reread.frame,NULL); /* 2nd pass */
          t=reread.frame*dt.cfg; }
        else {
          /* read playback instead of simulating */
          if (dt.plb<=0) ERROR(("Re-read mode (-m1) and dt.plb <= 0.\n\
*** Set dt.plb > 0 to the stride between PLBNAME.# files."))
          if (dt.cfg!=0 && dt.cfg<dt.plb)
            ERROR(("Re-read mode (-m1) and dt.cfg(%g) < dt.plb*reread.by(%g).\n\
*** Set dt.cfg, preferably to an integer multiple of dt.plb*reread.by.",
                   dt.cfg,dt.plb*reread.by))
          readplayback(reread.frame,1);
          t=(reread.frame-1)*dt.plb; }

        prt("frame %d of %d read, calculated time=%g",reread.frame,no,t);
        measure=1; }

      /* <<<<<<<<<< internal loop: one sweep = noint MD steps <<<<<<<<< */
      else loopto (iint,1,noint) {
#  ifdef POLAR
        /* hack for conductivity of polarizable models:
           store previous Drude positions */
        if (iint==noint) {
          if (!lastrpols) allocarray(lastrpols,sizeof(vector)*No.s);
          copy(lastrpols,polarrof(molec,cfg[0]->rp),sizeof(vector)*No.s); }
#  endif /*# POLAR */

        /* default time interval for statistics inside the noint-cycle */
        StaSet(DT/noint,lag.ierr,2,lag.in);

        if (option('s') && option('b')>3)
          fprintf(stderr,"%d/%d %d/%d t=%f   \r",iint,noint,icyc,no,t);

#  if defined(SLAB) && SLAB & 1
        /* hack for "zone melting", melting point in the slab geometry */
        if (tau.L>0)
          box.V=rescalecfg(cfg[0],rescale|RESCALE_L,0,Pscale);
#  endif /*# defined(SLAB) && SLAB & 1 */

        /* friction-like isobaric ensemble: Pscale kept during cycle */
        if (thermostat<T_NPT && tau.P>0)
          box.V=rescalecfg(cfg[0],rescale|RESCALE_L,0,Pscale);

        /* adjusting a cross size site-site parameter */
        if (tau.sig>0) {
          *sigvdWptr/=Pscale[0];
          setss(&sstab[tau.i][tau.j],tau.i,tau.j,LJcutoff,0);
          if (tau.i!=tau.j) sstab[tau.j][tau.i]=sstab[tau.i][tau.j]; }

        if (option('d')>2) distancecheck();

        if (sig==5) break;

        measure=iint==noint;

        if (gear.order==2) {
          /* Verlet + SHAKE */
          if (thermostat==T_LANGEVIN_CM || thermostat==T_LANGEVIN)
            ERROR(("Verlet+Langevin not implemented (yet)"))

          Shake(epsc,
                thermostat==T_MAXWELL||thermostat==T_MAXWELL_CM ? (double)justnow(tau.T,h/2) :
                thermostat==T_ANDERSEN||thermostat==T_ANDERSEN_CM ? h/tau.T : 0); }
        else {
          /* Gear integration */
#  ifdef POLAR
          Gear2pol(No.eq,No.s*DIM,(option('p')/10)%10,(option('p')/100)%10,cfg);
#  else /*# POLAR */
          Gear(No.eq,cfg);
#  endif /*#!POLAR */
          CPUtime("Gear"); }

#  if defined(COULOMB) && COULOMB<0
        /*
           Every-step Ewald-based dipole moment summation statistics
           incl. optional microwaves (Eext.f).  NB: Cannot move to
           ewald.c because could be called several times for POLAR.
        */
        if (Q) {
          loop (i,0,DIM) {
            /* StaAdd("Mx"), etc. */
            StaAdd(string("EwM%c",i+'x'),Q->M[i]);
            if (Eext.E[i] && Eext.f) {
              StaAdd(string("EwM%c*cos",i+'x'),Q->M[i]*cos(Eext.arg[i]));
              StaAdd(string("EwM%c*sin",i+'x'),Q->M[i]*sin(Eext.arg[i])); } }
          StaAdd("EwM^2",SQR(Q->M)); }
#  endif /*# defined(COULOMB) && COULOMB<=-1 */

        nsteps++;

        /* Gear-based thermostats, constraints, etc. */
        if (gear.order!=2) {
          int n;
          vector *rp=cfg[0]->rp;
          vector *rpvel=cfg[1]->rp;

          if (thermostat==T_ANDERSEN) Maxwell(-1,No.N,h/tau.T);
          if (thermostat==T_ANDERSEN_CM) MaxwellCM(h/tau.T);
          if (thermostat==T_MAXWELL && justnow(tau.T,h/2)) Maxwell(-1,No.N,1);
          if (thermostat==T_MAXWELL_CM && justnow(tau.T,h/2)) MaxwellCM(1);

          if (measure && option('c')&4) {
            En.r2=constrainterror(cfg[0],cfg[1]);
            En.v2=vconstrainterror; }

          if (!(option('c')&1)) {
            loop (n,FROM,No.N) if (molec[n].nc)
              constrit[molec[n].sp].nit[COR_BOTH] +=
                Scorrect(molec+n,rof(molec+n,rp),rpvel,epsc); }
          else
            Lcorrect(cfg[0],cfg[1]);
          CPUtime("corr"); }

        /* polishing and final constraint check */
        normalize(-1);

        if (measure && option('c')&8) {
          En.r3=constrainterror(cfg[0],cfg[1]);
          En.v3=vconstrainterror; }

        CPUtime("norm");

        if (option('y')) if (justnow(dt.plb,h/2)) writeplayback(); // for dt.plb!=0

#  ifdef POLAR
        if (scf.epsq>0) scf.eps = scf.eps*scf.epsq + (1-scf.epsq)*scf.eps0;
#  endif /*# POLAR */

        if (No.c) doconstrit(omegac);

        bounddr(drmax); /* every step since V3.6k */

        /* remove CoM/momentum/angM if requested */
        if ((drift&DRIFT_WHEN)==DRIFT_STEP) removedriftssta(option('v')&64);

        {
          /* for sure .. fool-proof test (remove later if OK) */
          double tt=No.t0+(noint*icyc+iint)*h;
          tt=(t-tt)/(t+tt);
          //          if (option('v')&4) prt("tt=%g",tt);
          /* should be |tt|< 2^-53, for sure I use 2^-51 */
          if (fabs(tt)>4.441e-16) WARNING(("t+=h error (%g) > 2^-51",tt))
        }

        /* to avoid cumulation of errors from t += h (as checked above) */
        t=No.t0+(noint*icyc+iint)*h;

        init_append|=1; } /* >>>>>>>>>>>>>>> iint >>>>>>>>>>>>>>> */

      /* re-set DT = 1 cycle */
      StaSet(DT,lag.err,2,lag.n);

      /* calculate forces etc. in the re-read mode */
      if (option('f')) {
#  ifdef GAUSSIANCHARGES
        static int pass=1;

        if (pass) {
          WARNING(("Virial pressure is not implemented in the re-read mode\n\
*** => pressure is wrong. See constrd.c how to fix it!"))
          pass=0; }
#  endif /*# GAUSSIANCHARGES */
        zeroEn();
        StaSet(0,2,2,0);
#  ifdef POLAR
        /* accurate (based on scf.epsx,scf.omegax) */
        memset(polarrof(molec,cfg[0]->rp),0,No.s*sizeof(vector));
        scforces(cfg[gear.order+1],cfg[0]); /* incl. En.pot += En.self + En.el; */

        /* taken from constrd.c: */
        StaAdd("Eel(excl. self-term)",En.el);
        StaAdd("Eel(incl. self-term)",En.el+En.self);
        StaAdd("Eself",En.self);
        // En.vir+=En.virc+En.el-En.self*2
        En.vir -= En.self*2;
#  else /*# POLAR */
        // ???? ecc(0);
        forces(cfg[gear.order+1],cfg[0]);
        En.pot += En.el; /* do not move to forces(), but is in scforces() */
#  endif /*#!POLAR */
        En.vir+=En.virc+En.el; /* En.virc zero for recalculation! */
        En.U=En.Unc=En.kin+En.pot;
#  ifndef FREEBC
        En.U+=En.corr/box.V;
#  endif /*# FREEBC */
        En.tot=En.U+En.ext; /* cf. constrd.c */
        En.kin=No.f*T/2;
        /* The kinetic pressure contains the nominal temperature T!

           Note that this is multipled by No.Pkinq later; with corr&4
           (kinetic pressure correction), No.Pkinq=No.f0/No.f, which means
           that the ideal gas EOS is exactly recovered.

           In NVT and without virial of constraint forces, the
           following correction (taken, e.g., from .cp) should be added:
              (Tkin-T)*Punit/box.V/3*No.Pkinq
           where No.Pkinq=kinetic pressure correction, box.V is in AA^3)

           The virtual volume change is calculated elsewhere */
      } /* option('f') */

      if (dt.cfg!=0 && gear.order<2) {
        static int lastn=-1,n;

        /* save during re-read mode (cf. saving during simulation below) */
        stoptime=mytime();
        if (tau.sig) sigvdW=*sigvdWptr; /* probably nonsense here */
        n=t/dt.cfg+0.5;
        if (n==lastn) n++;
        lastn=n;
        savecfg(n,(int4)stoptime,&sigvdW); }

      if (h!=0 && !isfinite(En.tot)) {
        sig=2;
        WARNING(("En.tot is infinite or not a number => stop as if Ctrl-C")) }

      /* statistics of constraint errors */
      if (No.c) {
        int i,sp;

        /* constraint iterations reported every cycle even if known every step */
        /* ... not for SHAKE */
        StaSet(0,2,2,0);
        loop (i,0,COR_LAST) loop (sp,0,nspec) if (constrit[sp].nit[i])
          StaAdd(string("constrit%d spec%d",i,sp),(double)constrit[sp].nit[i]/spec[sp]->N);
        if (En.v1) StaAdd("v constr err before",En.v1);
        if (En.v2) StaAdd("v constr err mid",En.v2);
        if (En.v3) StaAdd("v constr err after",En.v3);
        if (En.r1) StaAdd("constr err before",En.r1);
        if (En.r2) StaAdd("constr err mid",En.r2);
        if (En.r3) StaAdd("constr err after",En.r3); }

      //if (drift&DRIFT_DEPEND) /* new in 2.0b */
      depend_r(cfg[0],0); /* always in 3.6k (for sure, perhaps not needed) */

      box.V=PROD(box.L);
#  ifdef FREEBC
      box.V=1; /* ugly! */
#  endif /*# FREEBC */

      /*** pressure [Pa] ***/
      /* virial-based pressure w/o cutoff corrections (elst virial=energy) */
      En.Pevir.n=((En.kin*2*No.Pkinq+En.vir)/DIM)/box.V; // former En.Pnc

      /* as above, w. cutoff corrections added */
      En.Pevir.c=En.Pevir.n+En.corr/Sqr(box.V);
      // to move:      if (option('v')&4) StaAdd("En.P double check [Pa]",((En.kin*2*No.Pkinq+En.vir)/DIM+En.corr/box.V)/box.V*Punit);

#  ifdef ECC
    /* Since V3.4a, En.P becomes the ECC pressure,
       in statistics and CP, these are still marked by ECC */
      // TO BE FIXED AFTER CHANGES: should not *=Punit, use En.Pevir
      En.ECC_Pcorr*=Punit;
      En.Pnc+=En.ECC_Pcorr;
      En.ECC_Pvir=En.P; /* Pvir - see the paper */
      En.ECC_Pscaled=En.ECC_Pvir-En.ECC_P3ref; /* conventional w/o ECC */
      En.P+=En.ECC_Pcorr; /* true ECC pressure */
#  endif /*# ECC */

#  if (PRESSURETENSOR&PT_ANY) == PT_ANY
      /* pressure w. cutoff corrections based on the directly calculated pressure tensor */
      if (option('v')&4) StaAdd("En.trPt double check [Pa]",(SUM(En.Ptens)/3+En.corr/(box.V*box.V))*Punit);
#  endif /*# (PRESSURETENSOR&PT_ANY) == PT_ANY */

      /*** STOP control ***/
      if (nomax>0 && (issta?StaN("Tkin")+1:1)>=nomax) sig=3;

      if (stop>0 && t>stop-h/2) sig=3;

      if (Tstop) {
        static double Tlast=0;
        if (Tlast) {
          if ((En.T-Tstop)*(Tlast-Tstop)<0) sig=5; }
        Tlast=En.T; }

      /*** hack: check of remaining disk space ***/
      if (option('w')>1) expectedbuffer+=16;
      if (expectedbuffer>7500) {
        expectedbuffer=0;
        waitfordiskspace(64);
        fflush(out); }

      /*** print a line of results ***/
      if ( option('v')&1 && (icyc==0 || icyc+1==no || sig || justnow(dt.prt,halftcyc))) {

        if (option('w')>1) expectedbuffer+=88;

#  ifdef POLAR
        if (icyc==0) {
          if (tau.sig) header("   t/ps    Tkin/K   Epot/K    Etot/K   sigvdE/AA P/MPa  polerr polit  cerr1 ");
          else header("   t/ps    Tkin/K   Epot/K    Etot/K  rho/kgm^-3 P/MPa  polerr polit  cerr1 ");
          fflush(out); }                                                                    
        if (tau.sig)
          prt("%9.3f %6.1f %9.0f %12.2f %7.4f %8.3f%9.2e%2.0f%10.2e",
                t,    En.T,En.pot,En.tot,*sigvdWptr,
                                              En.Pref*(Punit/1e6),En.polstderr,En.polnit,En.r1);
        else
          prt("%9.3f %6.1f %9.0f %12.2f %7.2f %8.3f%9.2e%2.0f%10.2e",
               t,    En.T,En.pot,En.tot,No.mass/box.V*rhounit,
                                              En.Pref*(Punit/1e6),En.polstderr,En.polnit,En.r1);
#  else /*# POLAR */
        double pcer=En.r1+En.r2+En.r3;

        if (pcer>0) pcer=-log10(pcer); else pcer=9.99;

        if (icyc==0) {
          prt("\nNOTE: pcer = -log10(summary constraint error), recommended pcer>6");
          if (T) prt_("NOTE: the graph of Tkin is in the range of [0,2*T]=[0,%g] K",T*2);
          else prt_("NOTE: the graph of Tkin is in the range of [0,0.2] K");
          if (tau.sig) header("   t/ps    Tkin/K   Epot/K    Etot/K   sigvdE/AA P/MPa pcer0-------Tkin------2T");
          else         header("   t/ps    Tkin/K   Epot/K    Etot/K   rho/kgm-3 P/MPa pcer0-------Tkin------2T"); }

        /*** standard nonpolar prt ***/
        if (tau.sig)
          prt_("%9.3f %6.1f %9.0f %12.2f %7.4f %6.1f%5.2f",
                t,    En.T,En.pot,En.tot,*sigvdWptr,
                                               En.Pref*(Punit/1e6),pcer);
        else
          prt_("%9.3f %6.1f %9.0f %12.2f %7.2f %6.1f%5.2f",
                t,    En.T,En.pot,En.tot,No.mass/box.V*rhounit,
                                               En.Pref*(Punit/1e6),pcer);
        graph(En.T/(T+0.1*(T==0))*0.5,20);

        if (gear.order>2 && option('v')&4)
          if (En.r1>=0 || En.r2>=0 || En.r3>=0) prt(
           "                                                 %10.2e%10.2e%10.2e",
                                                             En.v1,En.v2,En.v3);
#    ifdef erminmax
        if (gear.order>2)
          prt("%8.5f %8.5f",sqrt(Erfc.min),sqrt(Erfc.max));
        Erfc.min=3e33; Erfc.max=0;
#    endif /*# erminmax */
#  endif /*#!POLAR */
      } /* print a line of results */

#  if !defined(FREEBC) && !defined(PERSUM)
      {
        double minLh=fmin(fmin(box.Lh[0],box.Lh[1]),box.Lh[2]),perc;
        static double oldperc;
        static int info=1;

        if (box.cutoff>minLh*1.005) {
          perc=(box.cutoff/minLh-1)*100;
          if (info)
            prt("\
Explanation: The box size is variable in NPT simulation while the cutoff is\n\
  constant. If the cutoff is inappropriately set (e.g., to minimum L/2 at\n\
  simulation start which is the default), it may happen that the cutoff\n\
  exceeds minL/2. Overshooting by less than 0.5 %% is acceped without warning.\n\
  More than a few %% is inaccurate and not acceptable in productive runs");
          if (perc>oldperc)
            WARNING(("cutoff=%f exceeding minL/2=%f by %.2f%%\n\
*** in particular cases (no charges, no rdf, no LINKCELL) this may be bening\n\
*** (more warnings for cutoff exceeding minL/2 by more than %.2f%%)",
                      box.cutoff,minLh,perc,oldperc=perc*1.3))
#    ifdef LINKCELL
            if (perc>70) ERROR(("...this is too much - stop"))
#    else /*# LINKCELL */
            if (perc>300) ERROR(("...this is too much - stop"))
#    endif /*#!LINKCELL */
            info=0; }
      }
#  endif /*# !defined(FREEBC) && !defined(PERSUM) */

      if (sig==5) break; /* 2nd attempt */

      depend_r(cfg[0],0); /* correct r for measurements (e.g., dipole moment) */

#  include "maincps.c"

      if (dt.cfg!=0 && gear.order>=2) {
        /* save during reread mode is elsewhere */
        int n=(t+halftcyc)/dt.cfg;

        if ( n > (int)((t-halftcyc)/dt.cfg) ) {
          stoptime=mytime();
          if (tau.sig) sigvdW=*sigvdWptr;
          savecfg(n,(int4)stoptime,&sigvdW); } }

      /* Berendsen-like simple NVE and NPH ensemble via velocity rescaling */
      if (tau.E!=0) {
        double Etot=En.tot;

        if (tau.P) Etot+=No.P*box.V; /* enthalpy in K; P is the parameter in Pa */
        rescalecfg(cfg[1],RESCALE_XYZ,exp(halftcyc*(E-Etot)/En.kin/tau.E),NULL); }

#  include "mainresc.c"

      /* V=L[0]*L[1]*L[2];  USED TO BE HERE - WHY? */

      init_append|=2;

#  ifdef FREEBC
      /* stop simulation if any atom farther than # from (0,0,0) */
      if (option('d')<0) {
        vector *r=cfg[0]->rp;
        double rr,maxrr=0;

        loop (i,0,No.s) {
          rr=SQR(r[i]);
          Max(maxrr,rr) }
        if (maxrr>Sqr((double)option('d'))) sig=6; }
#  endif /*# FREEBC */

#  ifdef SLAB
      /* molecule to be removed found */
      if (slab.out) {
        removemol.n=slab.torem;
        if (removemol.n>=0) {
          measure1drift(removemol.n);
          sig=7; } }
#  endif /*# SLAB */

      /* raises sig if SIMNAME.stp exists - batch mode */
      if (!option('s')) {
        /* batch mode SIMNAME.stp N-*/
        int checkstp=0;

        if (option('b')>0 && icyc%option('b')==option('b')-1) checkstp++;

        if (option('b')<0) {
          static double time0;
          if (time0) {
            double time=mytime();
            if (time-time0>-option('b')) {
              checkstp++;
              time0=time; } }
          else
            time0=mytime(); }

        if (checkstp) {
          FILE *stp=fopen(Fn("stp"),"rt");
          if (stp) {
            fprintf(stderr,"%s exists, interrupting\n",lastFn);
            sig=4;
            fclose(stp); } } }

#  ifdef POLAR
      /* scf.omega autoset */
      if (scfautoset(icyc,noint)) sig=99;
#  endif /*# POLAR */
    } /* >>>>>>>>>>>> icyc (# of cycles=no, or until ^C) >>>>>>>>>>>>> */

#  ifdef POLAR
    if (sig==99) {
      sig=0;
      if (scf.margin>0)
        WARNING(("scf.margin was positive, set to %g",scf.margin=-scf.margin))
        prt("SCF AUTOSET divergence for omega=%g => scf.omega:=%g\n\
scf.domega:=0 (no autoset in the next sweep)",
            scf.omega,scf.omega+scf.margin);
      fprintf(stderr,"SCF AUTOSET divergence for omega=%g => scf.omega:=%g\n",scf.omega,scf.omega+scf.margin);
      scf.omega+=scf.margin;
      scf.domega=0; }
#  endif /*# POLAR */

    if (option('v')&1) header("");
    if (option('d')) distancecheck();

#  ifdef POLAR
    if (scf.epsq>0) {
      put3(scf.epsq,scf.eps,scf.eps0)
      scf.eps=scf.eps0; scf.epsq=-fabs(scf.epsq); }
    prt("%g  %g %g  %g %g  %g %g  POLAR SCF:omega  iter err  1maxerr err  1stderr err  lastmaxerr err  laststderr err",
        scf.omega,
        StaMean("polar no of iter"),StaStdErr("polar no of iter"),
        StaMean("polar 1st iter stderr"),StaStdErr("polar 1st iter stderr"),
        StaMean("polar last iter stderr"),StaStdErr("polar last iter stderr"));
#  endif /*# POLAR */

    underline("sweep timing");
    if (gear.order>=2)
      prt("planned sweep (one ';' in the data) = %d cycles by %d steps by %g ps", no,noint,h);
    else
      prt("planned re-read mode: from %d to %d by %d",
          reread.from, reread.to, reread.by);

    stoptime=mytime();
    stop_start=stoptime-starttime;

    prt_("stopped at real time = %.4f s = %s",stoptime,myctime(stoptime));
    prt("duration %.4f s = %.3f min = %.5f h",stop_start,stop_start/60.,stop_start/3600.);
    prt("stopped at MD time t=%.5f ps",t);
    if (gear.order>=2)
      prt("sweep MD time = %.11g ps = %d MD steps (%d cycles by %d steps)", nsteps*h,nsteps,icyc,noint);
    if (nsteps && stop_start>0) {
      lastspeed=(double)nsteps/stop_start;
      prt("speed = %g cycles/min = %g cycles/hour = %g ns/day",
          60.0*lastspeed/noint,
          3600.0*lastspeed/noint,
          86.4*(h*lastspeed));
      lastspeed*=h; }

    if (gear.order>=2)
      fprintf(stderr,"%d of %d*%d steps (t=%.3f ps) stop at %s", nsteps,icyc,noint,t,myctime(stoptime));
    else
      fprintf(stderr,"t=%.3f ps stop at %s", t,myctime(stoptime));

    if (option('t')) {
#  ifdef CHEAPTIME
      if (nsteps) {
        i=(int)stop_start%nsteps;
        if (i>nsteps/2) i=nsteps-i;
        if (i<=(int)sqrt(nsteps)) prt("CHEAPTIME: possible interference"); }
#  endif /*# CHEAPTIME */
      printtime(); }

#  if PARALLEL
    printpartimes();
#  endif /*# PARALLEL */
    _n

    if (option('s') && option('b')>1) fputc('\a',stderr);

    /* results of this run (one ; in data)
       NB: mainres.c contains release(mark) */

#  include "mainres.c"

#endif  /*#!MARK */
    if (tau.sat<0) init=1; /* special: see mainres.c - some problems?! */
    else init=0; /* default in every pass */
  } /* >>> run >>> */

 TheEnd:

 checkranges(1);

  if (sig>0) {
    if (sig==3) prt("finished because nomax or stop reached");
    else if (sig==4) prt("quitted because of %s",Fn("stp"));
    else if (sig==5) prt("quitted because of Tstop=%g reached",Tstop);
    else if (sig==6) prt("quitted because of a molecule escaped");
    else if (sig==7) prt("quitted because of a molecule was removed");
    else prt("quitted with sig=%d",sig); }

  if (fclose(out)) Error("close out");
  fclose(in);

  if (option('s') && option('b')>0) fputc('\a',stderr);

  return 0;
}
