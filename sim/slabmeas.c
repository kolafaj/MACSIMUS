/* #included from simmeas.c */

/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        density profiles in z-direction
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

static struct dpr_s {
  int4 size;          /* SDS size */
  int n;              /* # of bins; negative value means [0,Lz] scaled to [0,1]
                         Should be the same for all density profiles! */
  int nmeas;          /* number of measurements */
  int xSP;            /* removed in V3.6l, see slab.sp instead */
  char name[8];       /* atom name or "mol0", "mol1",.. */
  double grid;        /* per 1 AA (if n>0) or per Lz (if n<0); all grids are the same */
  double Vbin;        /* sum of volumes of one bin (to be divided by nmeas) */
  double Lz;          /* sum of Lz (to be divided by nmeas) */
  double mass_charge; /* molecular mass (or charge - not used) */
  union {
    unsigned ihist[1]; /* for site (number) density profiles */
    double   dhist[1]; /* for mass and charge density profiles */
  } hist; /* [n] */
} **dpr,**dpr0,*dprz;
/* indices of dpr,dpr0:
  [0,nspec) : centers-of-mass of molecules, double dhist[]
  [nspec,2*nspec) : charges (by molecules), double dhist[]
  [2*nspec,2*nspec+nsites) : sites by types, unsigned ihist[] */

static long unsigned (*sumihist)[2]; /* 0: original, 1: extended */
static int extext; /* extext=1: extended dpr (the original = dpr0) */
static double surftenscorr[2];

void initdpr(double grid,double zmax,char *fn) /************* initdpr */
{
  int n,sp;

  if (grid==0) return;
  if (grid<0) ERROR(("initdpr: grid=%g<0 is invalid",grid))

  underline("density profile and walls");
  if (zmax==0) {
    n=(int)grid;
    if (n!=grid) ERROR(("integer number of grid points expected"))
    prt("slab.max=0 and slab.geom=%d: the density profile range is %s divided into %d bins",
        slab.geom,
        slab.geom%3==2?"Lz":slab.geom%3==1?"sqrt(Ly*Lz)":"cbrt(Lx*Ly*Lz)",n);
    /* negative n indicates dynamic range and for slab.geom also z-periodicity */
    n=-n; }
  else {
    n=zmax*grid+0.5;
    if (fabs(n/grid/zmax-1)>1e-13) ERROR(("slab.max is not an integer multiple of the bin width\n*** (prohibited in V2.9g because of problems with defining the charge profile)")) }
  prt("initdpr: slab.grid=%g  slab.max=%g  #_of_histogram_bins=%d",grid,zmax,n);
  if (abs(n)<2) WARNING(("only 1 histogram bin in the z-profile"))

  rallocarray(dpr,2*nspec+nsites);
  rallocarray(sumihist,nsites);

  loop (sp,0,2*nspec+nsites) {
    int elemsize=
      sp<2*nspec ? sizeof(dpr[sp]->hist.dhist[0])
                 : sizeof(dpr[sp]->hist.ihist[0]);
    sdsralloczero(dpr[sp],sizeof(dpr[sp][0])+elemsize*(abs(n)-1));
    dpr[sp]->n=n; /* must be the same for all density profiles */
    dpr[sp]->grid=grid; /* must be the same for all density profiles */
    if (sp<nspec) sprintf(dpr[sp]->name,"mol%d",sp);
    else if (sp<2*nspec) sprintf(dpr[sp]->name,"mol%d",sp-nspec);
    else strcpy(dpr[sp]->name,sitedef[sp-2*nspec].name); }

  loop (sp,0,nspec) {
    /* unnecessarily recalculated (with check) */
    int i;
    double m=0,q=0;

    loop (i,0,spec[sp]->ns) {
      siteinfo_t *si=spec[sp]->si+i;
      m += si->mass;
      q += si->charge;
#ifdef POLAR
      q += si->chargepol; /* bug fixed by T.Trnka */
#endif /*# POLAR */
    }
    dpr[sp]->mass_charge=m;
    dpr[sp+nspec]->mass_charge=q; /* not used */
    if (fabs(dpr[sp]->mass_charge-spec[sp]->mass)>1e-8) ERROR(("mass of molecule has changed"))
    if (fabs(dpr[sp+nspec]->mass_charge-spec[sp]->charge)>1e-8) ERROR(("charge of molecule has changed"))
  }

  if (slab.sp<nspec) {
    if (slab.sp>=0) prt("density profile autocenter using CoM of species %d",slab.sp);
    else prt("density profile autocenter using CoM of all species 0..%d (incl.)",-slab.sp); }
}

void loaddpr(double dprgrid) /************************************** loaddpr */
{
  double gr;
  int sp,nsp;

  if (dprgrid==0) return;
  if (!dpr) ERROR(("internal"))

  VarOpen(Fn("dpr"),"r");
  VarRead(&gr,sizeof(gr));
  if (gr!=dprgrid)
    ERROR(("%s: inconsistent grid (%g in file, %g in data)",lastFn,gr,dprgrid))
  VarRead(&nsp,sizeof(nsp));
  if (nsp!=nspec)
    ERROR(("%s: inconsistent # of species (%d in file, %d in data)",lastFn,nsp,nspec))
  loop (sp,0,2*nspec+nsites) VarReadSds(dpr[sp]);
  VarClose();
}

void savedpr(double dprgrid) /************************************** savedpr */
{
  int sp;

  if (dprgrid==0) return;
  if (!dpr) ERROR(("internal"))

  backup("dpr");
  VarOpen(Fn("dpr"),"w");
  VarPut(&dprgrid,sizeof(dprgrid));
  VarPut(&nspec,sizeof(nspec));
  loop (sp,0,2*nspec+nsites) VarPutSds(dpr[sp]);
  VarClose();
}

void measuredpr(int slabmode) /********************************** measuredpr */
{
  int n,nbins,nsym=0;
  double Grid;
  vector sym={0,0,0}; /* for export of asymmetry to En.sym */
  int i,k;
  int geom3=slab.geom%3;
  static int check=1;

#ifdef FREEBC
  ERROR(("implementation limitation: cannot measure density profiles in free b.c."))
#endif /*# FREEBC */

  /* autocenter for slab.sp=big but not 0x7fffffff, default if no measuredpr */
  VO(box.center,=0)

    if (check)
    switch (geom3) {
      case 2: // slab
        if (box.L[2]<=box.L[0] || box.L[2]<=box.L[1] || fabs(log(box.L[1]/box.L[0]))>0.01) {
          check=0;
          WARNING(("Slab geometry requested, but the box is not tailored for\n\
*** a z-perpendicular slab. Check slab.geom or conditions Lx=Ly<Lz.")) }
        break;
      case 1: // trickle or tunnel in x-direction
        if (box.L[0]>=box.L[1] || box.L[0]>=box.L[2] || fabs(log(box.L[1]/box.L[2]))>0.01) {
          check=0;
          WARNING(("Trickle/tunnel requested, but the box is not tailored for\n\
*** a x-trickle or x-tunnel geometry. Check slab.geom or conditions Lx<=Ly=Lz.")) }
        break;
      case 0: // droplet or cavity
        if (fabs(log(box.L[2]/box.L[0]))>0.01 || fabs(log(box.L[1]/box.L[0]))>0.01) {
          check=0;
          WARNING(("Droplet/cavity requested, but the box is not tailored for\n\
*** them. Check slab.geom or conditions Lx=Ly=Lz.")) }
        break;
      default:
        ERROR(("slab.geom=%d is not supported",slab.geom)) }

  if (!dpr) return;

  /* NB: all histogram sizes dpr[*]->n and grids dpr[*]->grid are the same */
  nbins=abs(dpr[0]->n);
  Grid=dpr[0]->grid;
  if (dpr[0]->n<0)
    switch (geom3) {
      case 2: Grid/=box.L[2]; break;
      case 1: Grid/=sqrt(box.L[1]*box.L[2]); break;
      case 0: Grid/=cbrt(PROD(box.L)); break;
      default: ERROR(("slab.geom=%d out of range",slab.geom)) }

  if (geom3==2 && slab.max!=0 && slab.max<box.L[2]*0.99999999999) {
    static int pass=1;
    if (pass) {
      pass=0;
      ERROR(("slab.max=%g < Lz=%g not allowed\n\
because the z-profiles would overlap",slab.max,box.L[2]))
      /* in fact, slab.max<Lz causes index overflow */ } }

  /* autocenter density profiles (for measurements only) */
  /* NB: for slab.sp>0 but <0x7fffffff, box.center = 0 (see above) */
  if (slab.sp==0x7fffffff)
    VV(box.center,=box.Lh) /* this is the default = box center */
  else if (slab.sp<nspec) {
    /* autocentering */
    vector s={0,0,0},c={0,0,0};

    /* droplet/cylinder/slab autocenter, based on molecular center of mass */
    loop (n,0,No.N) {
      vector cm={0,0,0}; /* position of the center of mass */
      molecule_t *mn=molec+n;
      int sp=mn->sp;
      siteinfo_t *si=spec[sp]->si;
      int ns=mn->ns;
      vector *r=rof(mn,cfg[0]->rp);

      if (slab.sp<0) { if (sp>abs(slab.sp)) continue; }
      else { if (sp!=slab.sp) continue; }

      /* slab.geom: 0=droplet, 1=trickle(cylinder), 2=slab(film),
                    3=cavity,  4=tunnel,            5=slab (inverted) */
      loop (i,0,ns) loop (k,geom3,3) cm[k]+=si[i].mass*r[i][k];

      loop (k,geom3,3) {
        cm[k]*=2*PI/dpr[sp]->mass_charge/box.L[k];
        s[k]+=sin(cm[k]);
        c[k]+=cos(cm[k]); }
    } /* n */

    // old: En.shift[2]=(atan2(c[2],s[2])/(2*PI)+1.25)*box.L[2];
    loop (k,geom3,3) {
      box.center[k]=(atan2(s[k],c[k])/(2*PI)+(slab.geom/3)*0.5)*box.L[k];
      if (box.center[k]>=box.L[k]) box.center[k]-=box.L[k];
      if (box.center[k]<0) box.center[k]+=box.L[k]; }
  } /* autocentering */

  slab.torem=-1; /* no molecule to be removed */

  loop (n,0,No.N) {
    vector cm={0,0,0}; /* z-position of the center of mass */
    vector rv;
    double rc; /* was only z for <V3.5h */
    molecule_t *mn=molec+n;
    int sp=mn->sp;
    siteinfo_t *si=spec[sp]->si;
    int ns=mn->ns;
    vector *r=rof(mn,cfg[0]->rp);
    int nh;
#ifdef POLAR
    vector *rpol=(vector*)((char*)r+polar_off); /* relative to r, chargepol */
    /* equivalent: rpol=r+No.s; or rpol=polaroff(mn,cfg[0]->rp); */
#endif /*# POLAR */

    loop (i,0,ns) {

      loop (k,geom3,3) {
        rv[k]=r[i][k]-box.center[k];
        cm[k]+=si[i].mass*rv[k];
        while (rv[k]>=box.Lh[k]) rv[k]-=box.L[k];
        while (rv[k]<-box.Lh[k]) rv[k]+=box.L[k]; }

      /* site density profiles */
      switch (geom3) {
        case 2: /* z-profile */
          rc=rv[2];
          /* periodic */
          nh=(int)(rc*Grid+nbins)%nbins;
          break;
        case 1: /* rc=from cylinder/tunnel center */
          rc=sqrt(Sqr(rv[2])+Sqr(rv[1]));
          goto nh1;
        case 0: /* rc = from droplet center */
          rc=sqrt(SQR(rv));
        nh1:
          nh=rc*Grid;
          /* the last bin contains everything further */
          Min(nh,nbins-1) }

      dpr[2*nspec+si[i].st]->hist.ihist[nh]++;

      /* charge density profile (slab only) */
      if (geom3==2) {
#if COULOMB<-2
#  include "gcprof.c"
#else /*# COULOMB<-2 */
#  include "cprof.c"
#endif /*#!COULOMB<-2 */
      }
    } /* i */

    loop (k,geom3,3) {
      rv[k]=cm[k]/dpr[sp]->mass_charge; /* cm[k] has already centered coordinates */
      while (rv[k]>=box.Lh[k]) rv[k]-=box.L[k];
      while (rv[k]<-box.Lh[k]) rv[k]+=box.L[k]; }

    if ( (slab.sym>=0 && sp==slab.sym) || (slab.sym<0 && sp<=abs(slab.sym)) ) {
      loop (k,geom3,3) sym[k]+=rv[k];
      nsym++; }

    /* molecular density profiles and detection of evaporating molecules */
    switch (geom3) {
      case 2: // slab, perpendicular to z
        rc=rv[2];
        /* periodic */
        nh=(int)(rc*Grid+nbins)%nbins;
        if (slab.out && fabs(rc)>slab.outx) slab.torem=n;
        break;
      case 1: // trickle (cylinder) or tunnel
        rc=sqrt(Sqr(rv[2])+Sqr(rv[1]));
        goto nh2;
      case 0: // droplet or spherical cavity
        rc=sqrt(SQR(rv));
      nh2:
        if (slab.out && rc>slab.outx) slab.torem=n;
        nh=rc*Grid;
        /* the last bin contains everything further */
        Min(nh,nbins-1) }

    dpr[sp]->hist.dhist[nh]+=dpr[sp]->mass_charge;
  } /* n */

  if (slab.torem>=0)
    fprintf(stderr,"molecule %d will be removed (slab.outx=%g)\n",slab.torem,slab.outx);

  slab.outx=slab.outx*0.75+slab.out*0.25;

  /* V,L statistics for dpr[0] only (is the same for all) */
  switch (geom3) {
    case 0:
      dpr[0]->Vbin+=1/Cub(Grid); // to be multiplied by 4*PI/3*(2*i*i+3*i+1)
      break;
    case 1:
      dpr[0]->Vbin+=box.L[0]/Sqr(Grid); break; // to be multiplied by PI*(2*i+1)
      break;
    case 2:
      dpr[0]->Vbin+=box.L[0]*box.L[1]/Grid;
      break; }

  dpr[0]->Lz+=box.L[2];
  dpr[0]->nmeas++;
  loop (k,geom3,3) En.sym[k]=fmod(sym[k]/nsym+box.L[k],box.L[k])-box.Lh[k];
}

void printdpr(int slabmode,int slabprt,double dV,char *PdVname) /* printdpr */
/*
  Print density profiles and report surface tension results
  slabmode = slab.mode&4:
    1 calculate surface tension and post-corrections
    2 test area version
    4 extended slab (bit), with slabprt&014 only
    8 (POLAR not Gauss): too short dipoles for charge z-profile made longer,
                         equivalent to triangle histogram bin
  slabprt = slab.prt&014:
    1 =  01 = .cm.AA-3.z
    2 =  02 = .cm.kgm-3.z
    4 =  04 = .site.AA-3.z
    8 = 010 = .site.kgm-3.z
   16 = 020 = .q.eAA-3.z
  dV, test area parameter, see slab.mode=2 and the manual
  PdVname = constrd.PdVname
    = rescale&RESCALE_CM ? "PdVmol [Pa]" : "PdVatom [Pa]"
*/
{
  int sp,i,nbins;
  FILE *f;
  double s,satom;
  int number; /* 0: use [kg.m-3] or [e.AA-3], 1: use [AA-3] */
  char *fmt[]={" %10.4f", " %10.7f", " %11.8f"};
            /*  [kg.m-3]   [AA-3]    [e.AA-3] */
  double DZ,Lzf;
  int ext=slabmode&4;
  int geom3=slab.geom%3;
  char *geom=geom3==2?"SLAB":geom3==1?"TRICKLE":"DROPLET";

  if (ext) ext=1; /* ext =0,1 expected because subscript of sumihist */
  extext=ext; /* static (external in these functions) */

  if (!dpr) return;

  slabmode&=3;
  dprz=NULL; /* first nonzero dpr[]: normally dprz=dpr[0]
                if called from extenddpr then dprz=dpr[nspec*2] */
  loop (sp,0,nspec*2+nsites) if (dpr[sp]) { dprz=dpr[sp]; break; }
  if (!dprz) ERROR(("internal: no data in dpr"))

  nbins=abs(dprz->n);

  if (dprz->n>0) DZ=1./dprz->grid;
  else DZ=dprz->Lz/(dprz->nmeas*dprz->grid);

  loop (number,0,2) {

    if (slabprt&"\2\1"[number]) {
      /*** molecules: center-of-mass density profile ***/
      f=fopen(Fn(number?"cm.AA-3.z":"cm.kgm-3.z"),"wt");
      fprintf(f,"# %s DENSITY PROFILES BY CENTER-OF-MASS OF MOLECULES\n\
# calculated from %s\n\
# z|r is in AA, density profiles are in %s\n",
              geom,
              Fn("dpr"),number?"AA-3":"kg.m-3");
      if (dprz->n<0)
        fprintf(f,"# accumulated in [0,Lz], using <Lz>=%.9g, dz=%.9g\n",dprz->Lz/dprz->nmeas,DZ);

      fprintf(f,"#   z|r  ");
      loop (sp,0,nspec) fprintf(f,"  %c=%-7s",sp+'B',dpr[sp]->name);
      fprintf(f,"    SUM\n");

      loop (i,0,nbins) {
        fprintf(f,"%6.3f ",(i+0.5)*DZ);

        s=0;
        loop (sp,0,nspec) {
          double x=dpr[sp]->hist.dhist[i]/dprz->Vbin;

          if (number) x/=dpr[sp]->mass_charge;
          else x*=rhounit;

          if (geom3==1) x/=PI*(2*i+1);
          if (geom3==0) x/=4*PI/3*(3*i*(i+1)+1);

          fprintf(f,fmt[number],x);
          s+=x; }
        fprintf(f,fmt[number],s);
        fprintf(f,"\n"); }

      fclose(f); }

    if (slabprt&"\10\4"[number]) {
      /*** sites-based density profile, also for ext ***/
      if (ext) f=fopen(Fn(number?"site.AA-3.extz":"site.kgm-3.extz"),"wt");
      else f=fopen(Fn(number?"site.AA-3.z":"site.kgm-3.z"),"wt");
      fprintf(f,"# %s DENSITY PROFILES BY SITES\n\
# calculated from %s\n\
# SUM(all)=sum over all sites, SUM(atoms)=only sites with nonzero mass\n\
# z|r is in AA, density profiles are in %s\n",
              geom,Fn("dpr"),number?"AA-3":"kg.m-3");
      if (dprz->n<0) fprintf(f,"# accumulated in [0,Lz], using <Lz>=%.9g, dz=%.9g\n",dprz->Lz/dprz->nmeas,DZ);

      fprintf(f,"#    z   ");
      loop (sp,nspec*2,nspec*2+nsites)
        fprintf(f,"  %c=%-7s",sp-2*nspec+'B',dpr[sp]->name);
      fprintf(f,"  SUM(all) SUM(atoms)\n");

      loop (sp,0,nsites) {
        sumihist[sp][ext]=0;
        loop (i,0,nbins) sumihist[sp][ext]+=dpr[2*nspec+sp]->hist.ihist[i]; }

      loop (i,0,nbins) {
        fprintf(f,"%6.3f ",(i+0.5)*DZ);

        satom=s=0;
        loop (sp,2*nspec,2*nspec+nsites) {
          double x=dpr[sp]->hist.ihist[i]/dprz->Vbin;
          double m=sitedef[sp-2*nspec].M/Munit;

          if (!number) x*=m*rhounit;

          if (geom3==1) x/=PI*(2*i+1);
          if (geom3==0) x/=4*PI/3*(3*i*(i+1)+1);

          fprintf(f,fmt[number],x);
          s+=x;
          if (m) satom+=x; }
        fprintf(f,fmt[number],s);
        fprintf(f,fmt[number],satom);
        fprintf(f,"\n"); }

      fclose(f); } }

  if (slabprt&16) {

    number=2; /* should be =2 after the above loop .. */

    if (geom3==2) {
      /** charge density profile (for slab only) */
      f=fopen(Fn("q.eAA-3.z"),"wt");
      fprintf(f,"# z-CHARGE DENSITY PROFILE BY MOLECULES\n\
# calculated from %s\n\
# z is in AA, density profile is in e.AA-3\n",Fn("dpr"));
    if (dprz->n<0) fprintf(f,"# accumulated in [0,Lz], using <Lz>=%.9g, dz=%.9g\n",dprz->Lz/dprz->nmeas,DZ);

    fprintf(f,"#    z   ");
    loop (sp,nspec,2*nspec) fprintf(f,"  %c=%-7s",sp-nspec+'B',dpr[sp]->name);
    fprintf(f,"    SUM\n");

    loop (i,0,nbins) {
      fprintf(f,"%6.3f ",(i+0.5)*DZ);
      s=0;
      loop (sp,nspec,2*nspec) {
        double x=dpr[sp]->hist.dhist[i]/dprz->Vbin;

        x/=electron*2; // see gcprof.c
        fprintf(f,fmt[2],x);
        s+=x; }
      fprintf(f,fmt[number],s);
      fprintf(f,"\n"); }

    fclose(f); } }

  if (geom3<2) return;

  Lzf=0.75e-10*dprz->Lz/dprz->nmeas;

  /* surface tension */
  if (slabmode&1) {
    double c[3]={slabcorr(0),slabcorr(1),slabcorr(2)};
    double Esta,Psta=StaMean("P(vir) [Pa]"),Ptsta=StaMean("Pt [Pa]");
    double sc[3],Pc,dP=StaStdErr("P(vir) [Pa]")*Lzf;
    double V=box.L[0]*box.L[1]*(dprz->Lz/dprz->nmeas);

    Esta=StaMean("Epot [J/mol]");
    if (StaError) Esta=StaMean("Epot [kcal/mol]")*kcal; /* not SI deprecated */

    underline(ext?"EXTENDED SLABS":"surface tension");

    if (ext) prt("\
Stacking error corrected using extended slabs (see files *.extz).\n\
The z-density profiles have been extended using the following parameters:\n\
  slab.sp=%d = %s\n\
  phase of lower conc. of slab.sp extended by slab.ext.zero=%d bins\n\
  phase of higher conc. of slab.sp extended by slab.ext.center=%d bins\n\
  slab.ext.span=%d bins used for extension\n\
%s\n\
The final stacking correction to be added\n\
  = the (*) value minus the previous (*) value",
                 slab.sp,slab.sp>=0?
                 "the species used for slab centering":
                 "species 0 .. |slab.sp| (incl.) used for slab centering",
                 slab.ext.zero,slab.ext.center,slab.ext.span,
                 option('v')&4?"VERBOSE FLAG -v4 USED, MANY NOT MEANINGFUL DATA ARE PRINTED BELOW":"(obtain technical info by flag -v4)");

    sc[0]=(c[0]-En.corr/V)*Eunit;
    sc[1]=(c[1]/V-En.corr/(V*V))*Punit;
    sc[2]=(c[2]/V-En.corr/(V*V))*Punit;

    if (slab.sp>=nspec) prt("WARNING: no autocentering, see the manual for slab.sp");

    if (!extext | option('v')&4) {
      prt("\
This module calculates the slab cutoff corrections to Epot, P, and\n\
Pt=(Pxx+Pyy-2*Pzz)/3 from accumulated z-density profiles (post-processing).\n\
WARNING: %s\n\
raw slab corrections: c[0]=%g  c[1]=%g  c[2]=%g\n\
raw homogeneous fluid: corr=%g  corr/V=%g  corr/V^2=%g",
          slab.K
          ?"The preferred Fourier-expansion-based cutoff corrections have been\n\
         already included. The following corrections are informative only and\n\
         should not be added."
          :"This may be less accurate, consider Fourier-expansion-based cutoff\n\
         corrections instead!",
          c[0],c[1],c[2],  En.corr,En.corr/V,En.corr/(V*V));
      if (En.corr) prt_("\
WARNING: The homogeneous fluid corrections have been by mistake (see corr)\n\
included in the calculations. To fix, they will be subtracted here and the\n\
slab-based corrections will be added instead.");

      header("                               Epot [J/mol]        P [Pa]         Pt [Pa]  ");
      prt("homogeneous fluid correction %14.9g  %14.9g  %14.9g",
          En.corr/V*Eunit,En.corr/(V*V)*Punit,En.corr/(V*V)*Punit);
      prt("uncorrected values           %14.9g  %14.9g  %14.9g",
          Esta-En.corr/V*Eunit,Psta-En.corr/(V*V)*Punit,Ptsta-En.corr/(V*V)*Punit);
      prt("slab-based cutoff correction %14.9g  %14.9g  %14.9g",
          c[0]*Eunit,c[1]*Punit/V,c[2]*Punit/V);
      prt("correction to be added       %14.9g  %14.9g  %14.9g",
          sc[0],sc[1],sc[2]);
      prt("final corrected value        %14.9g  %14.9g  %14.9g",
          Esta+sc[0],Psta+sc[1],Ptsta+sc[2]);
      Pc=Psta+sc[2]; /* for virial... */
      if (En.corr) prt("corr_slab/corr_homog.fluid   %12.5f %15.5f %15.5f\n\
thick slab limits                    1               1              -",
                       c[0]*V/En.corr,c[1]*V/En.corr,c[2]*V/En.corr);
      header(""); }

    prt("\n\
Surface tension calculation:\n\
gamma/[N/m] = -Pt/[Pa]*%g = -P/[Pa]*%g (for NPT in x,y)",Lzf,Lzf/1.5);
    prt("Post-processed z-profile-based correction (%s): %g N/m (*)",
        slab.K?"to be ignored":"as added below",surftenscorr[extext]=-c[2]/V*Punit*Lzf);
#if (COULOMB<=-1) && (PRESSURETENSOR&PT_VIR)
    prt("Yeh-Berkowitz correction to be added = %g N/m",En.YB[2]*Lzf/1.5);
#endif /*# (COULOMB<=-1) && (PRESSURETENSOR&PT_VIR)  */
    if (En.corr) {
      prt("wrong homogeneous correction:                      %13.8f N/m",
          En.corr/(V*V)*Punit*Lzf);
      prt("fix (to be added to homogeneously-corrected value):%13.8f N/m",
          -sc[2]*Lzf); }

    if (!extext | (option('v')&4)) {
      if (slab.K)
        prt("\n\
Fourier expansion cutoff corrections for the slab geometry (slab.K=%d):\n\
  - use the results marked \"FF-corrected\"\n\
  - the \"DO_NOT_USE\" column is doubly-corrected - ignore it",slab.K);

      header(slab.K?
             " surface tension gamma in N/m: DO_NOT_USE    stderr    FF-corrected":
             " surface tension gamma in N/m: corrected     stderr    uncorrected");
      prt("from P(vir) assuming Psat=0 %12.7f %12.7f %12.7f (a,b)",
          -Pc*Lzf, dP, -(Psta-En.corr/(V*V)*Punit)*Lzf);
      if (dV) {
        char shortname[8];

        memcpy(shortname,PdVname,7); shortname[7]=0;
        prt("%s%-7s) %12.7f %12.7f %12.7f%s",
            slabmode&2?"test-area (tangent ":
            "virtual V change  (",
            shortname,
            -(StaMean(PdVname)+sc[2])*Lzf,
            StaStdErr(PdVname)*Lzf,
            -(StaMean(PdVname)-En.corr/(V*V)*Punit)*Lzf,
            slabmode&2?"":" (a)"); }
#if PRESSURETENSOR&PT_VIR
      prt("from virial part of Pt      %12.7f %12.7f %12.7f (c)",
          -(StaMean("Pvir xy-zz [Pa]")+c[2]/V*Punit)*Lzf,
          StaStdErr("Pvir xy-zz [Pa]")*Lzf,
          -StaMean("Pvir xy-zz [Pa]")*Lzf);
#  if PRESSURETENSOR&PT_KIN
      prt("from full Pt (incl.kin.part)%12.7f %12.7f %12.7f (d)",
          -(StaMean("Pt [Pa]")+c[2]/V*Punit)*Lzf,
          StaStdErr("Pt [Pa]")*Lzf,
          -StaMean("Pt [Pa]")*Lzf);
#  endif /*# PRESSURETENSOR&PT_KIN */
#endif /*# PRESSURETENSOR&PT_VIR */
      header("");
      prt("(a) Assumes <Pzz>=0: subtract %g*(saturated pressure in Pa)",Lzf);
      prt("(b) Assumes elst.virial = -elst.energy, "
#if COULOMB<0
          "cf. =P(tens)-P(vir) [Pa]= for accuracy"
#else /*# COULOMB<0 */
          "not good with cut-off electrostatics"
#endif /*#!COULOMB<0 */
                                                   );
      prt("(c) Valid only for monoatomic models or vibrating bonds (cook -u9999),\n\
    incorrect for constrained bonds or rigid molecules");
#if (PRESSURETENSOR&PT_VIR) && !(PRESSURETENSOR&PT_KIN)
      prt("    (cook* compiled with PRESSURETENSOR&%d needed for this case)",PT_KIN);
#endif /*# (PRESSURETENSOR&PT_VIR) && !(PRESSURETENSOR&PT_KIN) */
      prt("(d) Preferred value (if available)");
      _n }

    if (extext)
      prt("\nThe final stacking correction: %.12g N/m\n\
  (to be added to the results of table \"surface tension\")\n",
          surftenscorr[1]-surftenscorr[0]); }

#if PRESSURETENSOR&PT_ANY
  else
    prt("surface tension = %10.7f %10.7f [N/m]",
                          -StaMean("Pt [Pa]")*Lzf,
                                  StaStdErr("Pt [Pa]")*Lzf);
#endif /*# PRESSURETENSOR&PT_ANY */
}

static
double slabpaircorr(int is,int js,int isP) /******************* slabpaircorr */
/***
  Cutoff corection calculation in the slab geometry from z-profiles for one
  pair of atom types:
    isP=0: energy correction returned,
    isP=1: isotropic pressure correction returned (not needed),
    isP=2: pressure-equivalent correction for surface tension returned.
  Cf. setss() in setss.c;
    Unlike setss(), slabpaircorr should be called after the density profile
    has been calculated.
  Called (for pairs) from slabcorr().
  V2.6c: extended to sum of periodic images in z upto 2*Lz (?)
***/
{
#ifdef NIBC
#  error SLAB + NIBC not supported
#endif /*# NIBC */

#ifdef SHARPCUTOFF
  return 0;
#else /*# SHARPCUTOFF */

  DECLARE_PARMS

  sitesite_t *ss=&sstab[is][js];
  int i,j,ig,n,nn,nbins,NBINS;
  double sumi,sumj,h,dz,xi,sum,x,rr,y,U=0,f,Urep,frep,z,xmax,hh,DZ;
  double *Iabsz;

  if (!dpr) return 0;
  nbins=abs(dprz->n);
  NBINS=nbins*3; /* will integrate over z-periodic images */

  if (ss->A==0 || ss->C1q==0) return 0;

  allocarrayzero(Iabsz,NBINS); /* to Lz*2 */

  if (dpr[is])
    if (dpr[is]->n!=dprz->n || dpr[js]->n!=dprz->n)
      ERROR(("z-profiles of different sizes not supported"))

  /* r-integrals (over Dz-distant slabs) in advance (because time-consuming):
    isP=0: Iabsz(Dz) = INT_0^infty 2 PI r dr Du(sqrt(r^2+Dz^2))
    isP=1: Iabsz(Dz) = INT_0^infty 2 PI r dr r*f_ij(sqrt(r^2+Dz^2))/3
    isP=2: Iabsz(Dz) = INT_0^infty 2 PI r dr (r^2-2z^2)/r*f(sqrt(r^2+Dz^2))/3
    method: Gauss 2 point formula (error ~ h^4)
  */
#  define NGauss4 64
#  define Gauss4 0.2113248654051871

  hh=1/sqrt(ss->C1q)/NGauss4;

  if (dprz->n>0) DZ=1./dprz->grid;
  else DZ=dprz->Lz/(dprz->nmeas*dprz->grid);

  //  put2(dprz->n,DZ) put3(dprz->Lz,dprz->nmeas,dprz->grid)

  loop (i,0,NBINS) {
    dz=i*DZ;
    xmax=1/fmax(sqrt(ss->C1q),dz);
    nn=(int)(xmax/hh+1.99);
    h=xmax/nn;
    sum=0;

    loop (n,0,nn) loop (ig,0,2) {
      if (ig) xi=(n+1-Gauss4)*h;
      else xi=(n+Gauss4)*h;
      rr=1/Sqr(xi);
      /* cf. LJM in internp.c */
      U=0; /* NB: we have U+= in SS_MEASURE but f= */
      SS_MEASURE_rep
      SS_MEASURE
      if (rr < ss->C2q) {
        x=rr-ss->C2q; U-=x*x*ss->A; f-=ss->A4*x; }
      switch (isP) {
        case 0: sum+=U*rr/xi; break;
        case 1: sum+=Sqr(rr)/xi*f/3; break;
        case 2: sum+=rr*(rr-3*Sqr(dz))/xi*f/3; } }
    Iabsz[i]=sum*PI*h; }

  /* optimize by calculating in advance ... */
  sumi=0; loop (i,0,nbins) sumi+=dpr[2*nspec+is]->hist.ihist[i];
  sumj=0; loop (j,0,nbins) sumj+=dpr[2*nspec+js]->hist.ihist[j];

  sum=0;
  loop (i,0,nbins)
    loop (j,0,nbins) {
      int k;
      double s=0;

      for (k=i-j-NBINS; k<NBINS; k+=nbins) if (abs(k)<NBINS) s+=Iabsz[abs(k)];

      if (extext) s*=(double)sumihist[is][1]*(double)(sumihist[js][1])
                   /((double)sumihist[is][0]*(double)(sumihist[js][0]));

      sum+=(double)dpr[2*nspec+is]->hist.ihist[i]
          *(double)dpr[2*nspec+js]->hist.ihist[j]*s; }

  // loop (i,0,NBINS) prt("%g %g Iabsz%d",(i+0.5)*DZ,Iabsz[i],isP);

  free(Iabsz);

  if (sumi*sumj==0) {
    if (sum!=0) ERROR(("slabpaircorr: wrong sums"))
    return 0; }

  return sum/(sumi*sumj);
#endif /*#!SHARPCUTOFF */
} /* slabpaircorr */

double slabcorr(int isP) /**************************************** slabcorr */
/***
  Counterpart of cutcorr() in siminit.c for slab geometry, obtained by
  integrating the z-density profile:
    isP=0: energy correction returned,
    isP=1: V*isotropic pressure correction returned (not needed),
    isP=2: V*pressure-equivalent correction for surface tension returned.
  However, this function should be called at the end of the simulation run!!!
  NOTE: gamma = -3/4*Lz*Pt, Pt=(Pxx+Pyy-Pzz*2)
    for isP=1,2: V*Delta P = -4/3* Lx Ly Delta gamma is returned
***/
{
  int n,m,i,j,np;
  double c,sumcorr=0;

#ifdef SHARPCUTOFF
  WARNING(("SLAB and SHARPCUTOFF: no slab cutoff correction"))
  return 0;
#endif /*# SHARPCUTOFF */

  loop (n,0,nspec)
    loopto (m,0,n) {
      if (m==n) np=spec[n]->N*(spec[n]->N-1)/2;
      else np=spec[n]->N*spec[m]->N;
      c=0;
      loop (i,0,spec[n]->ns) loop (j,0,spec[m]->ns)
        c += slabpaircorr(spec[n]->si[i].st,spec[m]->si[j].st,isP);
      sumcorr+=c*np; }

  return sumcorr/(box.L[0]*box.L[1]);
  /* pressure correction (isP=1,2) should be divided by V */
} /* slabcorr */

/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       NEW: Slab stacking corrected
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

static struct dpr_s **dpr0; /* original dpr stored here */

void extenddpr(int slabmode,int slabprt,double dV,char *PdVname,
               int nzero,int ncenter,int nspan) /***** extenddpr */
/*
  See printdpr!
  Example of bin replication scheme for nbins=26, nspan=2, nzero=4, ncenter=6
    (spaces for readability only)
       ABCDEFGHIJKLM        NOPQRSTUVWXYZ ->
  ZABY ABCDEFGHIJKLM MNOLMN NOPQRSTUVWXYZ
*/
{
  int n,nbins0,nbins,sp,ispan,i;

  if (!dpr) return;

  nbins=abs(dpr[0]->n); /* all sizes assumed the same */

  if (nbins&1) WARNING(("odd number of bins is not recommended, but accepted"))

  if (nspan<0 || nspan>nbins/6) {
    WARNING(("extenddpr: too long nspan"))
    return; }

  if (slab.sp>=nspec) {
    WARNING(("cannot extend the slab because not autocentered"))
    return; }

  nbins0=nbins;
  nbins=nbins+nzero+ncenter;
  n=sign(dpr[0]->n)*nbins;
  /* old dpr stored as dpr0 */
  dpr0=dpr;
  rallocarrayzero(dpr,2*nspec+nsites);

  loop (sp,2*nspec,2*nspec+nsites) {
    /* only site z-profiles allocated, code kept general as from sp=0 */
    int elemsize=
      sp<2*nspec ? sizeof(dpr[sp]->hist.dhist[0])
                 : sizeof(dpr[sp]->hist.ihist[0]);
    double q=(double)nbins/nbins0;

    sdsralloczero(dpr[sp],sizeof(dpr[sp][0])+elemsize*(abs(nbins)-1));
    dpr[sp]->n=n;
    dpr[sp]->grid=dpr0[0]->grid;
    if (dpr0[0]->n<0) dpr[sp]->grid*=q;
    dpr[sp]->nmeas=dpr0[0]->nmeas;
    dpr[sp]->Lz=dpr0[0]->Lz*q;
    dpr[sp]->Vbin=dpr0[0]->Vbin;
    dpr[sp]->mass_charge=dpr0[0]->mass_charge; /* not needed */
    if (sp<nspec) sprintf(dpr[sp]->name,"mol%d",sp);
    else if (sp<2*nspec) sprintf(dpr[sp]->name,"mol%d",sp-nspec);
    else strcpy(dpr[sp]->name,sitedef[sp-2*nspec].name);

    ispan=-nspan/2;
    loop (i,0,nzero) {
      /* ispan=[-nspan,nspan) */
      dpr[sp]->hist.ihist[i]=dpr0[sp]->hist.ihist[(nbins0+ispan)%nbins0];
      ispan++;
      if (ispan==nspan) ispan=-nspan; }

    loop (i,0,nbins0/2)
      dpr[sp]->hist.ihist[i+nzero]=dpr0[sp]->hist.ihist[i];

    loop (i,0,ncenter) {
      dpr[sp]->hist.ihist[i+nbins0/2+nzero]=dpr0[sp]->hist.ihist[nbins0/2+ispan];
      ispan++;
      if (ispan==nspan) ispan=-nspan; }

    loop (i,nbins0/2,nbins0)
      dpr[sp]->hist.ihist[i+nzero+ncenter]=dpr0[sp]->hist.ihist[i];
  } /* sp */

  printdpr(slabmode|4,slabprt&014,dV,PdVname);

  release(dpr);
  dpr=dpr0;
} /* extenddpr */
