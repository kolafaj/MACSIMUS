#include "ground.h"
#include "sds.h"
#include "simglob.h"
#include "simdef.h"
#include "interpot.h"
#ifdef COULOMB
#  include "elst.h"
#endif /*# COULOMB */
#include "setss.h"
#include "simils.h" /* brrr: only Fn while DIHHIST==-1 needed...*/

/* optimized water models: WATER passed from lj/sitesite.h */
#ifdef LINKCELL
#  undef WATER
#endif /*# LINKCELL */
#ifdef POLAR
#  undef WATER
#endif /*# POLAR */
#ifdef WATER
#  include "water.h"
#endif /*# WATER */

#include <time.h>
#include "simmeas.h"
#include "units.h"

#define TETRAHEDRAL 1.9106332362490184

/* global variables */

int nsites,nspec,ngroups,nnbfixes;
int ns; /* see checkatomn */
int n14=2; /* 1 if 1-4 are missing, 2 if 1-4 are present */

#ifdef SLAB
/*                     rho  z[0] z[1] */
struct wall_s wall = { 19320, {0,1} };
struct wsstab_s *wsstab;
#endif /*# SLAB */

/* the energy scaling factor in the ble-file, in p.u. = k_B * K
   added to the ble-file in V2.8a */
double eunit=kcal/Eunit; /* legacy default */

sitedef_t *sitedef;
vector **initsites;
nbfix_t *nbfix0;
double *sigvdWptr; /* =fixij(tau.i,tau.j)->onefour[0].sig */

nbfix_t *fixij(int i,int j) /**************************************** nbfixi */
/*
   returns pointer to the (nbfix structure corresponding to) sites i,j
   in the nbfixes table
*/
{
  nbfix_t *fix;

  for (fix=nbfix0; fix; fix=fix->next)
    if ( (fix->indx[0]==i && fix->indx[1]==j)
     ||  (fix->indx[0]==j && fix->indx[1]==i) )
      return fix;

  return NULL;
}

static char *line,*copyofline;

static char *realname(char *name)
/* strip leading ? or lowercase group name */
{
  if (name[0]=='?') return name+1; /* water recognizer support */
  if (islower(name[0]) || isdigit(name[0])) return name+1; /* rdf group support */

  return name;
}

int atomn(char *name,int spg,int ig) /******************************** atomn */
/***
    returns atom number if the name is given
    not very efficient implementation
    if name starts with ?, does not report error but returns -1
NEW:
    if ig>=0, checks if the atom (species spg, site ig) is a member of a group
    1st letter is the group ID, =1 char lowercase
***/
{
  int i,iempty,igroup,j;
  char *wname=realname(name);

  iempty=-1;
  loop (i,0,nsites) if (!strcmp(wname,realname(sitedef[i].name)) && !sitedef[i].g) {
    if (iempty>=0)
      ERROR(("duplicate default rdf group - check %s.s-s",simils.simname))
    iempty=i; }

  igroup=-1;

  if (ig>=0) loop (i,0,nsites) if (!strcmp(wname,realname(sitedef[i].name)) && sitedef[i].g)
    loop (j,0,sitedef[i].g->n) if (sitedef[i].g->stitem[j].sp==spg && sitedef[i].g->stitem[j].i==ig) {
      if (igroup>=0)
        ERROR(("duplicate atom in a rdf group - check %s.s-s",simils.simname))
      igroup=i; }

#if 0
  if (igroup>=0) prt("~g name=%s spg=%d ig=%d atomn=%d",name,spg,ig,igroup);
  else if (iempty>=0) prt("~e name=%s spg=%d ig=%d atomn=%d",name,spg,ig,iempty);
  sleep(1);
#endif /*# 0 */

  if (igroup>=0) return igroup;
  else if (iempty>=0) return iempty;

  if (name[0]!='?')
    ERROR(("unknown type %s (real=%s) of atom (site) %d",name,wname,ig))

  return -1;
}

static specinfo_t *specsp;
int isss=0;

static void checkatomn(char *name,int i) /*********************** checkatomn */
{
  char *nm;

  if (i<0 || i>=ns) ERROR(("atom number %d of %s out of range [0,%d)",i,name,nsites))
  nm=sitedef[specsp->si[i].st].name;

  if (strcmp(nm,name)) {
    if (isss) isss++;
    else WARNING(("wrong type %s of atom %d replaced by %s in the following line\n%s",
                  name,i,nm,copyofline)) }
}

static char *tokenline;

static char *mygetline(void) /**************************************** getline */
/***
    gets one `line' from file `in', skipping comments (lines beginning with !)
    and warning on ***
    NOTE: WARNINGs,ERRORs during reading getdata are not tested here!
***/
{
  char *c;
  int l,warning;
  static int warned;

  line[0]=0;
  tokenline=line;
  copyofline[0]=0;

  do {
    c=fgets(line,LINELEN,in);
    if (!c) return c;
    /* this is to allow DOS ble-files to be read by UNIX */
    while ( (l=strlen(line)) && line[l-1]<=' ') line[l-1]=0;
    memcpy(copyofline,line,LINELEN);
    if (option('v')&4) prt("B %s",line);
    warning=!memcmp(line,"***",3);
    if (warning) {
      prt("%s",line);
      if (!warned++) WARNING(("=== the blend-file contains warnings or errors\n\
*** === (more warnings suppressed)")) }
  } while (warning || line[0]=='!' );

  if (strstr(line,"WARNING") || strstr(line,"ERROR"))
    WARNING(("the ble-file contains ERROR or WARNING produced by blend"))

  return c;
}

char *foundtoken;

static int checkmodel(const char *wmodel)
{
  return !strcmp(foundtoken,wmodel);
}

static void find(char *st) /******************************************* find */
/***
    Searches for line with given keyword in file `in'
    then skips comment lines (beginning with !) and reads the first
    non-comment line into `line'
***/
{
  char *c;

  do {
    if (!mygetline()) goto err;
    c=strtok(line," \n\r\t");
  } while (!c || strcmp(c,st));

  if ( (c=strtok(NULL," \n\r\t")) ) strcpy(foundtoken,c);
  /*.....fprintf(stderr,"reading %s\n",st);*/

  if (mygetline()) return;

 err: ERROR(("table %s not found\n\
*** check blank lines after tables\n\
*** note that config needs to be ended by 2 blank lines",st))
}

int checktok=1;

static char *tok(int line) /******************************************** tok */
{
  char *tok=strtok(tokenline," \n\r\t");

  tokenline=NULL;
  if (!tok) {
    if (checktok)
      ERROR(("called from %s:%d\n\
*** missing data in the following line:\n%s",__FILE__,line,copyofline))
    tok="-2146495994"; }

  return tok;
}

static void readtorsion(char *token,int ntorsions,int sp,torsion_t **head)
                                                            /*** readtorsion */
{
  int i;
  torsion_t *torsion;
#ifdef DIHHIST
  char type[4][6];
#endif /*# DIHHIST */
  int indx[4];
  int n,j;
#ifdef DIHHIST
#  if DIHHIST==-1
  int count=1; /* for dihedral number of .ddh, numbered from 1 */
  int *ddh=NULL,nddh,NDDH,pass,nn,nd;
  char line[128];

  if (token[0]=='d' && sp==0) {
    FILE *f=fopen(Fn("ddh"),"rt");
    if (!f)
      prt("WARNING no %s file, all dihedrals included",lastFn);
    else {
      loop (pass,0,2) {
        nddh=0;
        while (fgets(line,128,f) && line[0]!='!' && strlen(line)>2) {
          if (sscanf(line,"%d%d",&nn,&nd)!=2) ERROR(("%s: : missing or bad format",lastFn))
          if (nn) {
            if (pass) {
              ddh[nddh]=nd;
              if (nd<=0 || nd>ntorsions) ERROR(("%s: dihedral %d not in expected [1,%d]",lastFn,nd,ntorsions))
              if (nddh>=NDDH) ERROR(("%s: file changed while reading?",lastFn)) }
            nddh++; } }
        if (!pass) {
          allocarray(ddh,NDDH=nddh);
          rewind(f); } }
      fclose (f); } }
#  endif /*# DIHHIST==-1 */
#endif /*# DIHHIST */

  find(token);
  loop (i,0,ntorsions) {

    loop (j,0,4) {
      indx[j]=atoi(tok(__LINE__));
#ifdef DIHHIST
      memset(type[j],0,6);  /* set all uninitialized bytes */
      strcpy(type[j],tok(__LINE__)); /* store site type */
      checkatomn(type[j],indx[j]);
#else /*# DIHHIST */
      checkatomn(tok(__LINE__),indx[j]); /* check and skip atom name */
#endif /*#!DIHHIST */
      }
    tok(__LINE__); /* skip nX */

    /* periodicity */
    n=atoi(tok(__LINE__));
    if (abs(n)>6) { ERROR(("periodicity n=%d too high",n)) n=6; }

    alloc(torsion,sizeof(torsion_t)+sizeof(ireal)*((n==0)-(n<0)*n));

    copy(torsion->indx,indx,sizeof(indx));
    torsion->parm.n=n;

    torsion->next=*head; *head=torsion;

    torsion->parm.K[0]=eunit*atof(tok(__LINE__));
    tok(__LINE__); tok(__LINE__);
    if (n==0)
      torsion->parm.K[1]=PI/180*atof(tok(__LINE__)); /* angle conversion */
    else if (n<0)
      loopto (j,1,abs(n)) torsion->parm.K[j]=eunit*atof(tok(__LINE__));

#ifdef DIHHIST
    torsion->parm.dihhist=NULL;
    if (token[0]=='d' && sp==0)
#  if DIHHIST==-1
      if (ddh) {
        loop (nddh,0,NDDH) if (count==ddh[nddh])
          includedihhist(torsion,type,sp); }
      else
#  endif /*# DIHHIST==-1 */
        includedihhist(torsion,type,sp);
#  if DIHHIST==-1
    count++;
#  endif /*# DIHHIST==-1 */
#endif /*# DIHHIST */
    mygetline(); }

#ifdef DIHHIST
#  if DIHHIST==-1
  if (ddh) free(ddh);
#  endif /*# DIHHIST==-1 */
#endif /*# DIHHIST */
}

static void zerofix(real *x,char *info) /*************************** zerofix */
{
  double sum=SUM(x);

  if (fabs(sum)>1e-5) WARNING(("sum of coefficients %s = %g (should be 0) - fixed",sum))
  sum/=3;
  VO(x,-=sum)
}

static void readdependants(int ndependants,depend_t **head)
                                                       /***** readdependants */
/*
  dependants table line format (example)
  !site  ndep dep1 w   dep2 w  dep3  w
  [M] 2 M     3   1 O     0.73612449  0 H     0.13193776  3 H     0.13193776  err=2e-17
  L  (as above without err) + X[3], Y[3], wz, fx[3], fy[3]
  R  (as above without err) + l

  the err=3e-9 field is info only
  M = "Middle" (old version: M is missing)
  L = "Lone" (general 3D)
  R = "Rowlinson"
*/
{
  int i,idep,n,k;
  int indx=0;
  depend_t *d;
  real sum;
  char *key;
  enum deptype_e deptype=0;

  find("dependants");
  *head=NULL;

  loop (i,1,DEP_N) No.depend[i]=0;
  loop (idep,0,ndependants) {
    key=tok(__LINE__);
    if (toupper(key[0])=='M') {
      indx=atoi(tok(__LINE__));
      deptype=DEP_M; }
    else if (toupper(key[0])=='R') {
      deptype=DEP_R;
      indx=atoi(tok(__LINE__)); }
    else if (toupper(key[0])=='L') {
      deptype=DEP_L;
      indx=atoi(tok(__LINE__)); }
    else if (!isdigit(key[0]))
      ERROR(("unknown key %s in the table of dependants",key))
    else {
      deptype=DEP_OLDM;
      prt("NOTE: Old dependant format, \"Middle\" (linear) assumed for compatibility");
      indx=atoi(key); }

    No.depend[deptype]++;

    checkatomn(tok(__LINE__),indx); /* check and skip atom name */
    n=atoi(tok(__LINE__));

    if (deptype>=DEP_R) {
      if (n!=3) ERROR(("out-of-plane dependant requires 3 parents (%d given)",n))
      alloc(d,sizeof(depend_t)); }
    else {
      depend_t *aux; /* to calculate size only */
      alloc(d,(char*)&aux->dep[0]-(char*)aux+n*sizeof(struct depitem_s)); }

    d->type=deptype;
    d->n=n;
    d->indx=indx;
    d->next=*head; *head=d;

    sum=0;
    loop (i,0,n) {
      d->dep[i].i=atoi(tok(__LINE__));
      checkatomn(tok(__LINE__),d->dep[i].i); /* check and skip atom name */
      sum += d->dep[i].w=atof(tok(__LINE__)); }

    /* normalize to sum of weights = 1 */
    if (fabs(sum-1)>1e-5) WARNING(("sum of dependant weights=%g (should be 1) - fixed",sum))
    loop (i,0,n) d->dep[i].w/=sum;

    if (deptype==DEP_R) {
      /* one extra parm d->wz==LO length (w.sign) */
      d->wz=atof(tok(__LINE__)); }
    else if (deptype==DEP_L) {
      loop (k,0,3) d->x[k]=atof(tok(__LINE__)); zerofix(d->x,"x");
      loop (k,0,3) d->y[k]=atof(tok(__LINE__)); zerofix(d->y,"y");
      d->wz=atof(tok(__LINE__));
      loop (k,0,3) d->tx[k]=atof(tok(__LINE__)); zerofix(d->tx,"tx");
      loop (k,0,3) d->ty[k]=atof(tok(__LINE__)); zerofix(d->ty,"ty"); }

    key=strtok(NULL," \n\r\t");
    if (!key || key[0]!='e')
      if (deptype!=DEP_OLDM) WARNING(("key 'e' (marking end of data) expected in the following line of dependants:\n%s",copyofline))

    if (option('v')&4) {
      char *depinfo[DEP_L+1]= {
        "Middle (linear), old format",
        "Middle (linear)",
        "Rowlinson (perpendicular to a flexible triangle)",
        "Lone (general 3D position to a rigid triangle)"};
      prt("\nDependant type %s: %d parents", depinfo[deptype],d->n);
      loop (i,0,d->n)
        prt("  i=%d w=%.8f",d->dep[i].i,d->dep[i].w);
      if (d->type==DEP_R) {
        prt("  LO=%11.8f",d->wz);  }
      if (d->n==DEP_L) {
        prt("  x =%11.8f %11.8f %11.8f",VARG(d->x));
        prt("  y =%11.8f %11.8f %11.8f",VARG(d->y));
        prt("  wz=%11.8f",d->wz);
        prt("  tx=%11.8f %11.8f %11.8f",VARG(d->tx));
        prt("  ty=%11.8f %11.8f %11.8f",VARG(d->ty)); } }

    mygetline(); }

 /* NB: in initNo(), No.depend[0]=total # if dependants */
  No.depend[1]+=No.depend[0];
}

int depmassmoved;

static void depmass(int i,int ns,depend_t *dep) /* ----------------- depmass */
{
  struct depitem_s *di;
  siteinfo_t *si=specsp->si;

  while (dep) {
    if (dep->indx==i) {
      for (di=dep->dep; di->i>=0; di++) {
        if (di->i<0 || di->i>=ns)
          ERROR(("depmass: wrong dependant table (check masses of sites)"))
        si[di->i].mass+=di->w*si[i].mass; }
      depmassmoved++;
      si[i].mass=0; }
    dep=dep->next; }
}

static struct {
  char *sites;
  int sp;
  int nc;
} waterchk = {NULL,-1,0};

static void zeroHT(void) /* ----------------------------------------- zeroHT */
/* for TIP3p: HT-HT and HT-OW terms are set zero */
{
  nbfix_t *fix;
  int j,HT=atomn("?HT",-1,-1),OT=atomn("?OW",-1,-1);
  static int done;

  if (done) return;
  done++;

  /* if OW found, OT not tested (alcohol O in charmm21) */
  if (OT<0) OT=atomn("?OT",-1,-1);

  if (OT>=0 && HT>=0) {
    prt("\nTIP3P: H-O and H-H Lennard-Jones terms turned off because of option -x%d",option('x'));
    looplist (fix,nbfix0)
      if ( (fix->indx[0]==HT && fix->indx[1]==HT)
        || (fix->indx[0]==HT && fix->indx[1]==OT)
        || (fix->indx[0]==OT && fix->indx[1]==HT) )
        loop (j,0,n14) fix->onefour[j].eps=0; }
  else if (OT>=0 || HT>=0)
    prt("WARNING: OT or HT but no both (TIP3P)");
}


void readblend(void) /******************************************** readblend */
/*
   this procedure reads in species definitions as generated by blend
   returns: minimum nonzero mass in g/mol
*/
{
  bond_t *bond;
  angle_t *angle;
  nbfix_t *fix;
  char name[4][8];
  int i,sp,nx14;
  int polar=-1, gnsites=0, nheavy;
  double charge,mass,chargeq;
  static int nparms=0;
  siteinfo_t *si;
  int sqrt_rule=(int)-2147483333; /* legacy */
  int comb_rule=0; /* default if missing in *.ble file */
  int distance14=4; /* distance14=0 will turn all 1--3, 1--4, 1--5 exceptions off */
#if SS_PARMS>0
  int k;
#endif /*# SS_PARMS>0 */
  double qscaling=abs(option('q'))/optionscaling;
  double EvdWscaling=1;
  double RvdWscaling=1;

  if (option('j')<0) RvdWscaling=-option('j')/optionscaling;
  else EvdWscaling=option('j')/optionscaling;

#ifdef QQTAB
  if (option('a')!=option('_')) ERROR(("(?) option -a%d not supported with QQTAB",option('a')))
  if (option('a')!=option('_')) ERROR(("(?) option -a%d not supported with QQTAB",option('a')))
#endif /*# QQTAB */

  alloc(line,LINELEN);
  alloc(foundtoken,LINELEN);
  alloc(copyofline,LINELEN);
  getdata
      /* auxiliary variables for calculations */
        static double
          aux=0, a=0, b=0, c=0, x=0, y=0, z=0;
        static int
          i=0,j=0,k=0,n=0;
    get(nsites) get(nspec) get(nnbfixes)
    get(factor14) get(distance14)
    get(sqrt_rule) get(comb_rule)
    get(polar)
    get(nparms)
    get(eunit)
        /* auxiliary variables */
        get(aux)
        get(a) get(b) get(c)
        get(x) get(y) get(z)
        get(i) get(j) get(k)
        get(n)
    checkdata
  enddata

  if (distance14) n14=2; else n14=1;

  if (sqrt_rule!=(int)-2147483333) {
    prt("\n\
NOTE: Your ble-file contains parameter sqrt_rule, which has been renamed\n\
      to comb_rule. I am assuming comb_rule=sqrt_rule.");
    comb_rule=sqrt_rule; }

  if (nparms!=SS_PARMS)
    ERROR(("ble-file has wrong number of extra parameters\n\
expected SS_PARMS=%d but nparms=%d",SS_PARMS,nparms))

  if (polar<0) prt("WARNING: obsolete ble-file");

#ifdef POLAR
  if (!polar)
    WARNING(("POLAR version of cook used with nonpolar force field"))
  else
#  if POLAR&64
    if (polar!=5) WARNING(("this is ADIM (POLAR&64) version of cook\n\
*** but polar=%d read from %s (expected polar=5)",polar,lastFn))
#  else /*# POLAR&64 */
    if (polar!=1) WARNING(("this is not ADIM (POLAR&64) version of cook\n\
*** but polar=%d read from %s (expected polar=1)",polar,lastFn))
#  endif /*#!POLAR&64 */
#else /*# POLAR */
  if (polar==1) WARNING(("nonpolar version of cook used on POLAR force field"))
#endif /*#!POLAR */

  factor14_1=factor14-1;

  /* determining gnsites (# of groups) */
  {
    FILE *f=fopen(Fn("s-s"),"rt");
    char line[LINELEN];
    char *tok;

    if (f) {
      prt("reading %s: determine # of groups",lastFn);

      while (fgets(line,LINELEN,f)) if (!strchr("!#",line[0])) {
        if (strlen(line)>=LINELEN-1) ERROR(("%s: too long line",lastFn))
        tok=strtok(line," \n\r\t");

        if (!tok) continue;
        if (strlen(tok)==1 && (islower(tok[0]) || isdigit(tok[0])))
          gnsites++; }
      if (gnsites) {
        prt("%d group(s) read from %s",gnsites,lastFn);
        isss=1; }
      fclose (f); }
  }

  underline("combining rules");
  combruleinfo(comb_rule);

  allocarray(sitedef,nsites+gnsites);
  alloc(initsites,nspec*sizeof(vector*));

  find(polar<0?"site-types":SS_TABLE);

#ifdef POLAR
  //  if (abs(option('^'))>999) prt("Note: Drude charge (shell) changed because of option -^%d",option('^'));
#endif /*# POLAR */

  nx14=0;

  loop (i,0,nsites) {
    sitedef_t *sd=sitedef+i;
    int j;
    char *t;

    sd->ncharges=0;

    tok(__LINE__); /* skip # */
    strcpy(sd->name,tok(__LINE__));
    sd->M=atof(tok(__LINE__));
    sd->g=NULL;

    loop (j,0,n14) {
      t=tok(__LINE__);
      if (!strcmp(t,"=") || !strcmp(t,"-2146495994")) { nx14++; break; }
      sd->LJ[1].alpha=sd->LJ[j].alpha
        =atof(t);
      sd->LJ[1].EvdW=sd->LJ[j].EvdW
        =EvdWscaling*pow(eunit,Epow[0].EvdW)*atof(tok(__LINE__));
      sd->LJ[1].RvdW=sd->LJ[j].RvdW
        =RvdWscaling*pow(eunit,Epow[0].RvdW)*atof(tok(__LINE__));
#if SS_PARMS>0
      loop (k,0,nparms)
        sd->LJ[1].parm[k]=sd->LJ[j].parm[k]=pow(eunit,Epow[0].parm[k])*atof(tok(__LINE__));
#endif /*# SS_PARMS>0 */
      checktok=0;
    }

    /* since V2.8a, initcombrule does not include kcal->K conversion */
    initcombrule(&sd->LJ[0],comb_rule);
    initcombrule(&sd->LJ[1],comb_rule);

#ifdef POLAR
    {
      double alphapol=atof(tok(__LINE__));
      double chargepol=atof(tok(__LINE__));
      double Esat=atof(tok(__LINE__)); /* *eunit later */
      double kappa=atof(tok(__LINE__));
#  if POLAR&1
      int rep=atoi(tok(__LINE__));
#  endif /*# POLAR&1 */

      /* NOTE: -2146495994. marks missing token */
      if (alphapol!=-2146495994.)
        sd->alphapol=alphapol;
      else {
        sd->alphapol=sd->LJ[0].alpha;
        prt("WARNING: %2i %4s : alphapol missing, combining rule value %g substituted",
            i,sd->name,sd->alphapol); }

      if (chargepol!=-2146495994.)
        sd->chargepol=chargepol;
      else {
        sd->chargepol=-1000;
        WARNING(("missing chargepol, %g e substituted\n\
*** (distance14 or other parameter?)",sd->chargepol)) }

      /* calculation of chargepol from the Kirkwood combining rule
         REMOVED in V2.8a *****
      else {
        if (sd->alphapol && option('a'))
          sd->chargepol = sqr(256/362.3376*sd->LJ[0].EvdW*Pow6(sd->LJ[0].RvdW))
            /Cub(sd->alphapol);
        else
          sd->chargepol=0;
          WARNING(("%2i %4s : chargepol missing, combining rule value %g substituted",
                 i,sd->name,sd->chargepol)) }
      */

      sd->alphapol*=option('a')/optionscaling;

      //      if (abs(option('^'))>999) sd->chargepol=option('^')/1000;
#  if POLAR&1
      sd->rep = rep?QTYPE_SHELL2:0; /* QTYPE_SHELL2=2 if active */
#  endif /*# POLAR&1 */

      if (Esat!=-2146495994.)
#  if POLAR&2
        sd->Esat=Esat*=eunit; else sd->Esat=0;
#  else /*# POLAR&2 */
        if (Esat!=0) WARNING(("Esat present but !(POLAR&2)"))
#  endif /*#!POLAR&2 */

      if (kappa!=-2146495994.)
#  if POLAR&1
      sd->kappa=kappa/=electron; else sd->kappa=0;
#  else /*# POLAR&1 */
      if (kappa!=0) WARNING(("kappa present but !(POLAR&1)"))
#  endif /*#!POLAR&1 */

      prt("ATOM %2i %4s: eps=%g R=%g pol=%g chargepol=%g Esat=%g kappa=%g",
          i,sd->name,sd->LJ[0].EvdW,sd->LJ[0].RvdW,
                                     sd->alphapol,sd->chargepol,Esat,kappa);

      sd->chargepol*=electron; /* scaling by -q removed again */
    }
#endif /*# POLAR */

#ifndef FREEBC
    /* optional scattering length */
    checktok=0;
    sd->sfweight=atof(tok(__LINE__));
    if (sd->sfweight==-2146495994.)
      sd->sfweight=sd->M;
    else
      if (sd->sfweight)
        prt("structure factor weight (scattering length) set to b=%g (instead of M)",sd->sfweight);
#endif   /*# FREEBC */

    checktok=1;
    mygetline();
  }

  if (nx14) prt("NOTE: %d undefined parms for 1-4 terms assigned from non-bonded terms",nx14);

  /* groups - continued */
  if (gnsites) {
    FILE *f=fopen(Fn("s-s"),"rt");
    char line[LINELEN];
    int n,j,an;
    stitem_t stitem[LINELEN/4];
    char *tok,prefix;

    i=nsites;

    prt("reading %s: analyze groups",lastFn);
    if (!f) ERROR(("%s disappeared",lastFn))

    while (fgets(line,LINELEN,f)) if (!strchr("!#",line[0])) {
      if (strlen(line)>=LINELEN-1) ERROR(("%s: too long line",lastFn))
      tok=strtok(line," \n\r\t");

      if (!tok) continue;
      prefix=tok[0];
      if (strlen(tok)==1 && (islower(prefix) || isdigit(prefix))) {
        if (i>=nsites+gnsites) ERROR(("%s has changed length",lastFn))

        tok=strtok(NULL," \n\r\t");
        if (!tok) ERROR(("%s: missing type after %c",lastFn,prefix))

        an=atomn(tok,-1,-1);
        memcpy(sitedef+i,sitedef+an,sizeof(sitedef[i]));
        sitedef[i].name[0]=prefix;

        strcpy(sitedef[i].name+1,sitedef[an].name);

        n=0;

        tok=strtok(NULL," \n\r\t");
        if (!tok) ERROR(("%s: empty group",lastFn))
        while (tok) {
          char *dot=strchr(tok,'.');
          if (!dot) ERROR(("%s: missing . in group member SPECIES.ATOM",lastFn))
          stitem[n].sp=atoi(tok);
          stitem[n].i=atoi(dot+1);
          n++;
          tok=strtok(NULL," \n\r\t"); }
        alloc(sitedef[i].g,sizeof(stgroup_t)+sizeof(int)*(n-1));
        sitedef[i].g->n=n;
        loop (j,0,n) sitedef[i].g->stitem[j]=stitem[j];
        i++; } }
    fclose (f); }

  /* backwards compatibility */
  find(polar<0?"NBFIX":"nbfixes");

  if (gnsites && nnbfixes) ERROR(("sorry, rdf groups + nbfixes not implemented (checked) yet"))

  nx14=0;
  loop (i,0,nnbfixes) {
    int j;
    char *t;

    alloc(fix,sizeof(nbfix_t));
    fix->next=nbfix0;
    nbfix0=fix;

    fix->indx[0]=atomn(tok(__LINE__),-1,-1);
    fix->indx[1]=atomn(tok(__LINE__),-1,-1);
    loop (j,0,n14) {
      t=tok(__LINE__);
      if (!strcmp(t,"=") || !strcmp(t,"-2146495994")) { nx14++; break; }
      fix->onefour[1].eps=fix->onefour[j].eps=pow(eunit,Epow[1].EvdW)*atof(t);
      fix->onefour[1].sig=fix->onefour[j].sig=pow(eunit,Epow[1].RvdW)*atof(tok(__LINE__));
#if SS_PARMS
      loop (k,0,nparms)
        fix->onefour[1].parm[k]=fix->onefour[j].parm[k]=pow(eunit,Epow[1].parm[k])*atof(tok(__LINE__));
#endif /*# SS_PARMS */
      checktok=0; }
    checktok=1;
    mygetline(); }

  if (nx14) prt("NOTE: %d undefined 1-4 nbfixes assigned from non-bonded terms",nx14);

  nsites+=gnsites;
  alloc(spec,nspec*sizeof(specinfo_t*));

#ifdef LINKCELL
  box.threebonds=0;
#endif /*# LINKCELL */

  loop (sp,0,nspec) {
#ifdef POLAR
    int naxials=0;
#  if POLAR&32
    int nfq4=0;
#  endif /*# POLAR&32  */
#endif /*# POLAR */
    int i=0,N=0;
    int water=0,config=0;
    int nc=0,nangles=0,ndihedrals=0,nimpropers=0,naromatics=0;
    double zero_energy=0;
    int ndependants=-1; /* compatibility: does not search for the table at all */
    int nvb=0;

    ns=0; /* see checkatomn */

    find("species");
    prts(copyofline);

    getdata
      get(i) get(N)
      get(ns) get(nc) get(nangles)
      get(ndihedrals) get(nimpropers) get(naromatics)
      get(ndependants)
      get(zero_energy)
      get(config) get(water)
#ifdef POLAR
      get(naxials)
#  if POLAR&32
      get(nfq4)
#  endif /*# POLAR&32 */
#endif /*# POLAR */
      checkdata
    enddata

    No.bonded+=nc+nangles+ndihedrals+nimpropers+naromatics;

    if (i!=sp) WARNING(("species index sp=%d (expected %d)",i,sp))
    alloc(specsp,sizeof(specinfo_t)-sizeof(siteinfo_t)+max(ns,nc)*sizeof(siteinfo_t));
    spec[sp]=specsp;
    specsp->name=strdup(foundtoken);

    if (config) {
      if (N==0) ERROR(("species sp=%d config: N=0",sp))
      alloc(initsites[sp],sizeof(vector)*ns*N); }
    else
      alloc(initsites[sp],sizeof(vector)*ns);

#ifdef METAL
    if (ns>1) ERROR(("sp=%d: METAL requires ns=1 but ns=%d found",sp,ns))
#endif /*# METAL */

    specsp->pot=0;

#ifndef WATER
    if (option('x')&1) {
      option('x')-=1;
      if (water) prt("\
WARNING: A registered rigid water model was passed to cook, but this version\n\
         of cook does not support direct optimized rigid water models.\n\
INFO: to export a rigid model to the ble-file, use 'blend -h WATER'");
      else
        prt("\noptimized water models not supported, option -x1 turned off");
    }
#endif /*# WATER */

    { char *c; for (c=foundtoken; *c; c++) *c=toupper(*c); }

    /* TIP3P: remove H-H and H-O LJ terms according to -x */
    if (checkmodel("TIP3P") || checkmodel("HOH"))
      if (option('x')&2) WARNING(("\
You are using the TIP3P (or HOH = nickname for TIP3P) water model with\n\
*** nonzero water-water H-H and H-O Lennard-Jones terms.  To clear these\n\
*** interaction as in the original TIP3P model, use -x%d (without flag 2)",option('x')-2))

    if ( (option('x')&2) == 0) zeroHT();

    waterchk.sites=NULL;

#ifdef PERSUM
    if (option('x')&1) WARNING(("Option -x1 with PERSUM is not recommended:\n\
*** wrong Ewald intramolecular correction\n\
*** may be correct if the r-space potential is ~zero at r=box size"))
#endif /*# PERSUM */

#ifdef WATER
    if (option('x')&1) {
#  ifdef WATERPLUS
      /* special potentials (not water...) - special project only */
#    include "xxxdef.c"
#  endif       /*# WATERPLUS  */
      if (checkmodel("TIP3P") || checkmodel("SPC") || checkmodel("SPCE") || checkmodel("HOH")) {
        specsp->pot=3;
        waterchk.sites="HHO";
        waterchk.nc=3;
        waters.iH=0; waters.iO=2; waters.iM=2; }
      if (checkmodel("TIP4P") || checkmodel("TIP4PEW") || checkmodel("TIP4P05") || checkmodel("TIP4PICE")) {
        specsp->pot=4;
        waterchk.sites="HOMH";
        waterchk.nc=3;
        waters.iH=0; waters.iO=1; waters.iM=2; }
      if (checkmodel("ST2")) {
        specsp->pot=5;
        waterchk.sites="OLHLH";
        // ?        waterchk.nc=9;
        waterchk.nc=6; /* not good for lone dependants */
        waters.iH=2; waters.iO=0; waters.iM=1; } }
#endif /* WATER */ /*# WATER */

    /* for waters, specsp->pot==ns */
    if (specsp->pot) {
      prt("\nspecies %s: optimized code for rigid water-water potential will be used",foundtoken);
      if (nangles)
        ERROR(("species %s: rigid model recognized by name (because -x1),\n\
*** but flexible angle(s) found:\n\
*** to use the flexible model irrespective of name, run cook* -x0\n\
*** to prepare a rigid model, use blend -h to constrain angles",foundtoken))
#ifdef WATERPLUS
#  if WATERPLUS>1
      if (specsp->pot<6)
#  endif /*# WATERPLUS>1 */
#endif /*# WATERPLUS */
      if (specsp->pot!=ns)
        ERROR(("species %s: %d sites, %d expected\n\
*** (if incorrectly recognized as water, rename species in %s.ble)",foundtoken,ns,specsp->pot,simils.sysname))
      if (waterchk.nc!=nc)
        ERROR(("species %s: %d bond constraints, %d expected\n\
*** (if incorrectly recognized as water, rename species in %s.ble)", foundtoken,nc,waterchk.nc,simils.sysname))
    }

    specsp->ns=ns;
    specsp->nc=0; /* the real count depends on -u = vibrating bonds */
    specsp->N=N;
    specsp->config=config;
    specsp->bond=NULL;
    specsp->angle=NULL;
    specsp->dihedral=specsp->improper=specsp->aromatics=NULL;

    find("sites");
    mass=charge=chargeq=0;
    nheavy=0;

    nx14=0;

    loop (i,0,ns) {
      real *v;
      int j,n;
      char *token;

      si=specsp->si+i;
      v=initsites[sp][i];
      j=atoi(tok(__LINE__));
      if (i!=j) ERROR(("sites: site number %d (expected %d)",j,i))
      si->st=atomn(token=tok(__LINE__),sp,i);
      if (waterchk.sites) if (waterchk.sites[i]!=token[0])
        ERROR(("sites: %s: bad order of sites (%s expected)",foundtoken,waterchk.sites))

      /* GLOBAL scaling of charges by option -q
         since V2.7r also Drude, removed again in 3.1c */
      si->charge=electron*qscaling*atof(tok(__LINE__));
      charge += si->charge;
      chargeq += Sqr(si->charge);

#if COULOMB<-2 && defined(SLAB)
      /* sqrt(2)*sigma, for charge density profile of GAUSSIANCHARGES */
      si->esig=sqrt(2.0)*sitedef[si->st].LJ[0].parm[SS_PARMS-1];
#endif

#ifdef POLAR
      {
        /* both scaled by -q */
        double qcenter=si->charge; /* total now */
        double qpol=sitedef[si->st].chargepol;

        if (qcenter) si->qtype |= QTYPE_CHARGE; /* nonzero atom (total) charge */
        if (qpol) si->qtype |= QTYPE_DRUDE; /* nonzero Drude charge */
#  if POLAR&2
        si->Esat=sitedef[si->st].Esat;
#  endif /*# POLAR&2 */
#  if POLAR&1
        si->qkappa=sitedef[si->st].kappa;
        // ?: si->qkappa=q*sitedef[si->st].kappa;
        si->qtype |= (si->qkappa!=0) | sitedef[si->st].rep; /* shell-core: 1|2 */
#  endif /*# POLAR&1 */
        si->chargepol=qpol;

        /* final central charge (total-Drude); small (rounding error) -> zero */
        si->charge -= qpol;
        if (si->charge!=0 && fabs(si->charge)/(fabs(qpol)+fabs(qcenter)+fabs(si->charge))<1e-15) {
          prt("WARNING: central Drude charge %g p.u. changed to exact zero",si->charge);
          si->charge=0; }
        if (si->charge) {
          si->qtype |= QTYPE_CENTER; /* nonzero central charge */
#  if defined(QQTAB) && QQTAB==0
          ERROR(("nonzero central charge not supported by QQTAB=0 (COS models)\n\
*** in simopt.h: #define QQTAB 1 (central+Drude)\n\
***              #define QQTAB 2 (central+Drude mixed with COS models)"))
#  endif /*# defined(QQTAB) && QQTAB==0 */
        }
        si->qqpol = qpol = qpol*qpol;
        if (qpol) qpol=1/qpol;
        si->alpha_qq = sitedef[si->st].alphapol*qpol;
        if (option('v')&4)
          prt("POLAR final charges in e for site[%d] : site=%g Drude=%g",i,si->charge/electron,si->chargepol/electron);
      }
#endif /*# POLAR */
      loop (j,0,DIM) v[j]=atof(tok(__LINE__)); /* not used if config */

      switch (sitedef[si->st].ncharges) {
        case 0: sitedef[si->st].charge=si->charge; sitedef[si->st].ncharges++; break;
        case 1: if (fabs(sitedef[si->st].charge-si->charge)>1e-13) sitedef[si->st].ncharges++; }

      mass += si->mass = sitedef[si->st].M/Munit;
      nheavy += sitedef[si->st].M>1.5;
      /* si->imass moved because of dependants */
#ifndef FREEBC
      si->sfweight=sitedef[si->st].sfweight;
#endif /*# FREEBC */

      /* exceptions */
      n=atof(tok(__LINE__));

      if (n<0) ERROR(("negative # of exceptions or wrong format of ble-file"))
      alloc(si->exc,(n+1)*sizeof(exception_t));
      loop (j,0,n) {
        token=tok(__LINE__);
        if (token[0]=='*') {
          nx14++;
          si->exc[j].indx=atoi(token+1);
          si->exc[j].type=ONEFOUR; }
        else {
          si->exc[j].indx=atoi(token); si->exc[j].type=EXCL; } }
      si->exc[n].indx=i; /* this is the actual algorithm sentinel */
      si->exc[n].type=END; /* this sentinel used by MPI broadcast */
      mygetline(); }

    if (distance14 && !nx14) prt("species %d: no 1--4 interaction found, 1-4 terms turned off",sp);
    if (!distance14 && nx14) ERROR(("1-4 terms missing in the ble-file (because distance14=0),\n\
*** but %d such terms found in the species",nx14))

    if (nheavy<=1 && nangles)
      WARNING(("Only %d heavy atom(s) and %d flexible angle(s) in species %d,\n\
*** consider a rigid molecule instead (use blend -h SPECIES).",
               nheavy,nangles,sp))

    /* small charge := 0 */
    if (Sqr(charge)/chargeq<Sqr(No.eps)*(ns+1)) {
      prt("INFO: small charge %d e of species %d is replaced by zero",charge,sp);
      charge=0; }
    specsp->charge=charge;
    specsp->mass=mass; /* molecule mass: for equalize.cfg will be changed in initNo() */

    si=specsp->si;

    if (waterchk.sites) {
      if (waterchk.sp>=0) ERROR(("%s: cannot combine two waters",foundtoken))
      waterchk.sp=sp;
#ifdef WATER
      waters.qq.HH=si[waters.iH].charge*si[waters.iH].charge;
      waters.qq.MH=si[waters.iM].charge*si[waters.iH].charge;
      waters.qq.MM=si[waters.iM].charge*si[waters.iM].charge;

      /* DOUBLECHECK of water models */
      {
        double qM=si[waters.iM].charge/electron;
        double qO=si[waters.iO].charge/electron;
        double qH=si[waters.iH].charge/electron;
        sitesite_t dummy;
        double Emin=sitedef[si[waters.iO].st].LJ[0].EvdW;
        double RvdW=sitedef[si[waters.iO].st].LJ[0].RvdW;
        int iw,iwiw;
        static struct {
          int wtype; char *wm; double qO,qH,Emin,RvdW; } w[] = {
            { 3,"TIP3P",-0.834,.417,-76.525816,1.768246},
/* SPC,SPCE are according to Berendsen etal, JPC 91, 6269 (1987): */
            { 3,"SPC", -0.82,0.41, -78.19734,1.7766093},
            { 3,"SPCE",-.8476,.4238,-78.19734,1.7766093},
/* SPC/GROMOS96 version: for combining rules only:
            { 3,"SPC-G96",-.82,.41,-86.575405,1.747353},
            { 3,"SPC-G96",-.82,.41,-86.4037,1.747353}, */
            { 4,"TIP4P",0,0.52, -78.02,1.7698858},
            { 4,"TIP4PEW",0,0.52422,-81.898887, 1.7759314},
            { 4,"TIP4P05",0,0.5564, -93.2, 1.77287268},
            { 4,"TIP4PICE",0,0.5897,-106.1,1.77730641},
            { 5,"ST2",  0,.2357,-38.101,  1.7398162},
/*  ??? - TIP5P does  not use the switch function, old version anyway
            { 5,"TIP5P",.241, -80.515, 1.7510408}, */
            { 0,NULL,0,0,0 } };

#ifdef SLAB
        dummy.Skk=NULL;
#endif
        setss(&dummy,si[waters.iO].st,si[waters.iO].st,15.0,2); /* dummy call, C2 is arbitrary */
        /* WARNING: this requires LJ (guaranteed by #define WATER) */
        if (fabs(Emin-dummy.a.E4/-4)+fabs(RvdW-sqrt(dummy.a.Sq)*0.5612310241546865)>1e-10) {
          prt("\
NOTE: water Oxygen parameters changed by nbfixes\n\
      Emin=%.8f RvdW=%.8f (generic O values)\n\
      Emin=%.8f RvdW=%.8f (nbfix values)",
             Emin,RvdW,
             dummy.a.E4/-4,sqrt(dummy.a.Sq)*0.5612310241546865);
          Emin=dummy.a.E4/-4;
          RvdW=sqrt(dummy.a.Sq)*0.5612310241546865; }

        prt("this model: qH=%6.4f Emin=%8.4f RvdW=%6.4f qM=%6.4f qO=%6.4f ",qH,Emin,RvdW,qM,qO);

        if (
            fabs(qO+qM*(waters.iM!=waters.iO)+qH*2)>1e-6 // TIP3,4P, SPC*
            &&
            fabs(qO+qM*2+qH*2)>1e-6 // ST2
            ) {
          if (option('x')&4)
            ERROR(("WATER doublecheck: molecule not neutral\n\
*** (if you are sure, turn off the -x4 flag to suppress this message)"))
          else
            WARNING(("WATER doublecheck: molecule not neutral")) }


        iwiw=-1;
        for (iw=0; w[iw].wtype; iw++) {
          int match=fabs(qH-w[iw].qH)<5e-5
                 && fabs(qO-w[iw].qO)<5e-5
                 && fabs(Emin-w[iw].Emin)<0.03
                 && fabs(RvdW-w[iw].RvdW)<1e-3;

          prt("%-12sqH=%6.4f Emin=%8.4f RvdW=%6.4f %s",
	      w[iw].wm, w[iw].qH,w[iw].Emin,w[iw].RvdW, match?"<- matches":"");
          if (match) iwiw=iw; }

	if (iwiw>=0) {
          if (w[iwiw].wtype==specsp->pot)
            prt("doublecheck: the parameters match the %s water model",w[iwiw].wm);
          else {
            if (option('x')&4)
              WARNING(("parameters do not match the water model"))
            else
              ERROR(("parameters do not match the water model\n\
*** (if you are sure, turn off the -x4 flag to suppress this message)")) } }
	else {
          if (option('x')&4)
            ERROR(("parameters do not correspond to any registered water model\n\
*** (if you are sure, turn off the -x4 flag to suppress this message)"))
          else
            WARNING(("parameters do not correspond to any registered water model")) }
      }

#  if defined(WATERPLUS) && WATERPLUS==2
      waters.qq.OO=si[waters.iO].charge*si[waters.iO].charge;
      waters.qq.OH=si[waters.iO].charge*si[waters.iH].charge;
      waters.qq.OM=si[waters.iO].charge*si[waters.iM].charge;
#  endif /*# defined(WATERPLUS) && WATERPLUS==2 */
#endif /*# WATER */
    }

    if (config) {
      int n,i,j;
      real *v;
#ifdef GOLD
      double minz=9e9;
#endif /*# GOLD */

      find("config");

      loop (n,0,N) {
        loop (i,0,ns) {
          v=initsites[sp][n*ns+i];
          loop (j,0,DIM) v[j]=atof(tok(__LINE__));
#ifdef GOLD
          Min(minz,v[2])
#endif /*# GOLD */
          mygetline(); }
        mygetline(); }

#ifdef GOLD
      if (wall.minz>0) {
        loop (n,0,N)
          loop (i,0,ns) {
            v=initsites[sp][n*ns+i];
            v[2]+=wall.minz-minz; }
        prt("config shifted to minz=%f",wall.minz); }
      else {
        prt("not shifted, minz=%f",minz);
        if (minz<=0) WARNING(("minz<=0")) }
#endif /*# GOLD */
      }
#ifdef GOLD
    else
      if (wall.minz<0) WARNING(("GOLD: not config and wall.minz=%f<=0",wall.minz))
#endif /*# GOLD */


    find("bonds");

    loop (i,0,nc) {
      real K,length;
      int i,j,vb;
      char why[64];

#ifdef MORSE
      real a;
      if (7!=sscanf(line, "%d%s%d%s" realfmt realfmt realfmt,
                    &i,name[0],&j,name[1],&K,&length,&a))
#else /*# missing #if */
      if (6!=sscanf(line, "%d%s%d%s" realfmt realfmt,
                    &i,name[0],&j,name[1],&K,&length))
#endif /*#!missing #if */
        ERROR(("bonds: missing or bad format"))

      checkatomn(name[0],i); checkatomn(name[1],j);

      /* vibrating or rigid bond ? */
      if (option('u')>=0) {
        vb=K<option('u')+0.5;
        sprintf(why,"K=%g kcal/mol/AA2 %c option -u",K,"><"[vb]); }
      else {
        /* new in V3.4d, negative -u: use vibrational wavenumber [cm-1]
           units: here K is in kcal/mol/AA2, formula U=K(r-r0)^2
           M is in g/mol */
        int ii=atomn(name[0],-1,-1);
        int jj=atomn(name[1],-1,-1);
        double Mi=sitedef[ii].M,Mj=sitedef[jj].M; /* in g/mol */
        /* constant: evu "to cm-1; \(1[kcal/mol AA2]*2/1[g/mol])/2/pi/c" */
        double nu=153.57137207786*sqrt(K*(1/Mi+1/Mj)); /* evu: \(1[kcal/mol AA2]*2/1[g/mol])/2/pi/c; to cm-1 */

        vb=nu<abs(option('u'))+0.5;
        sprintf(why,"nu=%g [cm-1] %c |option -u|",nu,"><"[vb]); }

      if (vb && K==0) {
        vb=0;
        strcpy(why,"K==0"); }

      if (vb && specsp->pot) {
        WARNING(("%s\nspecies %s (rigid model): cannot vibrate bonds\n\
*** (clear flag 1 of -x to use a generic nonoptimized code for water)",line,foundtoken))
        vb=0;
        strcpy(why,"the molecule is rigid (specsp->pot exits)"); }

      if (option('v')&4)
        prt("BOND %s[%d]-%s[%d] is %s because %s",
            name[0],i,name[1],j,vb?"VIBRATING":"RIGID",why);

      nvb+=vb;

      if (vb) {
        /* vibrating bond */
        alloc(bond,sizeof(bond_t));
        bond->next=specsp->bond;
        specsp->bond=bond;
        bond->indx[0]=i; bond->indx[1]=j;
        bond->parm.length=length;
        K*=(ireal)eunit;
        bond->parm.K=K;
#ifdef MORSE
        bond->parm.a=a;
        bond->parm.Ka=K*a;
        bond->parm.Ka2=7./12*K*Sqr(a);
        bond->parm.fa1=1.5*a;
        bond->parm.fa2=7./6.*Sqr(a);
#endif /*# MORSE */
        bond->parm.Ki2=-2*K; }
      else {
        /* constraint bond: K>option -u or K==0 */
        si->pair[0]=i; si->pair[1]=j;
        si->bondq=Sqr(length);
        specsp->nc++;
        si++; }
#ifdef LINKCELL
      /* since V2.4k, this serves for the initial check only - see update14() */
      Max(box.threebonds,length)
#endif /*# LINKCELL */
      mygetline(); }

    prt("species %d BONDs: %d RIGID, %d VIBRATING (option -v4 for details)",sp,nc-nvb,nvb);

    find("angles");

    loop (i,0,nangles) {
      ireal dummy;
      int items;

      alloc(angle,sizeof(angle_t));
      angle->next=specsp->angle;
      specsp->angle=angle;

      items=sscanf(line,"%d%s%d%s%d%s"
                   irealfmt irealfmt irealfmt irealfmt irealfmt irealfmt,
                   &angle->indx[0],name[0],
                   &angle->indx[1],name[1],
                   &angle->indx[2],name[2],
                   &angle->parm.K, &angle->parm.angle,
                   &dummy,&dummy,
                   &angle->parm.Kub,&angle->parm.length);
      if (items>=8 && items<=10) {
        angle->parm.Kub=angle->parm.length=0;
        items=12; }
      checkatomn(name[0],angle->indx[0]);
      checkatomn(name[1],angle->indx[1]);
      checkatomn(name[2],angle->indx[2]);
      if (items!=12) ERROR(("angles: missing or bad format"))
      angle->parm.angle*=(ireal)(PI/180);
      if (option('h')==4)
        if (fabs(angle->parm.angle-TETRAHEDRAL)<1e-4)
          /* angle 109.47 replaced by exact tetrahedral (109.5 not) */
          angle->parm.angle=TETRAHEDRAL;
      angle->parm.K*=(ireal)eunit;
      angle->parm.K2=angle->parm.K*2;
      angle->parm.Kub*=(ireal)eunit;
      angle->parm.Kubi2=-2*angle->parm.Kub;
      mygetline(); }

    readtorsion("dihedrals",ndihedrals,sp,&specsp->dihedral);
    readtorsion("impropers",nimpropers,sp,&specsp->improper);
    readtorsion("aromatics",naromatics,sp,&specsp->aromatics);

    specsp->intra=specsp->bond || specsp->angle
      || specsp->dihedral || specsp->improper || specsp->aromatics;

    if (ndependants>=0) readdependants(ndependants,&specsp->dependants);
    else specsp->dependants=NULL;

    /* New: what if some dependants have nonzero mass? redistribute! */
    depmassmoved=0;
    si=specsp->si; /* (bug found by Brian Kinnear) */
    loop (i,0,ns)
      if (si[i].mass!=0) depmass(i,ns,specsp->dependants);
    if (depmassmoved)
      WARNING(("Mass of %d dependants in species %d redistributed to neighbors.\n\
*** Kinetic properties will be affected!",depmassmoved,sp))

    /* molecule-based equalization */
    if (equalize.mol && (sp==equalize.sp || equalize.sp<0 && sp<=-equalize.sp) && ns>1) {
      /* equalize.mol masses */
      double M=0;
      int nm=0;

      loop (i,0,ns) {
        M+=si[i].mass;
        if (si[i].mass) nm++; }
      if (nm>1) {
        prt("\
Equalizing masses of atoms in species %d by factor equalize.mol=%g,\n\
mass of the molecule is unchanged.\n\
WARNING: Kinetic quantities are affected!",sp,equalize.mol);
      if (equalize.mol<0 || equalize.mol>1)
        WARNING(("equalize.mol = %g out of range [0,1]",equalize.mol))
        M/=nm;
        loop (i,0,ns) if (si[i].mass) {
          if (option('v')&2) prt_("%2d: %.4f", i, si[i].mass*Munit);
          si[i].mass=equalize.mol*M+(1-equalize.mol)*si[i].mass;
          if (option('v')&2) prt(" -> %.4f g/mol", si[i].mass*Munit); } } }

    specsp->zero_energy=eunit*zero_energy;

#ifdef POLAR
#  if POLAR&8
    /* read axials */
    loop (i,0,ns) specsp->si[i].tozz=-1;
    find("axials");
    loop (i,0,naxials) {
      int ii,tozz;
      real kappa,alpha,alphazz,Esat;

      if (sscanf(line,"%d%s%d%s"realfmt realfmt realfmt realfmt,
                 &ii,name[0],&tozz,name[1],&kappa,&alpha,&alphazz,&Esat)!=8)
        ERROR(("axials: missing or bad format"))
      checkatomn(name[0],ii);
      checkatomn(name[1],tozz);
      if (ii<0 || ii>=ns || tozz<0 || tozz>=ns || tozz==ii)
         ERROR(("axials: %d %d are bad atoms",ii,tozz))

      alpha*=option('a')/optionscaling;
      alphazz*=option('a')/optionscaling;

      si=specsp->si+ii;
      si->tozz=tozz;
      si->qkappa=kappa/electron*sitedef[si->st].chargepol;
      if (option('a')) {
        si->alpha_qq=alpha/si->qqpol;
        si->alphazz_qq=alphazz/si->qqpol; }
      else
        si->alpha_qq=si->alphazz_qq=0;
#    if POLAR&2
      si->Esat=Esat;
#    endif /*# POLAR&2 */
      mygetline(); }
#  else /*# POLAR&8 */
    if (naxials) ERROR(("naxials!=0 and !(POLAR&8)"))
#  endif /*#!POLAR&8 */

#  if POLAR&32

#    if defined(POLAR) && POLAR&32
    {
      fq4_t *fq4;
      int k;

      looplist (fq4,specsp->fq4)
        loop (k,0,4) {
          si=specsp->si+fq4->indx[k];
          if (si->qtype&QTYPE_DRUDE)
            ERROR(("spec=%d site=%d: both Drude and fluctuating charge not allowed",sp,fq4->indx[k]))
          si->qtype |= QTYPE_FQ; }
    }
#    endif /*# defined(POLAR) && POLAR&32 */


    find("fq4");
    loop (i,0,nfq4) {
      fq4_t *fq4;
      int k;
      double iden;

      alloconezero(fq4);
      fq4->next=specsp->fq4;
      specsp->fq4=fq4;
      loop (k,0,4) {
        fq4->indx[k]=atoi(tok(__LINE__));
        checkatomn(tok(__LINE__),fq4->indx[k]);
        if (specsp->si[fq4->indx[k]].qtype&QTYPE_DRUDE) ERROR(("both FQ and Drude for atom %d in species %d",fq4->indx[k],sp))
        specsp->si[fq4->indx[k]].qtype |= QTYPE_FQ; }

      WARNING(("fq4: this code should be checked"))
      /* if (option('a')) ? */
      iden=optionscaling/option('a');
      fq4->AHH=atof(tok(__LINE__))*iden;
      fq4->ALL=atof(tok(__LINE__))*iden;
      fq4->AHL=atof(tok(__LINE__))*iden;
      fq4->AHL2=fq4->AHL*2;

      iden=0.5/(fq4->AHH*(fq4->AHH+fq4->ALL-4*fq4->AHL));
      fq4->H1H1=(4*fq4->AHL-2*fq4->AHH-fq4->ALL)*iden;
      fq4->H1H2=(fq4->ALL-4*fq4->AHL)*iden;
      fq4->HL=fq4->AHH*iden;

      iden=0.5/(fq4->ALL*(fq4->AHH+fq4->ALL-4*fq4->AHL));
      fq4->L1L1=(4*fq4->AHL-2*fq4->ALL-fq4->AHH)*iden;
      fq4->L1L2=(fq4->AHH-4*fq4->AHL)*iden;
      fq4->LH=fq4->ALL*iden; } }
#  endif /*# POLAR&32 */
#endif /*# POLAR */

  } /* loop sp */

  free(copyofline);
  free(foundtoken);
  free(line);
  fprintf(stderr,"%i species read\n",nspec);
  _n

#ifdef LINKCELL
  prt("max bond=%f  pessimistic exclimit=%f",box.threebonds,box.threebonds*3);
  box.threebonds*=3; /* 3*max bond = limit for 1-4 interactions */
  /* this is pessimistic */
#endif /*# LINKCELL */

  if (isss>1)
    prt("\
rdf group support: %d atom types in bonds, angles, torsions\n\
                   replaced by their group copies",isss-1);

  if (!nspec) DISASTER(("no species - nothing to do"))

}

vector *initr(int sp) /*********************************************** initr */
/***
  Array[ns] of vectors containing the initial configuration for species sp
  is allocated and returned.  It should be freed when no longer needed!
  The initial configuration need not satisfy exactly the constraints.
  It is corrected later (by Scorrect in initcfg)
  COOK patch: if species is given by abs config, no array is allocated and
  only reference is returned : must NOT be freed!
***/
{
  vector *r;

  if (spec[sp]->config) return initsites[sp];

  alloc(r,sizeof(vector)*spec[sp]->ns);
  copy(r,initsites[sp],sizeof(vector)*spec[sp]->ns);

  return r;
} /* initr */


void initss(double LJcutoff) /*************************************** initss */
/*
   initialize tables of site-site potentials
   negative LJcutoff means cutoff in units of sigma
   part of repeatable initializations: will be cleared by global release
*/
{
  int i,j,sp,ntot,nLJ;

  /* zero required by Skk */
  ralloc2Darrayzero(sstab,nsites,nsites);
  if (n14==2) ralloc2Darrayzero(sstab14,nsites,nsites);

  loop (i,0,nsites) {
    if (sitedef[i].ncharges>1)
#ifdef QQTAB
      ERROR(("Site type %d = %s appears with different charges, which is not supported\n\
*** by the QQTAB version of cook\n\
*** turn off QQTAB in simopt.h (not for GAUSSIANCHARGES), duplicate site type",i,sitedef[i].name))
#else /*# QQTAB */
      if (option('v'))
        prt("Site type %d = %s appears with different charges (this is OK unless QQTAB)", i,sitedef[i].name);
#endif /*#!QQTAB */
  }

#ifdef WATER
  if (waterchk.sp>=0) {
    siteinfo_t *si=spec[waterchk.sp]->si;

#  ifndef WATERPLUS
    if (rdf) if (waters.iM!=waters.iO && waters.iM!=waters.iH) {
      /* do not consider rdf of M (L) sites at all */
      /* NOTE: rdf[][] is used by setss and copied to ss */
      j=si[waters.iM].st;
      loop (i,0,nsites) rdf[i][j]=rdf[j][i]=NULL; }
#  endif /*# WATERPLUS */

    setss(&waters.ss.HH, si[waters.iH].st,si[waters.iH].st, LJcutoff,0);
    setss(&waters.ss.OH, si[waters.iO].st,si[waters.iH].st, LJcutoff,0);
    setss(&waters.ss.OO, si[waters.iO].st,si[waters.iO].st, LJcutoff,0);
#  ifdef WATERPLUS
    setss(&waters.ss.MH, si[waters.iM].st,si[waters.iH].st, LJcutoff,0);
    setss(&waters.ss.OM, si[waters.iO].st,si[waters.iM].st, LJcutoff,0);
    setss(&waters.ss.MM, si[waters.iM].st,si[waters.iM].st, LJcutoff,0);
#  endif /*# WATERPLUS */
    waters.ss.dummy.rdf=NULL; }
#endif /*# WATER */

  loop (i,0,nsites) loopto (j,0,i) {
      setss(&sstab[i][j], i,j, LJcutoff,0);
      if (n14==2) setss(&sstab14[i][j], i,j, LJcutoff,1);
      if (i!=j) {
        sstab[j][i]=sstab[i][j];
        if (n14==2) sstab14[j][i]=sstab14[i][j]; } }

  /* set isLJ to 1 if there is any other atom (incl. nbfixes) to interact
     by LJ (does not apply to 1-4) */
  loop (sp,0,nspec) {
    int ns=spec[sp]->ns;
    siteinfo_t *si=spec[sp]->si;

    ntot=nLJ=0;

    loop (i,0,ns) {
      si[i].isLJ=0;
      ntot++;
      loop (j,0,nsites)
        if (sstab[si[i].st][j].A) si[i].isLJ=1;
      if (si[i].isLJ) nLJ++; }
    prt("species %d: %d (of %d) site(s) with (at least one) nonzero LJ interaction",sp,nLJ,ntot);
  }

#ifdef QQTAB
  loop (i,0,nsites)
    loop (j,0,nsites) {
#  ifdef POLAR
      /* NB: cannot use sstab[j][i]=sstab[i][j] because of qqtab[1][0],qqtab[0][1] */
      sstab[i][j].  qqtab[0][0].qq=sstab[i][j].  qqtabx[0][0].qq=sitedef[i].charge*sitedef[j].charge;
      sstab[i][j].  qqtab[0][1].qq=sstab[i][j].  qqtabx[0][1].qq=sitedef[i].charge*sitedef[j].chargepol;
      sstab[i][j].  qqtab[1][0].qq=sstab[i][j].  qqtabx[1][0].qq=sitedef[i].chargepol*sitedef[j].charge;
      sstab[i][j].  qqtab[1][1].qq=sstab[i][j].  qqtabx[1][1].qq=sitedef[i].chargepol*sitedef[j].chargepol;
      if (n14==2) {
        sstab14[i][j].qqtab[0][0].qq=sstab14[i][j].qqtabx[0][0].qq=sitedef[i].charge*sitedef[j].charge*factor14;
        sstab14[i][j].qqtab[0][1].qq=sstab14[i][j].qqtabx[0][1].qq=sitedef[i].charge*sitedef[j].chargepol*factor14;
        sstab14[i][j].qqtab[1][0].qq=sstab14[i][j].qqtabx[1][0].qq=sitedef[i].chargepol*sitedef[j].charge*factor14;
        sstab14[i][j].qqtab[1][1].qq=sstab14[i][j].qqtabx[1][1].qq=sitedef[i].chargepol*sitedef[j].chargepol*factor14; }
#  else /*# POLAR */
      sstab[i][j].qq=sitedef[i].charge*sitedef[j].charge;
      if (n14==2) sstab14[i][j].qq=sstab[i][j].qq*factor14;
#  endif /*#!POLAR */
    }
#endif /*# QQTAB */

} /* initss */

#ifdef QQTAB
/* separate charge-charge splines */

double setqq(ertab_p *tab,double qq, /******************************** setqq */
#  if COULOMB<-2
             double sigma1,double sigma2,
#  endif /*# COULOMB<-2 */
             double factor14)
/* interface to initerfc (splines); tab is allocated */
{
  return initerfc(tab,qq,
#  if COULOMB<-2
                  sigma1,sigma2,
#  endif /*# COULOMB<-2 */
                  factor14,
                  option('v')&1 ? el.grid : -el.grid,
#  ifndef GAUSSIANCHARGES
                  el.minqq,
#  endif /*# GAUSSIANCHARGES */
                  box.cutoff+el.rplus,
                  box.cutoff,
                  option('v')&4?-el.alpha:el.alpha,
                  el.rshift);
}

#  if COULOMB<-2
/* GAUSSIANCHARGES: GCPARM=sigma_i,sigma_j */
#    define GCPARM sitedef[i].LJ[0].parm[SS_PARMS-1],sitedef[j].LJ[0].parm[SS_PARMS-1],
#  else /*# COULOMB<-2 */
#    define GCPARM /* empty */
#  endif /*#!COULOMB<-2 */

void initrspace(void) /****************************************** initrspace */
/*
   QQTAB hack: initerfc using global parms
   called after setss
*/
{
  int i,j;

  /* this code follows the structure of initss() */
#  ifdef WATER
#    if COULOMB<-2
#      error "WATER+COULOMB<-2 not implemented"
#    endif /*# COULOMB<-2 */
  if (waterchk.sp>=0) {
    siteinfo_t *si=spec[waterchk.sp]->si;

    waters.ss.HH.sgrid=setqq(&waters.ss.HH.tab,waters.ss.HH.qq,1);
    waters.ss.OH.sgrid=setqq(&waters.ss.OH.tab,waters.ss.OH.qq,1);
    waters.ss.OO.sgrid=setqq(&waters.ss.OO.tab,waters.ss.OO.qq,1);
#    ifdef WATERPLUS
    waters.ss.MH.sgrid=setqq(&waters.ss.MH.tab,waters.ss.MH.qq,1);
    waters.ss.OM.sgrid=setqq(&waters.ss.OM.tab,waters.ss.OM.qq,1);
    waters.ss.MM.sgrid=setqq(&waters.ss.MM.tab,waters.ss.MM.qq,1);
#    endif /*# WATERPLUS */
    waters.ss.dummy.rdf=NULL; }
#  endif /*# WATER */

#  ifdef POLAR
  loop (i,0,nsites) loopto (j,0,i) {
    int p1,p2;

    loop (p1,0,2) loop (p2,0,2) {
      sstab  [i][j].qqtab[p1][p2].sgrid=setqq(&sstab[i][j].qqtab[p1][p2].tab,  sstab[i][j].qqtab[p1][p2].qq,GCPARM 1);
      sstab  [j][i].qqtab[p2][p1].sgrid=sstab[i][j].qqtab[p1][p2].sgrid;
      sstab  [j][i].qqtab[p2][p1].tab=sstab[i][j].qqtab[p1][p2].tab;
      /* exclusions: */
      sstab  [i][j].qqtabx[p1][p2].sgrid=setqq(&sstab[i][j].qqtabx[p1][p2].tab,sstab[i][j].qqtab[p1][p2].qq,GCPARM 0);
      sstab  [j][i].qqtabx[p2][p1].sgrid=sstab[i][j].qqtabx[p1][p2].sgrid;
      sstab  [j][i].qqtabx[p2][p1].tab=sstab[i][j].qqtabx[p1][p2].tab;
#    if COULOMB<-2
#      undef GCPARM
/* GAUSSIANCHARGES: GCPARM=sigma_i,sigma_j */
#      define GCPARM sitedef[i].LJ[1].parm[SS_PARMS-1],sitedef[j].LJ[1].parm[SS_PARMS-1],
#    endif /*# COULOMB<-2 */
      if (n14==2) {
        sstab14[i][j].qqtab[p1][p2].sgrid=setqq(&sstab14[i][j].qqtab[p1][p2].tab,sstab14[i][j].qqtab[p1][p2].qq,GCPARM factor14);
        sstab14[j][i].qqtab[p2][p1].sgrid=sstab14[i][j].qqtab[p1][p2].sgrid;
        sstab14[j][i].qqtab[p2][p1].tab=sstab14[i][j].qqtab[p1][p2].tab; }
    }
  }

  /* DEBUG */
  if (option('v')&4) loop (i,0,nsites) loop (j,0,nsites) {
    int p1,p2;
    loop (p1,0,2) loop (p2,0,2) {
      prt("%d-%c %d-%c qq=%g sgrid=%g [0]:Au=%g Bu=%g Cu=%g"
#    if SPLINES==3
" Du=%g"
#    endif /*# SPLINES==3 */
          ,i,"cp"[p1],j,"cp"[p2]
          ,sstab[i][j].qqtab[p1][p2].qq, sstab[i][j].qqtab[p1][p2].sgrid
          ,sstab[i][j].qqtab[p1][p2].tab[0].Au
          ,sstab[i][j].qqtab[p1][p2].tab[0].Bu
          ,sstab[i][j].qqtab[p1][p2].tab[0].Cu
#    if SPLINES==3
          ,sstab[i][j].qqtab[p1][p2].tab[0].Du
#    endif /*# SPLINES==3 */
); } }

#  else /*# POLAR */
  loop (i,0,nsites) loopto (j,0,i) {
    sstab  [i][j].sgrid=setqq(&sstab[i][j].tab,sstab[i][j].qq,GCPARM 1);
    sstab  [j][i].sgrid=sstab[i][j].sgrid;
    sstab  [j][i].tab=sstab[i][j].tab;
    /* no exclusions - they are still implemented by exact_eru - no splines */
    if (n14==2) {
      sstab14[i][j].sgrid=setqq(&sstab14[i][j].tab,sstab14[i][j].qq,GCPARM factor14);
      sstab14[j][i].sgrid=sstab14[i][j].sgrid;
      sstab14[j][i].tab=sstab14[i][j].tab; }
  }
#  endif /*#!POLAR */
}

#else /*# QQTAB */

/* one spline all r-space terms */
void initrspace(void)
{
#  if defined(COULOMB)
  initerfc(option('v')&1 ? el.grid : -el.grid,
           el.minqq,
           box.cutoff+el.rplus,
           box.cutoff,
           option('v')&4?-el.alpha:el.alpha,
           el.rshift);
#  endif /*# defined(COULOMB) && COULOMB!=0 */
}
#endif /*#!QQTAB */

#ifdef SLAB
void initwall(void) /********************************************** initwall */
/* called if wall.is */
{
  int n,sp,ns,i,last=0;
  siteinfo_t *si;
  molecule_t *mn;

  underline("initialization of walls");
  prt("gravity = %g AA/ps2 (in the z-direction)",wall.g);
  loop (i,0,2)
    if (abs(wall.n)&(1<<i))
      prt("%sive wall %d at z = %.9g*Lz = %.9g AA",
          wall.n>0?"repuls":"attract", i,wall.z[i],box.L[2]*wall.z[i]);
    else {
      /* for the initializer */
      prt("no wall %i, z=%g set",i,i-1);
      wall.z[i]=i-1; }

  if (!wall.n) prt("no wall present");

  loop (n,0,No.N) {
    mn=molec+n;
    sp=mn->sp;
    ns=mn->ns;
    si=spec[sp]->si;

    loop (i,0,ns) if (si[i].st==nsites-1) last++; }

  if (last) WARNING(("The last atom describing the wall material\n\
*** appears %d-times in the configuration.",last))

  wall.numden=wall.rho/rhounit/sitedef[nsites-1].M; /* M is before equalization */
  prt("density of wall atoms %.9g kg/m3",wall.rho);
  prt("number density of wall atoms=%.9g/AA^3",wall.numden);

  rallocarrayzero(wsstab,nsites);
  loop (i,0,nsites) {
    sitesite_t *ss=sstab[nsites-1]+i;
    double sig6,eps4;

    if (wall.LJ[i].sig>-2e33) {
      if (option('v')&2) prt("The calculated wall-(%d) LJ sigma is replaced by %g from the data.",
                             i,wall.LJ[i].sig);
      ss->a.Sq=Sqr(wall.LJ[i].sig); }
    if (wall.LJ[i].neps>-2e33) {
      if (option('v')&2) prt("The calculated wall-(%d) LJ epsilon is replaced by %g from the data.",
                             i,wall.LJ[i].neps/wall.numden);
      ss->a.E4=4*wall.LJ[i].neps/wall.numden; }

    sig6=Cub(ss->a.Sq);
    eps4=wall.numden*ss->a.E4*sig6;

    if (eps4<0) ERROR(("implementation error: eps4=%g should not be negative",eps4))

    wsstab[i].A=eps4*sig6*(PI/45);
    wsstab[i].AA=eps4*sig6*(PI/5);
    if (wall.n<0) {
      wsstab[i].B=eps4*(PI/6);
      wsstab[i].BB=eps4*(PI/2); }
    else {
      wsstab[i].B=0;
      wsstab[i].BB=0; } }

  wall.minz=wall.z[0]*box.L[2]+1;
  wall.scalez=wall.z[1]-wall.z[0]-2/box.L[2];

  if ( (abs(wall.n)&2)==0 && wall.g>=0 )
    WARNING(("Missing top wall and no or wrong gravity.\n\
*** A molecule may escape from the box (followed by simulation crash)"))
  if ( (abs(wall.n)&1)==0 && wall.g<=0 )
    WARNING(("Missing bottom wall and no or wrong gravity.\n\
*** A molecule may escape from the box (followed by simulation crash)"))
  if ((abs(wall.n)&3)==3 && wall.z[0]>=wall.z[1])
    ERROR(("bad positions of walls (wall.z[0]>=wall.z[1])"))
}
#endif /*# SLAB */
