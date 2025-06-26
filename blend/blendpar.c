/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLENDPAR.C

This module defines the internal representation of force field parameters
(for atoms, bonds, angles, dihedral and improper torsions) and performs
reading and writing the parameter files.
See also intrapot.h !

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#include "ground.h"
#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "options.h"

int NATOMS;
struct atom_s *atom; /* [NATOMS] */

sitesite_t (**sstab)[2]; /* [site1][site2][onefour] */

int nbonds;
bond_t *bond0;

int nangles;
angle_t *angle0;

int ndihedrals,nimpropers,ncisdihedrals;
torsion_t *dihedral0,*improper0,*cisdihedral0;

int nnbfixes;
nbfix_t *nbfix0;
/* ... bond0, angle0, etc. are heads of linked lists */

int nwaters;
struct water_s *water0;

int bbdata[9][NBBTYPE];
/* wasting memory: BBN=1, BBCA=2, BBC=3, BBCO=8 */

#ifdef POLAR
int npolarbonds;
struct polarbond_s *polarbond0;

int npolarangles;
struct polarangle_s *polarangle0;

void scalepolarizabilities(double scalealpha,double scaleEsat)
{
  struct polarbond_s *pb;
  struct polarangle_s *pa;
  int i;

  loop (i,0,NATOMS) {
    atom[i].isotropicparm.alpha*=scalealpha;
    atom[i].isotropicparm.Esat*=scaleEsat; }

  for (pb=polarbond0; pb; pb=pb->next) {
    pb->parm.alpha*=scalealpha;
    pb->parm.alphazz*=scalealpha;
    pb->parm.Esat*=scaleEsat; }

  for (pa=polarangle0; pa; pa=pa->next) {
    pa->parm.alpha*=scalealpha;
    pa->parm.alphazz*=scalealpha;
    pa->parm.Esat*=scaleEsat; }

  if (scalealpha!=1 || scaleEsat!=1)
    prt("! polarizabilities rescaled %gx, saturation energies %gx",
	scalealpha,scaleEsat);
}
#endif /*# POLAR */

int all_dihedrals=1;
int ar_dih_limit=-1;  /* incl. limit for dihedrals in aromatic rings; -1=off */
int all_impropers=0;
int column_X=0;       /* disable/enable column X */
int comb_rule=0;      /* the combining rule, see XXX/sitesite.h */
double factor14=0.5;  /* to multiply 1-4 forces */
int distance14=4;     /* distance of 1-4 interactions
                         (value distance_14=5 means thus 1-5 interactions */
int polar=0;          /* polar=1 expected in the par-file for POLAR */
int ar_14_limit=0;    /* 3 for GROMOS (no 1-4 in or at aromatic rings) */

char *getnoncommentline(FILE *file, char *line) /********* getnoncommentline */
/***
    gets one `line' from `file', skipping comments (lines beginning with !)
    returns line on success, otherwise NULL
    unlike stupid function fgets(), line length is tested (max option('l'))
    and no LF is appended
    comment lines may be longer than option('l')
***/
{
  char *c=line, *cmax=line+option('l')-1;
  int i;

  line[0]=0;

 nextline:

  /* 1st character in the line */
  if ( (i=fgetc(file))<0 ) return NULL;

  if (i=='!') {
    /* skip comment */
    do {
      if ( (i=fgetc(file))<0 ) return NULL;
      } while (i!='\n');
    goto nextline; }

  while (i>=0 && i!='\n') {
    if (c>=cmax) {
      ERROR(("line too long: increase option -l"))
      return NULL; }
    if (i!='\r')
      *c++=i;
    i=fgetc(file); }

  *c=0;

  // fprintf(stderr,"%s\n",line);

  return line;
}

/* bit more efficient token parsing */
int ntoks,endtok; /* # of tokens found so far, end mark */
char *tok,*toksep=" \t\n";

char *firsttok(char *l) /****************************************** firsttok */
{
  char *t=strtok(l,toksep);

  ntoks=endtok=0;
  if (t) ntoks++; else endtok++;

  return tok=t;
}

char *nexttok(void)
{
  if (endtok) return NULL;
  {
    char *t=strtok(NULL,toksep);
    if (t) {
      if (t[0]=='!') { endtok++; return tok=NULL; }
      ntoks++; }
    else endtok++;

    return tok=t;  }
}

/* file,line,mygetline are local here */

static FILE *file;
static char *line;

static char *mygetline(void) /************************************** getline */
{
  return getnoncommentline(file,line);
}

static int find(char *st) /******************************************** find */
/***
    Searches for a line with given keyword in the parameter file `file'
    then skips comment lines (beginning with !) and reads the first
    non-comment line into `line'
    returns 1 on success
***/
{
  char *c,*aux,*auxx;
  int pass;
  int atoms;
  static int warn=1;

  loop (pass,0,2) {
    c="top";
    atoms=0;
    do {
      aux=c;
      if (!mygetline()) goto err;
      c=strtok(line," \t\n");
      if (c && !strcmp(c,"atoms")) atoms++;
      if (c && warn && atoms && !aux) {
        aux=auxx=NULL;
        aux=strtok(NULL," \t\n");
        auxx=strtok(NULL," \t\n");
        if (auxx) {
          prt("! WARNING: the par-file may contain data lines which are not part of any table\n\
!          the first such line starts with items: %s %s %s\n\
! EXPLANATION: since any blank line ends a table, an accidental blank line\n\
!              in the par-file may cause ignoring some data",c,aux,auxx);
          warn=0; } }
      } while (!c || strcmp(c,st));
    if (option('v')&2) fprintf(stderr,"reading %s\n",st);
    if (mygetline()) return 1;
 err:
    if (!pass) rewind(file); }
  prt("! table \"%s\" not found",st);

  return 0;
}

int atomn_ERROR=1;

int atomn(char *name) /******************************************** atomn */
/***
    returns atom number if the name is given
    not very efficient implementation
***/
{
  int i;

  if (!name) ERROR(("internal"))
  loop (i,0,NATOMS) if (!strcmp(name,atom[i].name)) return i;
  if (atomn_ERROR) ERROR(("%s is illegal atom name",name))

  return -1;
}

static char nm[4][8];

static void prtnm(int n) /******************************************** prtnm */
{
  int i;

  loop (i,0,n-1) prt_("%s - ",nm[i]);
  prt("%s",nm[n-1]);
}

static int inlist(char *list,char *item)
{
/*
   takes items from a comma-separated list, returns 1 on success, 0 if empty
   read items are returned in item and deleted from list
*/
  char *l,*i;

  if (!*list) return 0;

  for (l=list,i=item; *l && *l!=','; l++,i++) *i=*l;
  *i=0;
  i=list;
  if (*l) {
    l++;
    while (*l) *i++=*l++; }
  *i=0;

  return 1;
}

static int readtorsions(char *name,torsion_t **torsion0)
/* read dihedrals and impropers ******************************* readtorsions */
{
  int nt=0,n,i;
  torsion_t *torsion,*t;
  ireal K,angle;
  char list0[128],list1[128],list2[128],list3[128];

  if (find(name)) for (;;) {

    i=sscanf(line,"%s%s%s%s"irealfmt"%d"irealfmt,
                   list0,list1,list2,list3,
                            &K,&n,&angle);
    if (i<=0) break;
    if (i!=7) ERROR(("%s\n%s table: format",name,line))

    if (K!=0 || option('v')&4) {

      if (n<0 || n>MAXPER || (n>0 && fabs(angle-90)!=90)) {
        ERROR(("%s\n%s table: invalid n or angle",name,line)) }

      if (n>0) if (angle==180) K=-K;

      while (inlist(list0,nm[0])) {
        char l1[128];

        strcpy(l1,list1);
        while (inlist(l1,nm[1])) {
          char l2[128];

          strcpy(l2,list2);
          while (inlist(l2,nm[2])) {
            char l3[128];

            strcpy(l3,list3);
            while (inlist(l3,nm[3])) {

              /* for n==0, K[1]==angle needed */
              alloc(torsion,sizeof(torsion_t)+(n==0)*sizeof(ireal));

              torsion->parm.K[0]=K;
              if (n==0) torsion->parm.K[1]=angle;
              torsion->parm.n=n;
              loop (i,0,4) torsion->indx[i]=atomn(nm[i]);

#define MATCH(I,J) (torsion->indx[I]==t->indx[J])

              if (option('v')&4) for (t=*torsion0; t; t=t->next)

                if (name[0]=='i') {
                  /* improper checked */
                  if (MATCH(0,0) && MATCH(3,3))
                    if ( (MATCH(1,1) && MATCH(2,2))
                      || (MATCH(2,1) && MATCH(1,2)) )
                      if (torsion->parm.n==t->parm.n) {
                        prts_("! improper duplicated: "); prtnm(4); } }

                else {
                  /* dihedral checked */
                  if ( (MATCH(0,0) && MATCH(1,1) && MATCH(2,2) && MATCH(3,3))
                    || (MATCH(3,0) && MATCH(2,1) && MATCH(1,2) && MATCH(0,3)) )
                    if (torsion->parm.n==t->parm.n) {
                      prts_("! dihedral duplicated: "); prtnm(4); } }

#undef MATCH

              torsion->next=*torsion0;
              *torsion0=torsion;
              nt++; } } } }
      }
    mygetline(); }

  return nt;
}

char *PARfn;
char *blendpath;
static char *fullfn;

static char *fullname(char *name) /******************************* fullname */
{
  char *s;

  alloc(s,strlen(name)+strlen(blendpath)+2);
  strcpy(s,blendpath);
  if (*s)
    if (s[strlen(s)-1] != '/') strcat(s,"/");
  strcat(s,name);

  return s;
}

static void newPARfn(char *fn) /********************************** newPARfn */
/* allocates PARfn, cuts extension (must have extension) */
{
  int l=strlen(fn);

  if (PARfn)
    ERROR(("%s->%s: misplaced or duplicated .par|.bin file",PARfn,fn))
  alloc(PARfn,l+1);
  strcpy(PARfn,fn);
  PARfn[l-4]=0;
  if (!(blendpath=getenv("BLENDPATH"))) blendpath="";
  fullfn=fullname(fn);
  prt("! parameter file %s",fullfn);
}

static char parinfoline[80]; /* must be array, not pointer! */
static void checkpolarforcefield(void)
{
#ifdef POLAR
  if (!polar) {
    prt("! WARNING: nonpolar force field is used with POLAR BLEND");
    option('~')=0; }
#else /*# POLAR */
  if (polar)
    ERROR(("polar force field and nonpolar blend"))
#endif /*#!POLAR */
}

void makesstab(void) /******************************************* makesstab */
{
  int i,j,k;
  nbfix_t *fix;
  sitesite_t *ss;
  pairparm_t pp[2];

  /* prepare data for the combining rule */
  loop (j,0,NATOMS) loop (k,0,2) if (atom[j].LJ[k].RvdW!=0)
    initcombrule(&atom[j].LJ[k],comb_rule);

  /* allocate and build site-site pair table */
  alloc2Darray(sstab,NATOMS,NATOMS);
  loop (i,0,NATOMS) if (atom[i].name[0]) if (strcmp(atom[i].name,"X")) {
    loop (j,0,NATOMS)  if (atom[j].name[0]) if (strcmp(atom[j].name,"X")) {
      ss=sstab[i][j];

      /* first try to find the pair in the nbfixes table */
      for (fix=nbfix0; fix; fix=fix->next)
        if ( (fix->indx[0]==i && fix->indx[1]==j)
         ||  (fix->indx[0]==j && fix->indx[1]==i) ) {
          copyarray(pp,fix->onefour,2);
          goto fixed; }

      loop (k,0,2)
        combrule(&pp[k],
                 &atom[i].LJ[k],&atom[j].LJ[k],
                 option('v')&256?string("pot%s %s-%s",k?"14":"NB",atom[i].name,atom[j].name):NULL);

    fixed:
      loop (k,0,2) initssaux(&ss[k].a,&pp[k]);
  } }
}

static void checkpar(void) /************************************** checkpar */
{
  if (distance14<3 || distance14>6)
    ERROR(("distance14=%d is invalid",distance14))
  if (ar_14_limit) {
    if (distance14!=4)
      ERROR(("ar_14_limit not implemented for distance14=%d",distance14))
    if (ar_14_limit<2 || ar_14_limit>4)
      ERROR(("ar_14_limit invalid")) }
}

void readPAR(char *fn) /******************************************* readPAR */
/***
    Reads the parameter file fn (ASCII version)
    fn is the name with extension .par
    BLENDPATH (environment variable) prepended
***/
{
  int i,j,A,X,Z,natoms,always14,LJsigma=0;
  int sqrt_rule=(int)-2147483333; /* legacy */
  bond_t *bond,*b;
  char chem[2];
  angle_t *angle,*a;
  nbfix_t *nbfix,*fix;
  struct water_s *water;
#ifdef POLAR
  struct polarbond_s *polarbond,*pb;
  struct polarangle_s *polarangle,*pa;
#endif /*# POLAR */
  FILE *f;
  double Escaling=1,Rscaling=1;

  if (option('x')) {
    if (option('x')>0) {
      Escaling=(option('x')%1000)/100.;
      if (option('x')/1000) Rscaling=(option('x')/1000)/1000.; }
    else {
      Escaling=log(-option('x')/1000.);
      Rscaling=exp(Escaling/3);
      Escaling=exp(-Escaling); }
    if (Escaling!=1 || Rscaling!=1)
      prt("! Lennard-Jones parameters rescaled: Emin scaled %g times, Rmin %g times",Escaling,Rscaling); }

  newPARfn(fn);
  if ( !(file=fopen(fullfn,"rt")))
    ERROR(("file %s not found (check BLENDPATH)",fullfn))
  free(fullfn);

  alloc(line,option('l'));

/* first pass over atoms to get the size */
/*  fgets(parinfoline,80,file); */
/* read longer and truncate in 2nd pass */
  fgets(line,option('l'),file);

  if (!find("atoms")) ERROR(("no \"atoms\" table"))
  NATOMS=0;

  for (;;) {
    if (sscanf(line,"%d",&j)!=1) break;
    Max(NATOMS,j)
    mygetline(); }

  prt("! max atom number = %d",NATOMS);
  NATOMS++;

  alloczero(atom,NATOMS*sizeof(struct atom_s));
  strcpy(atom[0].name,"X"); /* wildcard */
  natoms=0;

  rewind(file);
/*.....  fgets(parinfoline,80,file);*/
  fgets(line,option('l'),file);
  memcpy(parinfoline,line,80);
  parinfoline[79]=0;
  prt(parinfoline);

  f=in; in=file;

  if (distance14<3 || distance14>6)
    ERROR(("distance14=%d is invalid",distance14))

  j=_Id.echo; _Id.echo=0;
  getdata
    get(all_dihedrals) get(ar_dih_limit) get(all_impropers) get(ar_14_limit)
    get(factor14) get(distance14) get(column_X)
    get(sqrt_rule)
    get(comb_rule)
    get(LJsigma) get(polar)
    checkdata 
  enddata
  _Id.echo=j;
  in=f;

  if (sqrt_rule!=(int)-2147483333) {
    prt("\
! NOTE: the force field uses parameter sqrt_rule, which has been\n\
!       renamed to comb_rule. I am assuming comb_rule=sqrt_rule.");
    comb_rule=sqrt_rule; }

  if (option('|')!=-9) polar=option('|'); /* see #define Q in blend.c */

  checkpolarforcefield();

  if (!find("atoms")) ERROR(("no \"atoms\" table"))

  for (;;) {
    if (sscanf(line,"%d",&j)!=1) break;
    if (j>=NATOMS || j<=0)
      ERROR(("atom number %d is out of range",j))
    if (atom[j].name[0])
      ERROR(("%s\natom %d %s duplicated in the table",line,j,atom[j].name))
    if (sscanf(line,"%d%6s%d%d%lf%d%1s",
                     &j,
                       atom[j].name,
                         &X,
                           &A,
                             &atom[j].mass,
                               &Z,
                                chem) !=7 )
                                  ERROR(("%s\natoms: bad format",line))

    atom[j].Z=Z;
    if (Z<0) ERROR(("%s: Z<0",atom[j].name))
    atom[j].chem=chem[0];
    atom[j].X=X; atom[j].A=A;
    if (!column_X) atom[j].X=1; /* always match if option('x')=column_X=0 */
    mygetline();
    natoms++; }

  /* SS_TABLE is #defined in sim/XXX/sitesite.h */
  if (!find(SS_TABLE)) ERROR(("no \""SS_TABLE"\" table"))
  /* wildcards only of type X* */

  for (;;) {
    siteparm_t LJ[2];
    int no14;

    if (firsttok(line)) strcpy(nm[0],tok);
    loop (i,0,2) {
      if (!nexttok()) goto exss; LJ[i].alpha=atof(tok);
      if (!nexttok()) goto exss; LJ[i].EvdW=atof(tok)*Escaling;
      if (!nexttok()) goto exss; LJ[i].RvdW=atof(tok)*Rscaling;
#if SS_PARMS
      { 
        int j; 
        loop (j,0,SS_PARMS) {
          if (!nexttok()) goto exss; LJ[i].parm[j]=atof(tok); }
      }
#endif /*# SS_PARMS */
    }

  exss:

#if 0 /*OLD*/
    n=sscanf(line,"%s%lf%lf%lf%lf%lf%lf",
                   nm[0],
                      &LJ.alpha,&LJ.EvdW,&LJ.RvdW,
                               &LJ14.alpha,&LJ14.EvdW,&LJ14.RvdW);
#endif /*# 0 */
    if (ntoks<=0) break;
    if (ntoks!=4+SS_PARMS && ntoks!=7+2*SS_PARMS)
      ERROR(("%s\n"SS_TABLE" table: format",line))
    no14=ntoks!=7+2*SS_PARMS;
    if (no14) {
      LJ[1]=LJ[0];
      if (factor14!=1.0) {
        LJ[1].EvdW*=factor14; no14=0; } }

    if (LJsigma) loop (i,0,2) LJ[i].RvdW*=0.5612310241546865;

    /* new in 2.0p */
    loop (i,0,2) if (LJ[i].EvdW==0 && LJ[i].RvdW==0) {
      prt("! WARNING: %s 14=%d: empty LJ changed into Emin=0, Rmin=1",nm,i);
      LJ[i].RvdW=1; }

    if (nm[0][1]=='*') {
      /* fill all nm[0][0]*** atoms */
      loop (i,0,NATOMS) if (atom[i].name[0]==nm[0][0]) {
        atom[i].LJ[0]=LJ[0]; atom[i].LJ[1]=LJ[1]; atom[i].no14=no14; } }
    else {
      i=atomn(nm[0]);
      atom[i].LJ[0]=LJ[0]; atom[i].LJ[1]=LJ[1]; atom[i].no14=no14; }
    mygetline(); }

#ifdef POLAR /* ///////////////////////////////////////////////// POLAR PAR */

  if (find("polaratoms")) for (;;) {
    /* wildcards only of type X* */
    isotropicparm_t ip;
    int n;

    ip.rep=0; /* cf. table "shellrep" */
    n=sscanf(line,"%s"irealfmt irealfmt irealfmt irealfmt,
                   nm[0],&ip.alpha,&ip.shell,&ip.Esat,&ip.kappa);

    if (n<=0) break;

    if (n!=5)
      ERROR(("%s\npolaratoms table: format",line))

    if (nm[0][1]=='*') {
      /* fill all nm[0][0]*** atoms */
      loop (i,0,NATOMS) if (atom[i].name[0]==nm[0][0])
        atom[i].isotropicparm=ip; }
    else {
      i=atomn(nm[0]);
      atom[i].isotropicparm=ip; }
    mygetline(); }

  npolarbonds=0;
  if (find("polarbonds")) for (;;) {
    alloc(polarbond,sizeof(struct polarbond_s));
    i=sscanf(line,"%s%s"irealfmt irealfmt irealfmt,
                   nm[0],nm[1],
                        &polarbond->parm.alpha,
                                 &polarbond->parm.alphazz,
                                          &polarbond->parm.Esat);
    if (i<=0) { free(polarbond); break; }
    if (i!=5) ERROR(("table polarbonds: %s",line))
    loop (i,0,2) polarbond->indx[i]=atomn(nm[i]);

#  define MATCH(I,J) (polarbond->indx[I]==pb->indx[J])

    for (pb=polarbond0; pb; pb=pb->next)
      if (MATCH(0,0) && MATCH(1,1)) {
        prtnm(3);
        ERROR(("polarbond duplicated in the table")) }

#  undef MATCH

    polarbond->next=polarbond0;
    polarbond0=polarbond;
    npolarbonds++;
    mygetline(); }

  npolarangles=0;
  if (find("polarangles")) for (;;) {
    alloc(polarangle,sizeof(struct polarangle_s));
    i=sscanf(line,"%s%s%s"irealfmt irealfmt irealfmt,
                   nm[0],nm[1],nm[2],
                          &polarangle->parm.alpha,
                                   &polarangle->parm.alphazz,
                                            &polarangle->parm.Esat);
    if (i<=0) { free(polarangle); break; }
    if (i!=6) ERROR(("table polarangles: %s",line))
    loop (i,0,3) polarangle->indx[i]=atomn(nm[i]);

#  define MATCH(I,J) (polarangle->indx[I]==pa->indx[J])

    for (pa=polarangle0; pa; pa=pa->next)
      if (MATCH(0,0) && MATCH(1,1) && MATCH(2,2)) {
        prtnm(3);
        ERROR(("polarangle duplicated in the table")) }

#  undef MATCH

    polarangle->next=polarangle0;
    polarangle0=polarangle;
    npolarangles++;
    mygetline(); }

  if (find("shellrep")) {
    int n=0;
    for (;;) {
      i=sscanf(line,"%s",nm[0]);
      if (i<=0) break;
      if (i!=1) ERROR(("table shellrep: %s",line))
      i=atomn(nm[0]);
      atom[i].isotropicparm.rep=1;
      n++;
      mygetline(); }
    prt("! table \"shellrep\": %d repulsive counterparts (cations) of the shell-core terms",n); }

#endif /* ///////////////////////////////////////////////// end of POLAR PAR */ /*# POLAR */


  /* not given 1-4 term in the NBFIX table means factor14 times normal NBFIX */
  /* zero sigma 1-4 NBFIX means using combining rule (for 1-4) */
  nnbfixes=0;
  always14=0; /* set to 1 if 1-4 NBfixes ==> 1-4 interactions are
                 ALWAYS included (though it may happen that they are
               identical to the normal interactions) */

  if (find("nbfixes") || find("NBFIX")) for (;;) {
    /* in (most) LJ force fields:
         epsvdW = min energy, sigvdW = van der Waals diameter
       other force field generally:
         epsvdW = any (usu. energy-like) parameter
         sigvdW = any (usu. size-related) parameter
    */
    ireal epsvdW[2],sigvdW[2];
#if SS_PARMS
    ireal parm[2][SS_PARMS];
#endif /*# SS_PARMS */

    if (firsttok(line)) strcpy(nm[0],tok);
    if (nexttok()) strcpy(nm[1],tok);
    loop (i,0,2) {
      if (!nexttok()) goto ex; epsvdW[i]=atof(tok);
      if (!nexttok()) goto ex; sigvdW[i]=atof(tok);
#if SS_PARMS
      { 
        int j; 
        loop (j,0,SS_PARMS) {
          if (!nexttok()) goto ex; parm[i][j]=atof(tok); }
      }
#endif /*# SS_PARMS */
    }

  ex:

    if (ntoks<=0) break;
    if (ntoks!=4+SS_PARMS && ntoks!=6+2*SS_PARMS)
      ERROR(("%s\nNBFIX table: format",line))

    if (ntoks==6+2*SS_PARMS)
      always14=1; /* this is just for warning */
    else {
      /* missing 1-4 params ==> factor14*normal params */
      epsvdW[1]=epsvdW[0]*factor14;
      sigvdW[1]=sigvdW[0];
#if SS_PARMS
      { int j; loop (j,0,SS_PARMS) parm[1][j]=parm[0][j]; }
#endif /*# SS_PARMS */
    }

    alloc(nbfix,sizeof(nbfix_t));
    loop (i,0,2) {
      nbfix->indx[i]=atomn(nm[i]);
      if (LJsigma) sigvdW[i]*=1.122462048309373;
      nbfix->onefour[i].sig=sigvdW[i];
      nbfix->onefour[i].eps=epsvdW[i];
#if SS_PARMS
      { 
        int j; 
        loop (j,0,SS_PARMS) nbfix->onefour[i].parm[j]=parm[i][j]; 
      }
#endif /*# SS_PARMS */
    }

    /* auxiliary values to be used in sspot */
    loop (i,0,2)
      if (nbfix->onefour[i].sig<=0) {
        /* zero RvdW : use normal combining rule */
        WARNING(("pair %d-%d: nbfix contains RvdW<=0, combining rule used",
                 nbfix->indx[0],nbfix->indx[1]))
        combrule(&nbfix->onefour[i],
                 &atom[nbfix->indx[0]].LJ[i],
                 &atom[nbfix->indx[1]].LJ[i],
                 option('v')&256?string("fix%s %s-%s",i?"14":"NB",atom[nbfix->indx[0]].name,atom[nbfix->indx[1]].name):NULL); }

#define MATCH(I,J) (nbfix->indx[I]==fix->indx[J])
    for (fix=nbfix0; fix; fix=fix->next)
      if ( (MATCH(0,0) && MATCH(1,1)) || (MATCH(0,1) && MATCH(1,0)) ) {
        prtnm(2);
        ERROR(("NBFIX duplicated in the table")) }
#undef MATCH

    nbfix->next=nbfix0;
    nbfix0=nbfix;
    nnbfixes++;
    mygetline(); }

  nbonds=0;
  if (find("bonds")) for (;;) {
    bondparm_t parm;
    char list0[128],list1[128];

    i=sscanf(line,"%s%s"irealfmt irealfmt,
                   list0,list1,
                        &parm.K,&parm.length);
    parm.Ki2=-2*parm.K;
    if (i<=0) break;
    if (i!=4) ERROR(("%s\nbonds table: format",line))

    while (inlist(list0,nm[0])) {
      char l1[128];

      strcpy(l1,list1);
      while (inlist(l1,nm[1])) {

        alloc(bond,sizeof(bond_t));
        bond->parm=parm;
        loop (i,0,2) bond->indx[i]=atomn(nm[i]);

#define MATCH(I,J) (bond->indx[I]==b->indx[J])
        for (b=bond0; b; b=b->next)
          if ( (MATCH(0,0) && MATCH(1,1)) || (MATCH(0,1) && MATCH(1,0)) ) {
            prtnm(2);
            ERROR(("bond duplicated in the table")) }
#undef MATCH

        bond->next=bond0;
        bond0=bond;
        nbonds++; } }
    mygetline(); }


  nangles=0;
  if (find("angles")) for (;;) {
    angleparm_t parm;
    char list0[128],list1[128],list2[128];

#ifdef UREY_BRADLEY
    i=sscanf(line,"%s%s%s"irealfmt irealfmt irealfmt irealfmt,
                   list0,list1,list2,
                          &parm.K, &parm.angle,
                          &parm.Kub, &parm.length);
    parm.K2=parm.K*2;
    if (i<=0) break;
    if (i==5) { parm.Kub=parm.length=0; i=7; }

    if (i!=7) ERROR(("%s\nangles table: format",line))
#else /*# UREY_BRADLEY */
    i=sscanf(line,"%s%s%s"irealfmt irealfmt,
                   list0,list1,list2,
                          &parm.K, &parm.angle);
    if (i<=0) break;
    if (i!=5) ERROR(("angles: %s",line))
#endif /*#!UREY_BRADLEY */

    while (inlist(list0,nm[0])) {
      char l1[128];

      strcpy(l1,list1);
      while (inlist(l1,nm[1])) {
        char l2[128];

        strcpy(l2,list2);
        while (inlist(l2,nm[2])) {

          alloc(angle,sizeof(angle_t));
          angle->parm=parm;
          loop (i,0,3) angle->indx[i]=atomn(nm[i]);

#define MATCH(I,J) (angle->indx[I]==a->indx[J])

          for (a=angle0; a; a=a->next)
            if MATCH(1,1)
              if ( (MATCH(0,0) && MATCH(2,2))
                || (MATCH(0,2) && MATCH(2,0)) ) {
                prtnm(3);
                ERROR(("angle duplicated in the table")) }

#undef MATCH

          angle->next=angle0;
          angle0=angle;
          nangles++; } } }
    mygetline(); }

  ndihedrals=readtorsions("dihedrals",&dihedral0);

  if (0) {// DEBUG
    torsion_t *d;
    looplist (d,dihedral0)
      prt("%d-%d-%d-%d %d %g",
          d->indx[0],d->indx[1],d->indx[2],d->indx[3],
          d->parm.n,d->parm.K[0]);
  }
  
  nimpropers=readtorsions("impropers",&improper0);
  ncisdihedrals=readtorsions("cisdihedrals",&cisdihedral0);

  nwaters=0;
  if (find("waters")) for (;;) {
    char id[MAXWATERSITES][8];

    alloc(water,sizeof(struct water_s));
#if MAXWATERSITES!=4
#  error MAXWATERSITES changed, change the following statement:
#endif /*# MAXWATERSITES!=4 */
    /* this is NOT fool proof! short names (<=7 char) expected... */
    i=sscanf(line,"%s%d%s%lf%s%lf%s%lf%s%lf",
                   water->name,&water->ns,
                   id[0],&water->charge[0],
                   id[1],&water->charge[1],
                   id[2],&water->charge[2],
                   id[3],&water->charge[3]);
    if (i<=0) { free(water); break; }
    if (water->ns>MAXWATERSITES)
      ERROR(("%s\nwaters table: too many sites",line))
    if (i<2+2*water->ns) ERROR(("%s\nwaters table: format",line))
    loop (i,0,water->ns) water->type[i]=atomn(id[i]);

    water->next=water0;
    water0=water;
    nwaters++;
    mygetline(); }

  if (find("backbone")) loop (i,1,5) {
    char *tok=strtok(line," \t\n");
    int j=0;

    if (i==4) i=8; /* =BBCO (unknown here, which is mess) */

    if (!tok) ERROR(("backbone table: empty line"))

    while (tok) {
      if (j>=NBBTYPE) ERROR(("backbone table %d: too many atom types",i))
      bbdata[i][j++]=atomn(tok);

      tok=strtok(NULL," \t\n"); }

    while (++j<NBBTYPE) bbdata[i][j]=bbdata[i][j-1];

    mygetline(); }

  fclose(file);
  free(line);

  prt("! %d atoms   %d nbfixes\n\
! %d bonds  %d angles  %d dihedrals  %d impropers  %d cisdihedrals\n",
         natoms,    nnbfixes,
  nbonds,   nangles,   ndihedrals,   nimpropers,   ncisdihedrals);
  if (always14) {
    prt("! 1-%d interactions always included (may be inefficient)",distance14);
    loop (j,0,NATOMS) atom[j].no14=0; }

  if (PCHname) readPCH(PCHname);

  checkpar();
} /* readPAR */

static FILE *bin;
static size_t (*rw)(void*,size_t,size_t,FILE*);
static char *msg;
static int READ;

static void rwtab(int *n,any_t **a0,int len) /************************ rwtab */
{
  int i;
  any_t *a;

  if (rw(n,sizeof(int),1,bin)!=1 || *n<0 || *n>16000) ERROR((msg))
  if (READ) loop (i,0,*n) {
    alloc(a,len);
    a->next=*a0;
   *a0=a; }
  len-=sizeof(any_t*); /* WARNING this is not 100% portable ! */
  for (a=*a0; a; a=a->next)
    if (rw(a->indx,len,1,bin)!=1) ERROR((msg))
}

int rwBIN(char *fn,char *ext) /************************************** rwBIN */
/***
    read/write binary image of parameter file
    ext==NULL : read
    ext== position of `.' of extension in fn : write
    read: returns 1 on success
***/
{
  char version[4];
  memcpy(version,VERSION,4);

  READ=!ext;

  if (READ) {
    rw=fread;
    newPARfn(fn);
    bin=fopen(fullfn,"rb");
    free(fullfn);
    if (bin==NULL) return 0;
    msg="read .bin file"; }
  else {
    char bfn[128];

    strcpy(bfn,fn);
    strcpy(bfn+(ext-fn),".bin");
    prt("! writing %s",bfn);
    if (fopen(bfn,"rb")) WARNING(("file %s already exists",bfn))
    bin=fopen(bfn,"wb");
    if (bin==NULL) ERROR(("cannot write to %s",bfn))
    /* casting needed because fwrite has const void* 1st arg */
    rw=(size_t (*)(void*,size_t,size_t,FILE*))fwrite;
    msg="write .bin file"; }

#define RW(X) if (rw(&X,sizeof(X),1,bin)!=1) ERROR((msg))
  RW(version)
  if (memcmp(version,VERSION,3))
    ERROR(("%s has bad version %.4s (expected %c.%c+)\n\
*** generate the new bin-files by `blend -b' (or erase the bin-files)",
           fn,version,VERSION[0],VERSION[2]))
  RW(parinfoline)
  if (READ) prt(parinfoline);
  RW(all_dihedrals)
  RW(ar_dih_limit)
  RW(all_impropers)
  RW(ar_14_limit)
  RW(factor14)
  RW(distance14)
  RW(column_X)
  RW(comb_rule)
  RW(polar)
  checkpolarforcefield();
  RW(NATOMS)
  /* only some tests ... */
  if (column_X<0 || column_X>1)
    ERROR(("bad column_X in %s (damaged binary file or wrong architecture?)",fn))
  if (NATOMS<1 || NATOMS>1024)
    ERROR(("bad NATOMS in %s (damaged binary file or wrong architecture?)",fn))

  if (READ) alloc(atom,NATOMS*sizeof(struct atom_s));
  if (rw(atom,sizeof(struct atom_s),NATOMS,bin)!=NATOMS) ERROR((msg))
  rwtab(&nnbfixes,(any_t**)&nbfix0,sizeof(nbfix_t));

#ifdef POLAR
  rwtab(&npolarbonds,(any_t**)&polarbond0,sizeof(struct polarbond_s));
  rwtab(&npolarangles,(any_t**)&polarangle0,sizeof(struct polarangle_s));
#endif /*# POLAR */

  rwtab(&nbonds,(any_t**)&bond0,sizeof(bond_t));
  rwtab(&nangles,(any_t**)&angle0,sizeof(angle_t));
/***
!!! pessimistic length for dihedrals - in most cases will be shorter by ireal
***/
  rwtab(&ndihedrals,(any_t**)&dihedral0,sizeof(torsion_t)+sizeof(ireal));
  rwtab(&nimpropers,(any_t**)&improper0,sizeof(torsion_t)+sizeof(ireal));
  rwtab(&ncisdihedrals,(any_t**)&cisdihedral0,sizeof(torsion_t)+sizeof(ireal));
  rwtab(&nwaters,(any_t**)&water0,sizeof(struct water_s));

  RW(bbdata);

  fclose(bin);
  prt("! atom#<%d  %d nbfixes\n\
!  %d bonds  %d angles  %d dihedrals  %d impropers %d cisdihedrals\n",
         NATOMS,   nnbfixes,
      nbonds,nangles,   ndihedrals,   nimpropers,  ncisdihedrals);

  if (READ) {
    if (PCHname) readPCH(PCHname);
    /* makenbfixtab(); */ }

  return 1;
}

void readPARBIN(char *fn) /************************************* readPARBIN */
/***
    reads fn.bin; if it does not exist, then fn.par
    re-read checked
***/
{
  char *c,*rfn;

  alloc(rfn,strlen(fn)+5);
  strcpy(rfn,fn);
  c=rfn+strlen(fn);

  if (PARfn) {
    /* param file already read in */
    if (strcmp(rfn,PARfn))
/*.....    ERROR(("parameter set change (%s->%s) is illegal",rfn,PARfn))*/
      prt("! WARNING: parameter set change to %s not accepted, %s used",rfn,PARfn);
    goto freerfn; }

  strcat(rfn,".bin");

  if (rwBIN(rfn,NULL)) goto freerfn;

  free(PARfn);
  strcpy(c,".par");
  readPAR(rfn);
 freerfn: free(rfn);

  checkpar();

  makesstab();
}

char *atom_f;
char *NBFIXexc_f;
char *dihedral_f;
char *torsion_f;
char *z_f;
char *bond_f;
char *angle_f;
char *angleUB_f;
char *depend_f;
char *outtext_f;
char *charge_f;
char *infoxyz_f;
char *sites_f;
char *infor_f;
char *infoa_f;
char *infod_f;
#ifdef POLAR
char *polar_f;
char *axial_f;
#endif /*# POLAR */

char *fill1prec(char *template) /******************************** fill1prec */
/*
  adds precision (option -_ in abs. value) to format template
  in template, %n.mF (n,m=digits) is replaced by %N.Mf, where
    M=m+precision-2
    N=n+precision-2
  example: %4.1F -> %10.7f for -_6
*/
{
  char fmt[80];
  char *t,*f;

  for (t=template,f=fmt; *t; ) {
    *f++=*t++;
    if (t[-1]=='%') {
      char *F=strchr(t,'F');

      if (F-t==3 && t[1]=='.') {
        int n=t[0]+abs(option('_'))-2;
        int m=t[2]+abs(option('_'))-2;

        if (n>'9') *f++='1',n-=10; *f++=n;
        *f++='.';
        if (m>'9') *f++='1',m-=10; *f++=m;
        *f++='f';
        t+=4; } } }

  *f=0;

  return strdup(fmt);
}

void fillprec(void) /********************************************* fillprec */
{
  atom_f=fill1prec(" %6.2F %7.3F %5.2F");
  NBFIXexc_f=fill1prec(" %7.3F %5.2F");
  dihedral_f=fill1prec("%1d %2d %5.2F %7.2f %7.4f");
  torsion_f=fill1prec("%1d %2d %3.1F %7.2f %7.4f");
  z_f=fill1prec(" %6.2F");
  bond_f=fill1prec(" %4.1F  %5.2F  %6.3f  %8.3f");
  angle_f=fill1prec( "%5.1F %8.1F   %7.2f  %8.3f");
  angleUB_f=fill1prec( "%5.1F  %8.1F   %7.2f  %8.3f %5.1F %5.2F");
  depend_f=fill1prec("%7.4F");
  outtext_f=fill1prec("%0.3F");
  charge_f=fill1prec(" %5.2F ");
  infoxyz_f=fill1prec("%d %s %s %.4fe (%0.2F,%0.2F,%0.2F)");
  sites_f=fill1prec("%6.2F %6.2F %6.2F ");
  infor_f=fill1prec(" r[%d]=%0.2F");
  infoa_f=fill1prec(" a[%d]=%0.0F");
  infod_f=fill1prec(" d[%d]=%0.0F");
#ifdef POLAR
  polar_f=fill1prec(" %5.2F %2.0F %5.2F %5.2F %d");
  axial_f=fill1prec(" %6.3F %5.2F %5.2F %4.1F");
#endif /*# POLAR */
}

/*** partial charges ***/

char *PCHname,*REAname;

struct partial_s *partial0;

void readPCH(char *fn) /******************************************** readPCH */
{
  FILE *f;
  char *line,*tok,*eq;
  struct partial_s *partial,*plast;
  char *pathfn=NULL;

  if (!(f=fopen(fn,"rt")))
    if (blendpath) {
      prt("! partial charge file %s not found, trying BLENDPATH",fn);
      fn=pathfn=fullname(fn);
      f=fopen(fn,"rt"); }

  if (!f) {
    ERROR(("partial charge file %s not found",fn))
    goto ret; }

  prt("! partial charge file %s",fn);

  if (partial0) ERROR(("internal"))

  alloc(line,option('l'));

  while (fgets(line,option('l'),f)) if (line[0]!='!') {
    if (!(tok=strtok(line," \t\n"))) continue;

    if (!(eq=strchr(tok,'='))) {
      ERROR(("%s: bad syntax in %s",tok,fn))
      continue; }
    *eq=0;

    alloc(partial,sizeof(struct partial_s));
    if (partial0) plast->next=partial; else partial0=partial;
    plast=partial; partial->next=NULL;

    partial->center=atomn(tok);
    partial->centerq=atof(eq+1);

    partial->nnbr=0;
    while ( (tok=strtok(NULL," \t\n")) ) {
      if (tok[0]=='!') break;
      if (partial->nnbr>=MAXVAL) {
	ERROR(("%s: too many neighbors in %s",tok,fn))
	break; }
      if (!(eq=strchr(tok,'='))) {
	ERROR(("%s: bad syntax in %s",tok,fn))
        continue; }
      *eq=0;
      partial->nbr[partial->nnbr]=atomn(tok);
      partial->q[partial->nnbr]=atof(eq+1);
      partial->nnbr++; }

#if 0
    {
      int KO,n;

      do {
	loop (n,1,partial->nnbr) if (partial->nbr[n-1]>partial->nbr[n]) {
	  int sw=partial->nbr[n-1];
	  partial->nbr[n-1]=partial->nbr[n];
	  partial->nbr[n]=sw;
	  KO++; }
        } while (KO);
      }
#endif /*# 0 */

    }

  free(line);
  fclose(f);
 ret: if (pathfn) free(pathfn);
}

void scalebonds(void) /****************************************** scalebonds */
/*
  if option -b < 0, sets all bond constants to minus the option
  if option -b > 0, multiplies all bond constants by the option %
  NB: U=K*(r-r0)^2
*/
{
  bond_t *b;

  /* NB: -b1 means binary image .bin from .par (deprecated) */
  if (option('b')==100 || option('b')==1) return;

  if (option('b')<=0) looplist (b,bond0) b->parm.K=-option('b');
  else looplist (b,bond0) b->parm.K*=option('b')*0.01;

  looplist (b,bond0) b->parm.Ki2=-b->parm.K*2;
}

void scaleangles(void) /**************************************** scaleangles */
/*
  if option -a < 0, sets all angle constants to minus the option
  if option -a > 0, multiplies all angle constants by the option %
  NB: U=K*(alpha-alpha0)^2
*/
{
  angle_t *a;

  if (option('a')==100) return;

  if (option('a')<=0) looplist (a,angle0) a->parm.K=-option('a');
  else looplist (a,angle0) a->parm.K*=option('a')*0.01;

  looplist (a,angle0) a->parm.K2=a->parm.K*2;
}
