/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% blendgen.c

  This module assigns parameters of atoms, bonds, angles, dihedrals and
  impropers to the molecules and generates input file for `cook' and `cooks'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#include "ground.h"
#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "blendmed.h"
#include "blendgen.h"
#include "blendimp.h"
#include "blendmin.h"
#include "blendedt.h"
#include "blenddep.h"
#include "options.h"
#include "rndgen.h"

extern int measure;

int nsites,clust;
static int ns;       /* = spec->ns */
static char *specfn; /* = spec->fn */
site_t *site;    /* = spec->site ==> blendmin */

void listofsites(void) /************************************ listofsites */
/***
  list of all used atom (site) types
***/
{
  species_t *spec;
  int i,j,split=0,onefour;
  char *used;
  nbfix_t *fix;

  nsites=nnbfixes=0;

  looplist (spec,spec0) {
    split+=spec->N<0;

    /* mark used atoms and calculate nsites */
    loop (i,0,spec->ns) {
      used=&(atom[spec->site[i].type].used);
      if (!(*used)) {
        nsites++; *used=1; } } }

  /* NBFIX exceptions */
  loop (i,0,NATOMS) if (atom[i].used)
    loopto (j,0,i) if (atom[j].used)
      looplist (fix,nbfix0)
        if ( (fix->indx[0]==i && fix->indx[1]==j)
          || (fix->indx[1]==i && fix->indx[0]==j) ) { nnbfixes++; break; }

  /* warning: if molecule is split, field nspec= will be rewritten! */
  if (distance14!=4) prt("! WARNING 1-%d interactions instead of 1-4",distance14);
  prt("nspec=%d%s\n\
  nsites=%d  nnbfixes=%d  factor14=%f  distance14=%d  comb_rule=%d\n\
  polar=%d  nparms=%d;",
      nspec,split?" !?":"",nsites,nnbfixes,factor14,distance14,comb_rule,
      polar,SS_PARMS);

  combruleinfo(comb_rule);

  if (option('v')&1) {

    prts_("\n"SS_TABLE"\n\
! i atom    mass      alpha     EvdW     RvdW");
    loop (i,0,SS_PARMS) prts_(" parm");
    prt_(" alpha[1-%d] EvdW[1-%d] RvdW[1-%d]", distance14,distance14,distance14);
    loop (i,0,SS_PARMS) prts_(" parm");
    prt(" alphapol shell Esat arep rep");

    loop (i,0,NATOMS) if (atom[i].used) {
      prt_("%3d %-4s %9.5f", i,atom[i].name,atom[i].mass);
      loop (j,0,2) {
        prt_(atom_f,
             atom[i].LJ[j].alpha*!(comb_rule&1), /* =0 if comb_rule for eps */
           atom[i].LJ[j].EvdW, atom[i].LJ[j].RvdW);
#if SS_PARMS
        { int k; loop (k,0,SS_PARMS) prt_(z_f,atom[i].LJ[j].parm[k]); }
#endif
      } /* j */
#ifdef POLAR
      prt_(polar_f,
           atom[i].isotropicparm.alpha,
           atom[i].isotropicparm.shell*(option('@')/100.),
           atom[i].isotropicparm.Esat,
           atom[i].isotropicparm.kappa,
         atom[i].isotropicparm.rep);
#endif
      _n }

    /* old name: NBFIX */
    prt("\nnbfixes\n\
!atom  atom    Emin     diamvdW     Emin[1-4]     diamvdW[1-4]");

    /* NBFIX exceptions */
    loop (i,0,NATOMS) if (atom[i].used)
      loopto (j,0,i) if (atom[j].used)
        looplist (fix,nbfix0)
          if ( (fix->indx[0]==i && fix->indx[1]==j)
               || (fix->indx[1]==i && fix->indx[0]==j) ) {
            int l;

            prt_("%-4s %-4s", atom[i].name,atom[j].name);
            loop (l,0,2) {
              prt_(NBFIXexc_f, fix->onefour[l].eps, fix->onefour[l].sig);
#if SS_PARMS
              { int k; loop (k,0,SS_PARMS) prt_(z_f,fix->onefour[l].parm[k]); }
#endif
            }
            _n
            break; }
  }

  /* old version re-implemented by T.Trnka */
  if (option('v')&8) loop (onefour,0,2) {
    /* print combined LJ energy terms */
    prt("\n! %s site-site parameters (LJ or similar):\n\
! f denotes a nbfix replacing the previous combining rule value\n\
!  atom : atom   Emin[kcal/mol]  diamvdW  [parm..]    LJeps[kJ/mol]  LJsigma",
        onefour?"1-4":"normal");
    loop (i,0,NATOMS) if (atom[i].used)
      loopto (j,0,i) if (atom[j].used) {
        siteparm_t *lj0=&atom[i].LJ[onefour];
        siteparm_t *lj1=&atom[j].LJ[onefour];
        pairparm_t pp;
        char key=' ';
        nbfix_t *fix;

        combrule(&pp, lj0, lj1, NULL);

        do {
          prt_("!%c %4s : %-4s", key, atom[i].name, atom[j].name);
          prt_(NBFIXexc_f, pp.eps, pp.sig);
#if SS_PARMS
          { int i; loop (i,0,SS_PARMS) prt_(z_f,pp.parm[i]); }
#endif
          prt_("    ");
          prt_(NBFIXexc_f, -pp.eps*4.184, pp.sig*0.8908987181403393);
          _n

          looplist (fix,nbfix0)
            if ( (fix->indx[0]==i && fix->indx[1]==j)
              || (fix->indx[1]==i && fix->indx[0]==j) ) {
            /* special energy term for the pair */
              pp=fix->onefour[onefour];
              key+='f'-' '; }

        } while (key=='f');
      }
  } /* onefour */
}

static int nfdih,nfimp,nfaro;  /* # of not defined or zero dihedrals,
                                  impropers, aromatic dihedrals */
static int nexc[EXCDIM];       /* exceptions statistics */

/* ndihedrals... defined elsewhere which is mess... fix later! */
static int naromatics;

#define T(I) site[I].type

int prtatomoffset; /* global, offset of a cluster (if option -n-1) */

void prtatom(int i) /********************************************* prtatom */
{
  prt_("%3i %-4s ",i-prtatomoffset,atom[T(i)].name);
}

static void prt4atoms(int i,int j,int k,int l) /**************** prt4atoms */
{
  prtatom(i); prtatom(j); prtatom(k); prtatom(l);
}

static int findbond(int i,int j) /******************************* findbond */
/***
  finds parameters for bond i-j in the bond table and adds to the bond list
  i,j are indices of atoms in molecule
***/
{
  bond_t *bond,*b;

  looplist (bond,bond0)
    if ( (bond->indx[0]==T(i) && bond->indx[1]==T(j))
      || (bond->indx[1]==T(i) && bond->indx[0]==T(j)) ) {
      /***
        adds the bond to the linked list that will be used for calculating
        energy and forces.  Note that now indx[] contains atom indices
        in the molecule, not atom types
      ***/
      ralloc(b,sizeof(bond_t));
      *b=*bond;
      b->indx[0]=i; b->indx[1]=j;
      b->next=b0; b0=b;
      return 1; }

  /* not found */
  prtatom(i); prtatom(j);

  if (option('f')&1) {
    prt("! WARNING missing bond parameter, using K=50 kcal/mol r0=1.5 A");
    ralloc(b,sizeof(bond_t));
    b->indx[0]=i; b->indx[1]=j;
    b->parm.length=1.5;
    b->parm.K=50; b->parm.Ki2=-2*b->parm.K;
    b->next=b0; b0=b;
    return 1; }
  else
    if (Xopt.W&1) WARNING(("missing bond parameter"))

  return 0;
}

static vector dummyforce; /* static to be initialized to (0,0,0) */
static double Hlimit;
static int allH;
static int lightangleomit;

static int lightangle(angle_t *a) /****************************** lightangle */
/*
  returns 1 if site[i] or site[k] is light so that angle i--j--k
  can be constrained
  In the situation like C-NH-C, if option -h is positive, only one
  angle from both is constrained because otherwise the mechanical system
  is ill-conditioned
  if j is massless, the angle will by removed
  NEW: if water is detected, Hlimit is 1.5 (solved elsewhere)
*/
{
  int i=a->indx[0],j=a->indx[1],k=a->indx[2],l,*kk;
  double mi=atom[site[i].type].mass;
  double mk=atom[site[k].type].mass;
  double ml;
  int iH=mi<Hlimit;
  int kH=mk<Hlimit;

  if (atom[site[j].type].mass==0) return 1;

  if (!iH && !kH) return 0;

  /* at least one of i,k is below mass limit (hydrogen) */

  if (allH || site[j].nnbr!=3)
    return 1; /* ALL angles with H */
  else {
    /* 3 neigbors special : remove planar singularity */
    l=-1;
    loopnbr (kk,j) if (*kk!=i && *kk!=k) l=*kk;
    if (l<0) ERROR(("internal"))
    /* ... now we have i k l around j and either i or k is hydrogen */
    if (iH && kH) {
      ml=atom[site[l].type].mass;
      if (ml<Hlimit) {
        if (mi*mk*ml==0) return 1;
        else WARNING(("\
Constraints for 3 light atoms around one heavy generally unresolvable.\n\
- For a simple pyramid as NH3, use -h-1\n\
NB: negative -h means \"do not remove overdetermined constraints\"")) }
      /* both H: in cases like -NH2, H-N-H angle is free */
      return 0; }
    lightangleomit++;

    if (iH) return k<l; else return i<l; }
}

static int planarangleomit;
static int planarangleunsolved;
static int planarangles;

static int planarangle(angle_t *a) /**************************** planarangle */
/*
  returns 1 if the real angle i--j--k is close to 180 deg and the equilibrium
  angle parameter is less than 91 deg
  bug fixed 6/99: angle par.<0 means 180 deg!
  NOTE: uses coordinates, may give incorrect results if coordinates are wrong
  or are missing.
  To be used for removing extra angles in the planar cases like:

     B
     |
  A--X--C
     |
     D

  where A-X-C and B-X-D should not be included if the equilibrium
  angle is 90 deg
*/
{
  int i=a->indx[0],j=a->indx[1],k=a->indx[2];
  site_t *si=site+i,*sj=site+j,*sk=site+k;

  if (a->parm.angle>91*PI/180 || a->parm.angle<0) return 0;

  planarangles++;

  if (si->keep==WANTED || sj->keep==WANTED || sk->keep==WANTED) {
    planarangleunsolved++;
    return 0; }

  anglepot(si->r,sj->r,sk->r,
           dummyforce,dummyforce,dummyforce,&a->parm);

  if (phi>120*PI/180) {
    planarangleomit++;
    if (option('v')&4) {
      prts_("! "); prtatom(i); prtatom(j); prtatom(k);
      prts(" planar angle omitted"); }
    return 1; }

  return 0;
}


static angle_t *lastbonded,*toremove;
static int validangle;

static int a_bonded(int i,int j) /********************************* a_bonded */
/***
  returns 1 if sites i,j are already bonded in the bond list
  i,j are indices of atoms in molecule
  prepares some angle stuff
***/
{
  bond_t *b;
  angle_t *a;

  looplist (b,b0)
    if ( (b->indx[0]==i && b->indx[1]==j)
      || (b->indx[1]==i && b->indx[0]==j) ) return 1;

  lastbonded=NULL;

  looplist (a,a0) {
    if ( (a->indx[0]==i && a->indx[2]==j)
      || (a->indx[2]==i && a->indx[0]==j) ) if (lightangle(a)) {
        validangle++;
        toremove=a;
        return 1; }
    lastbonded=a; }

  return 0;
}

static int findangle(int i,int j,int k) /*********************** findangle */
/***
  finds parameters for bond angle i-j-k in the angles table and adds to the
  angles list
  i,j,k are indices of atoms in molecule
***/
{
  angle_t A,*angle,*a;

  looplist (angle,angle0)
    if (angle->indx[1]==T(j))
      if ( (angle->indx[0]==T(i) && angle->indx[2]==T(k))
        || (angle->indx[2]==T(i) && angle->indx[0]==T(k)) ) {
      A=*angle;
      A.indx[0]=i; A.indx[1]=j; A.indx[2]=k;
      if (A.parm.angle==180) A.parm.angle=-1; /* PI treated separately */
      A.parm.angle*=PI/180;
      A.parm.Kubi2=-2*A.parm.Kub;
      /* patch: omit planar angles */
      if (planarangle(&A)) return 0;
      ralloc(a,sizeof(angle_t));
      *a=A;
      a->next=a0; a0=a;
      return 1; }

  prtatom(i); prtatom(j); prtatom(k);

  if (option('f')&2) {
    prt("! WARNING missing angle parameter, using 120deg 20kcal/mol");
    A.indx[0]=i; A.indx[1]=j; A.indx[2]=k;
    A.parm.angle=PI*2/3;
    A.parm.K=20; A.parm.K2=A.parm.K*2;
#ifdef UREY_BRADLEY
    A.parm.Kub=0;
    A.parm.Kubi2=0;
    A.parm.length=0;
#endif
    ralloc(a,sizeof(angle_t));
    *a=A;
    a->next=a0; a0=a;
    return 1; }
  else
    if (Xopt.W&2) WARNING(("missing angle parameter"))

  return 0;
}

static int match; /* count exact matches */

static int Match(int tabtype,int type) /**************************** Match */
/***
  One atom matching. tabtype is a table entry.
  Any type but Hydrogen matches tabtype==0 (i.e. wildcard atom "X").
***/
{
  if (option('v')&8) prt("matching: %s %s",atom[type].name,atom[tabtype].name);
  if (tabtype==type) { match++; return 1; } /* exact match */
  if (tabtype==0) return atom[type].X; /* match if in the table */

  return 0; /* no match */
}

#define MATCH(I,i) Match(t->indx[I],T(i))


int finddihedral(int i,int j,int k,int l,torsion_t *torsion0)
/***                                                     **** finddihedral ***
  Finds parameters for dihedral

  i      l
   \    /
    j--k

  using torsion_t *torsion0 as the head of the look-up-table.

  torsion0==dihedral0 for normal dihedrals, torsion0==improper0 for forces
  U=K*phi^2 (phi=dihedral angle) that are used in aromatic rings to keep the
  ring planar.  The corresponding dihedral potential parameters have
  periodicity parameter n=0 and appear, for unknown reasons, in the tables
  of impropers and are also sometimes called `improper torsion' by CHARMM.

  NEW: also torsion0=cisdihedral0 for special term for cis dihedrals

  Wildcard atom X (type=0) in the table stands for any atom marked by 1 in
  column X in the table of atoms.

  If there are more matching entries for the same value of periodicity
  number n, that one with less wildcards is considered; if this selection is
  not unique, this is error.

  If there are more matching entries for different values of n, they are
  summed; this is marked by negative n (|n| is the maximum periodicity), and
  the length of array parm.K[] is n+1 ireals.

  It is not possible to combine n==0 ( U=K*(phi-angle)^2 ) and n!=0.

  Note: n=5 is special and marks `cisdihedrals'

  Returns 1 on success, 0 if the required dihedral is not found.
***/
{
  torsion_t *t,*torsions[MAXPER+1];
  int
    maxmatch[MAXPER+1], /* max # of exact atom matches (i.e. not matches to
                           wildcards) for given n */
    nmaxmatch[MAXPER+1];/* # of cases for which the max match was reached */
  int nn=0,m,M,size;
  int incltoaro=torsion0==improper0;
  char *torsionname;

  /* NEW: for ar_dih_limit<0, we search aromatics always in the table of
     dihedrals 
     -- WHY?  Removed in V3.4d (peptides of gromos96 not affected)
  if (ar_dih_limit<0 && incltoaro) torsion0=dihedral0;
  */
  
  torsionname=
    torsion0==improper0?"aromatic"
    :torsion0==cisdihedral0?"cisdihedral":"dihedral";

  /* 3-cycles skipped */
  if (i==l) return 0;
  loopto (m,0,MAXPER) maxmatch[m]=-1,nmaxmatch[m]=0;

  /* looks for the maximum match number (less X atoms) */
  looplist (t,torsion0) {
    match=0;
    if (MATCH(0,i) && MATCH(1,j) && MATCH(2,k) && MATCH(3,l))
      Max(maxmatch[t->parm.n],match)
    match=0;
    if (MATCH(3,i) && MATCH(2,j) && MATCH(1,k) && MATCH(0,l))
      Max(maxmatch[t->parm.n],match) }

  /* scans again to check whether the max match is unique */
  looplist (t,torsion0) {
    match=0;
    if (MATCH(0,i) && MATCH(1,j) && MATCH(2,k) && MATCH(3,l)) {
      if (maxmatch[t->parm.n]==match) {
        nmaxmatch[t->parm.n]++;
        torsions[t->parm.n]=t; } }
    else { /* to prevent double-match of symmetric cases */
      match=0;
      if (MATCH(3,i) && MATCH(2,j) && MATCH(1,k) && MATCH(0,l))
       if (maxmatch[t->parm.n]==match) {
         nmaxmatch[t->parm.n]++;
         torsions[t->parm.n]=t; } } }

  /***
  nn:=how many dihedral terms with possibly different periodicities n
       have been found
  M:=the largest periodicity found
  ***/
  loopto (m,0,MAXPER) if (nmaxmatch[m]) {
    nn++;
    M=m;
    if (nmaxmatch[m]>1) {
      prt4atoms(i,j,k,l);
      WARNING(("%s selection is not unique (%d times)",torsionname,nn)) } }

  /* some tests about that stupid aromatic case */
  if (torsion0==improper0) if (nn>1 || (nn==1 && M!=0)) {
    prt4atoms(i,j,k,l);
    ERROR(("bad aromatic-ring entry/ies in the table of impropers")) }

  /* no dihedral found */
  if (!nn) {
    if (option('v')&4) if (torsion0!=cisdihedral0) {
      prts_("! ");
      prt4atoms(i,j,k,l);
      prt("%s not found",torsionname); }
    return 0; }

  /***
  calculating size of torsion_t because the length of array parm.K[] is
  variable
  ***/
  size=sizeof(torsion_t);

  /* one term: if M=parm.n=0 then parm.K[1]=angle */
  if (nn==1) if (M==0) size+=sizeof(ireal);

  /* more terms: must allocate space for parm.K[0..M] */
  if (nn>1) {
    size+=M*sizeof(ireal);
    if (nmaxmatch[0]) {
      prt4atoms(i,j,k,l);
      ERROR(("%s: cannot combine n=0 and n!=0",torsionname)) } }

  ralloc(t,size);

  if (nn==1) {
    /* one term: if M=parm.n=0 then parm.K[1]=angle */
    copy(t,torsions[M],size);
    t->X=4-maxmatch[M];
    if (M==0) t->parm.K[1]*=PI/180; }

  else {
    /* more terms: must convert SUM{i} cos(i phi) --> SUM{i} cos^i phi */

    static int C[7][7]={ /* cos(n*x) = SUM{i=0..n} C[n][i] * cos^i(x) */
      { 1},
      { 0, 1},
      {-1, 0, 2},
      { 0,-3, 0, 4},
      { 1, 0,-8, 0,  8},
      { 0, 5, 0,-20, 0,16},
      {-1, 0,18, 0,-48, 0,32} };
    int a;

    t->X=0;

    loopto (a,0,M) t->parm.K[a]=0;
    loopto (m,1,M) if (nmaxmatch[m]) {
      t->parm.K[0]+=fabs(torsions[m]->parm.K[0]);
      t->X+=4-maxmatch[m]; }
    loopto (a,0,M)
      loopto (m,1,M) if (nmaxmatch[m])
        t->parm.K[a]+=torsions[m]->parm.K[0]*C[m][a];

    if (all_dihedrals&2) {
      /* re-set K[0] to reach Emin:=0 */
      ireal *K=t->parm.K;
      double x,E;
      K[0]=3e33;
      for (x=-1; x<=1; x+=1./1024) {
        E=0;
        for (a=M; a; a--) E=(E+K[a])*x;
        Min(K[0],E) }
      K[0]=-K[0]; }

    t->parm.n=-M;  }

  t->indx[0]=i; t->indx[1]=j; t->indx[2]=k; t->indx[3]=l;

  /* adding to the list */
  if (torsion0==cisdihedral0) {
    t->next=cis0; cis0=t;
    if (option('v')&4) {
      prt_("! %s dihedral: ",t->parm.K[0]>0?"cis":"trans");
      prt4atoms(i,j,k,l); _n } }
  else
    /* something problematic here... */
#if 0
    if (torsion0==dihedral0) { t->next=d0; d0=t; }
    else { t->next=ar0; ar0=t; }
#else
    if (incltoaro) { t->next=ar0; ar0=t; }
    else { t->next=d0; d0=t; }
#endif

  return 1;
}

void findanydihedral(int i,int j,int k,int l) /************* findanydihedral */
/* tries dihedral0 and cisdihedral0 */
{
  if (finddihedral(i,j,k,l,dihedral0)) ndihedrals++;
  else nfdih++;
  if (finddihedral(i,j,k,l,cisdihedral0)) ncisdihedrals++;
}

static int chirsg;

int findimproper(int i,int j,int k,int l) /******************** findimproper */
/***
  Finds parameters for the improper torsion

    j
     \
      i->l
     /
    k

  Wildcard atom X (type=0) in the table stands for any atom marked by 1 in
  column X in the table of atoms.

  If there are more matching entries, the one with less wildcards is considered.
  Returns 1 on success, 0 if the required improper is not found.
  The dihedral-like `impropers' in aromatic rings are looked for by procedure
  finddihedral.

  To have correct sign for chiral atoms, the order of atoms j-k-l must be a
  cyclic permutation of the order used to define chirality!

  NEW (V1.7b): all_impropers&2 enables a possibility to use the improper
      j                 j
       \                 \
  for   l->i instead of   i->l
       /                 /
      k                 k
  V1.8h: bug fixed -- active for impropers with zero angle only
***/

{
  torsion_t *torsion,*t;
  int m,maxmatch=-1,nmaxmatch=0;

  looplist (t,improper0) {
    match=0;
    if (MATCH(0,i) && MATCH(3,l)) {
      m=match;
      if (MATCH(1,j) && MATCH(2,k)) Max(maxmatch,match)
      else {
        match=m;
        if (MATCH(1,k) && MATCH(2,j)) Max(maxmatch,match) } }

    /* trying match with swapped i:l */
    /* only terms with angle = t->parm.K[1]==0 */
    if (all_impropers&2) if (t->indx[0]!=t->indx[3]) if (t->parm.K[1]==0) {
      match=0;
      if (MATCH(0,l) && MATCH(3,i)) {
        m=match;
        if (MATCH(1,j) && MATCH(2,k)) Max(maxmatch,match)
        else {
          match=m;
          if (MATCH(1,k) && MATCH(2,j)) Max(maxmatch,match) } } }
    } /* t */

  /* scans again to check whether the max match is unique */
  looplist (t,improper0) {
    match=0;
    if (MATCH(0,i) && MATCH(3,l)) {
      m=match;
      if (MATCH(1,j) && MATCH(2,k)) {
        if (maxmatch==match) { nmaxmatch++; torsion=t; } }
      else {
        match=m;
        if (MATCH(1,k) && MATCH(2,j))
          if (maxmatch==match) { nmaxmatch++; torsion=t; } } }

    /* trying match with swapped i:l */
    if (all_impropers&2) if (t->indx[0]!=t->indx[3]) {
      match=0;
      if (MATCH(0,l) && MATCH(3,i)) {
        m=match;
        if (MATCH(1,j) && MATCH(2,k)) {
          if (maxmatch==match) { nmaxmatch++; torsion=t; } }
        else {
          match=m;
          if (MATCH(1,k) && MATCH(2,j))
            if (maxmatch==match) { nmaxmatch++; torsion=t; } }
        } }
    } /* t */

  if (nmaxmatch>1) {
    prt4atoms(i,j,k,l);
    WARNING(("improper torsion selection is not unique (%d times)",nmaxmatch)) }

  if (!nmaxmatch) {
    /* no match found */
    if (option('v')&4) {
      prts_("! ");
      prt4atoms(i,j,k,l);
      if (!nmaxmatch) prt("improper torsion not found"); }
    return 0; }

  else {

    ralloc(t,sizeof(torsion_t)+sizeof(ireal));
    if (torsion->parm.n) ERROR(("n = %d != 0 for improper",torsion->parm.n))
    copy(t,torsion,sizeof(torsion_t)+sizeof(ireal));

    t->indx[0]=i; t->indx[1]=j; t->indx[2]=k; t->indx[3]=l;
    t->X=4-maxmatch;

    if (t->parm.K[1]) {
      if (!site[i].chir) {
        prtatom(i);
        WARNING(("%s: undefined chirality or improper sign\n\
  edit mol-file or try option -c\n\
  if -c then check input config (options -e -r)", specfn))
        site[i].chir=1; } }
    else
      if (site[i].chir) {
        site[i].chir=0;
        if (option('v')&4) {
          prtc('!'); prtatom(i); prt(" chir=0"); } }

    t->parm.K[1]*=PI/180*site[i].chir*chirsg;

    t->next=i0; i0=t;
    return 1; }
}

static int counttorsion(torsion_t *t0) /********************** counttorsion */
{
  torsion_t *t;
  int n=0;

  looplist (t,t0) n += site[t->indx[0]].clust==clust;

  return n;
}

static int nextcluster(int i) /******************************** nextcluster */
/* returns -1 if there is no following cluster, otherwise site # */
{
  int clust=site[i].clust;

  do
    if (++i>=ns) return -1;
   while (site[i].clust<0 || site[i].clust==clust);

  return i;
}

static void prttorsion(char *nm,torsionpot_t *pot,torsion_t *t0)
                                                          /****** prttorsion */
/***
  to print tables of dihedrals and impropers
***/
{
  int i,j,k,l,m;
  torsion_t *t;
  double U;

  _n
  prts(nm);
  prts("! i atom   i atom   i atom   i atom nX  n    K|K[0] calc.angle Upot  angle|K[1]");

  looplist (t,t0) if (site[t->indx[0]].clust==clust) {
    i=t->indx[0]; j=t->indx[1]; k=t->indx[2]; l=t->indx[3];

    U=pot(site[i].r,site[j].r,site[k].r,site[l].r,
          dummyforce,dummyforce,dummyforce,dummyforce,&t->parm);

    if (nm[0]=='d') {
      double phi0=phi;

      measure=3;
      (void)improperpot(site[i].r,site[j].r,site[k].r,site[l].r,
                        dummyforce,dummyforce,dummyforce,dummyforce,&t->parm);
      measure=2;
    if (fabs(fabs(phi)-phi0)>1e-6) ERROR(("internal")) }

    prt4atoms(i,j,k,l);

    prt_(nm[0]=='d'?dihedral_f:torsion_f,
         t->X,t->parm.n, t->parm.K[0], phi*(180/PI),U);

    if (t->parm.n==0)
     prt_(z_f,t->parm.K[1]*(180/PI)); /* angle */
    else if (t->parm.n<0)
      loopto (m,1,-t->parm.n) prt_(z_f,t->parm.K[m]);
    _n }
}


static double dist(int i,int j) /*********************************** dist */
/***
  bond length between sites i and j
  uses pre-compiled bond table b0->
***/
{
  bond_t *b;

  looplist (b,b00)
    if ( (b->indx[0]==i && b->indx[1]==j)
      || (b->indx[0]==j && b->indx[1]==i) ) return b->parm.length;
  WARNING(("unknown %d-%d bond (dist:=1)",i,j))

  return 1;
}

static int bonded(int i,int j) /************************************ bonded */
/* returns 1 if there is a bond connecting atoms i--j */
{
  int *k;

  loopnbr (k,i) if (*k==j) return 1;

  return 0;
}

static int aromatic4(int i,int j,int k,int l) /******************* aromatic4 */
/***
  Returns 1 if the chain i--j--k--l is a member of an aromatic ring.
  A 5-ring is called aromatic if it contains at least 3 aromatic atoms
  A 6-ring is called aromatic if it contains at least 4 aromatic atoms

  NOTE: in new versions of CHARMM, all atoms in aromatic rings are of
  aromatic type so that simpler algorithm could be used.  The presented
  algorithm is more general and allows some non-aromatic (general) atom
  types in aromatic rings if there are not enough true aromatic atom types.
***/
{
  int aroma = atom[site[i].type].A + atom[site[j].type].A
            + atom[site[k].type].A + atom[site[l].type].A;

  int *inbr,*lnbr;

  if (i==l) return 0; /* 3-ring */

  loopnbr (inbr,i) if (*inbr!=j) if (*inbr!=k) if (*inbr!=l)
    loopnbr (lnbr,l) if (*lnbr!=k) if (*lnbr!=j) if (*lnbr!=i) {
      /* chain *inbr--i--j--k--l--*lnbr */

      if (*inbr==*lnbr)
        /* it is a 5-ring */
        return aroma + atom[site[*inbr].type].A >= 3;
      else if (bonded(*inbr,*lnbr))
        return
          aroma + atom[site[*inbr].type].A + atom[site[*lnbr].type].A >= 4; }

  return 0;
}

static int aromatic3(int i,int j,int k) /************************ aromatic3 */
/***
  Returns 1 if the chain i--j--k is a member of an aromatic ring.
  See aromatic4 for a definition of an aromatic ring.
***/
{
  int aroma = atom[site[i].type].A + atom[site[j].type].A
            + atom[site[k].type].A;

  int *l,*inbr,*lnbr;

  loopnbr (l,k) if (*l!=i) if (*l!=j)
    loopnbr (inbr,i) if (*inbr!=j) if (*inbr!=k) if (*inbr!=*l)
      loopnbr (lnbr,*l) if (*lnbr!=k) if (*lnbr!=j) if (*lnbr!=i) {
        /* chain *inbr--i--j--k--*l--*lnbr */

        if (*inbr==*lnbr)
          /* it is a 5-ring */
          return aroma + atom[site[*l].type].A + atom[site[*inbr].type].A >= 3;
        else if (bonded(*inbr,*lnbr))
          return
            aroma
              + atom[site[*l].type].A
              + atom[site[*inbr].type].A + atom[site[*lnbr].type].A >= 4; }

  return 0;
}

static int aromatic2(int j,int k) /***************************** aromatic2 */
/***
  Returns 1 if the bond j--k is a member of an aromatic ring.
  See aromatic4 for a definition of an aromatic ring.
***/
{
  int aroma = atom[site[j].type].A + atom[site[k].type].A;

  int *i,*l,*inbr,*lnbr;

  loopnbr (i,j) if (*i!=k)
    loopnbr (l,k) if (*l!=*i) if (*l!=j)
      loopnbr (inbr,*i) if (*inbr!=j) if (*inbr!=k) if (*inbr!=*l)
        loopnbr (lnbr,*l) if (*lnbr!=k) if (*lnbr!=j) if (*lnbr!=*i) {
          /* chain *inbr--*i--j--k--*l--*lnbr */

          if (*inbr==*lnbr)
            /* it is a 5-ring */
            return
              aroma + atom[site[*i].type].A + atom[site[*l].type].A
                + atom[site[*inbr].type].A >= 3;
          else if (bonded(*inbr,*lnbr))
            return
              aroma + atom[site[*i].type].A + atom[site[*l].type].A
                + atom[site[*inbr].type].A + atom[site[*lnbr].type].A >= 4; }

  return 0;
}

static int includedih(int i,int j,int k,int l) /***************** includedih */
/* returns 1 if dihedral i-j-k-l is to be included */
{
  switch (ar_dih_limit) {
    case 0:
    case 1: if (aromatic2(j,k)) return 0;
    case 2: if (aromatic3(i,j,k) | aromatic3(l,k,j)) return 0;
    case 3: if (aromatic4(i,j,k,l)) return 0;
    case -1:
    case 4: return 1; /* all included */
    default: ERROR(("ar_dih_limit=%d is invalid",ar_dih_limit)) }

  return 0;
}

static void inclexc(int i,int nbr,enum exc_e type) /**************** inclexc */
/***
  Include nbr to the list of exceptions for site i if nbr<1
  If already in the list, chooses the lowest type (that is, 1-2 exception
  overrides 1-4)
  i should not change until all exceptions for i are found
  The list is ordered by index.
***/
{
  exception_t *exc,**excp;

  if (type==ONEFOUR) if (atom[site[i].type].no14 & atom[site[nbr].type].no14)
    /* no special LJ terms for 1-4 interaction */
    return;

  if (nbr>=i) return;
  looplist (exc,site[i].exc) if (exc->indx==nbr) {
    Min(exc->type,type) return; }

  ralloc(exc,sizeof(exception_t));

  for (excp=&site[i].exc; *excp && (*excp)->indx<nbr; excp=&(*excp)->next);
  exc->next=*excp;
  *excp=exc;
  exc->indx=nbr;
  exc->type=type;
}


static int comp(double x,double y) /********************************* comp */
{
  if (x==y) return 0;
  else if (x<y) return -1;
  else return 1;
}

static void prtchir(site_t *s,int which) /************************ prtchir */
/*
  prints one of the three characters: -?+
  which=0: just chir
  which=1: ordered by atom names (alphabetically)
  which=2: ordered by atomic masses
  for 1 and 2, only the three neighbors are considered; if the chirality
  cannot be determined, ? is printed
*/
{
  if (s->nnbr!=3) {
    prtc('?'); return; }
  else {
    int chir=s->chir;
    int i=site[s->nbr[0]].type,
        j=site[s->nbr[1]].type,
        k=site[s->nbr[2]].type;

    switch (which) {
      case 2:
        chir*=comp(atom[i].mass,atom[j].mass)
             *comp(atom[j].mass,atom[k].mass)
             *comp(atom[k].mass,atom[i].mass);
        goto doit;
      case 1:
        chir*=strcmp(atom[i].name,atom[j].name)
             *strcmp(atom[j].name,atom[k].name)
             *strcmp(atom[k].name,atom[i].name);
      doit:
        if (chir) chir = chir>0 ? 1 : -1;
      case 0:
        prtc("-?+"[chir+1]); } }
}

static struct {
  double Kfactor;
  double sidec;
} a2b;

void angle2bond(double a,double b,double gamma) /**************** angle2bond */
/*
  side c=AB of triangle ABC is returned in a2b.sidec
  multiplicative factor for recalculating K is returned in a2b.Kfactor
  Angle gamma=ACB is in radians
  Negative gamma means gamma=180 deg, Kfactor is replaced by 100 (would be infty)
*/
{
  double c2;
  int singular=gamma<0;

  if (singular) {
    /* FOOL PROOF TEST of angle=-1 meaning 180... */
    if (fabs(gamma+0.017453293)>1e-7) ERROR(("internal"))
    gamma=PI; }

  c2=a*a+b*b-2*a*b*cos(gamma);
  a2b.sidec=sqrt(c2);
  if (singular) a2b.Kfactor=100;
  else a2b.Kfactor=c2/sqr(a*b*sin(gamma));
}

/* WARNING: this function repeated in blendedt.c */
static int clustersize(int i) /********************************* clustersize */
/*
  returns # of sites in cluster started by site[i]
*/
{
  int cl=site[i].clust,n=0;

  while (i<ns && site[i].clust==cl) i++,n++;

  return n;
}

static char *watername="";

static void detectwater(int i) /******************************** detectwater */
/*
  detect water model written in the .par file as `waters'
  - order of sites must match
  - i is the 1st site of the cluster which is tested
  - the cluster must be contiguous (here)
  - watername is returned
*/
{
  int cl=site[i].clust,n,i0=i,isw;
  struct water_s *w;

  watername="";

  looplist (w,water0) {
    i=i0; n=0; isw=1;
    while (i<ns && site[i].clust==cl) {
      if (n>=w->ns
	  || w->type[n]!=site[i].type
	  || fabs(w->charge[n]-site[i].charge)>1e-6) isw=0;
      i++,n++; }
    if (isw && n==w->ns) {
      watername=w->name;
      prt("! %s water detected, rigid water model passed to cook",watername); }
  }
}

static int equalclusters(int i,int j) /*********************** equalclusters */
/*
  returns 1 if clusters are the same and 0 otherwise
  site[i] and site[j] are the two sites starting the clusters
WARNING:
  clusters having the same numbers and types of atoms and bonds on
  each atom are considered identical
  thus, clusters written with different order of atoms
  are not recognized as identical
  on the other hand, the following two clusters will be incorrectly
  considered as identical:
    A--B       A--B
    |  |  and   \/
    |  |        /\
    D--C       D--C
  (cluster must be contiguous here!)
*/
{
  int n=clustersize(i);

  if (n!=clustersize(j)) return 0;

  for ( ; n; n--,i++,j++) {
    if (site[i].type!=site[j].type
     || site[i].nnbr!=site[j].nnbr) return 0; }

  return 1;
}

void buildforcefield(species_t *spec) /********************* buildforcefield */
/***
assigns the force field to molecule of species spec
***/
{
  int nc; /* # of constraints (bonds) */
  int i,ii,tryagain,wasaromatic;
  int *j,*inbr,*jnbr;
  exception_t *exc;

  memset(nexc,0,sizeof(nexc));
  nfdih=nfimp=nfaro=0;
  _n _n
  *(spec->ext)=0;
  ii=strlen(spec->fn)+8;

  /* 11/2004 (detect water even if first), 4/2012 (one instance here) */
  detectwater(0);

  loop (i,0,ii) prtc('!');
  _n
  if (*watername)
    prt("species %s %s.0",watername,spec->fn);
  else
    prt("species %s",spec->fn);
  loop (i,0,ii) prtc('!');
  _n
  prt("! %s",spec->info);

  site=spec->site;

  imprison(spec);

  prt("! mass=%.4f g/mol",masscenter(spec,1));

  nc=nangles=ndihedrals=nimpropers=naromatics=ncisdihedrals=ndependants=0;

  /*********************************************************/
  /*  bonds, angles, dihedrals, and non-bonded exceptions  */
  /*********************************************************/

  planarangleomit=planarangleunsolved=planarangles=0;

  loop (i,0,ns) {

    /*
    initialize list of exceptions for site [i]
    for atom i, a list of exceptions contains those atoms with lower
    indices that either do not contribute to bonded forces at all
    (1-2 and 1-3 terms), or have different values of parameters (1-4)
    */
    site[i].exc=NULL;

    loopnbr (j,i) {

      measure=2; /* needed for planarangle */
      /*** angles ***/
      loopnbr (inbr,i) if (*j>*inbr)
        nangles += findangle(*inbr,i,*j);
      measure=1;

      /*** exceptions ***/

      /* 1-2 (bond) exceptions */
      inclexc(i,*j,EXCL);

      switch (distance14) {

        case 4: /* NORMAL case 1-4 exceptional, 1-2 and 1-3 removed */
          loopnbr (jnbr,*j) if (*jnbr!=i) {
            /* 1-3 exceptions */
            inclexc(i,*jnbr,EXCL);
            loopnbr (inbr,*jnbr) if (*inbr!=*j) {
              /* 1-4 exceptions */
              enum exc_e exc=ONEFOUR;

              switch (ar_14_limit) {
                case 4: if (aromatic4(i,*j,*jnbr,*inbr)) exc=EXCL;
                        break;
                case 3: if ( aromatic3(i,*j,*jnbr) || aromatic3(*j,*jnbr,*inbr) )
                        exc=EXCL;
		        break;
                case 2: if ( aromatic2(*j,*jnbr) ) exc=EXCL; }
              inclexc(i,*inbr,exc); } }
          break;

        case 3: /* STRANGE case 1-3 exceptional, 1-2 removed */
          loopnbr (jnbr,*j) if (*jnbr!=i) {
            /* 1-3 exceptions */
            inclexc(i,*jnbr,ONEFOUR); }
          break;

        case 5: /* OPLS case 1-5 exceptional, 1-2, 1-3, and 1-4 removed */
          loopnbr (jnbr,*j) if (*jnbr!=i) {
            /* 1-3 exceptions */
            inclexc(i,*jnbr,EXCL);
            loopnbr (inbr,*jnbr) if (*inbr!=*j) {
              int *fifth;
              /* 1-4 exceptions */
              inclexc(i,*inbr,EXCL);
              loopnbr (fifth,*inbr) if (*fifth!=*jnbr)
                /* 1-5 exception */
                /* i--*j--*jnbr--*inbr-*fifth */
                inclexc(i,*fifth,ONEFOUR); } }
          break;

        case 6: /* SPECIAL CASE 1-6 exceptional, 1-2, 1-3, 1-4, 1-5 removed */
          loopnbr (jnbr,*j) if (*jnbr!=i) {
            /* 1-3 exceptions */
            inclexc(i,*jnbr,EXCL);
            loopnbr (inbr,*jnbr) if (*inbr!=*j) {
              int *fifth;
              /* 1-4 exceptions */
              inclexc(i,*inbr,EXCL);
              loopnbr (fifth,*inbr) if (*fifth!=*jnbr) {
                int *sixth;
                /* 1-5 exception */
                inclexc(i,*fifth,EXCL);
                loopnbr (sixth,*fifth) if (*fifth!=*inbr)
                  /* 1--6 exception */
                  /* i--*j--*jnbr--*inbr-*fifth--*sixth */
                  inclexc(i,*sixth,ONEFOUR); } } }
          break;

        default:
          ERROR(("internal")) }

      /* dihedrals, incl. aromatic ones */

      if (*j<i) {
        /*** now (i-*j) runs over all bonds ***/

        /*  Note that bonds must be calculated before impropers
            because the table of bonds is used by function dist */
        nc+=findbond(i,*j);

        if (site[i].nnbr>1 && site[*j].nnbr>1 && all_dihedrals>=0) {

          if (all_dihedrals&1) {
            /* all possible dihedrals around the bond */

            loopnbr (inbr,i) if (*inbr!=*j)
              loopnbr (jnbr,*j) if (*jnbr!=i) {
                /* now *inbr--i--*j--*jnbr is a candidate for a dihedral */

	        /* to distinguish aromatic rings? */
                if (ar_dih_limit>=0 && aromatic4(*inbr,i,*j,*jnbr)) {
		  /* look for aromatics (n=0 dihedrals) in impropers */
                  if (finddihedral(*inbr,i,*j,*jnbr,improper0)) naromatics++;
		  else nfaro++;

		  /* in addition, look for normal dihedrals (but cis/trans) */
                  if (includedih(*inbr,i,*j,*jnbr)) {
		    if (finddihedral(*inbr,i,*j,*jnbr,dihedral0)) ndihedrals++;
                    else nfdih++; } }

                else
		  /* not aromatic - normal search or ar_dih_limit<0 */
                  findanydihedral(*inbr,i,*j,*jnbr); } }

          else { /* !(all_dihedrals&1) */
            /* one dihedral per a bond:              i--*j           */
            /* prefers that dihedral that           /     \          */
            /* has continued chain :         --*inbr       *jnbr--   */
            /* if not unique, `random' selection is made             */
            /* no cis-dihedrals ! */
            tryagain=wasaromatic=0;
            loop (ii,0,3) { /* ii=# of `free' ends */

              loopnbr (inbr,i) if (*inbr!=*j)
                loopnbr (jnbr,*j) if (*jnbr!=i)
                  if ( (site[*inbr].nnbr==1) + (site[*jnbr].nnbr==1) == ii) {
                    /* now *inbr--i--*j--*jnbr is a candidate for a dihedral
                       with ii free ending atoms */

                    if (tryagain) prt("! .. trying again");

                    if (includedih(*inbr,i,*j,*jnbr)) {
                      if (finddihedral(*inbr,i,*j,*jnbr,dihedral0)) {
                        if (tryagain) prt("! .. found equivalent (%d free ends)",ii);
                        ndihedrals++; goto dihfound; }
                      tryagain=option('v')&4; }
                    else wasaromatic++; }

              } /*ii*/

            if (!wasaromatic) nfdih++;

	  dihfound:

	    /* look for aromatics (n=0 dihedrals) in impropers */

            if (ar_dih_limit>=0)
              loopnbr (inbr,i) if (*inbr!=*j)
                loopnbr (jnbr,*j) if (*jnbr!=i) {
                  /* now *inbr--i--*j--*jnbr is a candidate for a dihedral */

                  if (aromatic4(*inbr,i,*j,*jnbr)) {
                    if (finddihedral(*inbr,i,*j,*jnbr,improper0)) naromatics++;
                    else nfaro++; } }

          } /*!all_dihedrals*/
        }
      } /* bond (i-*j) */
    } /* *j */

    /* count all exceptions */
    looplist (exc,site[i].exc) nexc[exc->type]++;

#ifdef POLAR
    site[i].polar=0;

    if (option('~')) {
      isotropicparm_t *ip=&atom[site[i].type].isotropicparm;
      ireal kappa=ip->kappa;
      struct shellrep_s *sr;

      /* shell-core - the repulsive atom */
      site[i].polar = 0;
      if (ip->rep) site[i].polar = POL_REP;

      if (ip->alpha) {
        site[i].polar |= POL_ISO | POL_SAT*(ip->Esat!=0) | POL_SHL*(kappa!=0);
        site[i].pol=(void*)ip;

        loopnbr (j,i) {
          struct polarbond_s *polarbond;

          looplist (polarbond,polarbond0)
            if (polarbond->indx[0]==T(i) && polarbond->indx[1]==T(*j)) {
              if (site[i].polar) {
                prts_("! "); prtatom(i); prtatom(*j);
                if (site[i].polar&POL_ISO)
                  prt("! isotropic polarizability overridden by bond axial");
                else
                  ERROR(("polarbond: multiple axial polarizibility")) }
              site[i].polar=POL_AXI | POL_SAT*(polarbond->parm.Esat!=0)
                          | POL_REP*(((isotropicparm_t*)site[i].pol)->kappa!=0);
              ralloc(site[i].pol,sizeof(axial_t));
              ((axial_t*)site[i].pol)->kappa=kappa;
              if (((isotropicparm_t*)site[i].pol)->kappa!=kappa)
                ERROR(("internal"))
              ((axial_t*)site[i].pol)->tozz=*j;
              ((axial_t*)site[i].pol)->parm=polarbond->parm; }

          loopnbr (jnbr,*j) {
            struct polarangle_s *polarangle;

            looplist (polarangle,polarangle0)
              if (polarangle->indx[0]==T(i)
                  && polarangle->indx[1]==T(*j)
                  && polarangle->indx[2]==T(*jnbr)) {
                if (site[i].polar) {
                  prts_("! "); prtatom(i); prtatom(*j); prtatom(*jnbr);
                  if (site[i].polar&POL_ISO)
                    prt("! isotropic polarizability overridden by angle axial");
                  else
                    ERROR(("polarangle: multiple axial polarizibility")) }
                site[i].polar=POL_AXI | POL_SAT*(polarangle->parm.Esat!=0)
                           | POL_REP*(((isotropicparm_t*)site[i].pol)->kappa!=0);
                ralloc(site[i].pol,sizeof(axial_t));
                ((axial_t*)site[i].pol)->kappa=kappa;
                if (((isotropicparm_t*)site[i].pol)->kappa!=kappa)
                  ERROR(("internal"))
                ((axial_t*)site[i].pol)->tozz=*j;
                ((axial_t*)site[i].pol)->parm=polarangle->parm; }
          } /* jnbr */
        } /* j */
      } } /* alpha,option '~' */
#endif

  } /* i */

  b00=b0; /* to enable function distance() */

  /* removed planar angles report */
  if (planarangleomit)
    prt("! %d planar angles from %d omitted",planarangleomit,planarangles);
  else if (planarangles)
    prt("! WARNING: %d planar angles and none omitted",planarangles);

  if (planarangleunsolved)
    WARNING(("%d unresolved planar angles (no coord)",planarangleunsolved))

  /***************/
  /*  impropers  */
  /***************/
  /*
  real impropers, not those keeping aromatic rings planar
  bonds must be calculated in advance
  */

  if (all_impropers>=0) loop (i,0,ns) if (site[i].nnbr==3) {
    int sw;
    int j=site[i].nbr[0],k=site[i].nbr[1],l=site[i].nbr[2];
    /* three atoms (j,k,l) around a central atom (i):
    j
     \
      i---l
     /
    k      */

#define SWAP(I,J) { sw=I; I=J; J=sw; chirsg=-chirsg; }

    chirsg=1; /* flipped by SWAP */

    if (all_impropers&1) {
      /*** all three possible improper torsion terms are chosen ***/
      if (findimproper(i,j,k,l)) nimpropers++; else nfimp++;
      if (findimproper(i,k,l,j)) nimpropers++; else nfimp++;
      if (findimproper(i,l,j,k)) nimpropers++; else nfimp++; }

    else {
      /***
      Only one of the three possible improper torsion terms is chosen
      The criteria are based on (1) valence of the three atoms around the
      central atom, (2) lengths of bonds - see the cases below
      ***/

      switch ((site[j].nnbr!=1) + (site[k].nnbr!=1) + (site[l].nnbr!=1)) {

        case 0: case 3:
          /*
          j             --j
           \               \           rearranges (j,k,l) so that
            i---l  or       i---l--    |ij| <= |ik| <= |il|
           /               /           (i--l is the longest bond)
          k             --k            */

          again:
            if (dist(j,i)>dist(k,i)) SWAP(j,k)
            if (dist(k,i)>dist(l,i)) { SWAP(k,l); goto again; }
          break;

        case 1:
          /*
          --j              j             j
             \              \             \        rearranges (j,k,l) so
              i--l  or       i--l  or      i--l--  that l is bonded
             /              /             /        (the last symmetric case)
            k            --k             k         */

          if (site[k].nnbr!=1) SWAP(l,k)
          else if (site[j].nnbr!=1) SWAP(l,j)
          break;

        case 2:
          /*
          --j               j            --j
             \               \              \      rearranges (j,k,l) so
              i--l--  or      i--l--  or     i--l  that l is not bonded
             /               /              /      (the last symmetric case)
            k             --k            --k       */

          if (site[k].nnbr==1) SWAP(l,k)
          else if (site[j].nnbr==1) SWAP(l,j)
          break;

        default: ERROR(("internal")) }

      if (!findimproper(i,j,k,l)) {
        /* OOPS! the improper was not found... try other arrangements */
        if (option('v')&4) prt("! .. trying again");
        if (!findimproper(i,k,l,j)) {
          /* if cases 0,3, points to second shortest |ij| */
          if (option('v')&4) prt("! .. trying again");
          if (!findimproper(i,l,j,k)) {
            /* definitively not found */
            nfimp++; goto impropernotfound; } } }

      nimpropers++;
      impropernotfound:; } }

}


static void printforcefield(species_t *spec) /************* printforcefield */
/*
  print tables to *.ble file, according to option -v
  if N<0, tries to split the molecule into submolecules (clusters)
*/
{
  int i,ii,j,iclust,*k,*l,*ll;
  int nc,nsclust,Nclust,cbonds;
  double U;
  bond_t *b;
  angle_t *a;
  exception_t *exc;
  struct dependant_s *d;

  if (!(option('v')&3)) {
    /***** no force field+sites output *******/
    newnspec=nspec; /* pretend number of species: does NOT calculate the
                       correct number from clusters */
    return; }

  measure=2; /* to force functions in intrapot to return angles */

  /* generating dependants */
  dep0=NULL;
  removeb=removea=removei=ii=0;

  if ( !(ii=readdependants(spec)) )
  loop (i,0,ns) if (atom[site[i].type].mass==0) switch (site[i].nnbr) {
    case 0:
      ERROR(("site %d massless and free",i))
      break;
    case 1:
      j=site[i].nbr[0];
      if (site[j].nnbr<2) {
        ERROR(("%d-%d massless atom in diatomisc",i,j))
        break; }
      goto adddep;
    default:
      j=i;
    adddep:
      ii++;
      alloc(d,sizeof(struct dependant_s));
      d->next=dep0; dep0=d;
      d->ndep=1;
      d->dep[0]=i;
      d->nnbr=site[j].nnbr;
      copy(d->nbr,site[j].nbr,sizeof(d->nbr));
      loop (iclust,0,d->nnbr) if (d->nbr[iclust]==i) d->nbr[iclust]=j;

      anymass=0;
      loopnbr (k,j) {
        removebond(*k,j);
        loopnbr (l,j) {
          removeangle(*k,j,*l);
          loopnbr (ll,j) removetorsion(&i0,j,*k,*l,*ll); } } }

  prt("! %d dependants: %d bonds, %d angles, %d impropers removed",
    ii,removeb,removea,removei);

  /* some angles will be constrained, i.e. printed as bonds */
  ii=0;
  lightangleomit=0;
  looplist (a,a0) ii+=lightangle(a);
  /* note that lightangleomit is # of accepted + # of omitted = 2*# of omitted */
  if (ii) prt("! %d angles constrained, %d singular not",ii,lightangleomit/2);

  prt_(
    "\n! %d dihedrals, %d impropers, and %d aromatics zero or not found\n",
            nfdih,        nfimp,            nfaro);

  if (nfdih+nfimp+nfaro) fprintf(stderr,
    "\n! %d dihedrals, %d impropers, and %d aromatics zero or not found\n",
         nfdih,        nfimp,            nfaro);

  prt("! %d pairs excluded (1-2 etc.)   %d interactions 1-%d",
         nexc[EXCL],                nexc[ONEFOUR],     distance14);

  if (spec->N>=0) {
    /* molecule will not be split: just pretend 1 cluster */
    loop (i,0,ns) site[i].clust=0;
    spec->nclust=1; }

  if (spec->nclust>1) prt("! split into %d submolecules",spec->nclust);

  for (;;) { /* (((((((((((((((((((((((((((((((((((((((((((((((((((((((((( */

    loop (iclust,0,ns) if ( (clust=site[iclust].clust) >= 0 ) goto clusterfound;

    break; /* all clusters done */

    clusterfound:

    /* recalculate numbers of sites for the cluster (and check for water models) */
    nsclust=clustersize(iclust);

    detectwater(iclust);

    if (clust) {
      spec->ext[0]=0;
      prt("\nspecies %s %s.%d",watername,spec->fn,clust); }

    /* try to find all occurences of the same submolecule */

    Nclust=0; /* # of the same clusters (submolecules) */

    for (ii=iclust; ii>=0; ii=nextcluster(ii))
      Nclust += equalclusters(iclust,ii);

    /* recalculate # of bonds nc, # of angles nangles etc. */
    nc=0;
    /* normal chemical bonds */
    looplist (b,b0)
      if (site[b->indx[0]].clust==clust) {
        nc++;
        /* test: are clusters really clusters ? */
        if (site[b->indx[1]].clust!=clust) ERROR(("internal")) }

    if (watername[0]) {
      if (!spec->opt_h) WARNING(("\
This registered water model (%s) is normally rigid, but you\n\
*** did not specify option -h so that it will be exported to cook with\n\
*** a flexible angle (flexible angles). IS THIS WHAT YOU WANT?\n\
*** For a rigid model: use -h in front of water name (-h- to cancel it later).\n\
*** For a flexible model: remove this warning from the ble-file.",watername))
      if (abs(spec->opt_h)>1)
        prt("! WARNING: nonstandard option -h%d used with water model %s.",spec->opt_h,watername); }

    Hlimit=abs(spec->opt_h)+0.5;

    if (option('v')&4) prt("! Hlimit=%g",Hlimit);

    nangles=0;
    /* angles + constraint angles */
    looplist (a,a0)
      if (site[a->indx[0]].clust==clust) {
        if (lightangle(a)) nc++;
        else nangles++; }

    if (spec->opt_h>0) {
      /* remove overdetermined constraints: tetrahedral H only */
      ii=0;
      loop (i,0,ns) if (site[i].nnbr==4) {
        int *k=site[i].nbr;

        validangle=0;
        if (a_bonded(k[0],k[1]) && a_bonded(k[0],k[2]) && a_bonded(k[0],k[3])
            && a_bonded(k[1],k[2]) && a_bonded(k[1],k[3]) && a_bonded(k[2],k[3]) ) {
          if (validangle) {
            if (option('v')&4) {
              prts_("! "); prtatom(i);
              prt("overdetermined angle removed (from %d)",validangle); }
            if (lastbonded) lastbonded->next=toremove->next;
            else a0=toremove->next;
            ii++; }
          else
            ERROR(("internal")) } }
      nc-=ii;
      prt("! %d overdetermined angle constraints removed",ii); }

    prt("\ni=%d  N=%d  config=%d  water=%d",
        newnspec++,spec->N<0?Nclust:spec->N,spec->N<0,watername[0]?1:0);
    prt("ns=%d  nc=%d  nangles=%d  ndihedrals=%d  nimpropers=%d  naromatics=%d",
        nsclust,nc,nangles,counttorsion(d0),counttorsion(i0),counttorsion(ar0));
    if (spec->zero_energy!=0) prt("zero_energy=%.12g",spec->zero_energy);
#ifdef POLAR
    ii=0;
    loop (i,0,ns) if (site[i].clust==clust) if (site[i].polar&POL_AXI) ii++;
    prt_("naxials=%d ",ii);
#endif
    prt("ndependants=%d ;",countdependants());

    if (ncisdihedrals) prt("! WARNING ncisdihedrals=%d not exported",counttorsion(cis0));
    prt("! charge = %.3f",totalcharge(spec,clust));

    /*** print sites in the molecule ***/
    prtatomoffset=iclust;
    prt("\nsites\n\
! i atom   charge %s #  excluded *1-%d  chir:nam",
        spec->N<0?"":"     x        y        z   ",distance14);
    loop (i,0,ns) if (site[i].clust==clust) {
      int nexc[EXCDIM],j;
      site_t *s=&site[i];

      prtatom(i);
      prt_(charge_f,s->charge);
      if (spec->N>=0) prt_(sites_f, s->r[0],s->r[1],s->r[2]);
      memset(nexc,0,sizeof nexc);
      looplist (exc,s->exc) nexc[exc->type]++;
      prt_(" %d  ",nexc[ONEFOUR]+nexc[EXCL]);
      looplist (exc,s->exc) {
        if (exc->type==ONEFOUR) prtc('*');
        prt_("%d ",exc->indx-prtatomoffset); }

      /* print normal chirality/alphabetical chirality/mass chirality */
      if (s->chir) loop (j,0,3) prtchir(s,j);

      /* print atom id - verbose mode only */
      if (option('v')&4) prt_(" =%s",s->id);

      _n }

    if (spec->N<0) {
      /* print config of all clusters */
      prts("\nconfig\n\
!    x        y        z");
      for (ii=iclust; ii>=0; ii=nextcluster(ii))
        if (equalclusters(iclust,ii)) {
          loop (i,ii,ii+nsclust) {
            site_t *s=site+i;
            prt(sites_f, s->r[0],s->r[1],s->r[2]); }
          _n }
      }

    if (option('v')&2) { /************/

      /*** print bonds ***/
      prts("\nbonds\n\
! i atom   i atom  K[kcal/mol/AA^2]  r[AA]   calc.     Upot");

      /* normal chemical bonds */
      looplist (b,b0) if (site[b->indx[0]].clust==clust) {
        int i=b->indx[0], j=b->indx[1];

        U=bondpot(site[i].r,site[j].r, dummyforce,dummyforce, &b->parm);
        prtatom(i); prtatom(j);
        prt(bond_f, b->parm.K,b->parm.length,phi,U); }

      /* bonds corresponding to constrained angles */
      cbonds=0;
      looplist (a,a0) if (site[a->indx[0]].clust==clust) {
        int i=a->indx[0], j=a->indx[1], k=a->indx[2];
        vector dr;

        if (lightangle(a)) {
          if (cbonds++==0) prt("! bonds equivalent to constrained angles:");
          U=anglepot(site[i].r,site[j].r,site[k].r,
                     dummyforce,dummyforce,dummyforce,&a->parm);
          prtatom(i); prtatom(k);
          VVV(dr,=site[i].r,-site[k].r)
          angle2bond(dist(i,j), dist(j,k), a->parm.angle);
          prt(bond_f, a->parm.K*a2b.Kfactor,a2b.sidec, sqrt(SQR(dr)),U); } }

      /*** print angles ***/
      prts("\nangles\n\
! i atom   i atom   i atom  K[kcal/mol] angle[deg]  calc.      Upot  (Urey-Bradley-bond:K   length)");

      looplist (a,a0) if (site[a->indx[0]].clust==clust) {
        int i=a->indx[0], j=a->indx[1], k=a->indx[2];

        if (lightangle(a)) continue; /* have been incl. as bonds */
        U=anglepot(site[i].r,site[j].r,site[k].r,
                   dummyforce,dummyforce,dummyforce,&a->parm);
        prtatom(i); prtatom(j); prtatom(k);
        if (a->parm.Kub)
          prt(angleUB_f,
              a->parm.K,a->parm.angle*(180/PI), phi*(180/PI), U,
              a->parm.Kub,a->parm.length);
        else
          prt(angle_f, a->parm.K,a->parm.angle*(180/PI), phi*(180/PI), U);
        }

      prttorsion("dihedrals",dihedralpot,d0);
      prttorsion("impropers",improperpot,i0);
      prttorsion("aromatics",dihedralpot,ar0);

      prtdependants();

#ifdef POLAR
      /* print axial polarizabilities */
      {
        int i;

        prt("axials\n\
! i atom-->i atom   arep    alpha  alphazz   Esat");

        loop (i,0,ns) if (site[i].clust==clust)
          if (site[i].polar&POL_AXI) {
            axial_t *ax=(axial_t *)site[i].pol;

            prtatom(i); prtatom(ax->tozz);
            prt(axial_f,ax->kappa,ax->parm.alpha,ax->parm.alphazz,ax->parm.Esat); }
        _n }
#endif

    } /*********** option('v')&2 */

    /* mark sites of the cluster as done */
    for (ii=nextcluster(iclust); ii>=0; ) {
      int ii0=ii;
      ii=nextcluster(ii); /* next cluster: before cluster ii is marked */
      if (equalclusters(iclust,ii0))
        loop (i,ii0,ii0+nsclust) site[i].clust=-1; }

    loop (i,iclust,iclust+nsclust) site[i].clust=-1;

  } /* for (;;) )))))))))))))))))))))))))))))))))))))))))))))))))))))))))) */

} /* printforcefield */


void build(species_t *spec) /***************************************** build */
/***
  assigns the force field to molecule of species spec
***/
{
  int *mark;
  ralloc(mark,sizeof(int)); /* ! like Mark(pointer) in Pascal */

  /* static variables */
  ns=spec->ns;
  site=spec->site;
  specfn=spec->fn;
  allH=spec->opt_h<0;

  buildforcefield(spec);
  if (spec->opt_j) readjet(spec,&b0,&fd0); /* readjet() is in blendedt.c */
  if (spec->Xopt.E) essential(spec);
  else if (spec->Xopt.G) inertiamatrix("gyration",spec);
  else if (spec->Xopt.I) inertiamatrix("inertia",spec);
  else if (spec->Xopt.A||spec->Xopt.D) anglemsd(spec);
  else if (spec->Xopt.V) virial(spec0,spec0->next,spec);
  else {
    if (spec->probe.ns) {
      int i,j,k;
      vector minr,maxr;
      int irg[3];
      double *ri;
      FILE *f;
      site_t *newsite;
      int newns=ns+spec->probe.ns*spec->probe.show; /* LONGER by probe.ns */
      int source=ns-spec->probe.ns;
      float *E;

      if (source!=spec->probe.i) ERROR(("wrong edit data for probing"))

      if (spec->probe.show) {
        alloczero(newsite,sizeof(site_t)*newns);
        copy(newsite,site,sizeof(site_t)*ns);
        alloc(E,spec->probe.show*sizeof(E[0]));
        loop (i,0,spec->probe.show) E[i]=4e33;
        /* site[].r ... temporarily shifted by 1 added molecule */
        loop (k,0,spec->probe.show) {
          static char fmt[]="probe%0?d";
          int offset=source+k*spec->probe.ns,l=8;
          char *id;

          i=spec->probe.show-1; fmt[7]='1';
          while (i>=10) i/=10,fmt[7]++,l++;
          alloc(id,l);
          sprintf(id,fmt,k);

          loop (i,0,spec->probe.ns) {
            site_t *si=&site[source+i];
            site_t *sn=&newsite[offset+i];

            if (k) {
              copy(sn,si,sizeof(site_t));
              loop (j,0,si->nnbr) sn->nbr[j]=si->nbr[j]+k*spec->probe.ns; }
            sn->id=id; } }
        free(site);
        spec->site=site=newsite; /* NOTE: spec->ns updated later! */
      }

      ri=site[spec->probe.i].r;
      strcpy(spec->ext,".pro");
      f=fopen(spec->fn,"wt");
      if (f==NULL) ERROR(("cannot write to %s",spec->fn))
      prt("! writing %s",spec->fn);

      VVV(minr,=maxr,=site[0].r)
      loop (i,1,source) loop (k,0,3) {
        Min(minr[k],site[i].r[k])
        Max(maxr[k],site[i].r[k]) }

      VO(minr,-=spec->probe.shell)
      VO(maxr,+=spec->probe.shell+spec->probe.d/2)

      site[spec->probe.i].keep=POLICE;

      fprintf(f,"# grid=%g shell=%g index of probe=%d -m%d\n",
              spec->probe.d, spec->probe.shell, spec->probe.i, spec->opt_m);
      fprintf(f,"# from x y z = %g %g %g\n",minr[0],minr[1],minr[2]);
      fprintf(f,"#  to  x y z = %g %g %g\n",maxr[0],maxr[1],maxr[2]);

      loop (k,0,3) irg[k]=(maxr[k]-minr[k])/spec->probe.d+1;
      fprintf(f,"# size = %d x %d x %d = %d\n",
              irg[0],irg[1],irg[2], irg[0]*irg[1]*irg[2]);

      for (ri[0]=minr[0]; ri[0]<maxr[0]; ri[0]+=spec->probe.d) {
        fprintf(stderr,"%.0f%%\r",(ri[0]-minr[0])*(100/(maxr[0]-minr[0])));
        for (ri[1]=minr[1]; ri[1]<maxr[1]; ri[1]+=spec->probe.d)
          for (ri[2]=minr[2]; ri[2]<maxr[2]; ri[2]+=spec->probe.d) {
            /* !!! there is a special patch for water inside minimize !!! */
            double m=minimize(spec,0,INJAIL);

            fprintf(f,"%9.5f %9.5f %9.5f %12.5g\n", ri[0],ri[1],ri[2],m);
            loop (j,0,spec->probe.show) if (m<E[j]) {
              for (i=spec->probe.show-1; i>j; i--) {
                loop (k,0,spec->probe.ns)
                  copy(site[ns+i*spec->probe.ns+k].r,
                       site[ns+(i-1)*spec->probe.ns+k].r, sizeof(site[0].r));
                E[i]=E[i-1]; }
              E[j]=m;
              loop (k,0,spec->probe.ns)
                copy(site[ns+j*spec->probe.ns+k].r,
                     site[source+k].r, sizeof(site[0].r));
              break; }

#if 0
            loop (i,0,spec->probe.show) prt("%2d %9.5f %9.5f %9.5f %g",i,
                                            site[3*i+ns].r[0],
                                            site[3*i+ns].r[1],
                                            site[3*i+ns].r[2],E[i]);
            _n
#endif
          } }
      fclose(f);

      if (spec->probe.show) {
        loop (i,ns,newns)
          copy(site[i-spec->probe.ns].r,site[i].r,sizeof(site[0].r));

        strcpy(spec->ext,".pwm");
        f=fopen(spec->fn,"wt");
        if (f==NULL) ERROR(("cannot write to %s",spec->fn))
        prt("! writing %s",spec->fn);
        strcpy(spec->ext,".pro");
        fprintf(f,"# %d minimum energies from %s -m%d\n",
                spec->probe.show,spec->fn,spec->opt_m);
        fprintf(f,"#     x         y         z       E [kcal/mol]\n");
        loop (i,0,spec->probe.show) {
          ri=site[i*spec->probe.ns+source].r;
          fprintf(f,"%9.5f %9.5f %9.5f %12.5g\n", ri[0],ri[1],ri[2],E[i]); }
        fclose(f);

        strcpy(spec->ext,".probe.mol");
        spec->ext=spec->fn+strlen(spec->fn)-4;
        spec->ns=newns-spec->probe.ns;
        if (!spec->opt_p) spec->opt_p=1;
        writeMOL(spec);
        write3D(spec,3);
        spec->opt_p=0;
        free(E); } }
    else
      minimize(spec,1,INJAIL);
  /*.....       minimize(spec,1,FREE);*/
    if (spec->Xopt.N) normalmodes(spec);
  }

  masscenter(spec,2);
  printforcefield(spec);

  release(mark); /* ! like in Pascal: release lists a0-> b0-> d0->... */
  a0=NULL;
  b0=NULL;
  d0=i0=ar0=cis0=NULL;
}
