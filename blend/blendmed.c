/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLENDMED.C

This module contains modules for reading and writing molecules (specied) in
various formats (MOL,CHE,PDB,ATM) and related stuff.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#include "ground.h"
#include "sds.h"
#include "varfile.h"
#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "blendmed.h"
#include "options.h"
#include <time.h>
#include <stdarg.h>
#include <unistd.h>

#ifdef POLAR
#  ifndef SS_MEASURE_rep
#    define SS_MEASURE_rep /*void*/
#  endif /*# SS_MEASURE_rep */
#endif /*# POLAR */

#define BADATOM 1.1e10 /* bad coordinate */
#define ZERODIST 0.01  /* dist. from the plane treated as zero to calc chir. */

/* #define RENUMBER : not supported anymore */

int nspec=0,newnspec=0,chargewarned=0;

double Uanglemax;
int ianglemax;
int PDBstyle=4; /* H style: 0=HCB1 1=HB1 2=1HCB 3=1HB 4=omit */

           /* dr  Jeps  amplitude T   core  phi   psi omega ropt toframe byframe */
Xopt_t Xopt={1e-5,1e-11,1.0,      300,2.0,  {-57, -47, 180}, 0,   -1,     1,
/* -Options: A, C,D,E,F,G,N,I,P,S  ,V,W,fn */
             0,-1,0,0,9,0,0,0,0,200,0,3,NULL};

species_t *spec0;
static species_t *spec;

/* file,line,mygetline are local here */
/* getnoncommentline is in blendpar */

static FILE *file;
static char *line;

static char *mygetline(void) /**************************************** getline */
{
  return getnoncommentline(file,line);
}


static site_t *site, *sitei;

static int diff(vector dr,int i1,int i0) /**************************** diff */
{
  if (site[i0].keep==WANTED || site[i1].keep==WANTED) return 0;
  VVV(dr,=site[i1].r,-site[i0].r)

  return 1;
}

void chiralities(species_t *spec,int option_c) /**************** chiralities */
/***
    Calculates chiralities from configuration (for atoms with 3 neighbors).
    Definition of the chirality of a three-bonded asymmetric atom i:
    ([0] stands for spec->site[i].nbr[0], etc.)

    [2]
       \
        i--[1]
       /
    [0]

    Positive if i is below the plane [0]-[1]-[2]
    Negative if i is above the plane [0]-[1]-[2]

    NOTE: "chiralities" may be (due to deviation from planarity) assigned
      also to atoms which are not chiral.  This is fixed in build.
    NOTE: When the force field is constructed, atoms are distinguishable.
      Thus the `chirality' is treated even if [0]=[1] etc. when the
      atom is not chiral from the chemical point of view
    option_c:
      0: check is made on bad signs; undefined chiralities are trapped in build
      1: undefined chiralities are assigned, defined ones are checked
      2: all chiralities are re-assigned
***/
{
  int ns=spec->ns,i;

  site=spec->site;

  loop (i,0,ns) {
    sitei=site+i;
    if (sitei->nnbr!=3) {
      if (option_c) sitei->chir=0; }
    else {
      int j=sitei->nbr[0];
      int k=sitei->nbr[1];
      int l=sitei->nbr[2];
      int sitei_chir=sitei->chir;
      int chir; /* calculated from cfg */
      vector k_j,l_j,vp;
      double t;

      if (diff(k_j,k,j) && diff(l_j,l,j)) {
	VECT(vp,k_j,l_j)
        /* note: vp points up plane (if [0]-[1]-[2]==j k l are as above) */
	t=sqrt(SQR(vp));
	if (t==0) chir=0;
	else {
	  VO(vp,/=t)
          diff(k_j,j,i);
          t=SCAL(vp,k_j); /* k_j = j-i */
          chir = fabs(t)<ZERODIST ? 0 : t<0 ? -1 : 1;
          if (option('v')&8)
            prt("%6s dist=%.3f chir=%3i",atom[sitei->type].name,t,chir); } }
      else
	chir=0;

      switch (option_c) {
        case 1: if (sitei_chir) break; /* chirality from mol-file left */
        case 2:
          spec->wrmol++; /* mark as to [re]write *.mol file */
          sitei->chir=chir; /* calculated assigned */
        case 0: break;
        default: ERROR(("-c%d is invalid",option_c)) }

      if (sitei->chir && ((option('v')&4 || option_c<2) && sitei->chir!=chir))
	prt("! %3i %-4s chir=%d != calc.chir=%d",
	    i,atom[sitei->type].name,sitei->chir,chir); } }

  prt("! %s chiralities %s",spec->fn,option_c ? "assigned" : "checked");
}

static void addtolist(void) /********************************** addtolist */
/***
    spec is added to the END of the list (to have output in order)
***/
{
  species_t *s;
  int t=option('t');

  if (spec0) {
    /* finding end of list */
    for (s=spec0; s->next; s=s->next);
    /* appending spec */
    s->next=spec;  }
  else
    /* new list */
    spec0=spec;

  spec->next=NULL;

  /* to the beginning of list:
     spec->next=spec0; spec0=spec;
  */

  spec->opt_m=option('m');
  spec->opt_h=option('h');
  spec->opt_w=option('w');
  spec->opt_k=option('k');
  spec->opt_j=option('j');
  spec->opt_y=option('y');
  spec->opt_p=option('p');
  spec->reordw=option('[');
  spec->Xopt=Xopt;

  /* cutoff */
  spec->C2=t+1;
  spec->C1=t-1;
  t=-t;
  if (t>0) {
    spec->C1=(t/100)*.1;
    spec->C2=spec->C1+(t%100)*.1; }

  if (spec->C1>1 && spec->C1<10) prt("!*****WARNING cutoff C1<10");
  spec->N=option('n');

  if ( (spec->opt_u=option('u')) )
    allocarray(spec->maxpair,abs(option('u')));
  else
    spec->maxpair=NULL;

  spec->wrmol=0;
  spec->chargewarned=chargewarned;
  spec->edit=option('e');
}

void cluster(void) /*********************************************** cluster */
/***
    Number of clusters (sub-molecules bound together by chemical bonds)
    is printed and returned.
    Cluster numbers are assigned.
    Normally, there should be one cluster per molecule.
    Possible exceptions: ion pairs, moleculed bound together by an H-bonds.

    Rewritten from project HTH, (c) J.Kolafa 1987 (PL/1).
***/
{
  char *tt;
  int *stack;
  int ns=spec->ns,ncltot=0,n,nt,ncl,nbl,i,l,*newindex;
  int linebreak=16;

  site=spec->site;
  spec->nclust=0;

  ralloc(tt,ns); memset(tt,1,ns);
  ralloc(stack,ns*sizeof(int));

  prts_("! cluster sizes:");

  loop (n,0,ns) if (tt[n]) {
    ncl=1; nt=n; i=0;

  NextMol:
    tt[nt]=0;
    site[nt].clust=spec->nclust;
    loop (l,0,site[nt].nnbr) {
      nbl=site[nt].nbr[l];
      if (nbl>0 && tt[nbl]) {
        i++; stack[i]=nt; nt=nbl; ncl++; goto NextMol; } }

    if (i>0) {
      nt=stack[i]; i--; goto NextMol; }

    /* ncl is the size of the cluster */
    if (ncl>ns) ERROR(("internal"))
    ncltot+=ncl;
    spec->nclust++;
    if (linebreak>77) { linebreak=3; prts_("\n! "); }
    linebreak+=prt_(" %d",ncl); }

  if (ncltot!=ns) ERROR(("internal"))
  _n
  release(tt);
  prt("\n! %d cluster(s) found",spec->nclust);

  if (option('n')>=0) return; /* no molecule split -- ordering irrelevant */

  /* is it necessary to renumber sites so that clusters are contiguous? */
  if (spec->reordw) goto renumber;
  loop (i,1,ns) if (site[i].clust<site[i-1].clust) goto renumber;

/* already properly ordered */
  return;

 renumber: /* RENUMBER: cf. also function fixorder() above */

  spec->wrmol=-32000; /* marks renumbering has occured */

  /* making renumbering table */
  alloc(newindex,ns*sizeof(int));
  l=0;
  loop (ncl,0,spec->nclust)
    loop (i,0,ns) if (spec->site[i].clust==ncl) newindex[i]=l++;

  /* again, with reordering ns=3 to HHO */
  if (spec->reordw) {
    l=0;
    loop (ncl,0,spec->nclust) {
      int cns=0;

      loop (i,0,ns) if (spec->site[i].clust==ncl) cns++;
      if (cns==spec->reordw) {
	loop (cns,0,2)
	  loop (i,0,ns)
            if (spec->site[i].clust==ncl
		&& ((atom[spec->site[i].type].name[0]=='H') ^ cns))
	      newindex[i]=l++; }
      else
	loop (i,0,ns) if (spec->site[i].clust==ncl) newindex[i]=l++; } }

  if (l!=spec->ns) ERROR(("internal"))

  loop (i,0,ns) if (newindex[i]!=i) goto contreord;

  prt("!? no reordering necessary");
  free(newindex);
  return;

 contreord:

  /* physical reordering sites according to cluster index:
     sites in individual molecules may be reordered
     faster, memory consuming  (cf. the old TINY version) */
  {
    site_t *newsite;
    alloc(newsite,ns*sizeof(site_t));

    loop (i,0,ns) newsite[newindex[i]]=site[i];
    free(spec->site);
    site=spec->site=newsite;
  }

  /* renumbering all references */
  loop (i,0,ns)
    loop (n,0,spec->site[i].nnbr)
      spec->site[i].nbr[n]=newindex[spec->site[i].nbr[n]];

  spec->newindex=newindex;

  prt("!*****WARNING INDEXes renumbered because clusters not contiguous and -n-1");
}

void perturb(species_t *spec,vector *r,int wave,double ran) /****** perturb */
/***
    wave>0: adds a wave perturbation (`wave' waves) to the molecule
  [REMOVED:    waves=0: adds a small random perturbation]
    ran: adds a random perturbation
    if r==NULL then spec->site[].r is used else r[]
***/
{
  int i,j,x,z;
  vector base;
  double A=0.5,f=0;
  double *ri;

  if (wave) {
    VO(base,=0.2)
    loop (j,0,wave) {
      x=j%3; z=(x+2)%3; /* to start waving z */
      base[x]/=-1.618;
      A*=1.2;
      loop (i,0,spec->ns) if (spec->site[i].keep==FREE) {
        ri = r ? r[i] : spec->site[i].r;
        f=SCAL(base,ri);
        ri[x]+=A*sin(f);
        ri[z]+=A*cos(f); } }
    prt("! waved %d times",wave); }

  if (ran) {
    loop (i,0,spec->ns) if (spec->site[i].keep==FREE) {
      ri = r ? r[i] : spec->site[i].r;
      VO(ri,+=ran*rndgauss()) }
    prt("! perturbed stdev=%g",ran); }
}

int checkneighbors(species_t *spec,int ns) /**************** checkneighbors */
/***
    consistency check for neighbors bound by bonds
    accepts non-contiguous numbering of sites upto ns-1 (?)
***/
{
  int i,nbr,n,m,k,ret=0;
  site_t *s,*site=spec->site;

 again:

  loop (i,0,ns) if ( (s=&site[i])->type>=0 ) {
    loop (n,0,s->nnbr) {
      nbr=s->nbr[n];
/*.....    prt("%s - %s",atom[s->type].name,atom[spec->site[nbr].type].name);*/

      if (nbr<0 || nbr>=ns || spec->site[nbr].type<0) {
	ERROR(("%s[%d]-->?[%d]: bonded neighbor does not exist",
	       atom[s->type].name,i,nbr));
        return 2; }
      else {
	loop (m,0,spec->site[nbr].nnbr)
	  if (i==spec->site[nbr].nbr[m]) goto OK;
	prt("! %s[%d] --> %s[%d] : no <-- bond - adding",
	    atom[s->type].name,i,
	    atom[spec->site[nbr].type].name,nbr);
	ret=1;
	if (spec->site[nbr].nnbr>=MAXVAL) {
	  ERROR(("cannot add - too many bonds"))
	return 2; }

	spec->site[nbr].nbr[spec->site[nbr].nnbr++]=i;
      OK:; }

      loop (m,0,n) if (s->nbr[m]==nbr) {
	/* there is multiple bond i-nbr: will be removed */

	prt("! %s[%d]===%s[%d] multiple bonded - fixed",
	    atom[s->type].name,i,
	    atom[spec->site[nbr].type].name,nbr);
	ret=1;

	s->nnbr--;
	loop (k,m,s->nnbr) s->nbr[k]=s->nbr[k+1];
	goto again; } } }

  return ret;
}

#ifdef RENUMBER
static void fixorder(void) /************************************** fixorder */
/***
    sites in a molecule of species spec are renumbered
    to go contiguously from 0 by 1 to spec->ns-1
***/
{
  int *newindex,i,j,n;

  prt("!*****WARNING INDEXes renumbered because they were not contiguously numbered");

  /* making renumbering table */
  alloc(newindex,spec->maxns*sizeof(int));
  j=0;
  loop (i,0,spec->maxns) if (spec->site[i].type>=0) newindex[i]=j++;
  if (j!=spec->ns) ERROR(("internal"))

  /* renumbering */
  j=0;
  loop (i,0,spec->maxns) if (spec->site[i].type>=0) {
    if (j<i) spec->site[j]=spec->site[i];
    loop (n,0,spec->site[j].nnbr)
      spec->site[j].nbr[n]=newindex[spec->site[j].nbr[n]];
    j++; }

  free(newindex);
}
#endif /*# RENUMBER */

static int indx;

static char *next(void) /********************************************* next */
{
  char *c=strtok(NULL," \t\n");

  if (!c) {
    put(indx) ERROR(("bad line")) }

  return c;
}

static int newspec(char *fn) /************************************* newspec */
/* fn is always appended: .che or .mol */
{
  char *c;

  file=fopen(fn,"rt");
  prt("\n! file %s",fn);
  if (file==NULL) return 0;

  nspec++;
  alloczero(spec,sizeof(species_t));
  spec->newindex=NULL;
  alloc(line,option('l'));

  /* longer (instead of .ext): .probe.mol, essmatch.mol., .nmf%04d.plb, +2 for sure */
  alloc(spec->fn,strlen(fn)+12);
  strcpy(spec->fn,fn);
  spec->ext=NULL;
  for (c=spec->fn; *c; c++) if (*c=='.') spec->ext=c;
  if (!spec->ext) ERROR(("newspec: no . found in %s",fn))
  *(spec->ext)=0;
  if (!mygetline()) ERROR(("bad format"))
  alloc(spec->info,strlen(line)+1);
  strcpy(spec->info,line);

  return 1;
}

static void newsite(void) /**************************************** newsite */
{
  int i,ns=spec->ns;
  int4 size=(int4)ns*sizeof(site_t);

  alloc(site,size);
  spec->site=site;
  memset(site,0,size);
  loop (i,0,ns) {
    site[i].type=-1;
    site[i].keep=FREE; }
}

void newemptyspec(void) /************************************** newemptyspec */
/* add empty species - to be filled, e.g., by appendspec */
{
  nspec++;
  alloczero(spec,sizeof(species_t));
  spec->info="created by blend";
  spec->fn=strdup("aux\0aux"); /* should be writable */
  spec->ext=spec->fn+3;
  spec->newindex=NULL;
  addtolist();
  spec->wrmol=0;
  spec->opt_w=spec->opt_p=0;
}

void appendspec(species_t *to,species_t *app) /****************** appendspec */
/* species *app is appended to species *to */
{
  int ii,i,oldns=to->ns;
  site_t *old=to->site;

  spec=to;
  spec->ns+=app->ns;

  newsite();

  site=spec->site;
  if (old) {
    memcpy(site,old,oldns*sizeof(site_t));
    free(old); }
  loop (ii,oldns,spec->ns) {
    site[ii]=app->site[ii-oldns];
    loop (i,0,site[ii].nnbr) site[ii].nbr[i]+=oldns; }
}

species_t *readMOL(char *fn) /************************************* readMOL */
/***
    reads a file in the format originally devised for the "Molekyle Editor"
    fn must have extension .mol
    returns the pointer on success
***/
{
  char *tok,*slash;
  int ns=0,foundns=0,i,j,wrmol=0;
  site_t *s;

  site=NULL;

  if (!newspec(fn)) return NULL;
  spec->zero_energy=0;

 try_again:

  if (foundns) spec->ns=ns=foundns;

  while (mygetline()) if ( (tok=strtok(line," =\t\n")) ) {

    if (!strcmp(tok,"parameter_set"))
      readPARBIN(strtok(NULL," =\t\n"));

    else if (!strcmp(tok,"zero_energy"))
      spec->zero_energy=atof((strtok(NULL," =\t\n")));

    else if (!strcmp(tok,"number_of_atoms")) {
      if (!foundns) {
	tok=strtok(NULL," =\t\n");
	spec->ns=ns=foundns?foundns:atoi(tok);
	prt("! number_of_atoms = %d",spec->ns); } }

    else if (!strcmp(tok,"atoms")) {
      if (ns<=0) {
	prt("! undefined or bad number_of_atoms - trying to count them");
	if (foundns) ERROR(("%s: cannot determine number_of_atoms",spec->fn))
        while (mygetline()) {
	  tok=line;
	  if (*tok=='*') tok++;
	  indx=atoi(strtok(tok," \t"));
	  if (foundns!=indx)
	    ERROR(("%s: atoms are not contiguously numbered %d %d",spec->fn,foundns,indx))
          foundns++; }
        if (!foundns) ERROR(("%s: no atom found",spec->fn))
        rewind(file);
	goto try_again; }

      if (!PARfn) readPARBIN(DEFAULTPAR);

      newsite();

      loop (i,0,ns) {
	s=site+i;
	if (!mygetline()) ERROR(("%s: missing lines",spec->fn))
        tok=line;
	if (*tok=='*') { tok++; s->keep=INJAIL; }
	indx=atoi(strtok(tok," \t"));
	if (i!=indx) ERROR(("%s: atoms are not contiguously numbered %d=%d",
			    spec->fn,i,indx))
        tok=next();
        slash=strchr(tok,'/'); /* Crad before / */
        if (slash) tok=slash+1;
        else wrmol++;
        alloc(s->id,strlen(tok)+1);
        strcpy(s->id,tok);
	s->type=atomn(next());
	s->charge=atof(next())*(option('q')/100.0);
	s->chir=atoi(next());
	if ( (s->nnbr=atoi(next())) > MAXVAL )
	  ERROR(("%s: %d %s too many bonds",spec->fn,i,atom[s->type].name))
        loop (j,0,s->nnbr) s->nbr[j]=atoi(next()); } }

    else
      WARNING(("%s in unknown keyword in %s",tok,fn)) }

  if (!site) ERROR(("%s: no atom",spec->fn))
  checkneighbors(spec,ns);
#ifdef RENUMBER
  if (badorder) {
      fixorder();
      checkneighbors(spec,ns); }
  if (badorder) ERROR(("%s: bad order"))
#endif /*# RENUMBER */

  addtolist();
  cluster();
  free(line);
  fclose(file);

  spec->wrmol=wrmol;

  return spec;
}

/*****************************************************************************/
/*                             chemical format                               */
/*****************************************************************************/

/* see also code in sim/clusters.c */

/* 'x' no longer means crossing bonds: use '+' */
const char CHESEP[]=" -/\\=|+\t\n";
static struct coor_s { int x0,x1,y; } *coor;
static int ns;

static int siten(int x,int y) /*************************************** siten */
{
  int i;

  loop (i,0,ns)
    if (y==coor[i].y)
      if (coor[i].x0<=x && x<=coor[i].x1) return i;

  ERROR(("%s: bad bond in rel.line=%d, col.=%d",spec->fn,y,x))

  return 0;
}

static char **screen;
static int si;

static void findnbr(int x,int y,int c,int dx,int dy) /************** findnbr */
{
  if (screen[y][x]!=c) return;

  do {
    x+=dx; y+=dy;
  } while (screen[y][x]==c || screen[y][x]=='+'); /* 'x' removed */

  if (spec->site[si].nnbr>=MAXVAL) {
    prts(atom[spec->site[si].type].name);
    ERROR(("%s: too many bonds",spec->fn)) }

  spec->site[si].nbr[spec->site[si].nnbr++]=siten(x,y);
}

static int atominfo(double *z,double *q,char *id,char *c) /******** atominfo */
/***
  the atom field *c is analyzed, atom number is returned, and
  *z=z-coordinate and *q=charge q, and id=identifier (after `:', if any)
  Examples of valid atom fields:  CH1E^.3:C1  Ovvn.55  CH3E:C2  MKp1
  ^ = atom is up by cca 1A from the plane
  vv = atom is down by cca 2A
  .3 = partial charge = 0.3 e
  p1 = partial charge = 1 e (cannot write MK1)
  n.55 = partial charge = -0.55 e (also m.55)
  NOTE: CH1E^.3C1 is the same as CH1E^.3:C1 (in this case `:' is not necessary)
***/
{
  char *e;
  char A[8],*a=A;
  int sg=1,n;

  *id=0;

  /* atom type */
  for (e=c; (*e>='A' && *e<='Z') || (*e>='0' && *e<='9'); e++,a++) *a=*e;
  *a=0;
  n=atomn(A);

  /* z-coordinate */
  *z=*q=0;
  for ( ; *e=='^' || *e=='v'; e++) if (*e=='^') *z+=1; else *z-=1;

  /* partial charge */
  if (*e=='n' || *e=='m') { sg=-1; e++; }
  else if (*e=='p') e++;
  if (!*e) return n;
  if (*e=='.' || (*e>='0' && *e<='9')) {
    *q=sg*atof(e)*(option('q')/100.0);
    while (*e=='.' || (*e>='0' && *e<='9')) e++; }

  /* identifier */
  if (*e==':') { e++; fprintf(stderr,"%s\n",e); }
  for ( ; *e; e++) *id++=*e;
  *id=0;

  /*.....else { prts(c); Error("bad atom field"); }*/

  return n;
}

static double scale2D(void) /*************************************** scale2D */
{
  static int warned;

   if (option('e')) {
     if (abs(option('e'))<2) {
       if (!warned)
         WARNING(("2D scaling only 1%% - check option -e"))
         warned++; }
     return abs(option('e'))*0.01; }

   return 1.0;
}

static FILE *file0;
static char *connect;
static int lconnect;

static void includeCHE(int pass) /******************************* includeCHE */
/*
  includeCHE() and nextline() handle the #include statement
  It cannot be nested (but easily extendable..)
*/
{
  char *c=strchr(line,'!');

  if (c) *c=0; /* erase from ! to eol; note that line[0]!='!' */

  if (line[0]=='#') {
    c=strtok(line+1,"\t\n \"");
    if (!strcmp("connect",c)) {
      c=strtok(NULL,""); /* until EOL - is this sufficiently portable? */
      lconnect+=strlen(c)+1;
      if (pass) {
        strcat(connect,c); strcat(connect," "); } }
    else if (!strcmp("include",c)) {
      /* !include statement */
      c=strtok(NULL,"\t\n \"");
      if (file0) {
        ERROR(("#include %s: nesting not supported",c))
        line[0]=0; return; }
      file0=file;
      file=fopen(c,"rt");
      if (!file) {
        ERROR(("#include %s: file not found",c))
        file=file0; file0=NULL; line[0]=0; return; }
      else
        mygetline(); }
    else {
      WARNING(("unknown #%s directive in che-file",c))
      line[0]=0; return; } }
}

static char *nextline(char *line) /******************************** nextline */
{
  char *g=line?line:mygetline();

  if (!g && file0) {
    fclose(file);
    file=file0; file0=NULL;
    g=mygetline(); }

  return g;
}

int connectsites(site_t *site, int n1,int n2) /**************** connectsites */
/*
  connects site[n1] to site[n1]
  returns 0 on error, 1 on success
  Inefficiency: see also command "ab" in blendedt
  Bug: does not detect double bonding: may lead to `too many atoms to
       bond', but if passes, is fixed later
*/
{
  site_t *s=&site[n1];

  if (s->nnbr>=MAXVAL) {
    ERROR(("%s: too many atoms to bond",s->id))
    return 0; }
  s->nbr[s->nnbr++]=n2;

  return 1;
}

static int findsiteCHE(site_t *site, char *id) /**************** findsiteCHE */
/* code repeats: see also findsite() in blendedt.c */
{
  int i;

  loop (i,0,ns) if (!strcmp(id,site[i].id)) return i;
  ERROR(("#connect %s - unknown atom id",id))

  return -1;
}

species_t *readCHE(char *fn) /************************************** readCHE */
/***
  reads *.che file in the chemical formula format.
  One file must contain just one molecule.
  Two passes through the file needed!
  If option -e100 (default) 1 Angstrom is represented by  1 line (y-direction)
    or 2.5 char in x-direction, otherwise scaled by -e%
***/
{
  int nx,ny,i,is,x,par=0;
  char *c,id[24];
  double *r,z,q,scale=scale2D();

  ns=0;

  if (!newspec(fn)) return NULL;

  /* now line=1st info line of the che-file*/
  do mygetline(); while (line[0]==0);
  /* now line=next non-comment and not empty line */
  /* cannot use strtok just now - writes to line ! */
  c=strstr(line,"parameter_set");
  if (!c) c=strstr(line,"parset");
  if (c) {
    readPARBIN(strtok(c+14," =\t\n"));
    par++;
    mygetline(); }

  if (!PARfn) readPARBIN(DEFAULTPAR);

/***
  first pass:
    # of sites calculated
    # of lines & columns calculated
***/
  nx=ny=0;
  lconnect=1;

  do {
    includeCHE(0);
    Max(nx,strlen(line));

    c=strtok(line,CHESEP);
    while (c) {
      if (c[0]=='*') c++;
      if (c[0]>='A' && c[0]<='Z') {
        ns++;
        atominfo(&z,&q,id,c); /* to check spelling only */ }
      else if (c[0]!='#') ERROR(("reading %s:\n\
*** atom (site) type in UPPERCASE expected, but \"%s\" found",fn,c))
      c=strtok(NULL,CHESEP); }

    ny++; }
    while (nextline(NULL));

  prt("! %d sites  %d lines  %d columns",ns,ny,nx);
  spec->ns=ns;

  nx+=2;

  ralloc(screen,(ny+2)*sizeof(*screen));
  ralloc(screen[0],nx); memset(screen[0],' ',nx);
  ralloc(coor,sizeof(struct coor_s)*ns);
  ralloc(connect,lconnect); connect[0]=0;

  newsite();

/***
  second pass:
***/

  rewind(file);
  mygetline();
  do mygetline(); while (line[0]==0);
  if (par) mygetline();

  is=0;
  loopto (i,1,ny) {
    ralloc(screen[i],nx);
    memset(screen[i],' ',nx);
    if (i==1) nextline(line);
    else if (!nextline(NULL)) ERROR(("internal"))
    includeCHE(1);
    if (strlen(line)>=nx) ERROR(("internal"))
    strcpy(screen[i]+1,line);

    c=strtok(line,CHESEP);
    while (c) {
      enum keep_e keep=FREE;

      if (c[0]=='*') c++,keep=INJAIL;
      if (c[0]>='A' && c[0]<='Z') {
        if (is>=ns) ERROR(("internal"))
        r=site[is].r;
        site[is].type=atominfo(r+2,&site[is].charge,id,c);
        site[is].keep=keep;
        if (id[0]==0) sprintf(id,"%d-%s",is,atom[site[is].type].name);
        alloc(site[is].id,strlen(id)+1); strcpy(site[is].id,id);
        coor[is].x0=(int)(c-line)+!keep;
        coor[is].x1=(int)(c-line)+strlen(c);
        coor[is].y=i;
        r[0]=(coor[is].x0+coor[is].x1-nx)*(0.2*scale);
        r[1] = (ny*0.5-coor[is].y)*scale;
        is++; }
      c=strtok(NULL,CHESEP); } }

  ralloc(screen[ny+1],nx); memset(screen[ny+1],' ',nx);

  /*** finding bonds ***/
  loop (si,0,ns) {

    /* shorter sides */
    findnbr(coor[si].x0-1,coor[si].y,'-',-1,0);
    findnbr(coor[si].x0-1,coor[si].y,'=',-1,0);

    findnbr(coor[si].x1+1,coor[si].y,'-',1,0);
    findnbr(coor[si].x1+1,coor[si].y,'=',1,0);

    /* longer sides + corners */
    loopto (x,coor[si].x0,coor[si].x1) {
      findnbr(x-1,coor[si].y-1,'\\',-1,-1);
      findnbr(x,coor[si].y-1,'|',0,-1);
      findnbr(x+1,coor[si].y-1,'/',1,-1);

      findnbr(x+1,coor[si].y+1,'\\',1,1);
      findnbr(x,coor[si].y+1,'|',0,1);
      findnbr(x-1,coor[si].y+1,'/',-1,1); } }

  /* #connect statements */
  {
    char *tok1,*tok2;
    int n1,n2;

    tok1=strtok(connect," ,\n\t");
    while (tok1) {
      if (!(tok2=strtok(NULL," ,\n\t")) )
        ERROR(("odd number of atoms to #connect in che-file"))
      n1=findsiteCHE(site,tok1);
      n2=findsiteCHE(site,tok2);
      if (n1>=0 && n2>=0) {
        connectsites(site,n1,n2);
        connectsites(site,n2,n1); }
      tok1=strtok(NULL," ,\n\t"); }
  }

  checkneighbors(spec,ns);
  cluster();
  chiralities(spec,2);
  perturb(spec,NULL,3*(option('e')>=0),0);
  addtolist();
  spec->wrmol++; /* marked as to generate *.mol file */
  spec->edit=0; /* cannot .che input AND edit mol-file */
  release(screen);
  free(line);

  return spec;
} /* readCHE */

/*****************************************************************************/
/*                      read/write configutation in binary                   */
/*****************************************************************************/

static void endian(void *c) /**************************************** endian */
/* changes endianess of 4 bytes pointed to by c */
{
  char x;

  x=((char*)c)[0]; ((char*)c)[0]=((char*)c)[3]; ((char*)c)[3]=x;
  x=((char*)c)[1]; ((char*)c)[1]=((char*)c)[2]; ((char*)c)[2]=x;
}

static int quieterr;
static int setframe(int *frame,int nframes,species_t *spec) /****** setframe */
/*
   input numbering: 0=1st frame, ...; -1=last frame
   output numbering: 0=1st frame,..., nframes-1=last frame
   returns 1 if read3D should return 0
*/
{
  if (*frame<0) *frame+=nframes;
  if (*frame>=nframes || *frame<0) {
    if (quieterr) return 1;
    ERROR(("%s: no frame # %d",spec->fn,*frame+1))
    *frame=0; }

  return 0;
}

int read3D(species_t *spec,int key) /******************************** read3D */
/***
  Reads configuration:
  key=0: binary .3db, DEPRECATED from V2.1a
  key=1: playback format .plb (= .3db with header of 2 floats), spec->frame
  key=-1: as key=1, but quietly returns 0 even if no frame
  key=2: ascii .3dt, DEPRECATED from V2.1a
  key=3: ascii .pla (new in V2.1a)
  The extensions are:
  If option('r')<0 then endian is changed on input
  if key>=0 then quietly returns 0 on error/eof, 1 on success
***/
{
  int ii,i,j,n,ns=spec->ns,varL;
  static float header[2]={0,0},r[3];
  double *d;
  static char *exts[4]={".3db",".plb",".3dt",".pla"};
  char *fn;

  quieterr=key<0;
  key=abs(key);
  site=spec->site;

  if (key>3) ERROR(("internal: bad key"))
  if ( (fn=spec->Xopt.fn) ) {
    file=fopen(fn,key<2 ? "rb" : "rt");
    if (!file) return 0; }
  else {
    strcpy(spec->ext,exts[key]);

    file=fopen(fn=spec->fn,key<2 ? "rb" : "rt");
    if (!file) return 0; }

  if (key>=2) alloc(line,option('l'));

  prt_("! reading %s",fn);

  if (key==1 || key==3) {
    int4 filesize;
    int nframes;

    if (key==1)
      fread(header,sizeof(header),1,file);
    else {
      mygetline();
      if (option('r')<0)
        ERROR(("invalid option -r (endian change with ascii input"))
      sscanf(line,"%f%f",header,header+1); }

    if (option('r')<0) endian(header);
    if (header[0]!=ns)
      ERROR(("%f is bad # of sites in the playback file header (expected %d)",
        header[0],ns))

    varL=(header[1]==-3);

    if (key==1) {
      /* reading specified frame implemented for binary file only */
      fseek(file,0,2);
      filesize=ftell(file);
      nframes=(int)(filesize/((int4)sizeof(float)*3*(ns+varL)));
      if (filesize!=(int4)sizeof(float)*3*(ns+varL)*nframes+sizeof(header))
        WARNING(("fractional number of frames in %s (wrong or truncated file)",fn));
      if (setframe(&spec->frame,nframes,spec)) return 0;
      if (setframe(&spec->Xopt.toframe,nframes,spec)) return 0;
      if (spec->Xopt.byframe<=0) spec->Xopt.byframe=1;
      prt_(", frame # %d (to %d by %d) %s",
        spec->frame+1,spec->Xopt.toframe+1,spec->Xopt.byframe,varL?"variable L":"old format");
      fseek(file,(int4)sizeof(float)*3*(ns+varL)*spec->frame+sizeof(header),0); } }

  if (spec->newindex) prt_(", renumbering while reading");

  varL=header[1]==-3;
  if (varL) {
    if (key==1)
      fread(spec->L,sizeof(spec->L),1,file);
    else {
      mygetline();
      sscanf(line,"%f%f%f",spec->L,spec->L+1,spec->L+2); } }
  else
    VO(spec->L,=header[1]);

  loop (ii,0,ns) {

    /* if sites have been renumbered... */
    if (spec->newindex) i=spec->newindex[ii]; else i=ii;

    d=site[i].r;

    if (key<2) { /* key=0,1: binary */
      if (fread(r,sizeof(r),1,file)!=1) ERROR(("%s read error",fn))
      if (option('r')<0) loop (j,0,3) endian(r+j);
      VV(d,=r) /* float -> double conversion */ }

    else { /* key==2,3: ascii */
      if (!mygetline()) {
        ERROR(("%s: unexpected EOF",fn))
        d[0]=1e6; n=3; }
      else
        n=sscanf(line,"%lf%lf%lf",d,d+1,d+2);
      if (n!=3) {
        prt("! line %d bad or empty in .3dt - atom pos. unknown",i);
        site[i].keep=WANTED; } }

    if (option('v')&4) prt_("\n! %3i %12.5g %12.5g %12.5g",i,d[0],d[1],d[2]);

    loop (j,0,3) {
      if (fabs(d[j])>=UNDEFATOM) site[i].keep=WANTED;
      if (fabs(d[j])>BADATOM)
        ERROR(("%s: atom %d %c=%g out of range (endian?)",
               fn,i,'x'+j,d[j])) } }

  prt(", done");

  fclose(file);

  if (spec->newindex) free(spec->newindex);
  if (key==2) free(line);

  return 1;
}

               /*          111111*/
               /*0123456789012345*/
char colors[17]="xxxxxxx  NXCOMSH";
               /*BBGCRMBGGBGCRMYW*/
               /*klryearrrlryeaeh*/

int atomcolor(int type) /***************************************** atomcolor */
/***
  Returns (suggested) atom color as 16-color VGA-like number for given type
  The color is derived from the 1st letter of CHARMM atom name
  ! This is cheap solution - file charmm.acr is NOT read
  Used both by Turbo C graphics and for file MOLECULE.gol
***/
{
  char ch=toupper(atom[type].chem);
  int col;

  if (ch<'A' || ch>'Z') ch=toupper(atom[type].name[0]);

  loop (col,7,16) if (ch==colors[col]) break;
  if (col==16) col=11; /* LIGHTCYAN */
  if (ch=='P') col=14; /* YELLOW, like S */
  if (ch=='Z') col=10; /* GREEN, as X */

  return col;
}

static void backup(species_t *spec,char *bakext,char *ext,int warn) /* backup */
{
  char *name;

  *(spec->ext)=0;
  if (warn) fprintf(stderr,
                    "WARNING: %s%s rewritten, old version backed up as %s%s\n",
                    spec->fn,ext,spec->fn,bakext);
  alloc(name,(int)(spec->ext-spec->fn)+60);
  strcpy(name,spec->fn);
  strcpy(spec->ext,bakext);
  strcat(name,ext);
  if (remove(spec->fn)==0) prt("! %s deleted",spec->fn);
  if (rename(name,spec->fn))
    ERROR(("cannot rename %s to %s",name,spec->fn))
  if (rename(tempname(),name))
    ERROR(("cannot rename %s to %s",tempname(),name))
  prt("! %s renamed to %s, %s written",name,spec->fn,name);
  free(name);
}

#include "ms.h"

double calcRvdW(int type) /**************************************** calcRvdW */
/* calculate
   - the van der Waals radius (r of potential minimum) of given atom type
   - if the minimum does not exist, find 1 kcal/mol level
   - added 2/2011: check for exp-6 problem of negative U for small r (bug found by T.Trnka)
*/
{
  sitesite_t *ss=&sstab[type][type][0];
  double rr,x,y,z,Urep,U=0,f,frep,zz;
  double Ull=0,Ul=0;
  double rrmin=0,rrmax=0,rrkcal=0;
  static double rho1,rho2,rho12; /* for METAL only */

  for (rr=200; rr>0.001; rr*=.99) {
#if defined(POLAR)
    SS_MEASURE_rep
#endif /*# defined(POLAR) */
    SS_MEASURE
    if (Ul<0 && U>Ul && Ul<Ull && rrmin==0) rrmin=rr/.99;
    if (Ul>0 && U<Ul && Ul>Ull && rrmax==0) rrmax=rr/.99;
    if (U>1 && Ul<=1 && rrkcal==0) rrkcal=rr/.995;
    if (U>1e3) break;
    Ull=Ul,Ul=U; }

  if (rrmin && rrmin>rrmax) {
    /* minimum detected, finding more precisely */
    MS_BEGIN(rrmin,1e-5)
      rr=rrmin;
      if (MS_it>20) {
        WARNING(("cannot determine the van der Waals radius for type=%d (precise minimum), 1 AA assumed",type))
        return 1; }
#if defined(POLAR)
      SS_MEASURE_rep
#endif /*# defined(POLAR) */
      SS_MEASURE
      MS_f=f;
    MS_END(rrmin,1)

    return sqrt(rrmin)/2; }

  if (rrkcal) return sqrt(rrkcal)*0.6; /* bit longer than U(r) = 1 kcal/mol */

  prt("! unknown radius of atom type=%d, 1 AA assumed",type);

  return 1.0;
}

void write3D(species_t *spec,int playback) /************************ write3D */
/***
    Write configuration in various binary or text formats.
    WARNING: changed since V2.4a, reverse endian removed, .gol removed, -w160 removed
    playback=1: write .plb
                -p-# write also .gol (legacy: .gol no longer needed
                -p# is vdW radius scaling in % (but -p=-p1 = 70%)
    playback=0: -w1           write .3db (DEPRECATED from V2.1a)
                -w-1          write .3db w. reversed endian
                -w#, #=2..9   write .3dt (xyz ascii)
                -w10[:Hstyle] write .pdb from scratch
                -w20          write .pdb using ID info
                -w30          write -w10 or -w20 (what appropriate)
                -w40          write .atm
                -w80          write .cfg (cook format >= 2.7a)
                -w160         write .cfg (cook format < 2.7, REMOVED)
                can be combined: -w-121 = .3db(reversed) + .atm + .cfg(normal)
                (.3db + .3dt cannot be combined)
***/
{
  int i,chendian,type,ascii=0,opt_w=abs(spec->opt_w),ext_w,torewrite;
  double q=0.70; /* default atom sizes are 70% of vdW radii */
  float r[3];
  char *fmt=(char*)r,*fn; /* [12]: not needed at the same time */
  FILE *f;

  if (abs(spec->opt_p)!=1) q=abs(spec->opt_p)*0.01;

  if (playback) opt_w=1;

  if (spec->opt_p<0) {
    /* write the goal file */
    static const char *colorname[9]={"GRAY",
      "GRAY","BLUE","GREEN","CYAN","RED","MAGENTA","YELLOW","WHITE"};

    strcpy(spec->ext,".gol");
    f=fopen(spec->fn,"wt");

    fprintf(f,"!sphere#1\n! %s\n",spec->fn);
    fprintf(f,"%d\n",spec->ns);
#if 0
    if (spec->probe.ns) {
      loop (i,0,spec->ns) {
        type=spec->site[i].type;
        fprintf(f,"%-7s %5.3f\n",
                i<spec->probe.i?colorname[atomcolor(type)-7]:
                i<(spec->probe.i*3+spec->ns)/4?"WHITE":
                i<(spec->probe.i+spec->ns)/2?"YELLOW":
                i<(spec->probe.i+spec->ns*3)/4?"GREEN":"BLUE",
                atom[type].LJ[0].RvdW*q); } }
    else
#endif /*# 0 */
    loop (i,0,spec->ns) {
      type=spec->site[i].type;
      fprintf(f,"%-7s %5.3f\n",
              colorname[atomcolor(type)-7],
              calcRvdW(type)*q); }
    fclose(f);
    prt("! %s written",spec->fn); }

  ext_w=opt_w/10;
  if (ext_w&3) {
    if ( ((ext_w&3)==3 && isdigit(spec->site[0].id[0])) || (ext_w&3)==1) writePDB1(spec);
    else writePDB2(spec); }
  if (ext_w&4) writeATM(spec);
  if (ext_w&8) writeCFG(spec);
  // REMOVED  if (ext_w&16) writeCFGV26(spec);

  opt_w=opt_w%10;

  if (!opt_w) return;

  if (opt_w>1) ascii=opt_w;

  /*** .plb, .3dt or .3db  ***/
  strcpy(spec->ext,playback ? ".plb" : ascii ? ".3dt" : ".3db");

  torewrite=1;

  f=fopen(spec->fn,"rb");
  if (f) {
    if (fread(r,sizeof(r),1,f)==1) torewrite=0;
    fclose(f); }

  fn=torewrite ? spec->fn : tempname();
  f=fopen(fn,ascii ? "wt" : "wb");
  if (!f) ERROR(("cannot write to %s",fn))

  if (playback) {
    /* header (variable L format) */
    r[0]=spec->ns; r[1]=-3;
    chendian=spec->opt_p<0;
    if (chendian) endian(r);
    fwrite(r,sizeof(float),2,f);
    fwrite(spec->L,sizeof(spec->L),1,f); }
  else
    chendian=spec->opt_w<0;

  sprintf(fmt,"%%%d.%df ",ascii+5,ascii);

  loop (i,0,spec->ns)
    if (ascii) {
      loop (type,0,3) fprintf(f,fmt,spec->site[i].r[type]);
      fputc('\n',f); }
    else {
      VV(r,=spec->site[i].r) /* double -> float conversion */
      if (chendian) {
        endian(r); endian(r+1); endian(r+2); }
      fwrite(r,sizeof(r),1,f); }

  if (fclose(f)<0) ERROR(("cannot close %s",spec->fn))

  if (!torewrite) backup(spec,
    ascii ? ".3dt~" : playback ? ".plb~" : ".3db~",
    ascii ? ".3dt" : playback ? ".plb" : ".3db",0);
  else
    prt("! %s written",spec->fn);
}

double totalcharge(species_t *spec,int clust) /***************** totalcharge */
{
  double charge=0,dtc;
  int i;

  loop (i,0,spec->ns)
    if (clust<0 || site[i].clust==clust) charge+=spec->site[i].charge;
  dtc=fabs(fmod(fabs(charge+0.5),1)-0.5);

  if (dtc > 1e-13) {
    if (!spec->chargewarned++)
      WARNING(("charge=%.13f (of cluster %d) is fractional",charge,clust)) }

  return charge;
}

double qeps;

double roundedcharge(species_t *spec,int clust) /************* roundedcharge */
{
  double charge=0;
  int i;

  loop (i,0,spec->ns)
    if (clust<0 || site[i].clust==clust) {
      char chs[16],*c;

      sprintf(chs,charge_f,site[i].charge);
      c=chs;
      while (*c==' ') c++;
      charge+=atof(c); }
  while (charge>0.5) charge-=1;
  while (charge<-0.5) charge+=1;

  return charge;
}

void writeMOL(species_t *spec) /*********************************** writeMOL */
/*
  Writes the molecule file *.mol
*/
{
  FILE *f;
  int i,ns=spec->ns,j,torewrite;
  double q=70; /* default atom sizes are 70% of vdW radii */

  if (abs(spec->opt_p)!=1) q=abs(spec->opt_p);

  /* NOTE: charge_f must be %X.2F */
  qeps=pow(10.,option('_'));

  if (option('_')<0) {
    /* fix small rounding errors in the charges */
    double dq=roundedcharge(spec,-1),r=0.01/sqrt(ns),maxerr=qeps*sqrt(ns);

    prt("! dq=%g maxerr=%g",dq,maxerr); // ???

    if (fabs(dq)>maxerr)
      WARNING(("fractional charge error dq is too large to be fixed"))
    else if (fabs(dq)>1e-13) {
      while (fabs(dq)>1e-13) {
        loop (i,0,ns) spec->site[i].charge-=dq*(1+rnd()*r)/(2*ns);
        r*=1.01;
        if (fabs(r)>2)
          ERROR(("cannot fix fractional charge error (bad algorithm)"))
        dq=roundedcharge(spec,-1); }
      prt("! fractional charges fixed"); } }

  strcpy(spec->ext,".mol");

  f=fopen(spec->fn,"rt");
  torewrite=f==NULL;
  if (!torewrite) fclose(f);

  f=fopen(torewrite?spec->fn:tempname(),"wt");
  if (f==NULL) ERROR(("cannot write to %s",spec->fn))
  prt("! writing %s",spec->fn);

  fputs(spec->info,f);
  fputc('\n',f);

  fprintf(f,"! total charge =");
  fprintf(f,charge_f,totalcharge(spec,-1));

  fprintf(f,"\n\nparameter_set = %s\n",PARfn);
  if (spec->zero_energy!=0) fprintf(f,"zero_energy = %.12g\n",spec->zero_energy);
  fprintf(f,"number_of_atoms = %d\n\n",ns);

  fputs("atoms\n! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",f);

  loop (i,0,ns) {
    sitei=site+i;
    if (spec->opt_p<0)
      fprintf(f,"%3i %-10s %-4s", i,sitei->id,atom[sitei->type].name);
    else
      fprintf(f,"%3i %c%03d/%-10s %-4s", i,
              "OOBGCRMYW"[atomcolor(sitei->type)-7],
              (int)(calcRvdW(sitei->type)*q+0.5),
              sitei->id,
              atom[sitei->type].name);
    fprintf(f,charge_f, sitei->charge);
    fprintf(f,"%2i %4i ",sitei->chir,sitei->nnbr);
    loop (j,0,sitei->nnbr) fprintf(f," %d",sitei->nbr[j]);
    if (fputc('\n',f)<0) ERROR(("cannot write to %s",spec->fn)) }

  fputc('\n',f);
  if (fclose(f)) ERROR(("cannot write to %s",spec->fn))

  if (!torewrite) backup(spec,".mol~",".mol",1);
  else prt("! %s written",spec->fn);
}

void writeMARKED(species_t *spec,int key) /********************* writeMARKED */
/*
  Writes marked/kept sites (to be read by mol2mol)
*/
{
  FILE *f;
  int i,ns=spec->ns;

  strcpy(spec->ext,key?".keep":".mark");
  site=spec->site;

  f=fopen(spec->fn,"wt");
  prt("! writing %s",spec->fn);

  fprintf(f,"! %s sites (readable by mol2mol, option -@%s)\n",
    key?"kept":"marked",spec->fn);

  loop (i,0,ns) if (key?site[i].keep==1:site[i].count) fprintf(f,"%d\n",i);

  fclose(f);
}

void readMARKED(species_t *spec,int keep) /********************** readMARKED */
/*
  Read marked or (if keep) kept sites: see @M @K -k8 -k4
*/
{
  int i,ns=spec->ns;

  site=spec->site;
  strcpy(spec->ext,keep?".keep":".mark");

  file=fopen(spec->fn,"rt");
  if (!file) return;
  prt("! reading %s",spec->fn);

  alloc(line,option('l'));
  while (mygetline()) {
    sscanf(line,"%d",&i);
    if (i>=0 && i<ns)
      if (keep) site[i].keep=INJAIL;
      else site[i].count=1; }
  fclose(file);
  free(line);
}


#include "mendeleyev.c"

/********** PDB stuff *************/

static int isnbr(int j,site_t *si) /********************************** isnbr */
{
  int n;

  if (j<0) return 1;
  loop (n,0,si->nnbr) if (j==si->nbr[n]) return 1;

  return 0;
}

/* type of atom according to position in the peptide backbone */
enum bb_e {
  BBNONE=0, /* sidechain atom (or H) */
  BBN=1, /* backbone (peptide bond) nitrogen */
  BBCA=2, /* alpha-carbon */
  BBC=3, /* carbonyl carbon */
  BBNTERM=4, /* N-terminus atom (incl. patches) */
  BBCTERM=5, /* C-terminus atom (incl. patches) */
  BBCO=8, /* backbone carbonyl oxygen */
  BBCB=9 /* beta-carbon; used temporarily and later erased */
  } ;

static enum bb_e bbtype(int i) /************************************* bbtype */
{
 /* returns type of atom according to position in the backbone */

 int t;
 enum bb_e j;

 loopto (j,BBN,BBC) loop (t,0,NBBTYPE) if (site[i].type==bbdata[j][t]) return j;
 loop (t,0,NBBTYPE) if (site[i].type==bbdata[BBCO][t]) return BBCO;

 return BBNONE;
}

FILE *pdb;
static char Hnrid[6]="H",*nrid /*=Hnrid+1*/;

static void makeHnrid(int i) /************************************ makeHnrid */
/* generates PDB atom name as e.g. CG1, not for H */
{
 site_t *si=site+i;
 char *at=atom[si->type].name;
 signed char bb=si->backbone;

 if (si->nest>8) si->nest=9;

 sprintf(nrid,"%c%c%c",
         at[0],
         " ABGDEZHX"[si->nest-1],
         " 1234567?"[si->count]);

 if (bb && bb!=2) Hnrid[2]=Hnrid[3]=' ';
}

static void pdbdigit(char *id,char digit) /*********************** pdbdigit */
{
 if (!memcmp(nrid,id,2)) nrid[2]=digit;
}

static char pdbeol[10]="1$$$%4d\n";

static void pdbprintf(const char *format,...) /****************** pdbprintf */
{
 char line[128];
 static int ln=0;
 va_list args;

 va_start(args,format);

 memset(line,' ',80);
 vsprintf(line,format,args);
 *strchr(line,0)=' ';

 sprintf(line+72,pdbeol,++ln);
 fputs(line,pdb);
 va_end(args);
}

static int nonHnbrs(int i) /************************************** nonHnbrs */
/* returns # of non-H neighbors of given site */
{
 int j,n=0;

 loop (j,0,site[i].nnbr) n+=atom[site[site[i].nbr[j]].type].name[0]!='H';

 return n;
}

int globrsd;
static int pdbatomn;

static void putPDB(int i,char *XXX) /******************************* putPDB */
{
 site_t *si=site+i;
 int het=si->rsd==0;
 char *at=atom[si->type].name;
 static char rn[4]="???";
 char aaa[8],*a;
 char *resname;
 int j;

 static struct amino_s {
   char rsd[4]; int key; } amino[]= {
/* decimal digit with given number of atoms, ascending order,
   T=terminus (1=nter, 2=cter) */
/*     TSONC */
{"XXX",    0},
{"HOH",  100},
{"GLY",  112},
{"ALA",  113},
{"VAL",  115},
{"LEU",  116},
{"PHE",  119},
{"PEG",  120}, /* N-subst glycine */
{"LYS",  126},
{"TRP",  131},
{"HIS",  136},
{"ARG",  146},
{"SER",  213},
{"THR",  214},
{"TYR",  219},
{"ASN",  224},
{"GLN",  225},
{"ASP",  314},
{"GLU",  315},
{"HEM",  474}, /* ? */
{"HEM",  674}, /* ? */
{"CYS", 1113},
{"MET", 1115},
{"PRO", 2222}, /* special key, regular=115 */
{"ILE", 3333}, /* special key, regular=116 */

/* common caps - new in 2.0d */
{"ACE",10102},
{"AMI",20010}, /* -NH2  - also named CT2 */
{"CT3",20011}, /* -N-CH3 */
{"CT1",20101}, /* -O-CH3 */
/*     TSONC */
{"",0} };

 for (j=0; amino[j].rsd[0]; j++) if (si->key==amino[j].key) break;
 if (!amino[j].rsd[0]) j=0;
 resname=amino[j].rsd;

#if 0
 /* DEBUG PDB conversion */
 prt("%3d rsd=%-2d clust=%d chain=%c nest=%-2d key=%05d bb=%d %2d=%s %s",
     i,si->rsd,si->clust,si->chain+'`',si->nest,si->key,si->backbone,j,resname,at);
#endif /*# 0 */

 if (!resname) resname=rn;

 if (het) {
   copy(Hnrid,atom[si->type].name,4); Hnrid[3]=0;
   strcpy(rn,Hnrid);
   /* detect water here: */
   if (si->nnbr) {
     int O=i;
     site_t *siO;

     if (si->nnbr==1) O=si->nbr[0];
     if ( (siO=site+O)->nnbr==2
          && atom[siO->type].name[0]=='O'
          && atom[site[siO->nbr[0]].type].name[0]=='H'
          && atom[site[siO->nbr[1]].type].name[0]=='H' ) strcpy(rn,"HOH"); }
   while (strlen(rn)<3) strcat(rn,"X"); }

 if (at[0]!='H') makeHnrid(i);
 else if (si->nnbr) makeHnrid(si->nbr[0]);

 /* some fixes to have conventional numbering */
 if (!strcmp(resname,"HIS")) {
   pdbdigit("CD",'2');
   pdbdigit("ND",'1');
   pdbdigit("NE",'2');
   pdbdigit("CE",'1'); }
 if (!strcmp(resname,"ASN")) {
   pdbdigit("OD",'1');
   pdbdigit("ND",'2'); }
 if (!strcmp(resname,"GLN")) {
   pdbdigit("OE",'1');
   pdbdigit("NE",'2'); }
 if (!strcmp(resname,"THR")) {
   pdbdigit("OG",'1');
   pdbdigit("CG",'2'); }
 if (!strcmp(resname,"ILE")) {
   if (!memcmp(nrid,"CG",2)) nrid[2]='3'-nonHnbrs(i); }
 if (!strcmp(resname,"TRP")) {
   if (!memcmp(nrid,"CD",2)) nrid[2]='0'-1+nonHnbrs(i);
   pdbdigit("NE",'1');
   if (!memcmp(nrid,"CE",2)) nrid[2]='5'-nonHnbrs(i);
   if (!memcmp(nrid,"CZ",2)) {
     int j,n;
     loop (j,0,si->nnbr)
       if (si->nest-1==site[n=si->nbr[j]].nest) nrid[2]='5'-nonHnbrs(n); }
   pdbdigit("CH",'2'); }

 if (!strcmp(resname,"ACE")) {
   if (nrid[0]=='C' && nonHnbrs(i)==1) nrid[1]='A'; }

 if (!memcmp(resname,"CT",2) || !strcmp(resname,"AMI")) {
   if (!memcmp(nrid,"N ",2)) nrid[1]='T';
   if (nrid[0]=='O') nrid[1]='T',nrid[2]='2';
   if (nrid[0]=='C' && nonHnbrs(i)==1) nrid[1]='T'; }

/*
ATOM      1  N   THR     1      17.047  14.099   3.625  1.00 13.79      1CRN  70
111111222223334444445555556666666666667777777788888888999999000000-------===xxxx
ATOM    255 1HE2 GLN    16      10.127   0.799 -10.791  1.00  2.18      1FRC 352
*/
#if 0
 pdbprintf("%6s%5d %c%-4s%3s%6d%12.3f%8.3f%8.3f%6.2f%6.2f",
           het?"HETATM":"ATOM  ",++pdbatomn,
           at[0]=='H'?" 1234567?"[si->count]:' ',
           Hnrid+(at[0]!='H'),
           resname,globrsd,
           si->r[0],si->r[1],si->r[2],1.00,0.0);
#else /*# 0 */
 if (at[0]=='H') {
   if (PDBstyle&1) {
     /* strip off nbr but N */
     if (Hnrid[1]!='N' || Hnrid[2]!=' ') Hnrid[1]=Hnrid[2],Hnrid[2]=' '; }
   if (PDBstyle&2)
     sprintf(aaa,"%c%s"," 1234567?"[si->count],Hnrid);
   else {
     sprintf(aaa," %s",Hnrid);
     a=strchr(aaa+1,' ');
     if (a) *a=" 1234567?"[si->count]; } }
 else
   sprintf(aaa," %s",Hnrid+1);
 aaa[5]=0;

 if (XXX) strcpy(aaa,XXX);

 if (PDBstyle&4 && at[0]=='H') return; /* omit H */

 pdbprintf("%6s%5d %-5s%3s %c%4d%12.3f%8.3f%8.3f%6.2f%6.2f",
           het?"HETATM":"ATOM  ",++pdbatomn,
           aaa,
           resname,si->chain+'@',globrsd,
           si->r[0],si->r[1],si->r[2],1.00,0.0);
#endif /*#!0 */
}


static void marknbrx(site_t *si,int rsd,int chain,int nest) /**** marknbrx */
/* marks all yet unmarked neighbors of *si by rsd and chain (recursive) */
{
  int n;

  if (si->rsd) return; /* already marked */

  if (nest>12) ERROR(("marknbrx: nesting exceeds limit (not peptide)"))
  si->chain=chain;
  si->rsd=rsd;

  /* WARNING: check whether proline causes problems ... */

  /* sulphur-sulphur breaks the chain */
  if (atom[si->type].name[0]=='S')
    loop (n,0,si->nnbr)
      if (atom[spec->site[si->nbr[n]].type].name[0]=='S') return;

  loop (n,0,si->nnbr) marknbrx(spec->site+si->nbr[n],rsd,chain,nest+1);

  return;
}

static void markfrom(site_t *si,int rsd,int chain) /************** markfrom */
/* marks all yet unmarked neighbors of *si by rsd and chain (start here) */
{
  si->rsd=0;
  marknbrx(si,rsd,chain,1);
}

static void markcalpha(site_t *si, int nest, int maxnest)
/*
  recursive marking of the distance from Calpha
  NEW: from N, N=1, alpha=2, etc.
*/
{
  int n;

  if (nest>maxnest) return;
  if (si->rsd!=globrsd) return;

  if (si->nest==0) si->nest=nest;
  if (si->nest!=nest) return;

  loop (n,0,si->nnbr) markcalpha(&spec->site[si->nbr[n]],nest+1,maxnest);
}

static int thesameid(int i,int ns) /***************************** thesameid */
{
 int j;
 site_t *si=site+i,*sj;

 loop (j,0,ns) if (i!=j)
   if (globrsd==(sj=site+j)->rsd
    && sj->count==0
    && atom[si->type].name[0]==atom[sj->type].name[0]
    && si->nest==sj->nest
    && si->backbone==sj->backbone) {
     if (si->nest==2 && atom[si->type].name[0]=='H')
       if (si->nnbr==1 && sj->nnbr==1 && si->nbr[0]==sj->nbr[0]) return j;
       else return -1;
     else return j; }

 return -1;
}

static int torewrite;

static void openpdb(species_t *spec,int frame) /******************** openpdb */
{
 int i;

 if (frame<=0) {
   strcpy(spec->ext,".pdb");
   pdb=fopen(spec->fn,"rt");
   torewrite=pdb==NULL;
   if (!torewrite) fclose(pdb); }
 else {
   torewrite=1; /* ==> does not make backup */
   sprintf(spec->ext,".%04d.pdb",frame); }

 loop (i,0,3)
   pdbeol[i+1]=isalpha(spec->fn[i])?toupper(spec->fn[i]):'X';

 pdbatomn=0;

 pdb=fopen(torewrite?spec->fn:tempname(),"wt");
 if (pdb==NULL) ERROR(("cannot write to %s",spec->fn))
 prt("! writing %s",spec->fn);

 pdbprintf("HEADER %s",spec->info);
 pdbprintf("REMARK total charge = %.3f",totalcharge(spec,-1));
 pdbprintf("REMARK parameter_set = %s",PARfn);
 if (spec->zero_energy!=0)
   pdbprintf("REMARK zero_energy = %.12g",spec->zero_energy);
 pdbprintf("REMARK number_of_atoms = %d",spec->ns);
}

static void closepdb(void) /************************************** closepdb */
{
 pdbprintf("END");
 if (fclose(pdb)) ERROR(("cannot write to %s",spec->fn))

 if (!torewrite) backup(spec,".pdb~",".pdb",1);
 else prt("! %s written",spec->fn);
}

static int makersdinfo(species_t *spec) /********************** makersdinfo */
/*
 - breaks a peptide into chains (takes into account S-S bonds)
 - finds the backbone
 - numbers residues (incl. NTER,CTER patches), consecutively
   (new chain does not restart numbering)
 - using # of atoms + topology, assigns keys to residues (later, keys
   are translated to residue names)
 - peptide MUST be in the standard order and any N (of the same type
   as peptide bond N) MUST NOT precede the peptide N
*/
{
  site_t *si;
  int ns=spec->ns,i,n,rsd=0,rsd0,endkey=-1;
  site_t *nitrogen,*calpha,*oc;

  site=spec->site;

  ns=spec->ns;

  /*** backbone analysis ***/
  spec->nchains=0;

  loop (i,0,ns) {
    si=site+i;
    si->backbone=si->count=0;
    si->nest=0;
    si->rsd=0;
    si->chain=0; }

 Chain: /* starting new chain */

  rsd0=rsd;
  spec->nchains++;

  n=0;
  /* finding first triple N-CA-CO and nter */
  loop (i,0,ns) {
    si=site+i;

    if (si->chain>0) continue; /* ignore previous chains */

    if (bbtype(i)==BBCA) {
      int iN=-1,iC=-1,bb,nx,j;

      loop (j,0,si->nnbr) {
        bb=bbtype(nx=si->nbr[j]);
        if (bb==BBC) iC=nx;
        if (bb==BBN) iN=nx; }

      if (iC>=0 && iN>=0) {
        /* N-CA-CO found, finding nterm */
        nx=min(iC,iN); nx=min(nx,i);
        while (nx--) {
          si=site+nx; /* (si re-used) */

	  if (si->chain && si->chain!=spec->nchains) continue;

          if (atom[si->type].name[0]!='H') {
            si->backbone=BBNTERM;
            si->chain=spec->nchains;
            si->rsd=rsd=rsd0+1;
            n=1;
            loop (j,0,si->nnbr) {
              site_t *six=site+si->nbr[j];

              if (atom[six->type].name[0]=='H') {
                six->backbone=BBNTERM;
                six->chain=spec->nchains;
                six->rsd=rsd; } } } }
        i=iN;
        goto Nitrogen; }
      } }

  /* no BBN found ==> all chains have been found */

  loopto (globrsd,1,rsd) {
    int key=0, CB=-1,nCO=0;

    loop (i,0,ns) if (globrsd==(si=site+i)->rsd) {
      if (si->backbone==BBNTERM) if (key<3333) key+=10000;
      if (si->backbone==BBCTERM) if (key<3333) key+=20000;
      if (si->backbone==BBCO) nCO++;
      switch (atom[si->type].name[0]) {
        case 'H': break;
        case 'N': key+=10; break;
        case 'C': key++; break;
        case 'O': key+=100; break;
        case 'S': key+=1000; break;
        default: key=40000; }
      /* mark the sidechain first, CA=2 (this is because of PRO) */
      if (si->backbone==BBCA) loop (n,2,9) markcalpha(si,2,n);
      /* mark the N-sidechain (peptoids) */
      if (si->backbone==BBN) loop (n,1,9) markcalpha(si,1,n);
      if (si->backbone==BBCB) CB=i, si->backbone=BBNONE; }

    if (nCO>1) key-=100; /* -COOH changed into -CO */

    /* LEU-ILE and GLY,VAL-PRO problem
       CB is the index of Cbeta:
       if there is no Cbeta (GLY), CB=-1
       rings: any of Cbeta's */
    if (CB>=0) {
      if (key==116) if (nonHnbrs(CB)==3) key=3333; /* ILE, not LEU */
      if (key==115) if (nonHnbrs(CB)==2) key=2222; /* PRO, not VAL */
      }

    loop (i,0,ns) if ( globrsd==(si=site+i)->rsd) si->key=key; }

  /* RETURN */
  return rsd;

  /* ===================================================================== */

 Nitrogen: /* N (of peptide bond) found */
  endkey=-1;
  (si=site+i)->backbone=BBN;
  nitrogen=si;
  rsd++;
  si->rsd=rsd;
  si->chain=spec->nchains;
  si->nest=1; /* to prevent wrong PRO sidechain numbering */

  /* oops - mark hydrogens bound to this N */
  loop (n,0,si->nnbr) {
    site_t *s=&site[si->nbr[n]];

    if (atom[s->type].name[0]=='H') {
      s->rsd=rsd;
      s->backbone=BBNONE;
      s->chain=spec->nchains;
      /*s->nest=2;*/ } }

  loop (n,0,si->nnbr)
    if (site[i=si->nbr[n]].backbone==0 && bbtype(i)==BBCA) {
      /* problems with PRO: Calpha must have C neighbor */
      int j;

      loop (j,0,site[i].nnbr)
        if (site[site[i].nbr[j]].backbone==0 && bbtype(site[i].nbr[j])==BBC)
          goto Calpha; }

  /* Cterm, not full residue */
  markfrom(nitrogen,rsd,spec->nchains);
  goto End;

 Calpha:
  (si=site+i)->backbone=BBCA;
  calpha=si;

  si->rsd=rsd; si->chain=spec->nchains;

  /* here we mark all neighbors of N */
  markfrom(nitrogen,rsd,spec->nchains);

  /* finding Cbeta */
  loop (n,0,si->nnbr) {
    i=si->nbr[n];
    if (site[i].backbone==0 && bbtype(i)!=BBC && atom[site[i].type].name[0]=='C')
      site[i=si->nbr[n]].backbone=BBCB; }

  loop (n,0,si->nnbr)
    if (site[i=si->nbr[n]].backbone==0 && bbtype(i)==BBC) goto Carbonyl;

  /* Cterm, not full residue */
  markfrom(calpha,rsd,spec->nchains);
  goto End;

 Carbonyl:
  (si=site+i)->backbone=BBC;
  si->rsd=rsd; si->chain=spec->nchains;
  /*si->nest=1;*/

  /* marking sidechain (all neighbors of Calpha but S-S) */
  markfrom(calpha,rsd,spec->nchains);

  /* marking carbonyl O */
  endkey=0;
  loop (n,0,si->nnbr) {
    oc=&site[si->nbr[n]];
    if (bbtype(si->nbr[n])==BBCO) {
      oc->chain=spec->nchains;
      endkey++;
      oc->backbone=BBCO;
      /*oc->nest=2;*/
      oc->rsd=rsd; } }

  loop (n,0,si->nnbr)
    if (site[i=si->nbr[n]].backbone==0 && bbtype(i)==BBN)
      goto Nitrogen;

 End: /* C-terminus */

  /* endkey: -1: partial residue, >=0: # of carbonyl/carboxyl O */
  switch (endkey) {
    case -1: /* this (partly already analyzed) residue is a patch */
    case 0: /* this is without any carboxyl/carbonyl O, will not be
               treated as aminoacid residue */
      si->backbone=BBCTERM;
      markfrom(si,rsd,spec->nchains);
      break;
    case 1: /* normal residue without one O, trying aminoacid residue */
    case 2: /* -COOH perhaps */
      loop (n,0,si->nnbr)
        if (atom[site[si->nbr[n]].type].name[0]!='H'
            && !site[si->nbr[n]].backbone) {
          rsd++;
          si=&site[si->nbr[n]];
          si->backbone=BBCTERM;
          markfrom(si,rsd,spec->nchains);
          /* at least 1 non-H atom found ==> patched cter found ==> numbered as residue */
          goto CterPatch; }
      break;
    default: WARNING(("%d O at Cter",endkey)) }

 /* -COO- and -COOH */
      /* if (endkey==2) prt("-COO found"); */

 CterPatch:
   prt("! backbone %c: %d residues",'`'+spec->nchains,rsd-rsd0);
   if (spec->nchains>26) WARNING(("too many chains (backbones)"))

   goto Chain;
}

#include "blendhlx.c"

void writePDB1(species_t *spec) /******************************* writePDB1 */
/*
  Writes the PDB-like file *.pdb, generating it from scratch
*/
{
 site_t *si;
 int i,ns=spec->ns,j,maxrsd,ii;
 enum bb_e bb;

 nrid=Hnrid+1;

 maxrsd=makersdinfo(spec);

 for (;; spec->frame+=Xopt.byframe) {

   if (Xopt.ropt) {
     if (spec->frame>Xopt.toframe) break;
     if (!read3D(spec,-1)) break; }

   openpdb(spec,Xopt.ropt?spec->frame+1:-1);

   loop (i,0,ns) site[i].count=0;

   loopto (globrsd,1,maxrsd) {

     loop (i,0,ns) if ( globrsd==(si=site+i)->rsd ) {
       ii=1;
       while ( (j=thesameid(i,ns))>=0 ) {
         si->count=1;
         site[j].count=++ii; }
       }

#if 0
     pdbprintf("%d %s",key,rsdname);
     if (globrsd==maxrsd) loop (i,0,ns) {
       si=site+i;
       pdbprintf("%2d %4s rsd=%d bb=%d nest=%d count=%d",
         i,atom[si->type].name,si->rsd,si->backbone,si->nest,si->count); }
#endif /*# 0 */

     /* N CA C O first */
     loopto (bb,BBN,BBCO) {
       int n=0;

       loop (i,0,ns)
         if ( (globrsd==(si=site+i)->rsd) && si->backbone==bb) {
           n++;
           putPDB(i,n>1 && bb==BBCO ? " OXT":NULL); } }

     loop (i,0,ns) if ( (globrsd==(si=site+i)->rsd)
                        && (si->backbone==BBNONE || si->backbone>=BBCB) )
       putPDB(i,NULL);
     }

   loop (i,0,ns) if ( !(si=site+i)->rsd ) putPDB(i,NULL);

   closepdb();
   if (!Xopt.ropt) break; }
}

void writeATM(species_t *spec) /********************************** writeATM */
/*
  writes file in the following format:
  NS

  At   X             Y          Z       mass    charge
  At   X             Y          Z       mass    charge
  ...
*/
{
  site_t *si;
  FILE *f;
  int i;

  strcpy(spec->ext,".atm");
  f=fopen(spec->fn,"wt");

  site=spec->site;

  fprintf(f,"%d\n\n",spec->ns);

  loop (i,0,spec->ns) {
    si=site+i;
#if 0 /* old */
    fprintf(f,"%3s %11.5f %12.5f %12.5f %7.3f %8.5f\n",
	    Mendeleyev(atom[si->type].mass), si->r[0], si->r[1], si->r[2], atom[si->type].mass, si->charge);
#else /*# 0 */
    if (atom[si->type].mass>0.5)
      fprintf(f,"%-3s %11.5f %12.5f %12.5f\n",
	      Mendeleyev(atom[si->type].mass), VARG(si->r));
#endif /*#!0 */
  }
  fclose(f);
}

void writePDB2(species_t *spec) /******************************* writePDB2 */
/*
  Writes the PDB-like file *.pdb, uses id info
*/
{
  site_t *si;
  int i,ns=spec->ns,j=0/*needed?*/,resno,het;
  char id[24],*c,chain;

  site=spec->site;
  resno=site[0].id[3]-'0';
  if (resno<0 || resno>9)
    ERROR(("%s mol-file info not enough to generate pdb: try option -w10",spec->fn))

  for (;; spec->frame+=Xopt.byframe) {

    if (Xopt.ropt) {
      if (spec->frame>Xopt.toframe) break;
      if (!read3D(spec,-1)) break; }

    openpdb(spec,Xopt.ropt?spec->frame+1:-1);

    loop (i,0,ns) {
      //      tofix=0;
      c=(si=site+i)->id;
      if ( !isupper(c[0]) || !isupper(c[1]) || (!isupper(c[2]) && !isdigit(c[2])) )
	/* accepts AAA AA1 ..., but not A11 nor 111 ... */
	ERROR(("%s: %s: 3 letter code expected: try option -w10 to generate PDB",spec->fn,c))
        c[23]=0;
      strcpy(id,c); c=id;

#if 0 /* ??? */
    again:
      /* remove possible by-hand marking ...*/
      strcpy(id,c); c=id;
      if (strlen(id)>2) {
	while (*c<'A' || *c>'Z') c++;
	c=strend(id)-1;
	if (c>id)
	  while ( !((*c>='A' && *c<='Z') || (*c>='0' && *c<='9')) ) *c--=0; }

      if (strlen(id)<5 || id[3]<'0' || id[3]>'9') {
	prt("! bad id=%s ... trying to fix",id);
	sprintf(id,"XXX%d",resno);
	if (!tofix && si->nnbr) {
	  tofix++; c=site[si->nbr[0]].id; goto again; } }
#endif /*# 0 */

      c=id+3;
      if (*c<'0' || id[3]>'9') resno++;
      else {
	j=atoi(c);
	if (j>9 && j>=10*resno) j/=10;
	resno=j; }
      do *c++=0; while (j/=10); // ? undefined j here
      chain='~';
      if (*c=='~' || (*c>='a' && *c<='z')) chain=*c++;
      if (chain=='~') chain=' ';
      chain=toupper(chain);
      if (isdigit(*c)) j=*c++;
      else j=' ';

      het=!si->nnbr || !strcmp(id,"HOH") || !strcmp(id,"WAT");

/*
ATOM      1  N   THR     1      17.047  14.099   3.625  1.00 13.79      1CRN  70
111111222223334444445555556666666666667777777788888888999999000000-------===xxxx
*/

      if (!*c /*|| tofix*/) {
	copy(c,atom[si->type].name,3);
	c[3]=0; }

      pdbprintf("%6s%5d %c%-4s%3s %c%4d%12.3f%8.3f%8.3f%6.2f%6.2f",
		het?"HETATM":"ATOM  ",i+1,
		j,c,id,chain,resno,
		si->r[0],si->r[1],si->r[2],1.00,0.0); }

    closepdb();
    if (!Xopt.ropt) break; }
}

void averagecharges(species_t *spec, int *center) /********** averagecharges */
/*
  charges on marked atoms are averaged
  useful e.g. if quantum chemistry charges are loaded
*/
{
  int i,ns=spec->ns,nq=0;
  double sumq=0;

  loop (i,0,ns)
    if (spec->site[i].count) nq++,sumq+=spec->site[i].charge;
  if (nq>1) {
    sumq/=nq;
    spec->wrmol=1; /* mark rewrite mol-file */
    spec->chargewarned=chargewarned; /* mark verbose check fractional charges */
    prt("! %d charges averaged to %g e",nq,sumq);
    loop (i,0,ns) {
      if (spec->site[i].count) spec->site[i].charge=sumq;
      site[i].count=0; }
    *center=-1; }
}

#if 0
/* removed in V2.4a */
void writeCFGV26(species_t *spec) /***************************** writeCFGV26 */
/* write configuration in the cook format (<V2.7) */
{
  vector L,rmin;
  struct a_s {
    int  size;    /* whole struct in bytes */
    int  padd;    /* padded to 8 bytes */
    double logs;  /* log of the Nose variable s */
    vector rp[1]; /* contiguous array of all r[ns] (and p[ns]) */
  } *a;
  int i,k,cfgkey,norm=spec->opt_y&8;
  static double t,h=0.001,Entot,RvdW,zero;
  double RT=Xopt.T*0.8314472;

/* sim/optlist.c for cook */
#  include "optlist.c"

  strcpy(spec->ext,".cfg");

  if (!(spec->opt_y&8)) WARNING(("no box and output to %s: -y13 suggested",spec->fn))

  VarOpen(spec->fn,"w");
  prt("writing %s",spec->fn);
  cfgkey=8+16; /* version 2.4 */
  VarPut(&cfgkey,sizeof(cfgkey));
  optionlist['m'&31]=1+(RT!=0); /* 1=just positions,2=also velocities */
  VarPut(optionlist,sizeof(optionlist));
  VarPut(&spec->nclust,sizeof(spec->nclust)); /* should be number of molecules */
  put2(cfgkey,spec->nclust)
  VV(L,=spec->L)
  VarPut(L,sizeof(L[0]));
  VarPut(&t,sizeof(t));
  VarPut(&h,sizeof(h));
  VarPut(&Entot,sizeof(Entot));
  put3(t,h,Entot)

  i=12-2;
  /* added in 2.4b */
  VarPut(&RvdW,sizeof(RvdW)); i--;

  while (i--) VarPut(&zero,sizeof(zero));
  VarPut(L+1,sizeof(L[1]));
  VarPut(L+2,sizeof(L[2]));
  putv(L)

  sdsalloc(a,sizeof(*a)+(spec->ns-1)*sizeof(vector));

  loop (k,0,3) rmin[k]=0;

  loop (i,0,spec->ns)
    loop (k,0,3) {
      Min(rmin[k],spec->site[i].r[k])
      a->rp[i][k]=spec->site[i].r[k]; }

  loop (k,0,3) if (rmin[k]<0)
    prt("WARNING: negative %c-coordinate=%g: %s",
                       k+'x',    rmin[k],
                       norm?"cfg will be shifted":"some versions of cook may fail");

  if (norm) loop (i,0,spec->ns) loop (k,0,3) a->rp[i][k]-=rmin[k];

  VarPut(a,a->size);
  if (RT) {
    loop (i,0,spec->ns) loop (k,0,3)
      a->rp[i][k]=h*sqrt(RT/atom[spec->site[i].type].mass)*rndgauss();
    VarPut(a,a->size); }

  KeyClose(0);
  prt("cfg written");
}
#endif /*# 0 */

void writeCFG(species_t *spec) /********************************* writePDB2 */
/* write configuration in the cook format (V>=2.7a)*/
{
  vector rmin;
  struct a_s {
    int  size;    /* whole struct in bytes */
    int  padd;    /* padded to 8 bytes */
    double logs;  /* log of the Nose variable s */
    vector lambda;/* barostat variable */
    vector shape; /* (reserved) */
    vector rp[1]; /* contiguous array of all r[ns] (and p[ns]) */
  } *cfg;
  struct rec_s {
    int key;
    int intval;
    vector vecval;
  } rec;
  int i,k,cfgkey,norm=spec->opt_y&8;
  double RT=Xopt.T*0.8314472,h=0.001;

/* sim/optlist.c for cook */
#include "optlist.c"

  strcpy(spec->ext,".cfg");

  if (!(spec->opt_y&8)) WARNING(("no box and output to %s: -y13 suggested",spec->fn))

  VarOpen(spec->fn,"w");
  prt("writing %s",spec->fn);
  cfgkey=1; /* version 2.7a */
  VarPut(&cfgkey,sizeof(cfgkey));
  optionlist['m'&31]=1+(RT!=0); /* 1=just positions,2=also velocities */
  VarPut(optionlist,sizeof(optionlist));
  put2(cfgkey,spec->nclust)

  memset(&rec,0,sizeof(rec));

  rec.key=1;
  rec.intval=spec->nclust;
  rec.vecval[1]=h;
  VarPut(&rec,sizeof(rec));

  rec.key=2;
  rec.intval=spec->ns;
  VV(rec.vecval,=spec->L)
  VarPut(&rec,sizeof(rec));

  sdsalloczero(cfg,sizeof(*cfg)+(spec->ns-1)*sizeof(vector));

  loop (k,0,3) rmin[k]=0;

  loop (i,0,spec->ns)
    loop (k,0,3) {
      Min(rmin[k],spec->site[i].r[k])
      cfg->rp[i][k]=spec->site[i].r[k]; }

  loop (k,0,3) if (rmin[k]<0)
    prt("WARNING: negative %c-coordinate=%g: %s",
                       k+'x',    rmin[k],
                       norm?"cfg will be shifted":"some versions of cook may fail");

  if (norm) loop (i,0,spec->ns) loop (k,0,3) cfg->rp[i][k]-=rmin[k];

  VarPut(cfg,cfg->size);
  if (RT) {
    loop (i,0,spec->ns) loop (k,0,3)
      cfg->rp[i][k]=h*sqrt(RT/atom[spec->site[i].type].mass)*rndgauss();
    VarPut(cfg,cfg->size); }

  KeyClose(0);
  prt("cfg written");
}

void randomcfg(species_t *spec) /******************************** randomcfg */
{
  int i,j,ns=spec->ns,it;
  vector ri;
  double *rj;
  /* with this L0, there would be on average 1 site per 64 AA^3 */
  double L0=pow((double)ns,1.0/3)*(4*0.707);
  double minrr=2;

  loop (i,0,ns) {

    it=0;
    again:
      if (it++>20) { minrr*=0.99; it=0; }
      VO(ri,=(rnd()-rnd())*L0)
      loop (j,0,i) {
        rj=spec->site[j].r;
        if (Sqr(ri[0]-rj[0])+Sqr(ri[1]-rj[1])+Sqr(ri[2]-rj[2])<minrr) goto again;
        }
    copy(spec->site[i].r,ri,sizeof(vector)); }

  prt("! random initial configuration, min distance > %g",sqrt(minrr));
}

int read2D(species_t *spec,const char *mode) /********************** read2D */
/***
  Reads 2D configuration
  mode must be either "rt" or "rb"
***/
{
  FILE *f;
  long int r[2];
  int i;
  /* default scaling is 1A = 10 pixels */
  double scale=scale2D();

  strcpy(spec->ext,".2d"); strcat(spec->ext,mode+1);
  f=fopen(spec->fn,mode);
  prt("! reading %s",spec->fn);
  if (f==NULL) return 0;

  loop (i,0,spec->ns) {
    if (mode[1]=='t') {
      if (fscanf(f,"%li%li",r,r+1)!=2)
        ERROR(("%s: format or unexpected EOF",spec->fn)) }
    else
      if (fread(r,sizeof(long int),2,f)!=2)
         ERROR(("%s: unexpected EOF",spec->fn))
    spec->site[i].r[0]=scale*(double)r[0];
    spec->site[i].r[1]=scale*(double)r[1];
    spec->site[i].r[2]=0; }
  perturb(spec,NULL,3*(option('e')>=0),0);
  spec->edit=0; /* cannot .che input AND edit mol-file */

  return 1;
}

double masscenter(species_t *spec,int pass) /******************** masscenter */
/*
  returns the mass of the molecule
  if (pass&spec->opt_y) then the molecule is centered (center-of-mass:=0)
  pass=1: was called before minimizing
  pass=2: was called after minimizing
*/
{
  int i,k,ns=spec->ns;
  double m,M=0,sig;
  vector r0,minr,maxr;

  site=spec->site;
  memset(r0,0,sizeof(vector));
  loop (i,0,ns) {
    m=atom[site[i].type].mass;
    M+=m;
    VV(r0,+=m*spec->site[i].r) }

  if (pass&spec->opt_y) {
    /* both before and after minimizing */
    m=1/M;
    loop (i,0,ns) VV(site[i].r,-=m*r0) }

  if (pass==2) {
    if (spec->opt_y&2) VO(spec->L,=0)
    if (spec->opt_y&12) {
      VO(minr,=9e5)
      VO(maxr,=-9e5)
      loop (i,0,ns) {
        if (spec->opt_y&16) sig=0;
        else sig=atom[spec->site[i].type].LJ[0].RvdW*0.89;
        loop (k,0,3) {
          Min(minr[k],spec->site[i].r[k]-sig)
          Max(maxr[k],spec->site[i].r[k]+sig) } } }

    if (spec->opt_y&4) {
      if (spec->opt_y&2) WARNING(("-y2 and -y4 are exclusive (-y4 applies)"))
      loop (i,0,ns) VV(site[i].r,-=minr)
      VV(maxr,-=minr) }

    if (spec->opt_y&8) {
      if (spec->opt_y&2) WARNING(("-y2 and -y8 are exclusive (-y8 applies)"))
      VV(spec->L,=maxr) } }

  return M;
}

char *tempname(void) /********************************************* tempname */
/*
  creates unique file name, returns it again in subsequent calls
  UGLY AND OBSOLETE, needs #include <unistd.h>
*/
{
  static char name[]="#blendtmp.XXXXXX";
  static int pass=0;
  int i;

  if (!pass) {
    i=mkstemp(name);
    if (i<0) { sleep(2); i=mkstemp(name); }
    if (i<0) { sleep(6); i=mkstemp(name); }
    if (i<0) ERROR(("%s: cannot create temporary file name",name))
    if (option('v')&4) fprintf(stderr,"using temporary file %s\n",name);
    pass++; }

  return name;
}

void partialcharges(species_t *spec) /*********************** partialcharges */
{
  int i,m,n,ns;
  int map[MAXVAL];
  struct partial_s *p;
  double q;
  site_t *s;

  if (partial0) {
    site=spec->site;
    ns=spec->ns;
    loop (i,0,ns) site[i].pch=0;

    loop (i,0,ns) {
      sitei=&site[i];
      for (p=partial0; p; p=p->next) {
        if (p->center==sitei->type && p->nnbr<=sitei->nnbr) {
          loop (n,0,sitei->nnbr) map[n]=-1;
          loop (m,0,p->nnbr) {
            loop (n,0,sitei->nnbr) if (map[n]<0)
              if (site[sitei->nbr[n]].type==p->nbr[m]) {
                map[n]=m; goto found; }
            goto notfound;
            found:; }

          /* match found */
          q=p->centerq;
          if (sitei->pch) {
            prt("! site %4d %-9s %-4s charge assignment overlap #%d (%.4f + %.4f)",
                i,sitei->id,atom[sitei->type].name,sitei->pch,
                sitei->charge,q);
            sitei->charge+=q; }
          else
            sitei->charge=q;
          sitei->pch++;

          loop (n,0,sitei->nnbr) if (map[n]>=0) {
            s=&site[sitei->nbr[n]];
            q=p->q[map[n]];
            if (s->pch) {
              prt("! site %4d %-9s %-4s charge assignment overlap #%d (%.4f + %.4f)",
                  i,s->id,atom[s->type].name,s->pch,
                  s->charge,q);
              s->charge+=q; }
            else
              s->charge=q;
            s->pch++; }
          goto nextone; }
        notfound:; }
      nextone:; }

    n=m=0;
    loop (i,0,ns) {
      s=&site[i];
      if (s->pch) n++;
      else if (s->charge==0) m++;
      if (option('v')&4)
        prt("! %-4d %-10s %-4s %7.4f %s",
            i,s->id,atom[s->type].name,s->charge,
            s->pch ? "new" : s->charge==0 ? "": "!"); }
    prt("! %d charges assigned, %d not (of these %d zero, %d nonzero)",
        n,ns-n,m,ns-n-m);
    if (ns-n) prt("! WARNING there are unassigned charges");
    spec->wrmol++; }
}
