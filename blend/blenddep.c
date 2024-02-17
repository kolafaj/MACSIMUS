#define REMOVENESTED /* does not allow dependants dependent on another
                        dependants, new in 1.9c */

#include "ground.h"
#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "blendmed.h"
#include "blendmin.h"
#include "blendgen.h"
#include "blenddep.h"
#include "blendedt.h"
#include "options.h"
#include "rndgen.h"

int ndependants;

struct dependant_s *dep0;

int removeb,removea,removei;
bond_t *b00; /* b00 -> (removed bonds) -> ... -> b0 real bonds -> NULL */

/* real distance of sites site[i],site[j] */
#define DIST(I,J) sqrt(SQRD(site[I].r,site[J].r))

int readdependants(species_t *spec) /******************** readdependants */
/* missing: remove bonded terms with dependants */
{
  int ii=0,i,j;
  int *k,*l,*m;
  FILE *f;
  struct dependant_s *d;
  bond_t *bdep,*b,*blast=NULL;

  site=spec->site;

  strcpy(spec->ext,".dep");

  if ( (f=fopen(spec->fn,"rt"))) {
    char *line,*tok,*sep;
    int nbr[MAXVAL],*n,nn,nnbr,ast;

    prt("! reading %s",spec->fn);

    alloc(line,option('l'));
    alloc(n,option('l')/2); /* pessimistic */

    while (fgets(line,option('l'),f)) {
      if (option('v')&4) prt("! %s",line);
      if (strlen(line)>1 && !strchr("!#",line[0])) {
        if ( (sep=strchr(line,'!')) ) *sep=0;
        if ( !(sep=strchr(line,':')) )
          ERROR(("%s: missing `:'",spec->fn))
        else {
          *sep++=0;

          /* after `:' */
          nnbr=0;
          tok=strtok(sep," \t\n,");
          while (tok) {
            if (nnbr>=4) ERROR(("%s: too many parents",spec->fn))
            nbr[nnbr++] = findsite(spec,tok,1);
	    tok=strtok(NULL," \t\n,"); }

	  if (nnbr<2) {
	    ERROR(("%s: not enough parents",spec->fn))
            break; }

	  /* before `:' */
	  nn=0;
	  tok=strtok(line," \t\n,");
	  if ( !(ast=!strcmp(tok,"*")) ) while (tok) {
	    n[nn++] = nbr[nnbr++] = findsite(spec,tok,1);
	    tok=strtok(NULL," \t\n,"); }
	  else
	    nn=spec->ns-nnbr;

	  alloc(d,sizeof(struct dependant_s)-sizeof(int)+nn*sizeof(int));
	  d->next=dep0; dep0=d;
	  d->nnbr=nnbr;
	  copy(d->nbr,nbr,sizeof(d->nbr));
	  d->ndep=0;
	  if (ast) loop (i,0,spec->ns) {
	    loop (j,0,nnbr) if (i==nbr[j]) goto found;
	    if (d->ndep>=nn) {
	      ERROR(("%s: wrong parents",spec->fn))
              break; }
	    d->dep[d->ndep++]=i;
	  found:; }
	  else {
	    d->ndep=nn;
	    copy(d->dep,n,nn*sizeof(int)); }
	  if (d->ndep!=nn) ERROR((""))
	  ii+=nn; }
        } }

    fclose(f);
    free(n);
    free(line);

    anymass=1;
    bdep=NULL;
    for (d=dep0; d; d=d->next) {

      /* creating new constraints between parents,
	 to be merged with the list of bonds later */
      loop (i,0,d->nnbr) loop (j,0,i) {
	alloc(b,sizeof(bond_t));
	b->next=bdep; bdep=b;
	b->indx[0]=d->nbr[i];
	b->indx[1]=d->nbr[j];
	b->parm.K=b->parm.Ki2=0;
	b->parm.length=DIST(b->indx[0],b->indx[1]); }

      loop (j,0,d->ndep) {
	/* remove bonded terms containing dependants
	   the algorithm thus differs from the massless case 
	   (which is tailored to tip4p) */
	i=d->dep[j];
	loopnbr (k,i) {
	  removebond(*k,i);
	  loopnbr (l,*k) {
	    removeangle(*l,*k,i);
	    loopnbr (m,*l) {
	      removetorsion(&d0,*m,*l,*k,i);
	      removetorsion(&ar0,*m,*l,*k,i); }
	    loopnbr (m,*k)
	      removetorsion(&i0,*m,*l,i,*k); }
	  loopnbr (l,i) {
	    removeangle(*l,i,*k);
	    loopnbr (m,*l) {
	      removetorsion(&d0,*m,*l,i,*k);
	      removetorsion(&ar0,*m,*l,i,*k); }
	    loopnbr (m,i)
	      removetorsion(&i0,*m,*l,*k,i); }
          } } }

    /* the new constraint list bdep is appended to b00 */
    for (b=b00; b; b=b->next) blast=b;
    if (blast) blast->next=bdep; else b00=bdep;
    if (!b0) b0=bdep; /* ? WARNING: bugs possible here */
    }

  *spec->ext=0;

  return ii;
}

int anymass; /* if set then no massless test performed in functions
                removebond etc. and bond etc. removed unconditionally */

void removebond(int i,int j) /********************************** removebond */
/*
  only if at least one atom is massless
  i,j are indices of atoms in molecule
Note:
  removed bonds are moved to a list starting from b00 and continuing via b0
  (but the last: then b0:=NULL and the last is not moved.  Otherwise an
  infinite loop results.)
  see also removeangle
*/
{
  bond_t *b,**bptr;

  for (b=b0,bptr=&b0; b; b=b->next)
    if ( ((b->indx[0]==i && b->indx[1]==j) || (b->indx[1]==i && b->indx[0]==j))
	 && (anymass || atom[site[i].type].mass==0 || atom[site[j].type].mass==0) ) {
      /* removing bond b: must keep separately from b00 */
      if (option('v')&4) {
	prts_("! dependants:"); prtatom(i); prtatom(j);
	prts(" bond removed"); }
      removeb++;
      if (bptr==&b0)
	b0=NULL;
      else {
	*bptr=b->next;
	b->next=b00; b00=b; } }
    else
      bptr=&b->next;
}

void removeangle(int i,int j,int k) /************************* removeangle */
/*
  note: a general method for removing items from a list is used
  could remove any number of items given by the `if' condition
*/
{
  angle_t *a,**aptr;

  for (a=a0,aptr=&a0; a; a=a->next)
    if (a->indx[1]==j
	&& ( (a->indx[0]==i && a->indx[2]==k)
	     || (a->indx[2]==i && a->indx[0]==k) )
	&& (anymass || atom[site[i].type].mass==0 || atom[site[k].type].mass==0) ) {
      /* removing angle a */
      if (option('v')&4) {
	prts_("! dependants:"); prtatom(i); prtatom(j); prtatom(k);
	prts(" angle removed"); }
      removea++;
      *aptr=a->next; }
    else
      aptr=&a->next;
}

void removetorsion(torsion_t **t0,int i,int j,int k,int l) /*** removetorsion */
{
  torsion_t *t,**tptr;

  for (t=*t0,tptr=t0; t; t=t->next)
    if ( ((t->indx[0]==i && t->indx[3]==l) || (t->indx[3]==i && t->indx[0]==l))
	 && ((t->indx[1]==j && t->indx[2]==k) || (t->indx[2]==j && t->indx[1]==k))
	 && (anymass || atom[site[i].type].mass==0 || atom[site[l].type].mass==0) ) {
      /* removing torsion t */
      if (option('v')&4) {
	prts_("! dependants:"); prtatom(i); prtatom(j); prtatom(k); prtatom(l);
	prts(" torsion removed"); }
      removei++;
      *tptr=t->next; }
    else
      tptr=&t->next;
}

int countdependants(void) /******************************** countdependants */
{
  struct dependant_s *d;
  int i,n=0;

  for (d=dep0; d; d=d->next)
    loop (i,0,d->ndep) n += site[d->dep[i]].clust==clust;

  return n;
}

void prtdependants(void) /*********************************** prtdependants */
{
  struct dependant_s *d;
  double w[MAXVAL],ww[MAXVAL];
  vector r[MAXVAL],r0;
  double f,ff,rad;
  int i,nnbr,idep,iter,niter;

  _n
  prts("dependants");
  prts("!i atom   #   i atom  weight...");

  for (d=dep0; d; d=d->next) 
    loop (idep,0,d->ndep) 
      if (site[d->dep[idep]].clust==clust) { /* cluster support added 6/00 */
	nnbr=d->nnbr;

#ifdef REMOVENESTED
	loop (i,0,nnbr)
	  if (atom[site[d->nbr[i]].type].mass==0) {
	    int j;
	    prt("! %d: nested dependant %d removed",d->dep[idep],d->nbr[i]);
	    d->nnbr=--nnbr;
	    loop (j,i,nnbr) d->nbr[j]=d->nbr[j+1]; }
#endif

	VV(r0,=site[d->dep[idep]].r)
        loop (i,0,nnbr) {
	  VV(r[i],=site[d->nbr[i]].r)
          ww[i]=1./nnbr; }
	rad=0.03; ff=9e9;

/* tip4p
ww[2]=0.1280120648689;
ww[1]=0.1280120648689;
ww[0]=1-2*ww[1];rad=1e-20;
*/

#define NITER 10000000
        iter=0; niter=NITER;
        /* twice as many iterations as needed to reach precision 1e-8 */
	while (iter++<niter) {
	  double sum=0;
	  vector c;

	  loop (i,0,nnbr) sum+=w[i]=ww[i]+(rnd()-rnd())*rad;
	  loop (i,0,nnbr) w[i]/=sum;
	  VO(c,=0)
          loop (i,0,nnbr) VV(c,+=w[i]*r[i])
	  f=SQRD(c,r0);
	  if (f<ff) {
	    loop (i,0,nnbr) ww[i]=w[i];
	    ff=f;
	    rad*=1+0.01*nnbr; }
	  else {
	    rad*=0.995; 
	    if (niter==NITER && rad<1e-8) niter=iter,iter=0; } }
#undef NITER

        if (rad>1e-8) ERROR(("cannot calculate dependants (rad=%g)",rad))

	prtatom(d->dep[idep]);
	prt_(" %d ",nnbr);
	loop (i,0,nnbr) {
	  prtatom(d->nbr[i]);
	  prt_(depend_f,ww[i]); }
	prt("  err=%.2g",sqrt(ff)); }
  _n
}
