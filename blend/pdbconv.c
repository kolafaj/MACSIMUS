#include "ground.h"
#include "options.h"
#include "pdbbasic.h"
#include "pdbconv.h"

void connect(site_t *s1,site_t *s2) /******************************* connect */
/*
  creates a bond between sites s1 and s2
  if already connected and option('v')>0, a message is printed
*/
{
  int n1=s1->n,n2=s2->n,i,err=0,shift=0;

  loop (i,0,s1->nnbr) if (s1->nbr[i]==n2) goto msg; /* already connected */
  loop (i,0,s2->nnbr) if (s2->nbr[i]==n1) goto msg; /* already connected */

  if (s1->nnbr>=option('l') || s2->nnbr>=option('l'))
    err=1;
  else {
    s1->nbr[s1->nnbr++]=n2;
    s2->nbr[s2->nnbr++]=n1;
    /* connected OK */
    shift=2; }

 msg:
  if (option('v')>0 || err) {
    prt_("%4d %s[%s] -%4d %s[%s] ",s1->n,s1->ids->id,s1->type,s2->n,s2->ids->id,s2->type);
    if (err) {
      _n
      ERROR(("cannot connect - too many bonds"))
      return; }
    prt("reconnected"+shift); }
}

void disconnect(site_t *pat,site_t *site)
/*
  if there is a bond from a site from list site-> which points to a
  site from pat->, this bond is removed form site->
*/
{
  site_t *p,*s;
  struct reconnect_s *re;
  int i,j;

  looplist (p,pat)
    looplist (s,site)
      loop (i,0,s->nnbr) if (s->nbr[i]==p->n) {
        if (option('v')>1)
          prt("patch disconnect: %d %s [%s] - %d %s [%s]",
              p->n,p->ids->id,p->type,
              s->n,s->ids->id,s->type);
        allocone(re);
        re->site[0]=p;
        re->site[1]=s;
        re->next=reconnecthead;
        reconnecthead=re;
        loop (j,i+1,option('l')) s->nbr[j-1]=s->nbr[j];
        s->nnbr--; }
}

void renum(site_t *site)
/*
  renumber list site-> (incl. all bonds) from zero
*/
{
  site_t *s;
  int n,m=0,i;
  int *ren;

  looplist (s,site) Max(m,s->n)
  m++;

  allocarray(ren,m);
  loop (n,0,m) ren[n]=-1;
  n=0; looplist (s,site) ren[s->n]=n++;

  looplist (s,site) {
    s->n=ren[s->n];
    loop (i,0,s->nnbr) {
      s->nbr[i]=ren[s->nbr[i]];
      if (s->nbr[i]<0)
        ERROR(("pdb internal error: patch has bond outside patch or patched residue"))
      } }
  free(ren);
}

#if 0
/* WARNING: NEVER USED, NEVER DEBUGGED */
void mergeres(residue_t *r1,residue_t *r2) /*********************** mergeres */
/* merge residues, incl. renumbering sites, but not proper naming of rsd */
{
  int i,ns,*ren,n;
  site_t *s,**lasts;
  atom_t *a,**lasta;

  /* renumbering */
  ns=0;
  looplist (s,r2->site) Max(ns,s->n)
  ns++;
  allocarray(ren,ns);
  i=0;
  looplist (s,r2->site) ren[s->n]=i++;
  looplist (s,r2->site) {
    s->n=ren[s->n];
    loop (n,0,s->nnbr) s->nbr[n]=ren[s->nbr[n]]; }
  free(ren);

  /* merging sites */
  lasts=&r1->site;
  looplist (s,r1->site) lasts=&s->next;
  *lasts=r2->site;

  /* merging atoms */
  lasta=&r1->atom;
  looplist (a,r1->atom) lasta=&a->next;
  *lasta=r2->atom;

  r1->next=r2->next;
}

int nsites(site_t *site) /****************************************** nsites */
/* # of sites in the list is returned */
{
  int n;

  while (site) n++,site=site->next;

  return n;
}
#endif

site_t *patch(site_t *pat,site_t *site,int merge) /****************** patch */
/*
  applies patch starting with site `pat' to residue starting with site `site'
  merge=1: a merged list site-> is generated (to be held by one residue)
           NULL is returned
  merge=0: a patched list site-> is generated, but additional sites
           not used in patching (renumbered) are still kept in `pat'
           and should be held in a separate residue; in this case,
           list of residue-residue bonds to connect later is created
           new pat is returned
  turns off all unique flags - strict matching required (this is
    because I'm too lazy to recalculate these flags correctly)
  IMPLEMENTATION LIMITATION: pat must be in 1 contiguous array
*/
{
  site_t *p,*s,*tail=NULL;
  int i,j,k,n=0,pns;
  struct id_s *idl,*idp;

  if (option('v')>1) prt("patching");

  /* find max site index in site-> and the last site (the tail) */
  looplist (s,site) {
    s->unique=0;
    Max(n,s->n)
    tail=s; }

  if (!tail) {
    ERROR(("empty residue to patch"))
    return 0; }

  /* prepare renumbering table */
  for (p=pat,i=0; p; p=p->next,i++) {
    p->unique=0;
    if (p->n!=i) ERROR(("pdb internal error: patch numbering"))
    if (p->next) if (p->next-p!=1) ERROR(("pdb internal error: non-contig patch"))
    if (p->patch[0]) {
      /* to be replaced */
      looplist (s,site) looplist (idl,s->ids)
        if (!strcmp(p->patch,idl->id)) {
          p->n=s->n; /* number from site */
          goto found; }

      if (option('v')) {
        WARNING(("patch: cannot find %s which should replace %s",p->ids->id,p->patch))
        prt("atom has been added"); }
      p->n=++n;
     found:; }
    else
      p->n=++n; /* new atom */ }
  pns=i;

  /* renumber sites */
  for (p=pat,i=0; p; p=p->next,i++)
    loop (j,0,p->nnbr) {
      if (p->nbr[j]<0 || p->nbr[j]>=pns)
        ERROR(("wrong bonds in patch starting at %s (renum %s nbr[%d])",
               pat->ids->id,p->ids->id,j))
      p->nbr[j]=pat[p->nbr[j]].n; }

  /* replace atoms */
  looplist (p,pat) if (p->patch[0])
   looplist (s,site) looplist (idl,s->ids)
     if (!strcmp(p->patch,idl->id)) {
          
       /* replace atom (list) by atom (list) */

       freeid(s);
       looplist (idp,p->ids) appendid(s,idp->id);
/*.....      strcpy(s->id,p->id);*/
       strcpy(s->patch,p->patch);
       strcpy(s->type,p->type);
       s->charge=p->charge;
       s->chir=p->chir;

       /* add bonds */
       loop (j,0,p->nnbr) {
         loop (k,0,s->nnbr)
           if (s->nbr[k]==p->nbr[j]) goto alreadybonded;
         if (s->nnbr>=option('l')) {
           ERROR(("patch: too many bonds - check option -l!"))
           s->nnbr=option('l'); }
         else
           s->nbr[s->nnbr++]=p->nbr[j];
       alreadybonded:; }
       break; /* looplist (idl,s->ids) */
     }

  if (merge) {
    /*
      one merged site list: move atoms from pat to the end of site
      (sites are correctly numbered and all bonds have been created)
    */
    looplist (p,pat) if (!p->patch[0]) {
      n++;
      tail->next=p;
      tail=p; }
    tail->next=NULL;
    return NULL; }

  else {
     /* separate site lists:
          site=patched list
          pat=additional (not replacing) atoms
     */
    site_t **last=&pat;

    /* remove the patched sites leaving the extra sites */
    looplist (p,pat)
      if (p->patch[0]) *last=p->next;
      else n++,last=&p->next;

    /* remove bonds between pat-> and site-> (standard nbr assumed) */
    disconnect(pat,site);
    disconnect(site,pat);

    /* and guarantee standard numbering of the sites left in pat */
    renum(pat);
    return pat; }
}

/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%% PDB converter support %%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

int renumber(site_t *site,int offset) /***************************** renumber */
/*
  renumber sites by adding `offset'
  site[0] gets number `offset', etc.
  number of sites is returned
*/
{
  int j,ns=0;

  while (site) {
    ns++;
    site->n+=offset;
    loop (j,0,site->nnbr) site->nbr[j]+=offset;
    site=site->next; }

  return ns;
}

site_t *findid(residue_t *res,char *id1,char *id2) /***************** findid */
/*
  finds 1st site with given id1 (id2 if not found) starting from site res->from
*/
{
  site_t *from;
  struct id_s *idl;
  char line[128];

  sprintf(line,"can't make peptide bond because can't find atom %s [%s] in %s %4s %c %d %c",
          id1,id2?id2:"",
          res->moltype==AMINOACID?"residue":"patch",
          res->resnm,res->chain,res->resno,res->resins);

 again:
  looplist (from,res->site) looplist (idl,from->ids)
    if (!strcmp(idl->id,id1)) return from;
  if (id2) {
    if (option('v')) prt("trying to connect %s because %s not found",id1,id2);
    id1=id2; id2=NULL;
    goto again; }

  if (res->moltype==AMINOACID) ERROR((line))
  else if (option('v')) prt(line);

  return NULL;
}


int delete(site_t *del) /******************************************** delete */
/* deletes site *del and all its references */
{
  site_t *s,*head;
  residue_t *res;
  int i,j,ns=0;

  looplist (res,reshead)
    for (s=res->site,head=NULL; s; s=s->next) {
      if (ns!=s->n) {
        ERROR(("numbering (pdb internal error): site->n=%d ns=%d",s->n,ns))
        ns=s->n; }
      ns++;
      if (s==del) {
        /* remove from list */
        if (head) head->next = del->next;
        else res->site = del->next;
        if (option('v')) {
          prtres(res);
          prt(": atom %s deleted",s->ids->id); }
        /* cannot free(del) ! */
        goto deleted; }
    head=s; }

  ERROR(("cannot find site %s (temporary n=%d) to delete",del->ids->id,del->n))
  return 0;

 deleted:

  /* renumber */
  looplist (res,reshead)
    looplist (s,res->site) {
     if (s->n>del->n) s->n--;
      for (i=0; i<s->nnbr; ) {
        if (s->nbr[i]>del->n) s->nbr[i++]--;
        else if (s->nbr[i]==del->n) {
          s->nnbr--;
          loop (j,i,s->nnbr) s->nbr[j]=s->nbr[j+1]; }
        else i++; } }

  return 1;
}

#if 0
/* old version: locations must be consecutive */
int altlocations(void) /************************************** altlocations */
{
  residue_t *res;
  int altlocations=0;

  looplist (res,reshead) {
    atom_t *a;

    for (a=res->atom; a; ) {
      if (a->altloc>='A') {
        double sumoccup=0,maxoccup=0;
        char *idA=a->id;
        char altch=a->altloc;
        vector r;
        atom_t *aA=a;

        VO(r,=0)

        while (a && !strcmp(idA,a->id)) {
          if (altch!=a->altloc)
            ERROR(("residue %4s %c %d %c, atom %s:\n\
*** alternative locations not in order",
                   res->resnm,res->chain,res->resno,res->resins,a->id))
          sumoccup+=a->occup;

          switch (option('a')&3) {
            case 0: if (altch=='A') VV(r,=a->r) break;
            case 1: VV(r,+=a->occup*a->r) break;
            case 2: if (a->occup>maxoccup) {
              maxoccup=a->occup;
              VV(r,=a->r) break; } }

          altch=a->altloc+1; /* normally = altch++ */
          a=a->next; }

        altlocations++;
        if (fabs(sumoccup-1)>1e-4)
          if (option('a')&1) {
            if (option('v')) {
              WARNING(("residue %4s %c %d %c, atom %s:\n\
     sum of occupancies=%.4f (should be 1)",
                       res->resnm,res->chain,res->resno,res->resins,a->id,
                       sumoccup))
              prt("rescaled"); } }
          else {
            static int warn;
            
            if (!warn && option('v'))
              prt("WARNING: sum of occupancies not 1. Further messages suppressed.");
            warn=1; }
        if (option('a')&1) VO(r,/=sumoccup)
        VV(aA->r,=r)
        aA->next=a; }
      else
        a=a->next;
    }
  }
  
  return altlocations;
}
#else
/* new: locations may be ordered arbitrarily within a residue */
int altlocations(void) /************************************** altlocations */
{
  residue_t *res;
  int altlocations=0;

  looplist (res,reshead) {
    atom_t *a;

    looplist (a,res->atom)
      if (a->altloc>='A') {
        /* alternate locations */
        double sumoccup=0,maxoccup=0;
        vector r;
        atom_t *aa,*last;
        char check[32],expected[32];
        int noccup=0;
        int Afound=0;

        VO(r,=0)

        looplist (aa,a) if (!strcmp(aa->id,a->id)) {
          sumoccup+=aa->occup;
          if (noccup>30)
            ERROR(("residue %4s %c %d %c, atom %s:\n\
*** too many alternate locations",
               res->resnm,res->chain,res->resno,res->resins,aa->id))
          check[noccup]=aa->altloc;
          expected[noccup]='A'+noccup;
          noccup++;

          switch (option('a')&3) {
            case 0:
              if (aa->altloc=='A') {
                Afound++;
                VV(r,=aa->r) }
              break;
            case 1:
              VV(r,+=aa->occup*aa->r)
              break;
            case 2:
              if (aa->occup>maxoccup) {
                maxoccup=aa->occup;
                VV(r,=aa->r) }
              break;
            default:
              ERROR(("bad option -a3")) } }

          if ((option('a')&3) == 0 && Afound!=1)
            ERROR(("residue %4s %c %d %c, atom %s:\n\
*** %d locations A (1 expected)",
                   res->resnm,res->chain,res->resno,res->resins,aa->id,Afound))

          check[noccup]=expected[noccup]=0;
          if (strcmp(check,expected))
            WARNING(("residue %4s %c %d %c, atom %s:\n\
*** alternate locations not in order, missing or other error\n\
*** order found: %s, order expected: %s",
             res->resnm,res->chain,res->resno,res->resins,a->id,check,expected))

        altlocations++;

        if (fabs(sumoccup-1)>1e-4) {
          if (option('a')&1) {
            if (option('v')) {
              WARNING(("residue %4s %c %d %c, atom %s:\n\
       sum of occupancies=%.4f (should be 1)",
                res->resnm,res->chain,res->resno,res->resins,a->id,sumoccup))
              prt("rescaled"); } }
          else {
            static int warn;
            if (!warn && option('v'))
              prt("WARNING: sum of occupancies not 1. Further messages suppressed.");
            warn=1; } }
        if (option('a')&1) VO(r,/=sumoccup)
        VV(a->r,=r)
        a->altloc=' '; /* marked as done */

        /* delete unused locations */
        last=a;
        for (aa=a->next; aa; aa=last->next) {
          if (!strcmp(aa->id,a->id)) {
            last->next=aa->next;
            /* cannot free(aa) because allocated as a member of a longer buffer */
            }
          else
            last=aa; } } }

  return altlocations;
}
#endif

void findCYSCYS(void) /****************************************** findCYSCYS */
/*
  find CYS-CYS bonds by distances
*/
{
  residue_t *r1,*r2;
  atom_t *a1,*a2;
  ssbond_t *ss;
  float dist;
  vector dr;
  float SSlimit=3;

  if (option('c')>1) SSlimit=option('c');

  if (option('v'))
    prt("CYS-CYS bridges from coordinates, SG-SG distance limit=%.1f A",SSlimit);

  /* all pairs: */
  looplist (r1,reshead) if (!strcmp(r1->resnm,"CYS")) {
    looplist (a1,r1->atom) if (!strcmp(a1->id,"SG")) {
      for (r2=reshead; r2!=r1; r2=r2->next) if (!strcmp(r2->resnm,"CYS")) {
        looplist (a2,r2->atom) if (!strcmp(a2->id,"SG")) {
          VVV(dr,=a1->r,-a2->r)
          dist=sqrt(SQR(dr));
          if (dist<SSlimit+1) {
            if (option('v')) prt_("CYS %c %d - CYS %c %d = %.3f A : ",
              r1->chain,r1->resno, r2->chain,r2->resno, dist);
            if (dist>SSlimit) {
              if (option('v')) prts("WARNING just above SG-SG limit"); }
            else {
              if (option('v')) prts("to be connected");
              alloc(ss,sizeof(ssbond_t));
              ss->chain[0]=r1->chain; ss->chain[1]=r2->chain;
              ss->resno[0]=r1->resno; ss->resno[1]=r2->resno;
              includeSS(ss); } }
          } /* a2=SG: checked in a1 that SG always exists */
        } /* r2=CYS */
      goto SGfound;
      } /* a1=SG */
    ERROR(("residue %4s %c %d %c: atom SG not found",
           r1->resnm,r1->chain,r1->resno,r1->resins))
    SGfound:; }
}

int sitebyoption(char *id,residue_t *res,int opt) /*********** sitebyoption */
/*
  returns 1 if site is a candidate to ignoring or deleting acording to
    options -i and -d
*/
{
  if (opt>0 && res->moltype!=MOLECULE) return 0;

  switch (abs(opt)) {
    case 0:  return 0;
    case 1:  return removedigits(id)[0]=='H'; /* Hygrogens */
    case 2:  return removedigits(id)[0]!='H'; /* heavy atoms */
    case 3:  return 1;                        /* all atoms */
    default: ERROR(("bad option -i or -d")); }

  return 0;
}

int freeH(site_t *s) /************************************************ freeH */
/*
  returns 1 if s is free H
*/
{
  struct id_s *idl;
  int ret=-1,r;

  looplist (idl,s->ids) {
    r=s->nnbr==0 && removedigits(idl->id)[0]=='H';

    if (ret>=0) if (r!=ret)
      ERROR(("cannot determine `free H' condition for incompatible aliases %s,%s",
             s->ids->id,idl->id)) }

  return ret;
}

int removesites(remove_f *condition,char *msg) /*************** removesites */
/*
  removes all sites s for which condition condition(s) is true
*/
{
  residue_t *res;
  site_t *s;
  int no=0;

  if (option('v')) prt("removing sites because %s",msg);

  looplist (res,reshead)
    looplist (s,res->site)
      if (condition(s)) no+=delete(s);

  if (option('v')) prt("%d sites removed",no);

  return no;
}
