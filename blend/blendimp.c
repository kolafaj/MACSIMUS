#include "ground.h"
#include "options.h"
#include "rndgen.h"
#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "blendmed.h"
#include "blendgen.h"
#include "blendimp.h"

static int knownneighbors(int i) /************************** knownneighbors */
{
  int nn=0,*nbr;

  loopnbr (nbr,i) nn+=(site[*nbr].keep!=WANTED);

  return nn;
}

void imprison(species_t *spec) /********************************** imprison */
/***
If there are atoms with unknown positions (WANTED) their
approximate positions are calculated.
If option -k == 0: all atoms are made FREE
If option -k == 1: WANTED atoms are always made FREE,
                   atoms marked by * in .mol file are kept INJAIL
                   atoms not marked by * in .mol file are FREE
If option -k == 2: all atoms but WANTED are kept INJAIL
If option -k == 3: if there is any WANTED found, 2 applies,
                   otherwise 1
***/
{
  int i,j,k,l,*nbr,ns=spec->ns,nonHfilled=0,totnw=0;
  vector dr[MAXVAL];
  enum keep_e nbrstat[MAXVAL];
  int nnbr,nknown,firstwanted;
  int stat[3];
  double rr,dif,COS;
  double L0=pow((double)ns,1.0/3)*1.4;
  int fillnest;
  static int first=1;
  int opt_k=spec->opt_k;
  
  static enum keep_e convert[3][3]={
  /*  FREE   INJAIL WANTED */
    { FREE,  FREE,  FREE }, /* option -k == 0 */
    { FREE,  INJAIL,FREE }, /* option -k == 1 */
    { INJAIL,INJAIL,FREE }  /* option -k == 2 */
    };
  
  if (spec->opt_k<-3 || spec->opt_k>3) {
    spec->opt_k=(spec->opt_k>0?1:-1)*(abs(spec->opt_k)&3);
    prt("! blendimp: option -k reset to -k%d",spec->opt_k); }
  
  site=spec->site; /* site is global */
  prtatomoffset=0;
  
  /* "infinite" # of passes to fill long backbones... */
  for (fillnest=20; fillnest>=0; fillnest--)
  
  /* calculate unknown positions */
  loop (i,0,ns) if (site[i].keep==WANTED) {
    totnw++;
    nonHfilled += toupper(atom[site[i].type].name[0])!='H';
    if (option('v')&4) { prtc('!'); prtatom(i); }
  
    switch (knownneighbors(i)) {
  
      case 0:
        if (fillnest) {
  	if (option('v')&4) prts("no known neighbor - postponed");
  	continue; }
        /* WANTED atom (i) is not bonded to anything known => random position */
      randomly:
        if (option('v')&4) prts_(" randomly");
        VO(site[i].r,=rndgauss()*L0)
      prtfilled:
        site[i].keep=FILLED;
        if (option('v')&4) prt(" filled %d-val",site[i].nnbr);
        break;
  
      case 1:
        /* WANTED atom (i) has just 1 known neighbor */
        j=-1;
        loopnbr (nbr,i)
  	if (site[*nbr].keep != WANTED) {
  	  if (j>=0) ERROR(("%d %d internal",i,j))
            j=*nbr; }
        if (j<0) ERROR(("%d internal",i))
  
  /*.....      j=site[i].nbr[0];*/
  
        /*
                 /
          (i)--(j)--     [ atom (i) is WANTED, (j) is known ]
                 \
        */
  
        firstwanted=-1;
        nknown=nnbr=0;
  
        loopnbr (nbr,j) {
          if ( (nbrstat[nnbr]=site[*nbr].keep) != WANTED) {
            VVV(dr[nnbr],=site[*nbr].r,-site[j].r)
            nknown++; }
          else 
            if (firstwanted<0) firstwanted=*nbr;
          nnbr++; }
  
        /*
          atom (j) has  nknown  known neighbors 
          and  nnbr-nknown  unknown (WANTED), incl. (i)
          all WANTED neighbors of (j) will be filled at once 
        */
        
        if (i!=firstwanted) {
          if (option('v')&4) prts(" elsewhere");
          /* atoms will be/have been filled in another pass - not now */
          break; }
  
        if (option('v')&4)
          loopnbr (nbr,j) if (site[*nbr].keep==WANTED && *nbr!=i) prtatom(*nbr);
  
        if (nnbr<=1) {
          if (fillnest) continue;
          else goto randomly; }
                
        COS=-1./(nnbr-1); /* desired value of cos(bond angle) */
#define ANGLEPREC 0.03
        dif=ANGLEPREC*nnbr; /* initial precision of cos(bond angle) */
  
        /* normalize known vectors */
        loop (k,0,nnbr) 
          if (nbrstat[k]!=WANTED) {
            rr=sqrt(SQR(dr[k]));
            VO(dr[k],/=rr) }
  
      again:
        dif*=1.0000025;
  
#if 0 /* old simple slow */
  
        loop (k,0,nnbr) {
          if (nbrstat[k]==WANTED) rndsphere(dr[k]);
          loop (l,0,k) {
            rr=SCAL(dr[k],dr[l]);
            if (fabs(rr-COS)>dif) goto again; } }
  
#else /* faster: last vector created by `reflection' */
        {
          int inw=nnbr-nknown;
          vector reflect;
  
          loop (k,0,nnbr) {
            if (nbrstat[k]==WANTED) {
  
              if (inw==1) {
                int sitel;
                /* last atom to fill */
  
                VO(reflect,=0)
                loop (l,0,nnbr) if (l!=k) {
                  sitel=site[j].nbr[l];
                  VV(reflect,-=dr[l]) }
                rr=1./sqrt(SQR(reflect)); 
                VV(dr[k],=rr*reflect);
                if (nnbr==2) {
                  int *ll;
                  vector rdih;
                  int nrdih=0;
  
                  /* patch to solve unknown OH in -COOH or so:
                     H is beeing filled so that it is trans to the backbone */
                  VO(rdih,=0)
                  loopnbr (ll,sitel) if (*ll!=j)
                    if (site[*ll].nnbr>1) {
                      nrdih++;
  		    /* positin of *ll rather than *lll ! */
                      VVV(rdih,+=site[*ll].r,-site[sitel].r) }
  
                  if (nrdih) {
                    double s=SCAL(dr[k],rdih)/SQR(dr[k]);
                    VV(rdih,-=s*dr[k]) /* perpendicular to dr */
                    VV(dr[k],-=(0.5/nrdih)*rdih) }
  
                  /* to avoid singular angle 180 for (i)-(j)- : 
                     add small random almost perpendicular vector */
                  do rndsphere(rdih); while (fabs(SCAL(rdih,dr[k]))>0.2);
                  VV(dr[k],+=0.1*ANGLEPREC*rdih) }
                break; }
  
              rndsphere(dr[k]);
              inw--; }
  
            loop (l,0,k) {
              rr=SCAL(dr[k],dr[l]);
              if (fabs(rr-COS)>dif) goto again; } } 
        }
#endif
  
        if (option('v')&4) {
          prt_("(");
          prtatom(j);
          prt_(") dif=%.3f cos=",dif);
          loop (k,0,nnbr)
            loop (l,0,k) prt_(" %.3f",SCAL(dr[k],dr[l]));
        }
  
        if (dif>0.3 && first) {
          WARNING(("dif=%g: neighbor(s) of %d-bonded atom %d\n\
*** might be calculated inaccurately\n\
*** (more warnings suppressed)",dif,nnbr,j))
          first=0; }
        loop (k,0,nnbr)
          if (nbrstat[k]==WANTED)
            VVV(site[site[j].nbr[k]].r,=site[j].r,+dr[k])
  
        goto prtfilled;
  
      default:
        /*
          WANTED atom (i) is bonded by more than 1 bond:
          center of bonded atoms used if possible
          NEW: trying neighbors of neighbors to bend (to avoid zero angles)
        */
        {
  	int nnext;
  
  	VVO(dr[1],=dr[0],=0)
          nnext=nknown=0;
  	loopnbr (nbr,i) if (site[*nbr].keep!=WANTED) {
  	  int *nn;
  
  	  loopnbr (nn,*nbr) if (site[*nn].keep!=WANTED) {
  	    VV(dr[1],+=site[*nn].r)
              nnext++; }
            VV(dr[0],+=site[*nbr].r)
            nknown++; }
  	if (nknown<2) ERROR(("internal"))
  	if (nnext) {
  	  double x=-0.3/nnext;
  			  
  	  rr=1.3/nknown;
  	  VVV(site[i].r,=rr*dr[0],+x*dr[1])
  	    rr=1e-3; }
  	else {
  	  rr=1.0/nknown;
  	  VV(site[i].r,=rr*dr[0])
  	    rr=1e-6; }
  	VO(site[i].r,+=rndgauss()*rr) 
          goto prtfilled; }
        
      } /* switch (site[i].nnbr) */
    } /* if (site[i].keep=WANTED) */
  
  loop (i,0,ns) if (site[i].keep==FILLED) site[i].keep=WANTED;
  
  if (totnw) prt("! %d atom(s) filled",totnw);
  else if (abs(spec->opt_k)==2) {
    WARNING(("%s: -k2 and no atom to fill",spec->fn))
    spec->opt_k=1; }
  if (nonHfilled) prt("!*****WARNING %d non-H atom(s) filled",nonHfilled);
  
  if (spec->opt_k==3) spec->opt_k-=1+!totnw;
  if (spec->opt_k==-3) spec->opt_k+=1+!totnw;
  
  /* change keep status according to option -k */
  loop (i,0,ns) site[i].keep=convert[abs(spec->opt_k)][site[i].keep];
  
  /* negative spec->opt_k: make single atoms always free */
  if (spec->opt_k<0) {
    k=0;
    for (i=0; i<ns; i=j) {
      for (j=i; j<ns && site[i].clust==site[j].clust; j++);
      if (j-i==1) {
        site[i].keep=FREE; /* single atom */
        k++; } }
    prt("! %d single atoms made free",k); }
  
  memset(stat,0,sizeof(stat));
  loop (i,0,ns) stat[site[i].keep]++;
  prt("! %d atoms kept  %d free",stat[INJAIL],stat[FREE]);
  
  if (opt_k&8) readMARKED(spec,0);
  if (opt_k&4) {
    readMARKED(spec,1);
    memset(stat,0,sizeof(stat));
    loop (i,0,ns) stat[site[i].keep]++;
    prt("! replaced by: %d atoms kept  %d free",stat[INJAIL],stat[FREE]); }
}
