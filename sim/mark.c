#include "ground.h"
#include "simglob.h"
#include "units.h"
#include "siminit.h"
#include "norm.h"
#include "rhs.h"
#include "interpot.h"
#include "simdef.h"
#include "simils.h"
#include "pakcp.h"
#include "mark.h"

extern int sig;

struct mark_s *mark_tab;

double mark_el,U;
static int mark_n;

void mark_init(void)
{
  FILE *f=fopen(Fn(option('r')>0?"mkr":"mkp"),"rt");
  char line[128],id[32];
  int i,j,e,m;
  
  if (!f) ERROR((lastFn))
  do fgets(line,128,f); while (line[0]=='!');
  mark_n=atoi(line);
  if (mark_n>253) ERROR(("%d is too many groups",mark_n))
  alloczero(mark_tab,(mark_n+2)*sizeof(struct mark_s));
  loopto (i,0,mark_n) alloczero(mark_tab[i].set,No.s*sizeof(mark_tab[i].set[0]));
  
  loop (i,0,No.s) {
    do
      if (!fgets(line,128,f)) {
        prt("only %d sites in %s",i,lastFn);
        goto done; }
      while (line[0]=='!');
    e=sscanf(line,"%d%s%d",&j,id,&m);
    if (e==3) {
      if (m<0 || m>mark_n) ERROR((line))
      if (i!=j) ERROR(("index in %s is %d (expected %d)",lastFn,j,i))
      mark_tab[m].set[i]++;
      if (!*mark_tab[m].id1) strcpy(mark_tab[m].id1,id);
      strcpy(mark_tab[m].id2,id);
      } }
  
  done: fclose(f);
  
  header("grp from_atom   to_atom   # of at");
  loopto (m,0,mark_n) {
    j=0;
    loop (i,0,No.s) j+=mark_tab[m].set[i];
    prt("%3d %-11s %-11s %5d",m,mark_tab[m].id1,mark_tab[m].id2,j); }
  header("");
}

void mark_do(int no,int endian)
{
  int i,m;
  FILE *cp=fopen(Fn("cp"),"wb");
  float h=CPmark;
  int4 ncp;
  
  mark_init();
  
  ncp=mark_n;
  if (ncp<2) ncp=2;
  
  fwrite(&h,sizeof(h),1,cp);
  fwrite(&ncp,sizeof(ncp),1,cp);
  loopto (m,3,ncp) fwrite(m<=mark_n?mark_tab[m].id1:"-",4,1,cp);
  
  loopto (i,1,no) {
    double R;
    
    if (option('r')) {
      loadcfg(i,endian,1,&R); /* 1st pass */
      loadcfg(i,endian,1,NULL); /* 2nd pass */ }
    else 
      readplayback(0,endian);
    loopto (m,1,mark_n) 
      mark_tab[m].LJ=mark_tab[m].el/*=mark_tab[m].bond*/=0;
    measure=1;
    rhs(a[2],a[0],a[1]);
    if (i==1)
      header("#         LJ+elst             LJ           elst");
    loopto (m,1,mark_n) {
      h=mark_tab[m].LJ+mark_tab[m].el/*+mark_tab[m].bond*/;
      prt("%2d %14.6g %14.6g %14.6g",
        m,h,mark_tab[m].LJ,mark_tab[m].el/*,mark_tab[m].bond*/);
      fwrite(&h,4,1,cp); }
    h=0;
    loop (m,mark_n,ncp) fwrite(&h,4,1,cp);
    _n
    if (sig) break; }
  fclose(cp);
}
