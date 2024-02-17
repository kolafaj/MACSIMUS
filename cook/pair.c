/* 
   pair forces (LJ+r-space Ewald) by the direct all-pair method
   optional cache optimization
   #included by cook/forces.c #ifndef LINKCELL 
   #ifdef LOG: En.LJ, En.intra, ... recorded
   WARNING: switch LOG should be double checked -- may unnecessarily included if (!measure)
*/

static vector *Arp,*Brp;

/* pair forces between molecules m--n */
static void onepair(int n,int m) /* -------------------------------- onepair */
{
  molecule_t *mn=molec+n,*mm=molec+m;

#ifdef LOG
  double Emn = (*pot[mn->sp][mm->sp]) (rof(mn,Brp),rof(mm,Brp),mn,mm,Arp);
  if (m==n) En.intra += Emn;
  else En.pot += Emn;
  En.LJ += En.LJmn;
#else
  En.pot += (*pot[mn->sp][mm->sp]) (rof(mn,Brp),rof(mm,Brp),mn,mm,Arp);
#endif  
}

/* recursive */
static void subDivide(int x,int y, int l) /* --------------------- subDivide */
/*    #
    o o
  o o o
* o o % *=(x,y), %=(x+l,y) #=(x+l,y+l) */
{
  if (l==0)
    onepair(x,y);
  else if (abs(l)==1) {
    onepair(x,y);
    onepair(x+l,y);
    onepair(x+l,y+l); }
  else if (abs(l)==2) {
    /* unnecessary optimization */
    int lh=l/2;

    onepair(x,y);
    onepair(x+lh,y);
    onepair(x+lh,y+lh);
    onepair(x+l,y);
    onepair(x+l,y+lh);
    onepair(x+l,y+l); }
  else {
    int lh=l/2,sg=l>0?1:-1,lh1=lh+sg,lh2=l-lh1;

    subDivide(x,y,lh);
    subDivide(x+lh1,y,lh2);
    subDivide(x+lh1,y+lh1,lh2);
    if (l&1) subDivide(x+l-sg,y+lh,sg-lh2);
    else     subDivide(x+l,y+lh,-lh2); }
}

static void pairforces(ToIntPtr B, ToIntPtr A) /***************** pairforces */
{
  int n,m,sp;
  vector *f,*rp=A->rp;
  molecule_t *mn,*mm;
  double Emn;
  double Enel0=En.el; /* since V2.7c: with POLAR, En.el may be large, so we 
                         try to avoid rounding errors if small contributions
                         are added; further improved in V2.7d */
  double Enel=0;

  if (option('v')&4) prt("start pairforces: En.el=%.14g",Enel0);

  En.el=0;

#ifdef LOG  
  if (cache<=0 && FROM) {
    WARNING(("cache=0 (recursive triangulation) does not support:\n\
*** monitoring energy of the 1st molecule (No.first)\n\
*** fixed positions of several first molecules (option -j)\n\
*** => cache=64 set instead"))
      cache=64; }
#endif
  //  fprintf(stderr,"QQQ start pairforces En.el=%g Pvir=%g %g %g\n",En.el,VARG(En.Pvir));

  if (cache<=0) {
    Arp=A->rp;
    Brp=B->rp;
    /* recursive triangulation, not good for small molecules */
    subDivide(0,0,No.N-1); }
  else if (cache==1) loop (n,FROM,No.N) {
    /* simplest code, No.first accepted */
    /*** pair forces and real part of Ewald sums ***/
    /*** in the Hamilton formalism this loop may be merged to the next loop ***/
    mn=molec+n;
    sp=mn->sp;
    f=rof(mn,B->rp);
    loopto (m,0,n) {
      mm=molec+m;
#ifdef LOG
      Emn = (*pot[sp][mm->sp]) (f,rof(mm,B->rp),mn,mm,rp);
      if (m==n) En.intra += Emn;
      else En.pot += Emn;
      En.LJ += En.LJmn;
      if (No.first) {
        if (n<No.first) {
          En.LJ0+=En.LJmn; En.el0+=(En.el-Enel); En.pot0+=Emn; }
        else if (m<No.first) {
          En.LJX+=En.LJmn; En.elX+=(En.el-Enel); En.potX+=Emn; }
        Enel=En.el; }
#else
      En.pot += (*pot[sp][mm->sp]) (f,rof(mm,B->rp),mn,mm,rp);
#endif
    }
  } /* n */
  else {
    /* cache > 1: systematic blocking, no No.first */
    int n0,nf,nt,m0,mt;
    /*** pair forces and real part of Ewald sums ***/
    /*** in the Hamilton formalism this loop may be merged to the next loop ***/
    for (n0=0; n0<No.N; n0+=cache) {
      nf=max(FROM,n0);
      nt=min(No.N,n0+cache); /* not incl. */
      for (m0=0; m0<=n0; m0+=cache) {
        loop (n,nf,nt) {
          mt=min(m0+cache-1,n);
          mn=molec+n;
          sp=mn->sp;
          f=rof(mn,B->rp);
          loopto (m,m0,mt) {
            mm=molec+m;
#ifdef LOG
            Emn = (*pot[sp][mm->sp]) (f,rof(mm,B->rp),mn,mm,rp);
            if (m==n) En.intra += Emn;
            else En.pot += Emn;
            En.LJ += En.LJmn;
#else
            En.pot += (*pot[sp][mm->sp]) (f,rof(mm,B->rp),mn,mm,rp);
#endif
          } } /* m,n */
      } } /* m0,n0 */
  } /* cache */

#ifdef LOG
  if (option('v')&4)
    prt("after mol-mol: intramol=%g intermol=%g En.el=%g En.el(sum)=%.14g",
        En.intra,En.pot,En.el,En.el+Enel0);
#endif
  //  fprintf(stderr,"QQQ end pairforces ΔEel=%g Pvir=%g %g %g\n",En.el,VARG(En.Pvir));

  En.el += Enel0;
#ifdef LOG
  En.pot += En.intra;
#endif
}
