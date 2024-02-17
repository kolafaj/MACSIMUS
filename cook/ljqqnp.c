/***
  elementary site-site interactions for the linked-cell list method, nonpolar
  - see ljqqpol.c for POLAR
  - Lennard-Jones + Coulomb (r-space Ewald or short-range)
  - #included by lc.c twice (with and without #define MEASURE)

  calculates LJ+QQ forces and energy of all sites in list l-> with site ls
  atom numbers in the pair (in comparison with internp.c):
    ls         ... atom 1 (with sitesite table start ssi and charge qi)
    l,l->next  ... atom 2
    dr = r1-r2 ... ls->r-l->r

  #define MEASURE: measurement version, otherwise no-measurement
  #define TEST: see forces.c
***/

#undef DUMP0
#undef DUMP1
#define DUMP0 /*void*/
#define DUMP1(A,B,C,D) /*void*/

#ifdef MEASURE

#  ifdef DUMP
/* dump of pair energies: only partly implemented, cf. internp.c 
   MOLECULE numbers printed (SITES in internp.c)
*/
double Uelst0;
#    define DUMP0 Uelst0=ENEL;
#    define DUMP1(A,B,C,D) dumpsslc(A,B,C,D);
static void dumpsslc(char *info,linklist_t *l1,linklist_t *l2,double u)
{
  if (l1->n>0 || l2->n>0) return;
  prt("%-7s: %3d-%-3d LJ=%-12g QQ=%-12g  %.7f %.7f %.7f  %.7f %.7f %.7f",
      info,
      l1->n,l2->n,
      u,ENEL-Uelst0,
      VARG(l1->r),VARG(l2->r));
}
#  endif /*#!DUMP */

#  define XADDFORCES MADDFORCES
static void LJqqm                                              /****** LJqqm */
#else /*# MEASURE */
#  define XADDFORCES ADDFORCES
static void LJqq                                                /****** LJqq */
#endif /*#!MEASURE */
  (linklist_t *l, linklist_t *ls, sitesite_t *ssi,double qi
#if PARALLEL==1
                                                          ,pll_global_t *glob
#else /*# PARALLEL==1 */
#  ifdef MEASURE
                                                          ,double *Uptr
#  endif /*# MEASURE */
#endif /*#!PARALLEL==1 */
                                                                              )

{
  vector dr;
  int i,j,sp;
  real f,x,y,z,rr,qq,fdr;
#ifdef MEASURE
  real U;
  rdf_t *ssrdf;
#endif /*# MEASURE */
  exception_t *exc;
  sitesite_t *ss;
  ermacrodcl

  for (; l; l=l->next) {

    if (l->n<FROM && ls->n<FROM) continue;

    VVV(dr,=ls->r,-l->r)

#ifdef TEST
    if (box.cutoff<1 || box.cutoff>50) ERROR(("cutoff=%g suspicious",box.cutoff))
    loop (i,0,3) {
      if (dr[i]<box.cutoff-box.L[i])
        ERROR(("%c %f<%f %f",i+'x',ls->r[i],l->r[i],dr[i]/box.Lh[i]))
      if (dr[i]>box.L[i]-box.cutoff)
        ERROR(("%c %f>%f %f",i+'x',ls->r[i],l->r[i],dr[i]/box.Lh[i])) }
#endif /*# TEST */

    if ((rr=SQR(dr))>=box.cq) goto nothing;

#ifdef TEST
    if (rr<0.01)  {
#  define NNN(r1,r2) \
  sqr(fabs(fabs(r1[0]-r2[0])-box.Lh[0])-box.Lh[0]) +\
  sqr(fabs(fabs(r1[1]-r2[1])-box.Lh[1])-box.Lh[1]) +\
  sqr(fabs(fabs(r1[2]-r2[2])-box.Lh[2])-box.Lh[2])
      int i,j;
      loop (i,0,No.s) if (NNN(a[0]->rp[i],l->r)<1e-9) prt("l:  %d",i);
      loop (i,0,No.s) if (NNN(a[0]->rp[i],ls->r)<1e-9) prt("ls: %d",i);
      WARNING(("sqrt(rr)=%f",sqrt(rr))); }
#endif /*# TEST */

    ss=ssi+l->si->st;
#ifdef MEASURE
    U=0;
#endif /*# MEASURE */

    if (rr<excrrlimit && ls->n==l->n) { /************** possible exception ***/
      sp=ls->sp;
      i=ls->si-spec[sp]->si;
      j=l->si-spec[sp]->si;
#ifdef TEST
      if (l->sp!=sp || i==j) ERROR(("l->sp=%d sp=%d  i=%d j=%d",l->sp,sp,i,j))
#endif /*# TEST */

      if (i>j) { sp=i,i=j,j=sp; exc=ls->si->exc; }
      else exc=l->si->exc;

      /* now i<j and exc[*] is the ordered exception table of (j) */
      /* it is ended by exc[*]=j */

      while (exc->indx<=i) {
        if (exc->indx==i) {
          /* this pair is an exception */
          Max(MAX14,rr)
          f=0;

#if COULOMB<=-1
/* -------------------------------------------------------------------- EXCL */

          if (exc->type==EXCL) {
            /* pair excluded (no LJ, no qq); fix for q-q for Ewald necessary */
            if ( (qq=qi*l->si->charge) ) {
#  ifdef MEASURE
#    if COULOMB==-2
              ENEL+=exacterud_sqrt(rr,&f)*qq; f*=qq;
#    else /*# COULOMB==-2 */
              ENEL+=(erud(rr)-(y=1/x))*qq; f=(byerd-y/rr)*qq;
#    endif /*#!COULOMB==-2 */
#  else /*# MEASURE */
#    if COULOMB==-2
              exacterud_sqrt(rr,&f); f*=qq;
#    else /*# COULOMB==-2 */
              f=(erd(rr)-sqrt(rr)/Sqr(rr))*qq;
#    endif /*#!COULOMB==-2 */
#  endif /*#!MEASURE */
            }
          } /* 1-2, 1-3 exclusion */

/* --------------------------------------------------------------------- 1-4 */

          else {
            /* now exc->type==ONEFOUR: special for 1-4 interaction */
            ss=&sstab14[l->si->st][ls->si->st];
            DUMP0
#  ifdef MEASURE
            if (rr < ss->C2q) {
              if (rr < ss->C1q)
                SS_MEASURE
              else {
                x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
              ENVIR+=rr*f; }
            if ( (qq=qi*l->si->charge) ) {
              y=factor14_1/sqrt(rr);
              ENEL+=(erud(rr)+y)*qq; f+=(byerd+y/rr)*qq; }
#  else /*# MEASURE */
            if (rr < ss->C2q) {
              if (rr < ss->C1q)
                SS_NOMEASURE(=)
              else {
                x=rr-ss->C2q; f=ss->A4*x; } }
            if ( (qq=qi*l->si->charge) )
              f+=(erd(rr)+factor14_1/sqrt(rr)/rr)*qq;
#  endif /*#!MEASURE */
            DUMP1("LJQQM14",ls,l,U)
            }
#elif COULOMB==2 || COULOMB==0
          /* (no Ewald EXCL fix for COULOMB==0,2) - 1-4 only */

          if (exc->type==ONEFOUR) {
            ss=&sstab14[l->si->st][ls->si->st];
#  ifdef MEASURE
            if (rr < ss->C2q) {
              if (rr < ss->C1q)
                SS_MEASURE
              else {
                x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
              ENVIR+=rr*f; }
            if ( (qq=qi*l->si->charge) ) {
              qq*=factor14;
              ENEL+=erud(rr)*qq; f+=byerd*qq; }
#  else /*# MEASURE */
            if (rr < ss->C2q) {
              if (rr < ss->C1q)
                SS_NOMEASURE(=)
              else {
                x=rr-ss->C2q; f=ss->A4*x; } }
            if ( (qq=qi*l->si->charge) )
              f+=erd(rr)*factor14*qq;
#  endif /*#!MEASURE */
          } /* 1-4 */
#endif /*#!COULOMB<=-1!COULOMB==2 || COULOMB==0 */
          goto addUf; } /* is exception */
        exc++; } /* exception loop */
      } /* possible exception */

/* --------------------------------------------------------------- nonbonded */

    f=0;
    DUMP0
#ifdef MEASURE
    if (rr < ss->C2q) {
      if (rr < ss->C1q)
        SS_MEASURE
      else {
        x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
      ENVIR+=rr*f; }

    if ( (qq=qi*l->si->charge) ) {
      ENEL+=erud(rr)*qq; f+=byerd*qq; }

#  if PARALLEL==1
    ssrdf=ss->rdfs[glob->ith];
#  else /*# PARALLEL==1 */
    ssrdf=ss->rdf;
#  endif /*#!PARALLEL==1 */
    if (ssrdf) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
#else /*# MEASURE */
    if (rr < ss->C2q) {
      if (rr < ss->C1q)
        SS_NOMEASURE(=)
      else {
        x=rr-ss->C2q; f=ss->A4*x; } }

    if ( (qq=qi*l->si->charge) ) f+=erd(rr)*qq;
#endif /*#!MEASURE */

    DUMP1("LJQQM",ls,l,U)

   addUf:
#if PARALLEL==1
#  ifdef MEASURE
    glob->U+=U;
#  endif /*# MEASURE */
    XADDFORCES(ls->f1,l->f2,f,dr)
#else /*# PARALLEL==1 */
#  ifdef MEASURE
    *Uptr += U;
#  endif /*# MEASURE */
    XADDFORCES(ls->f,l->f,f,dr)
#endif /*#!PARALLEL==1 */

#ifdef TEST /* serial */
  loop (i,0,3) {
      if (fabs(ls->f[i])>1e10) ERROR(("ls->f[%d]=%g",i,ls->f[i]))
      if (fabs(l->f[i])>1e10) ERROR(("l->f[%d]=%g",i,l->f[i])) }
#endif /*# TEST */

 nothing:; }

} /* LJqq[m] */

#undef XADDFORCES
