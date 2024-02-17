/* C2c ljqqtabpol.C; cppnest ljqqtabpol.c

  #####################################################
  USE VERSION WITH EXTENSION .C, then run C2c + cppnest
  #####################################################

  !!! see ljqqnp.c for notation !!!

  QQTAB version of ljqqpol.c, compatible with GAUSSIANCHARGES
  compatible with POLAR=0

***/

#ifndef QQTAB
#  error "QQTAB required"
#endif /*# QQTAB */
#if POLAR&(1+32+64)
#  error "unsupported POLAR version"
#endif /*# POLAR&(1+32+64) */
#if !defined(COULOMB)
#  error "unsupported COULOMB version"
#endif /*# !defined(COULOMB) || COULOMB>=0 */

#ifdef MEASURE
#  define XADDFORCES MADDFORCES
static void polarLJqqm                                    /****** polarLJqqm */
#else /*# MEASURE */
#  define XADDFORCES ADDFORCES
static void polarLJqq                                      /****** polarLJqq */
#endif /*#!MEASURE */
  (linklist_t *l, linklist_t *ls, sitesite_t *ssi,double qi,double qipol
#if PARALLEL==1
                                                          ,pll_global_t *glob
#else /*# PARALLEL==1 */
#  ifdef MEASURE
	,double *Uptr
#  endif /*# MEASURE */
#endif /*#!PARALLEL==1 */
                                                                              )
{
  vector dr,dr0;
  int i,j,sp;
  real f,x,y,z,rr,qq,qqpp,fdr,U;
  exception_t *exc;
  sitesite_t *ss;
#ifdef MEASURE
	rdf_t *ssrdf;
#endif /*# MEASURE */
  ermacrodcl

  for (; l; l=l->next) {
    if (l->n<FROM && ls->n<FROM) continue;

    VVV(dr,=ls->r,-l->r)

#ifdef TEST
    loop (i,0,3) {
      if (dr[i]<box.cutoff-box.L[i])
        ERROR(("%c %f<%f %f",i+'x',ls->r[i],l->r[i],dr[i]/box.Lh[i]))
      if (dr[i]>box.L[i]-box.cutoff)
        ERROR(("%c %f>%f %f",i+'x',ls->r[i],l->r[i],dr[i]/box.Lh[i])) }
#endif /*# TEST */

    if ((rr=SQR(dr))>=box.cq) goto nothing;

    ss=ssi+l->si->st;
#ifdef MEASURE
	U=0;
#endif /*# MEASURE */
    VV(dr0,=dr)

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

/* -------------------------------------------------------------------- EXCL */

          if (exc->type==EXCL) {
            /* pair excluded (no LJ, no qq); fix for q-q for Ewald necessary */

            /* center-center */
#undef ERFC
#define ERFC ss->qqtabx[0][0]
#ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#else /*# MEASURE */
	f=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr)
#endif /*#!PARALLEL==1 */

            /* pol-center */
#undef ERFC
#define ERFC ss->qqtabx[1][0]
            VVV(dr,=dr0,+ls->rpol)
            rr=SQR(dr);
#ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#else /*# MEASURE */
	f=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->f,f,dr)
#endif /*#!PARALLEL==1 */

#undef ERFC
#define ERFC ss->qqtabx[0][1]
            VVV(dr,=dr0,-=l->rpol)
            rr=SQR(dr);
#ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#else /*# MEASURE */
	f=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2pol,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->fpol,f,dr)
#endif /*#!PARALLEL==1 */

#undef ERFC
#define ERFC ss->qqtabx[1][1]
            VV(dr0,+=ls->rpol)
            rr=SQR(dr0);
#ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#else /*# MEASURE */
	f=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr0)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr0)
#endif /*#!PARALLEL==1 */
          } /* EXCL */

/* --------------------------------------------------------------------- 1-4 */

          else {
            /* now exc->type==ONEFOUR: special for 1-4 interaction */
            /*
               WARNING: this part has not been checked / debugged
               it is identical to EXCL but LJ, f+=, and qqtab instead of qqtabx
            */
            ss=&sstab14[l->si->st][ls->si->st];
#ifdef MEASURE
            if (rr < ss->C2q) {
              SS_MEASURE_rep
              if (rr < ss->C1q)
                SS_MEASURE
              else {
                x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
              ENVIR+=rr*f; }
#else /*# MEASURE */
            if (rr < ss->C2q) {
              SS_NOMEASURE_rep
              if (rr < ss->C1q)
                SS_NOMEASURE(=)
              else {
                x=rr-ss->C2q; f=ss->A4*x; } }
#endif /*#!MEASURE */

            /* note f+= because center-center added to LJ */

            /* center-center */
#undef ERFC
#define ERFC ss->qqtab[0][0]
#ifdef MEASURE
	ENEL+=erud(rr); f+=byerd;
#else /*# MEASURE */
	f+=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr)
#endif /*#!PARALLEL==1 */

            /* pol-center */
#undef ERFC
#define ERFC ss->qqtab[1][0]
            VVV(dr,=dr0,+ls->rpol)
            rr=SQR(dr);
#ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#else /*# MEASURE */
	f=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->f,f,dr)
#endif /*#!PARALLEL==1 */

#undef ERFC
#define ERFC ss->qqtab[0][1]
            VVV(dr,=dr0,-=l->rpol)
            rr=SQR(dr);
#ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#else /*# MEASURE */
	f=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2pol,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->fpol,f,dr)
#endif /*#!PARALLEL==1 */

#undef ERFC
#define ERFC ss->qqtab[1][1]
            VV(dr0,+=ls->rpol)
            rr=SQR(dr0);
#ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#else /*# MEASURE */
	f=erd(rr);
#endif /*#!MEASURE */
#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr0)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr0)
#endif /*#!PARALLEL==1 */
          } /* 1-4 */
          goto addforce; } /* is exception */
        exc++; } /* exception loop */
    } /* possible exception */

#ifdef MEASURE
#  if PARALLEL==1
	ssrdf=ss->rdfs[glob->ith];
#  else /*# PARALLEL==1 */
	ssrdf=ss->rdf;
#  endif /*#!PARALLEL==1 */
    if (ssrdf) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
#endif /*# MEASURE */

/* --------------------------------------------------------------- nonbonded */

    f=0;
    if (rr < ss->C2q) {
#ifdef MEASURE
      SS_MEASURE_rep
      if (rr < ss->C1q)
        SS_MEASURE
      else {
        x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
      ENVIR+=rr*f;
#else /*# MEASURE */
      SS_NOMEASURE_rep
      if (rr < ss->C1q)
        SS_NOMEASURE(=)
      else {
        x=rr-ss->C2q; f=ss->A4*x; }
#endif /*#!MEASURE */
    }

#if POLAR&4
    if (ls->si->neutral | l->si->neutral) {
      XADDFORCES(ls->f,l->f,f,dr) /* LJ only */
#  ifdef MEASURE
	return U;
#  else /*# MEASURE */
	return;
#  endif /*#!MEASURE */
    }

  /* ...at least one site is charged now */
#endif /*# POLAR&4 */

#if QQTAB==0

      /* for models with no center charges (COS-COS) */
      /* center-center */
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr) /* LJ */
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr) /* LJ */
#  endif /*#!PARALLEL==1 */

      /* pol-pol */
#  undef ERFC
#  define ERFC ss->qqtab[1][1]
      VVV(dr,+=ls->rpol,-l->rpol)
      rr=SQR(dr);
#  ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#  else /*# MEASURE */
	f=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr)
#  endif /*#!PARALLEL==1 */

#elif QQTAB==1

    /* both center and pol charges, inefficient for COS models only */
    VV(dr0,=dr)

    /* center-center */
#  undef ERFC
#  define ERFC ss->qqtab[0][0]
#  ifdef MEASURE
	ENEL+=erud(rr); f+=byerd;
#  else /*# MEASURE */
	f+=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr) /* charge-charge + LJ */
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr) /* charge-charge + LJ */
#  endif /*#!PARALLEL==1 */

    /* pol-center */
#  undef ERFC
#  define ERFC ss->qqtab[1][0]
    VVV(dr,=dr0,+ls->rpol)
    rr=SQR(dr);
#  ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#  else /*# MEASURE */
	f=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2,f,dr)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->f,f,dr)
#  endif /*#!PARALLEL==1 */

      /* center-pol */
#  undef ERFC
#  define ERFC ss->qqtab[0][1]
      VVV(dr,=dr0,-=l->rpol)
      rr=SQR(dr);
#  ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#  else /*# MEASURE */
	f=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2pol,f,dr)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->fpol,f,dr)
#  endif /*#!PARALLEL==1 */

      /* pol-pol */
#  undef ERFC
#  define ERFC ss->qqtab[1][1]
      VV(dr0,+=ls->rpol)
      rr=SQR(dr0);
#  ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#  else /*# MEASURE */
	f=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr0)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr0)
#  endif /*#!PARALLEL==1 */

#elif QQTAB==2

    /* both center and pol charges, optimized for center+pol/only pol (COS) */
    VV(dr0,=dr)

    /* center-center */
#  undef ERFC
#  define ERFC ss->qqtab[0][0]
    if (ERFC.qq) {
#  ifdef MEASURE
	ENEL+=erud(rr); f+=byerd;
#  else /*# MEASURE */
	f+=erd(rr);
#  endif /*#!MEASURE */
    }
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr) /* charge-charge + LJ */
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr) /* charge-charge + LJ */
#  endif /*#!PARALLEL==1 */

    /* pol-center */
#  undef ERFC
#  define ERFC ss->qqtab[1][0]
    if (ERFC.qq) {
      VVV(dr,=dr0,+ls->rpol)
      rr=SQR(dr);
#  ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#  else /*# MEASURE */
	f=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2,f,dr)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->f,f,dr)
#  endif /*#!PARALLEL==1 */
    }

      /* center-pol */
#  undef ERFC
#  define ERFC ss->qqtab[0][1]
    if (ERFC.qq) {
      VVV(dr,=dr0,-l->rpol)
      rr=SQR(dr);
#  ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#  else /*# MEASURE */
	f=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2pol,f,dr)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->fpol,f,dr)
#  endif /*#!PARALLEL==1 */
    }

      /* pol-pol */
#  undef ERFC
#  define ERFC ss->qqtab[1][1]
    /* no if (ERFC.qq) test because it is assumed that 
       pol-pol (Drude-Drude) term is always present
       (i.e., there are no nonpolar charges) */
    VVV(dr0,+=ls->rpol,-l->rpol)
    rr=SQR(dr0);
#  ifdef MEASURE
	ENEL+=erud(rr); f=byerd;
#  else /*# MEASURE */
	f=erd(rr);
#  endif /*#!MEASURE */
#  if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr0)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr0)
#  endif /*#!PARALLEL==1 */

#else /*#!QQTAB==0!QQTAB==1!QQTAB==2 */
#  error "QQTAB must be one of 0,1,2"
#endif /*#!QQTAB==0!QQTAB==1!QQTAB==2 */

   addforce: /* going here from 1--4 */
#if PARALLEL==1
#  ifdef MEASURE
    glob->U+=U;
#  endif /*# MEASURE */
    //    XADDFORCES(ls->f1,l->f2,f,dr)
#else /*# PARALLEL==1 */
#  ifdef MEASURE
    *Uptr += U;
#  endif /*# MEASURE */
    //    XADDFORCES(ls->f,l->f,f,dr)
#endif /*#!PARALLEL==1 */

#ifdef TEST /* serial */
  loop (i,0,3) {
      if (fabs(ls->f[i])>1e10) ERROR(("ls->f[%d]=%g",i,ls->f[i]))
      if (fabs(l->f[i])>1e10) ERROR(("l->f[%d]=%g",i,l->f[i])) }
#endif /*# TEST */

 nothing:; }

} /* polarLJqq[m] */

#undef XADDFORCES
