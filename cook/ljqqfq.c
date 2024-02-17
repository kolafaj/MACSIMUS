/* C2c ljqqfq.C; cppnest ljqqfq.c

  !!! DO NOT EDIT ljqqfq.c !!!
      * edit ljqqfq.C and run C2c ljqqfq.C
      * optionally, run cppnest ljqqfq.c

  fluctuation charge version of ljqqpol.c
  or link-cell list version if interfq.c
***/

#ifndef POLAR
#  error "POLAR not #defined"
#endif /*# POLAR */

#if !(POLAR&32)
#  error "interfq: POLAR&32 expected"
#endif /*# !(POLAR&32) */

#if POLAR&64
#  error "interfq: POLAR&64 not supported for FQ (probably OK for Drude - to be checked)"
#endif /*# POLAR&64 */

#ifdef STARS
#  error "STARS not supported for POLAR"
#endif /*# STARS */

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
#if POLAR&1
    double frep=0,Urep=0;

    int anion1=ls->si->qtype&QTYPE_SHELL1 && l->si->qtype&QTYPE_SHELL2;
    int anion2=l->si->qtype&QTYPE_SHELL1 && ls->si->qtype&QTYPE_SHELL2;
#endif /*# POLAR&1 */

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

           if (ls->si->qtype&QTYPE_DRUDE && l->si->qtype&QTYPE_DRUDE) {
            /* center-center */
            if ( (qq=qi*l->si->charge) ) {
#ifdef MEASURE
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#else /*# MEASURE */
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	exacterud_sqrt(rr,&y); f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	f=(erd(rr)-sqrt(rr)/Sqr(rr))*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr)
#endif /*#!PARALLEL==1 */
            }

            /* pol-center */
            if ( (qq=qipol*l->si->charge) ) {
              VVV(dr,=dr0,+ls->rpol)
              rr=SQR(dr);
#ifdef MEASURE
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#else /*# MEASURE */
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	exacterud_sqrt(rr,&y); f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	f=(erd(rr)-sqrt(rr)/Sqr(rr))*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->f,f,dr)
#endif /*#!PARALLEL==1 */
            }

            if (l->si->chargepol) {

              VVV(dr,=dr0,-=l->rpol)

              /* center-pol */
              if (qi) {
                qq=qi*l->si->chargepol;
                rr=SQR(dr);
#ifdef MEASURE
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#else /*# MEASURE */
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	exacterud_sqrt(rr,&y); f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	f=(erd(rr)-sqrt(rr)/Sqr(rr))*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2pol,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->fpol,f,dr)
#endif /*#!PARALLEL==1 */
              }

              if (qipol) {
                qq=qipol*l->si->chargepol;
                VV(dr0,+=ls->rpol)
                rr=SQR(dr0);
#ifdef MEASURE
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	ENEL+=exacterud_sqrt(rr,&y)*qq; f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	ENEL+=(erud(rr)-(y=1/sqrt(rr)))*qq; f=(byerd-y/rr)*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#else /*# MEASURE */
#  if COULOMB>=0
	/* no exclusion */
#  elif COULOMB==-2
	exacterud_sqrt(rr,&y); f=y*qq;
#  else /*#!COULOMB>=0!COULOMB==-2 */
	f=(erd(rr)-sqrt(rr)/Sqr(rr))*qq;
#  endif /*#!COULOMB>=0!COULOMB==-2 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr0)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr0)
#endif /*#!PARALLEL==1 */
              } }
           }

           /*** FQ-FQ ***/
           else if (ls->si->qtype&QTYPE_FQ && l->si->qtype&QTYPE_FQ) {
#if defined(FREEBC) || defined(NIBC) || COULOMB>=0
#  error internal
             /* void */
             return 0;
#else /*# defined(FREEBC) || defined(NIBC) || COULOMB>=0 */
             qq=ls->rpol[0]*l->rpol[0];
#  if COULOMB==-2
	x=exacterud_sqrt(rr,&y); f=y*qq;
#  else /*# COULOMB==-2 */
	y=1/sqrt(rr); x=erud(rr)-y; f=(byerd-y/rr)*qq;
#  endif /*#!COULOMB==-2 */
             ENEL+=x*qq;
             ls->fpol[0]+=l->rpol[0]*x;
             l->fpol[0]+=ls->rpol[0]*x;
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr)
#  endif /*#!PARALLEL==1 */
#endif /*#!defined(FREEBC) || defined(NIBC) || COULOMB>=0 */
           }
           //#  if !(POLAR&4) // OPTIMIZE?!
           else if (ls->si->qtype+l->si->qtype==0) {
             /* uncharged nonpolarizable */
             continue;
           }
           //#  endif   // OPTIMIZE?!
           else {
             // this includes, e.g., 36 vs. 0 , which is zero interaction
             //      ERROR(("qtypes %d vs. %d not implemented",si1->qtype, si2->qtype))
           }

          } /* EXCL */

/* --------------------------------------------------------------------- 1-4 */

          else {
            /* now exc->type==ONEFOUR: special for 1-4 interaction */
            ERROR(("1-4 for FQ not implemented"))
            // INACTIVE:
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
            if ( (qq=qi*l->si->charge) ) {
#ifdef MEASURE
#  if COULOMB>=0
	ENEL+=(y=qq*factor14)*erud(rr); f+=y*byerd;
#  else /*# COULOMB>=0 */
	ENEL+=(erud(rr)+(y=factor14_1/sqrt(rr)))*qq; f+=(byerd+y/rr)*qq;
#  endif /*#!COULOMB>=0 */
#else /*# MEASURE */
#  if COULOMB>=0
	f+=qq*factor14*erd(rr);
#  else /*# COULOMB>=0 */
	f+=(erd(rr)+(factor14_1/(sqrt(rr)*rr)))*qq;
#  endif /*#!COULOMB>=0 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr)
#endif /*#!PARALLEL==1 */
            }

            /* pol-center */
            if ( (qq=qipol*l->si->charge) ) {
              VVV(dr,=dr0,+ls->rpol)
              rr=SQR(dr);
#ifdef MEASURE
#  if COULOMB>=0
	ENEL+=(y=qq*factor14)*erud(rr); f=y*byerd;
#  else /*# COULOMB>=0 */
	ENEL+=(erud(rr)+(y=factor14_1/sqrt(rr)))*qq; f=(byerd+y/rr)*qq;
#  endif /*#!COULOMB>=0 */
#else /*# MEASURE */
#  if COULOMB>=0
	f=qq*factor14*erd(rr);
#  else /*# COULOMB>=0 */
	f=(erd(rr)+(factor14_1/(sqrt(rr)*rr)))*qq;
#  endif /*#!COULOMB>=0 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->f,f,dr)
#endif /*#!PARALLEL==1 */
            }

            if (l->si->chargepol) {

              VVV(dr,=dr0,-=l->rpol)

              /* center-pol */
              if (qi) {
                qq=qi*l->si->chargepol;
                rr=SQR(dr);
#ifdef MEASURE
#  if COULOMB>=0
	ENEL+=(y=qq*factor14)*erud(rr); f=y*byerd;
#  else /*# COULOMB>=0 */
	ENEL+=(erud(rr)+(y=factor14_1/sqrt(rr)))*qq; f=(byerd+y/rr)*qq;
#  endif /*#!COULOMB>=0 */
#else /*# MEASURE */
#  if COULOMB>=0
	f=qq*factor14*erd(rr);
#  else /*# COULOMB>=0 */
	f=(erd(rr)+(factor14_1/(sqrt(rr)*rr)))*qq;
#  endif /*#!COULOMB>=0 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2pol,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->fpol,f,dr)
#endif /*#!PARALLEL==1 */
              }

              if (qipol) {
                qq=qipol*l->si->chargepol;
                VV(dr0,+=ls->rpol)
                rr=SQR(dr0);
#ifdef MEASURE
#  if COULOMB>=0
	ENEL+=(y=qq*factor14)*erud(rr); f=y*byerd;
#  else /*# COULOMB>=0 */
	ENEL+=(erud(rr)+(y=factor14_1/sqrt(rr)))*qq; f=(byerd+y/rr)*qq;
#  endif /*#!COULOMB>=0 */
#else /*# MEASURE */
#  if COULOMB>=0
	f=qq*factor14*erd(rr);
#  else /*# COULOMB>=0 */
	f=(erd(rr)+(factor14_1/(sqrt(rr)*rr)))*qq;
#  endif /*#!COULOMB>=0 */
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr0)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr0)
#endif /*#!PARALLEL==1 */
            } } }
          goto addforce; } /* 1-4 */

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
    if (ls->si->qtype+l->si->qtype==0) {
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr) /* LJ */
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr) /* LJ */
#  endif /*#!PARALLEL==1 */
#  ifdef MEASURE
	return U;
#  else /*# MEASURE */
	return;
#  endif /*#!MEASURE */
    }

  /* ...at least one site is charged now */
#endif /*# POLAR&4 */

   if (ls->si->qtype&QTYPE_DRUDE && l->si->qtype&QTYPE_DRUDE) {

    VV(dr0,=dr)

    /* center-center */
    if ( (qq=qi*l->si->charge) ) {
#if POLAR&1
#  ifdef MEASURE
      if (anion1)
        U+=ls->si->qkappa*Urep,ENVIR+=ls->si->qkappa*frep*rr,f+=ls->si->qkappa*frep;
      if (anion2)
        U+=l->si->qkappa*Urep,ENVIR+=l->si->qkappa*frep*rr,f+=l->si->qkappa*frep;
#  else /*# MEASURE */
      if (anion1) f+=ls->si->qkappa*frep;
      if (anion2) f+=l->si->qkappa*frep;
#  endif /*#!MEASURE */
#endif /*# POLAR&1 */

#ifdef MEASURE
	ENEL+=erud(rr)*qq; f+=byerd*qq;
#else /*# MEASURE */
	f+=erd(rr)*qq;
#endif /*#!MEASURE */
    }

#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr) /* charge-charge + LJ */
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr) /* charge-charge + LJ */
#endif /*#!PARALLEL==1 */

    /* pol-center */
    if ( (qq=ls->si->chargepol*l->si->charge) ) {
      VVV(dr,=dr0,+ls->rpol)
      rr=SQR(dr);

#ifdef MEASURE
	ENEL+=erud(rr)*qq; f=byerd*qq;
#else /*# MEASURE */
	f=erd(rr)*qq;
#endif /*#!MEASURE */

#if POLAR&1
      if (anion1 && rr < ss->C2q) {
#  define f fff
#  ifdef MEASURE
        double z,frep=0,Urep=0,fff,UUU;
#    define U UUU
        SS_MEASURE_rep
        if (rr < ss->C1q)
          SS_MEASURE
#    if 0 /* this applied when f==frep, now U,f not needed */
        else {
          x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
#    endif /*# 0 */
        frep*=ls->si->qkappa;
        ENVIR-=rr*frep;
#    undef U
        U-=ls->si->qkappa*Urep;
#  else /*# MEASURE */
        double z,frep=0,Urep=0,fff;
        SS_NOMEASURE_rep
        if (rr < ss->C1q)
          SS_NOMEASURE(=)
#    if 0 /* this applied when f==frep, now U,f not needed */
        else {
          x=rr-ss->C2q; f=ss->A4*x; }
#    endif /*# 0 */
        frep*=ls->si->qkappa;
        U-=ls->si->qkappa*Urep;
#  endif /*#!MEASURE */
#  undef f
        f-=frep; }
#endif /*# POLAR&1 */

#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->f,f,dr)
#endif /*#!PARALLEL==1 */
    }

    if (l->si->chargepol) {

      VVV(dr,=dr0,-=l->rpol)

      /* center-pol */
      if (qi) {
        qq=qi*l->si->chargepol;
        rr=SQR(dr);

#ifdef MEASURE
	ENEL+=erud(rr)*qq; f=byerd*qq;
#else /*# MEASURE */
	f=erd(rr)*qq;
#endif /*#!MEASURE */

#if POLAR&1
        if (anion2 && rr < ss->C2q) {
#  define f fff
#  ifdef MEASURE
          double fff,UUU;
#    define U UUU

          SS_MEASURE_rep
          if (rr < ss->C1q)
            SS_MEASURE
#    if 0 /* this applied when f==frep, now U,f not needed */
          else {
            x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
#    endif /*# 0 */
          frep*=l->si->qkappa;
          ENVIR-=rr*frep;
#    undef U
          U-=l->si->qkappa*Urep;
#  else /*# MEASURE */
          double fff;

          SS_NOMEASURE_rep
          if (rr < ss->C1q)
            SS_NOMEASURE(=)
#    if 0 /* this applied when f==frep, now U,f not needed */
          else {
            x=rr-ss->C2q; U=x*x*ss->A; f=ss->A4*x; }
#    endif /*# 0 */
          frep*=l->si->qkappa;
#  endif /*#!MEASURE */
#  undef f
          f-=frep; }
#endif /*# POLAR&1 */

#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2pol,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->fpol,f,dr)
#endif /*#!PARALLEL==1 */
      }

      /* pol-pol */
      if (qipol) {
        qq=qipol*l->si->chargepol;
        VV(dr0,+=ls->rpol)
        rr=SQR(dr0);
#ifdef MEASURE
	ENEL+=erud(rr)*qq; f=byerd*qq;
#else /*# MEASURE */
	f=erd(rr)*qq;
#endif /*#!MEASURE */

#if PARALLEL==1
	XADDFORCES(ls->f1pol,l->f2pol,f,dr0)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->fpol,l->fpol,f,dr0)
#endif /*#!PARALLEL==1 */
      } }
    }

   else if (ls->si->qtype&QTYPE_FQ && l->si->qtype&QTYPE_FQ) {
     qq=ls->rpol[0]*l->rpol[0];
#ifdef COULOMB
     y=erud(rr);
     ls->fpol[0]+=y*l->rpol[0];
     l->fpol[0]+=y*ls->rpol[0];
     ENEL+=y*qq; f+=byerd*qq;
#else /*# COULOMB */
     y=sqrt(rr);
     ls->fpol[0]+=l->rpol[0]/y;
     l->fpol[0]+=ls->rpol[0]/y;
     ENEL+=y=qq/y; f+=y/rr;
#endif /*#!COULOMB */
#if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr)
#else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr)
#endif /*#!PARALLEL==1 */
   }
#if !(POLAR&4) 
   else if (ls->si->qtype+l->si->qtype==0) {
#  if PARALLEL==1
	XADDFORCES(ls->f1,l->f2,f,dr)
#  else /*# PARALLEL==1 */
	XADDFORCES(ls->f,l->f,f,dr)
#  endif /*#!PARALLEL==1 */
   }
#endif /*? !(POLAR&4) */ /*# !(POLAR&4)  */
   else {
     // this includes, e.g., 36 vs. 0 , which is zero interaction
     //      ERROR(("qtypes %d vs. %d not implemented",si1->qtype, si2->qtype))
   }

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
