/* C2c internp.C

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! DO NOT EDIT internp.c                  !!!
  !!! * edit internp.C and run C2c internp.C !!!
  !!! * optionally, run cppnest internp.c    !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  with QQTAB support

  elementary site-site interactions: Lennard-Jones + r-space Ewald
  - nonpolar version
  - #included by interpot.c
  - used by the the all-pair version and the configuration initializer
    (the linked-cell list version uses ljqqnp.c instead, however, this
    module is still needed for the configuration initializer)

  WARNING: SS_NOMEASURE_rep, SS_MEASURE_rep not supported!

  NOTES:
    LJ cutoff must not be greater than the r-space Ewald cutoff
    virial is defined as  r*force = -r*u'(r)
    only non-electrostatic part of energy is returned by the functions, the
      r-space part of the electrostatic energy is summed to En.el
    virial of electrostatic forces is equal to the electrostatic energy and
       will be added to the total virial with the k-space part included
    only LJM and LJQQM needed with the link-cell method

  PERSUM: internp.c will be #included twice, once standard, second time
  periodic sum of the same molecule. This is controlled by switch PERSUM_MM.

  See intermac.h for macros:
  beginCellImage .. endCellImage runs over all required images (unless PERSUM
    just the nearest image)
  beginDistance .. endDistance just calculates the distance
*/

#ifdef SS_MEASURE_rep
#  error "SS_MEASURE_rep not supported in nonpolar version"
#endif /*# SS_MEASURE_rep */

#ifdef DUMP
/* to be turned on globally, affects internp, intrapol
   SITE numbers printed (unlike ljqqnp.c)
 */
double Uelst0;
#  define DUMP0 Uelst0=ENEL;
#  define DUMP1(A,B,C,D) dumpss(A,B,C,D);
void dumpss(char *info,vector r1,vector r2,double u)
{
  prt("%-7s: %3d-%-3d LJ=%-12g QQ=%-12g  %.7f %.7f %.7f  %.7f %.7f %.7f",
      info,
      (vector*)r1-cfg[0]->rp,(vector*)r2-cfg[0]->rp,
      u,ENEL-Uelst0,
      VARG(r1),VARG(r2));
}

#else /*# DUMP */
#  define DUMP0 /*void*/
#  define DUMP1(A,B,C,D) /*void*/
#endif /*#!DUMP */

#ifdef PERSUM_MM
/* Interaction of different images of the same site 
   Condition HALFBOXES guarantees that pairs are not calculated twice
   (it would be more efficient to modify loop ranges in beginCellImage, 
   but am too lazy).
   NB: I want to keep -Wparentheses for other purposes, that's why the
   stupid superfluous parentheses...  */
#  define HALFBOXES ix[0]>0 || (ix[0]==0 && ix[1]>0) || (ix[0]==0 && ix[1]==0 && ix[2]>0)
#endif  /*# PERSUM_MM */

#ifdef PERSUM_MM
	double LJM_MM
#else /*# PERSUM_MM */
	double LJM
#endif /*#!PERSUM_MM */
  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 GLOB) /******* LJM */
  /* LJ site-site interaction, measurements incl. */
{
  double U=0,rr,f,x,y,z;
  real fdr;
  rdf_t *ssrdf;
  vector dr;

  DUMP0

  beginCellImage(0,box.cq)
#ifdef PERSUM_MM
	if (HALFBOXES) {
#endif /*# PERSUM_MM */
    if ( (ssrdf=ss->rdf) ) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
#ifdef NIBC
    SS_MEASURE
    ENVIR+=rr*f;
    MADDFORCES(f1,f2,f,dr)
#else /*# NIBC */
    if (rr < ss->C2q) {
      if (rr < ss->C1q)
        SS_MEASURE
      else {
        x=rr-ss->C2q; U+=x*x*ss->A; f=ss->A4*x; }
      ENVIR+=rr*f;
      MADDFORCES(f1,f2,f,dr) }
#endif /*#!NIBC */
#ifdef PERSUM_MM
	}
#endif /*# PERSUM_MM */
  endCellImage

  DUMP1("LJM",r1,r2,U)
  MARK_PATCH
  return U;
} /* LJM */

#ifdef PERSUM_MM
	void QQM_MM
#else /*# PERSUM_MM */
	void QQM
#endif /*#!PERSUM_MM */
  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB)
/* real-space charge-charge interactions, measurements incl.       ***** QQM */
{
  double rr,f;
  real fdr;
  rdf_t *ssrdf;
  ermacrodcl
  vector dr;

  DUMP0
  beginCellImage(,box.cq)
#ifdef PERSUM_MM
	if (HALFBOXES) {
#endif /*# PERSUM_MM */
    SG
    if ( (ssrdf=ss->rdf) ) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
#ifdef COULOMB
#  ifdef QQTAB
	ENEL+=erud(rr); f=byerd;
#  else /*# QQTAB */
	ENEL+=erud(rr)*qq; f=byerd*qq;
#  endif /*#!QQTAB */
#else /*# COULOMB */
    ENEL+=f=qq/sqrt(rr); f/=rr;
#endif /*#!COULOMB */
    MADDFORCES(f1,f2,f,dr)
#ifdef PERSUM_MM
	}
#endif /*# PERSUM_MM */
  endCellImage

  DUMP1("QQM",r1,r2,0)
  MARK_PATCH
} /* QQM */

#ifdef PERSUM_MM
	void rdfM_MM
#else /*# PERSUM_MM */
	void rdfM
#endif /*#!PERSUM_MM */
  (sitesite_t *ss,vector r1,vector r2) /******************************* rdfM */
/* measuring rdf of sites not having q-q nor LJ (e.g. some water models) */
{
  double rr;
  rdf_t *ssrdf;
  vector dr;

  beginCellImage(,box.cq)
#ifdef PERSUM_MM
	if (HALFBOXES) {
#endif /*# PERSUM_MM */
    if ( (ssrdf=ss->rdf) ) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
#ifdef PERSUM_MM
	}
#endif /*# PERSUM_MM */
  endCellImage
} /* rdfM */

#ifdef PERSUM_MM
	double LJQQM_MM
#else /*# PERSUM_MM */
	double LJQQM
#endif /*#!PERSUM_MM */
  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB)
/* LJ + real-space charge-charge interactions, measurements incl. **** LJQQM */
{
  double U=0,rr,f,x,y,z;
  real fdr;
  rdf_t *ssrdf;
  ermacrodcl
  vector dr;

  DUMP0
  beginCellImage(0,box.cq)
#ifdef PERSUM_MM
	if (HALFBOXES) {
#endif /*# PERSUM_MM */
    SG
    if ( (ssrdf=ss->rdf) ) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }

#ifdef NIBC
    SS_MEASURE
    ENVIR+=rr*f;
#else /*# NIBC */
    f=0;
    if (rr < ss->C2q) {
      if (rr < ss->C1q)
        SS_MEASURE
      else {
        x=rr-ss->C2q; U+=x*x*ss->A; f=ss->A4*x; }
      ENVIR+=rr*f; }
#endif /*#!NIBC */
#ifdef COULOMB
#  ifdef QQTAB
	ENEL+=erud(rr); f+=byerd;
#  else /*# QQTAB */
	ENEL+=erud(rr)*qq; f+=byerd*qq;
#  endif /*#!QQTAB */
#else /*# COULOMB */
    ENEL+=x=qq/sqrt(rr); f+=x/rr;
#endif /*#!COULOMB */
    MADDFORCES(f1,f2,f,dr)
#ifdef PERSUM_MM
	}
#endif /*# PERSUM_MM */
  endCellImage

  DUMP1("LJQQM",r1,r2,U)
  MARK_PATCH
  return U;
} /* LJQQM */

#ifndef LINKCELL

#  ifdef PERSUM_MM
	void LJQQ_MM
#  else /*# PERSUM_MM */
	void LJQQ
#  endif /*#!PERSUM_MM */
  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ) /** LJQQ */
/* LJ + real-space charge-charge interactions, no measurements */
{
  double rr,f,x,y;
  real fdr;
  ermacrodcl /* warning: byerd not needed */
  vector dr;

  beginCellImage(,box.cq)
#  ifdef PERSUM_MM
	if (HALFBOXES) {
#  endif /*# PERSUM_MM */
   SG
#  ifdef COULOMB
#    ifdef QQTAB
	f=erd(rr);
#    else /*# QQTAB */
	f=erd(rr)*qq;
#    endif /*#!QQTAB */
#  else /*# COULOMB */
    f=qq/(rr*sqrt(rr));
#  endif /*#!COULOMB */
#  ifdef NIBC
    SS_NOMEASURE(+=)
#  else /*# NIBC */
    if (rr < ss->C2q) {
      if (rr < ss->C1q)
        SS_NOMEASURE(+=)
      else {
        x=rr-ss->C2q; f+=ss->A4*x; } }
#  endif /*#!NIBC */
    ADDFORCES(f1,f2,f,dr)
#  ifdef PERSUM_MM
	}
#  endif /*# PERSUM_MM */
  endCellImage

  return;
} /* LJQQ */

#  ifdef QQTAB
#    ifdef PERSUM_MM
	void QQ_MM
#    else /*# PERSUM_MM */
	void QQ
#    endif /*#!PERSUM_MM */
  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2) /************* QQ */
#  else /*# QQTAB */
#    ifdef PERSUM_MM
	void   QQ_MM
#    else /*# PERSUM_MM */
	void   QQ
#    endif /*#!PERSUM_MM */
  (vector r1,vector r2,vector f1,vector f2,double qq) /****************** QQ */
#  endif /*#!QQTAB */
{
  double rr,f;
  real fdr;
  ermacrodcl /* warning: byerd not needed */
  vector dr;

  beginCellImage(,box.cq)
#  ifdef PERSUM_MM
	if (HALFBOXES) {
#  endif /*# PERSUM_MM */
    SG
#  ifdef COULOMB
#    ifdef QQTAB
	f=erd(rr);
#    else /*# QQTAB */
	f=erd(rr)*qq;
#    endif /*#!QQTAB */
#  else /*# COULOMB */
    f=qq/(rr*sqrt(rr));
#  endif /*#!COULOMB */
    ADDFORCES(f1,f2,f,dr)
#  ifdef PERSUM_MM
	}
#  endif /*# PERSUM_MM */
  endCellImage

  return;
}

#  ifdef PERSUM_MM
	void LJ_MM
#  else /*# PERSUM_MM */
	void LJ
#  endif /*#!PERSUM_MM */
  (sitesite_t *ss,vector r1,vector r2,vector f1,vector f2) /************* LJ */
/* LJ site-site interaction, no measurements */
{
  double rr,f,x,y;
  real fdr;
  vector dr;

  beginCellImage(,ss->C2q)
#  ifdef PERSUM_MM
	if (HALFBOXES) {
#  endif /*# PERSUM_MM */
#  ifdef NIBC
    SS_NOMEASURE(=)
#  else /*# NIBC */
    if (rr < ss->C1q)
      SS_NOMEASURE(=)
    else {
      x=rr-ss->C2q; f=ss->A4*x; }
#  endif /*#!NIBC */
    ADDFORCES(f1,f2,f,dr)
#  ifdef PERSUM_MM
	}
#  endif /*# PERSUM_MM */
  endCellImage

  return;
} /* LJ */


#  ifndef PERSUM_MM

/*** exceptions 1-2, 1-3, 1-4 */

double LJQQ14M(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ GLOB)
                                                                   /* LJQQM */
/* scaled 1-4 charge-charge term (by factor14), 
   Ewald: correction (1-factor14)
   LJ: 1-4 term 
   measure version */
{
  double U=0,rr,f,x,y,z,sr;
  real fdr;
  rdf_t *ssrdf;
  ermacrodcl
  vector dr;

  DUMP0
  beginDistance(0,box.cq)
    SG
    f=0;
    sr=sqrt(rr);
    if ( (ssrdf=ss->rdf) ) {
      int i=(int)(sr*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
#    ifdef NIBC
    SS_MEASURE
    ENVIR+=rr*f;
#    else /*# NIBC */
    if (rr < ss->C2q) {
      if (rr < ss->C1q)
        SS_MEASURE
      else {
        x=rr-ss->C2q; U+=x*x*ss->A; f=ss->A4*x; }
       ENVIR+=rr*f; }
#    endif /*#!NIBC */

#    ifdef COULOMB
#      if COULOMB==2 || COULOMB==0
    qq*=factor14; ENEL+=erud(rr)*qq; f+=byerd*qq;
#      else /*# COULOMB==2 || COULOMB==0 */
#        ifdef QQTAB
	ENEL+=erud14(rr); f+=byerd;
#        else /*# QQTAB */
	ENEL+=(erud(rr)+(y=factor14_1/sr))*qq; f+=(byerd+y/rr)*qq;
#        endif /*#!QQTAB */
#      endif /*#!COULOMB==2 || COULOMB==0 */
#    else /*# COULOMB */
    ENEL+=x=factor14*qq/sr; f+=x/rr;
#    endif /*#!COULOMB */

    MADDFORCES(f1,f2,f,dr)
  endDistance

  DUMP1("LJQQM14",r1,r2,U)
  MARK_PATCH
  return U;
} /* LJQQ14M */

void LJQQ14(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 DOUBLEQQ)
                                                                   /* LJQQM */
/* as above, no measure version */
{
  double rr,f,x,y;
  real fdr;
  ermacrodcl /* warning: byerd not needed */
  vector dr;

  beginDistance(,box.cq)
    SG

#    ifdef COULOMB
#      if COULOMB==2 || COULOMB==0
    f=erd(rr)*factor14*qq;
#      else /*# COULOMB==2 || COULOMB==0 */
#        ifdef QQTAB
	f=erd14(rr);
#        else /*# QQTAB */
	f=(erd(rr)+factor14_1/(rr*sqrt(rr)))*qq;
#        endif /*#!QQTAB */
#      endif /*#!COULOMB==2 || COULOMB==0 */
#    else /*# COULOMB */
    f=factor14*qq/(rr*sqrt(rr));
#    endif /*#!COULOMB */

#    ifdef NIBC
    SS_NOMEASURE(+=)
#    else /*# NIBC */
    if (rr < ss->C2q) {
      if (rr < ss->C1q)
        SS_NOMEASURE(+=)
      else {
        x=rr-ss->C2q; f+=ss->A4*x; } }
#    endif /*#!NIBC */
    ADDFORCES(f1,f2,f,dr)
  endDistance
} /* LJQQ14 */

#    ifdef COULOMB

#      ifdef QQTAB
#        if COULOMB!=-2
#          error "Implementation limitation: Ewald corrrections require exacteru(d). Use COULOMB=-2"
#          define DOUBLEQQ , double qq
#        endif /*# COULOMB!=-2 */
#      endif /*# QQTAB */

double XQQM(vector r1,vector r2,vector f1,vector f2,double qq GLOB)  /* XQQM */
/* corrections to 1-2 and 1-3 interactions because of Ewald, measure */
{
#      if COULOMB==2 || COULOMB==0 /* cutoff electrostatics */
  return 0;
#      else /*# COULOMB==2 || COULOMB==0 */
  real fdr;
  double U=0,rr,f,x;
  /*rdf_t *ssrdf;*/
#        if COULOMB>=-1
	double y; ermacrodcl
#        endif /*# COULOMB>=-1 */
  vector dr;

  beginDistance(0,box.cq)
    SG
    x=sqrt(rr);
    /* ?! do not calculate rdf */
#        if COULOMB==-2
	ENEL+=exacterud_sqrt(rr,&f)*qq; f*=qq;
#        else /*# COULOMB==-2 */
	ENEL+=(erud(rr)-(y=1/x))*qq; f=(byerd-y/rr)*qq;
#        endif /*#!COULOMB==-2 */
    MADDFORCES(f1,f2,f,dr)
  endDistance
  MARK_PATCH
  return U;
#      endif /*#!COULOMB==2 || COULOMB==0 */
} /* XQQM */

void XQQ(vector r1,vector r2,vector f1,vector f2 DOUBLEQQ) /********** XQQ */
/* corrections to 1-2 and 1-3 interactions because of Ewald, no measure */
/* parameter ss never used ! */
{
#      if COULOMB==2 || COULOMB==0 /* cutoff electrostatics */
  return;
#      else /*# COULOMB==2 || COULOMB==0 */
  double rr,f;
  real fdr;
#        if COULOMB>=-1
	ermacrodcl /* warning: byerd not needed */
#        endif /*# COULOMB>=-1 */
  vector dr;

  beginDistance(,box.cq)
    SG
#        if   COULOMB==-2
	exacterud_sqrt(rr,&f); f*=qq;
#        else /*#   COULOMB==-2 */
	f=(erd(rr)-sqrt(rr)/Sqr(rr))*qq;
#        endif /*#!  COULOMB==-2 */
    ADDFORCES(f1,f2,f,dr)
  endDistance

  return;
#      endif /*#!COULOMB==2 || COULOMB==0 */
} /* XQQ */

#    endif /*# COULOMB */

#  endif /*? !PERSUM_MM */ /*# PERSUM_MM */

#endif /*# LINKCELL */
