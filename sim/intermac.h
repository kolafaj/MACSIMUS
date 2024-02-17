/* macros for boundary conditions and measurements in 3D
   DEFAULT: periodic b.c.
   VERSIONS: SLIT (former WALL), FREEBC, NIBC
   MODIFIERS: CUT, PERSUM, WORM
   see also old+misc/intermac-before-begin-end.h

   WARNING: watercut.c, waterion.c, watertuc.c NOT UPGRADED

   Usage:
     vector r1,r2; // input
     vector dr;    // output difference, with b.c.
     double rr;    // SQR(dr), with b.c.
     beginCellImage(RET,CQ)
       now: is in cutoff (n.a. for NIBC)
       dr[],rr available
       if measure: calculate energy and forces
       else:       calculate forces
     endCellImage

   Note: unless NIBC it is possible (see xforces.c)
     beginCellImage(RET,CQ)
       is in cutoff
     } else {
       is not in cutoff
     endCellImage

  beginDistance .. endDistance do not check cutoff and do not run
  over images (PERSUM only)

     RET = used as 'return RET;' for return if outside of cutoff
           typical values: 0 if energy returned
           empty arg in void functions: 'return ;'
     CQ  = squared cutoff,
           CQ=box.cq with a few exceptions:
             pure LJ, where CQ=ss->C2q
             cluster analysis
   ENPVIR = En.Pvir (serial) or glob->Pvir (parallel)

  optimized add forces macro: call by ADDFORCES(f1,f2,f*dr)
  f1=forces to atom 1 (vector),
  f2=forces to atom 2 (vector),
  f = value of force/|dr|,
  dr=r1-r2 (vector)
  real fdr  must be known
*/

#define ADDFORCES(F1,F2,F,DR) { \
  F1[0] += fdr = F*DR[0]; F2[0] -= fdr; \
  F1[1] += fdr = F*DR[1]; F2[1] -= fdr; \
  F1[2] += fdr = F*DR[2]; F2[2] -= fdr; }

#if PRESSURETENSOR&PT_VIR
/* NOTE: now unnecessarily duplicated ENVIR+- */
#  if PRESSURETENSOR&PT_OFF
#    define MADDFORCES(F1,F2,F,DR) { \
  F1[0] += fdr = F*DR[0]; F2[0] -= fdr; \
  ENPVIR[0]+=DR[0]*fdr; \
  ENPVIR[5]+=DR[1]*fdr; \
  F1[1] += fdr = F*DR[1]; F2[1] -= fdr; \
  ENPVIR[1]+=DR[1]*fdr; \
  ENPVIR[3]+=DR[2]*fdr; \
  F1[2] += fdr = F*DR[2]; F2[2] -= fdr; \
  ENPVIR[2]+=DR[2]*fdr; \
  ENPVIR[4]+=DR[0]*fdr; }
#  else /*# PRESSURETENSOR&PT_OFF */
#    define MADDFORCES(F1,F2,F,DR) { \
  F1[0] += fdr = F*DR[0]; F2[0] -= fdr; \
  ENPVIR[0]+=DR[0]*fdr; \
  F1[1] += fdr = F*DR[1]; F2[1] -= fdr; \
  ENPVIR[1]+=DR[1]*fdr; \
  F1[2] += fdr = F*DR[2]; F2[2] -= fdr; \
  ENPVIR[2]+=DR[2]*fdr; }
#  endif /*#!PRESSURETENSOR&PT_OFF */
#else /*# PRESSURETENSOR&PT_VIR */
#  define MADDFORCES(F1,F2,F,DR) ADDFORCES(F1,F2,F,DR)
#endif /*#!PRESSURETENSOR&PT_VIR */

#ifdef GOLD
/*** charge inversion (GOLD surface) - not tested recently */
#  include "goldmac.h"
#else  /*# GOLD  */

#  define NI_f(I)  dr[I]=r1[I]-r2[I];

#  ifdef SLIT /* plain SLIT, not GOLD */
#    define NI2 NI_f
#  else /*# SLIT */
#    define NI2 NI
#  endif /*#!SLIT */

/* NO sum, NO nearest image - for exceptions -12, 1-3, 1-4 only
   If you are not sure that the 1-2, 1-3, 1-4 atoms are within cutoff,
   use beginCellImage .. endCellImage instead (not for PERSUM) */

#  define beginDistance(RET,CQ) NI_f(0) NI_f(1) NI_f(2) rr=SQR(dr);
#  define endDistance /* void */

#  ifdef FREEBC
/*** free boundary conditions */

#    define beginCellImage(RET,CQ) NI_f(0) NI_f(1) NI_f(2) rr=SQR(dr); if (rr<CQ) {
#    define endCellImage }

#    ifdef NIBC
#      error cannot FREEBC && NIBC
#    endif /*# NIBC */

#  else
/*** periodic b.c. ***/

#    ifdef WORM
  /* modifyier:* molecules of length up to 3/2*L are allowed */
#      define NI(I) \
  dr[I]=r1[I]-r2[I]; \
  if      (dr[I]<-box.Lh[I]) { dr[I]+=box.L[I]; if (dr[I]<-box.Lh[I]) dr[I]+=box.L[I]; } \
  else if (dr[I]> box.Lh[I]) { dr[I]-=box.L[I]; if (dr[I]> box.Lh[I]) dr[I]-=box.L[I]; }
#    else /*# WORM */
  /* molecules must be shorter than L/2 */
#      define NI(I) \
   dr[I]=r1[I]-r2[I]; \
   if (dr[I]<-box.Lh[I]) dr[I]+=box.L[I]; \
   else if (dr[I]>box.Lh[I]) dr[I]-=box.L[I];
#    endif /*#!WORM */

#    ifdef CUT
/* efficient if cutoff is significantly less than half box size
   METAL: see also eldens.c - to be improved !!! */
#      define beginCellImage(RET,CQ) \
  NI(0) if ((rr =dr[0]*dr[0])>=CQ) return RET; \
  NI(1) if ((rr+=dr[1]*dr[1])>=CQ) return RET; \
  NI2(2)     rr+=dr[2]*dr[2]; \
  if (rr<CQ) {
#    define endCellImage }

#    elif defined(PERSUM) /* cutoff longer than half box */

#    if defined(NIBC) || defined(CUT)
#      error illegal combination of PERSUM with NIBC or CUT
#    endif

/* just a bit optimized */
#      define beginCellImage(RET,CQ) { \
    int ix[3]; \
    double rr0,rr1; \
    loopto (ix[0],-No.nimg[0],No.nimg[0]) { \
      dr[0]=r1[0]-r2[0]+ix[0]*box.L[0]; \
      rr0=dr[0]*dr[0]; \
      if (rr0<CQ) loopto (ix[1],-No.nimg[1],No.nimg[1]) { \
        dr[1]=r1[1]-r2[1]+ix[1]*box.L[1]; \
        rr1=rr0+dr[1]*dr[1]; \
        if (rr1<CQ) loopto (ix[2],-No.nimg[2],No.nimg[2]) { \
          dr[2]=r1[2]-r2[2]+ix[2]*box.L[2]; \
          rr=rr1+dr[2]*dr[2]; \
          if (rr<CQ) {
#      define endCellImage }}}}}

#    elif defined(NIBC) /* nearest-image b.c without spherical cutoff */
#      define beginCellImage(RET,CQ)  NI(0) NI(1) NI2(2) rr=SQR(dr);
#      define endCellImage /* no closing } */

#    else  /*#!CUT !defined(NIBC) */
/*** standard periodic boundary conditions */
#      define beginCellImage(RET,CQ)  NI(0) NI(1) NI2(2) rr=SQR(dr); if (rr<CQ) {
#      define endCellImage }
#    endif /*#!CUT !defined(NIBC) */

#  endif  /*#!FREEBC */
#endif /*#!GOLD  */
