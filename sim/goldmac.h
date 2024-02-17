/* WARNING: NOT CHECKED FOR A LONG TIME 
   supported: SLIT (former WALL), FREEBC, NIBC
   macros changed as in intermac.h
*/

#define NI2(I)  dr[I]=r1[I]-r2[I];
#define Reflect0 dr[2]=r1[2]+r2[2]; rr=xxyy+Sqr(dr[2]); qq=-qq;

#define ADDFORCES2(F1,F2,F) { \
  F1[0] += fdr = F[0]; F2[0] -= fdr; \
  F1[1] += fdr = F[1]; F2[1] -= fdr; \
  F1[2] += fdr = F[2]; F2[2] += fdr; }

#ifdef FREEBC

#  define NI(I)  dr[I]=r1[I]-r2[I];
#  define NearestImage(RET,CQ) NI(0) NI(1) NI2(2) \
  rr=(xxyy=Sqr(dr[0])+Sqr(dr[1]))+Sqr(dr[2]); if (rr<CQ)
#  define Reflect Reflect0

#  ifdef NIBC
#    error cannot FREEBC && NIBC
#  endif /*# NIBC */
      
#else /* periodic b.c. with Ewald */ /*# FREEBC */

#  ifdef WORM
  /*** molecules of length up to 3/2*L are allowed ***/
#    define NI(I) \
  dr[I]=r1[I]-r2[I]; \
  if (dr[I]<-Lh[I]) { dr[I]+=L[I]; if (dr[I]<-Lh[I]) dr[I]+=L[I]; } \
  else if (dr[I]>Lh[I]) { dr[I]-=L[I]; if (dr[I]>Lh[I]) dr[I]-=L[I]; }
#  else /*# WORM */
/*** molecules of length up to L/2 are allowed ***/
#    define NI(I) \
   dr[I]=r1[I]-r2[I]; if (dr[I]<-Lh[I]) dr[I]+=L[I]; else if (dr[I]>Lh[I]) dr[I]-=L[I];
#  endif /*#!WORM */

#  ifdef CUT /* efficient if cutoff is significantly less than half box size */
#    define NearestImage(RET,CQ) \
  NI(0) if ((xxyy=dr[0]*dr[0])>=CQ) return RET; \
  NI(1) if ((xxyy+=dr[1]*dr[1])>=CQ) return RET; \
  NI2(2) rr=xxyy+dr[2]*dr[2]; \
  if (rr<CQ)
#    define Reflect Reflect0 if (rr<box.cq)

#  elif defined(NIBC) /* nearest-image b.c without spherical cutoff */
#    define NearestImage(RET,CQ) NI(0) NI(1) NI2(2) \
  rr=(xxyy=Sqr(dr[0])+Sqr(dr[1]))+Sqr(dr[2]); 

#    define Reflect Reflect0

#  else 
#    define NearestImage(RET,CQ) NI(0) NI(1) NI2(2) \
  rr=(xxyy=Sqr(dr[0])+Sqr(dr[1]))+Sqr(dr[2]); if (rr<CQ)
#    define Reflect Reflect0 if (rr<box.cq)
#  endif 

#endif 
