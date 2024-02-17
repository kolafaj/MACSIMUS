/* this 2D part has not been used for a long time ... */
/* NIBC not implemented */

/*
  optimized add forces macro: call by ADDFORCES(f1,f2,f*dr)
  f1=forces to atom 1 (vector),
  f2=forces to atom 2 (vector),
  f = value of force/|dr|,
  dr=r1-r2 (vector)
  real fdr  must be known
*/
#define ADDFORCES(F1,F2,F) { \
  F1[0] += fdr = F[0]; F2[0] -= fdr; \
  F1[1] += fdr = F[1]; F2[1] -= fdr; }

#ifdef FREEBC

#define NI(I)  dr[I]=r1[I]-r2[I];
#define NearestImage0 NI(0) NI(1) rr=SQR(dr); if (rr<cq)
#define NearestImage NearestImage0
      
#else /* periodic b.c. with Ewald */

#ifdef WORM
  /*** molecules of length up to 3/2*L are allowed ***/
#define NI(I) \
  dr[I]=r1[I]-r2[I]; \
  if (dr[I]<-Lh) { dr[I]+=L; if (dr[I]<-Lh) dr[I]+=L; } \
  else if (dr[I]>Lh) { dr[I]-=L; if (dr[I]>Lh) dr[I]-=L; }
#else
/*** molecules of length up to L/2 are allowed ***/
#define NI(I) \
   dr[I]=r1[I]-r2[I]; if (dr[I]<-Lh) dr[I]+=L; else if (dr[I]>Lh) dr[I]-=L;
#endif

#ifdef CUT /* efficient if cutoff is significantly less than half box size */
#define NearestImage0 \
  NI(0) if ((rr=dr[0]*dr[0])>=cq) return 0.0; \
  NI(1) rr+=dr[1]*dr[1]; \
  if (rr<cq)
#define NearestImage \
  NI(0) if ((rr=dr[0]*dr[0])>=cq) return; \
  NI(1) rr+=dr[1]*dr[1]; \
  if (rr<cq)
#else
#define NearestImage0 NI(0) NI(1) rr=SQR(dr); if (rr<cq)
#define NearestImage NearestImage0
#endif

#endif /* !FREEBC */

