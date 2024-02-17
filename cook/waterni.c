#undef NI
static vector dni;

/* version with nearest image shift vector dni */
#ifdef FREEBC
#define NI(I) dr[I]=r1[I]-r2[I];
#else
#define NI(I) \
   dr[I]=r1[I]-r2[I]; dni[I]=0; \
   if (dr[I]<-box.Lh[I]) dr[I]+=(dni[I]=box.L[I]); else if (dr[I]>box.Lh[I]) dr[I]+=(dni[I]=-box.L[I]);
#endif
