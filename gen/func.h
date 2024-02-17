/* FN(xxx) -> xxxl, xxxf, xxx according to PRECISION */
REAL FN(sqr)(REAL x);
REAL FN(cub)(REAL x);
REAL FN(pow4)(REAL x);
REAL FN(pow5)(REAL x);
REAL FN(pow6)(REAL x);
REAL FN(powi)(REAL x,int i);
REAL FN(trunc)(REAL x);

#define sign(X) ((X)>0 ? 1 : (X<0) ? -1 : 0)

#ifndef C99
/* fmin,fmax are part of C99 standard and in some implementations may
   interfere with these declarations */
REAL FN(fmax)(REAL x,REAL y);
REAL FN(fmin)(REAL x,REAL y);
#endif
