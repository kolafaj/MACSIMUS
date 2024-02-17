/*** some ugly globals with sums of charges etc. ***/
extern struct charges_s {
  real sum;   /* sum_i q_i */
  real max;   /* max_i |q_i| */
  real sumq;  /* sum_i q_i^2 */
  real summq; /* sum_i (q_i/m_i)^2 */
  real bgE;   /* backgtound energy, to be divided by V */
#ifdef POLAR
  real fqSumq;
#else
  /* for ECC - EXPERIMENTAL */
  real molsumdq; /* sum_n mu_n^2: of uncharged molecules, rigid geometry assumed */
  real molsumq;  /* sum_n q_n^2 */
#endif /*# POLAR */
#ifdef GAUSSIANCHARGES
  real sumaq;
#endif
} charges;

#if PARALLEL
#  define CACHELINE_DBL (CACHELINE/8)
#  if PRESSURETENSOR
#    define Pvir_DIM (PT_DIM+2+CACHELINE_DBL)/CACHELINE_DBL*CACHELINE_DBL-3
#    if Pvir_DIM < PT_DIM
#      error wrong Pvir_DIM - fix me!
#    endif /*# Pvir_DIM < PT_DIM */
#  endif /*# PRESSURETENSOR */
extern struct pll_ewald_s {
  double E;
  double Ereturned;
  int4 t,padd;  /* this must have the same length as double */
#  if PRESSURETENSOR
  double Pvir[Pvir_DIM];
#  else /*# PRESSURETENSOR */
  double dpadd[CACHELINE_DBL-3];
#  endif /*#!PRESSURETENSOR */
} *pll_ewald;
#endif /*# PARALLEL */

typedef struct { real re,im; } complex;

#ifdef GAUSSIANCHARGES
typedef struct /*SDS*/ {
  int size,N;
  vector M;         /* dipole moment */
  struct Qk_s {
    complex Q,T,S; 
  } Qk[1];  /* array[Q_t.N] : sums Q */
} Q_t;
#else
typedef struct /*SDS*/ {
  int size,N;
  vector M;         /* dipole moment */
  complex Qk[1];    /* array[Q_t.N] : sums Q */
} Q_t;
#endif

extern Q_t *Q;      /* Ewald pass 1 -> pass 2 interface variable */

double Ewald(int pass,vector *frp,vector *rp);

#ifndef NIBC
int Ewaldtest(double *setL);
#endif

