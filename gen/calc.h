/* unit order: m,kg,s,K,mol,A */

#ifndef CALC
#  error "CALC undefined"
#endif /*# CALC */

#ifndef CALC_NUNITS

#  if CALC&1
#    define CALC_NUNITS 6
typedef struct unit_s {
  char *name;
  double val;              /* former number *nr */
  double pow[CALC_NUNITS]; /* powers: m kg s K mol A (not implemented: cd sr) */
  int pu;                  /* 0=normal, 1=MACSIMUS program unit */
} CALC_t;
extern char *unitserr;
void opunits(struct unit_s *a,struct unit_s *b,char assignopb);
extern double unitfactor;
char *unit(struct unit_s *res);
#  else /*# CALC&1 */
#    define CALC_NUNITS
typedef double CALC_t;
#  endif /*#!CALC&1 */

CALC_t Calc(char *buf,char **endscan);

#endif /*# CALC_NUNITS */
