#include "prec.h"

typedef struct preditem_s {
  REAL *var;    /* pointer of variable */
  REAL New,Old; /* two values of the linear predictor */
  } preditem_t;

typedef struct pred_s {
  struct pred_s *next;
  int na;
  int pass;
  preditem_t item[1]; /* variable length */
  } pred_t;

void pred_begin(REAL *var, ...);
void predict(REAL eps);
void pred_save(void);
int pred_end(void); /* returns 0 if no more to free */
