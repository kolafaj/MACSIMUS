/*
  linear predictor - tailored for the NSK project
  use by (example):

    REAL independent_var,var1,var2,...
    ...
    pred_begin(&independent_var,&var1,&var2,...,NULL);
    do {
      <change independent_var>
      predict(eps);// eps=precision, 1st call irrelevant
      calculate(); // don't change independent_var, calculate new var1,var2,...
      pred_save();
    }
    end_pred();

  may be nested

  if calculations crash, clean everything by:
    while (pred_end());
*/

#include "ground.h"
#include "predict.h"
#include <stdarg.h>

static pred_t *predhead;

void pred_begin(REAL *var, ...) /******************************** pred_begin */
{
  va_list va;
  int i,na;
  pred_t *pr;
  REAL *adr;

  na=0;
  va_start(va,var);
  do {
    if (na++>10) DISASTER(("pred_begin: too long or unterminated arg list"))
  } while (va_arg(va,char*)!=NULL);
  va_end(va);

  alloc(pr,sizeof(pred_t)+na*sizeof(preditem_t));
  pr->next=predhead;
  predhead=pr;

  pr->pass=0;
  pr->na=na;

  pr->item[0].var=var;

  i=0;
  va_start(va,var);
  while ( (adr=va_arg(va,REAL*)) )
    pr->item[++i].var=adr;
  va_end(va);

  pred_save();
}

void pred_save(void) /******************************************** pred_save */
{
int i;
preditem_t *item;

loop (i,0,predhead->na) {
  item=&predhead->item[i];
  item->Old=item->New;
  item->New=*item->var; }
}

void predict(REAL eps) /******************************************** predict */
{
  int i,sg;
  preditem_t *item;
  REAL Q;

  if (predhead->pass++<2) return;

  item=&predhead->item[0];
  Q=item->New-item->Old;
  if (Q==0) {
    ERROR(("predict[%d]",predhead->na))
    Q=3e33; }
  sg=Q<0?-1:1;
  Q=(*item->var-item->New)/sqrt(Q*Q+eps*eps)*sg;

  loop (i,1,predhead->na) {
    item=&predhead->item[i];
    if (item->New!=*item->var) WARNING(("predict[%d]: data changed",predhead->na))
    *item->var+=Q*(item->New-item->Old); }
}

int pred_end(void) /*********************************************** pred_end */
{
  pred_t *pr;

  if (predhead==NULL) return 0;
  pr=predhead->next;
  free(predhead);
  predhead=pr;
  
  return 1;
}
