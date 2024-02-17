/*
  linear regression Y=A+BX
  rewritten form PASCAL
  see also: linregr1.c

  add data by:
    void LRAdd(char *name,double weight,double X,double Y);
  
  print all results
    void LRPrint(char key)  (key='0' for Y=BX regression)
    
  selective results:
    double LRRes(char *name,char *key)

  where key= (case insensitive)
    "?"   only check whether name exists (doesn't print error message)
    "A"   abs term A
    "B"   lin term B
    "C"   corr coeff
    "X"   <X>
    "Y"   <Y>
    "N"   # of points recorded
    "W"   sum of weights
    "DA"  error of A
    "DB"  error of B
    "DY"  error of Y
    "GA"  error of A provided that weights were 1/y_error^2
    "GB"  error of B provided that weights were 1/y_error^2
*/

#include "ground.h"
#include "linregr.h"

#define LRNAME 24

typedef struct LR_s {
  struct LR_s *next;
  char name[LRNAME];
  long n;
  double sum,sumX,sumY,sumXY,sumXX,sumYY; } LR_t;

LR_t *LR_head;

void LRAdd(char *name,double weight,double X,double Y) /************** LRAdd */
{
  LR_t *LR;

  for (LR=LR_head; LR; LR=LR->next)
    if (!strcmp(LR->name,name)) goto doit;

  alloczero(LR,sizeof(LR_t));
  LR->next=LR_head;
  LR_head=LR;

  if (strlen(name)>LRNAME-1) {
    ERROR(("LR:%s too long name",name))
    return; }

  strcpy(LR->name,name);

doit:

  LR->n++;
  LR->sum+=weight;
  LR->sumX +=X*weight;
  LR->sumY +=Y*weight;
  LR->sumXY+=X*Y*weight;
  LR->sumXX+=Sqr(X)*weight;
  LR->sumYY+=Sqr(Y)*weight;
}

static int BX;
static LR_t *LR;

static double LR_Res(char key) /* ----------------------------------- LR_Res */
{
  double B=0,A,d=1,dy;

  if (LR->sum==0) {
    fprintf(stderr,"LR_Res: zero weight\n");
    exit(0); }

  if (BX) {
    if (LR->sumXX==0) fprintf(stderr,"LR_Res: B: division by 0\n");
    else B=LR->sumXY/LR->sumXX;
    A=0;
    dy=0;
    if (LR->n>1)
      dy=(LR->sumYY-LR->sumXY*B)/(double)(LR->n-1); }
  else {
    d=LR->sumXX-Sqr(LR->sumX)/LR->sum;
    if (d==0) {
      fprintf(stderr,"LR_Res: d: division by 0\n");
      return 0; }
    B=(LR->sumXY-LR->sumX*LR->sumY/LR->sum)/d;
    A=(LR->sumXX*LR->sumY-LR->sumXY*LR->sumX)/LR->sum/d;
    dy=0;
    if (LR->n>2)
      dy=((LR->sumX*LR->sumY/LR->sum-LR->sumXY)*B-Sqr(LR->sumY)/LR->sum+LR->sumYY)
	/(double)(LR->n-2); };

  switch (key) {
    case 'B': return B;
    case 'A': return A;
    case 'X': return LR->sumX/LR->sum;
    case 'Y': return LR->sumY/LR->sum;
    case 'C': {
      double s=d*(LR->sumYY-Sqr(LR->sumY)/LR->sum);
  
      if (s>=0)
        return (LR->sumXY-LR->sumX*LR->sumY/LR->sum)/sqrt(s);
      else {
        fprintf(stderr,"LR_Res: C undefined\n");
        return -9; } }
  /* D? */
    case 'y': 
      if (dy<=0) return 0;
      else return sqrt(dy);
    case 'a': 
      if (BX || dy<=0) return 0;
      else return sqrt(LR->sumXX/(LR->sumXX*LR->sum-Sqr(LR->sumX))*dy);
    case 'b':
      if (dy<=0) return 0;
      else if (BX) return sqrt(dy/LR->sumXX);
      else return sqrt(LR->sum*dy/(LR->sumXX*LR->sum-Sqr(LR->sumX)));
  
  /* G? */    
  
    case 'g': return sqrt(LR->sum*Sqr(LR->sumXX)
  			-2*Sqr(LR->sumX)*LR->sumXX+LR->sumXX*Sqr(LR->sumX))
                      /d/LR->sum;
    case 'h': return  sqrt(LR->sumXX*Sqr(LR->sum)
  			 -2*Sqr(LR->sumX)*LR->sum+LR->sum*Sqr(LR->sumX))
                      /d/LR->sum;
    case 'N': return (double)LR->n;
    case 'W': return LR->sum;
    default: ERROR(("%c is bad key",key)) }

  return 0;
}

double LRRes(char *name,char *what) /********************************* LRRes */
{
  for (LR=LR_head; LR; LR=LR->next)
    if (!strcmp(LR->name,name)) {
      if (*what=='?') return 1;
      BX=0;
      if (*what=='0') { BX=1; what++; }
      if (*what==' ') what++;
      if (what[1]==0) return LR_Res(toupper(what[0]));
      if (what[0]=='d' || what[0]=='D') return LR_Res(tolower(what[1]));
      if (what[0]=='g' || what[0]=='G')
	return LR_Res(tolower(what[1])+('a'-'g'));
      ERROR(("LR:%s is bad key",what)) }

  if (*what!='?') prt("LR: %s unknown",name);

  return 0;
}

void LRPrint(char what) /******************************************* LRPrint */
{
  BX=what=='0';

  /* Y = (%g +- %g) + (%g +- %g) X\n\ */

  for (LR=LR_head; LR; LR=LR->next) {
    prt("\nLINEAR REGRESSION \"%s\" %ld points (sum weight=%g)\n\
<X>=%g  <Y>=%g  Y=A+B*X\n\
%g %g # A stderr\n\
%g %g # B stderr\n\
corr coeff=%.7f  dy=%g",
	LR->name,LR->n,LR->sum,
	LR_Res('X'),LR_Res('Y'),
	LR_Res('A'), LR_Res('a'),
	LR_Res('B'), LR_Res('b'),
	LR_Res('C'), LR_Res('y')); }
}

void LRFree(void) /************************************************** LRFree */
{
  while (LR_head) {
    LR=LR_head; LR_head=LR_head->next; free(LR); }
}
