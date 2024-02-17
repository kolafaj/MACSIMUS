/*
   fitting data to a+bt+c/sqrt(t) using the same interface as linregr.[ch]
   WARNING: only selected keys available, others return 0

   MAPLE source:

     restart;
     eq1:=a*S   +b*St +c*Stih=Sy;
     eq2:=a*St  +b*St2+c*Sth =Syt;
     eq3:=a*Stih+b*Sth+c*Sti =Sytih;
     solve({eq1,eq2,eq3},{a,b,c});
     assign(%);
     d:=(-S*St2*Sti+S*Sth^2+St^2*Sti-2*Sth*St*Stih+St2*Stih^2);
     convert(simplify(a*d),string);
     convert(simplify(b*d),string);
     convert(simplify(c*d),string);

  add data by:
    void SDAdd(char *name,double weight,double X,double Y);

  print all results
    void SDPrint(char key)

  selective results:
    double SDRes(char *name,char *key)

  where key= (case insensitive)
    "?"   only check whether name exists (doesn't print error message)
    "A"   abs term A
    "B"   lin term B
    "C"   term C/sqrt(X)
    "X"   <X>
    "Y"   <Y>
    "N"   # of points recorded
    "W"   sum of weights

*/

#include "ground.h"
#include "alloc.h"
#include "sdfit.h"

#define SDNAME 24

typedef struct SD_s {
  struct SD_s *next;
  char name[SDNAME];
  long n;
  double S,St,Stih,Sy,St2,Sth,Syt,Sti,Sytih,Syy; } SD_t;

SD_t *SD_head;

void SDAdd(char *name,double weight,double X,double Y)
{
  SD_t *SD;

  for (SD=SD_head; SD; SD=SD->next)
    if (!strcmp(SD->name,name)) goto doit;

  alloconezero(SD);
  SD->next=SD_head;
  SD_head=SD;

  if (strlen(name)>SDNAME-1) {
    ERROR(("SD:%s too long name",name))
    return; }

  strcpy(SD->name,name);

doit:

  SD->n++;
  SD->S    += weight;
  SD->St   += weight*X;
  SD->Stih += weight/sqrt(X);
  SD->Sy   += weight*Y;
  SD->St2  += weight*Sqr(X);
  SD->Sth  += weight*sqrt(X);
  SD->Syt  += weight*Y*X;
  SD->Sti  += weight/X;
  SD->Sytih+= weight*Y/sqrt(X);
  SD->Syy  += weight*Y*Y;
}

static SD_t *SD;

static double SD_Res(char key)
{
  double A,B,C,d,ssq;

  if (SD->S==0) {
    fprintf(stderr,"SD_Res: zero weight\n");
    exit(0); }

  d=-SD->S*SD->St2*SD->Sti+SD->S*Sqr(SD->Sth)+Sqr(SD->St)*SD->Sti-2*SD->Sth*SD->St*SD->Stih+SD->St2*Sqr(SD->Stih);
  if (d==0) {
    fprintf(stderr,"SD->SD_Res: d: division by 0\n");
    return 0; }

  A=(-SD->St*SD->Sth*SD->Sytih+SD->St*SD->Syt*SD->Sti-SD->Stih*SD->Sth*SD->Syt+SD->Stih*SD->St2*SD->Sytih-SD->Sy*SD->St2*SD->Sti+SD->Sy*Sqr(SD->Sth))/d;
  B=(SD->S*SD->Sth*SD->Sytih-SD->S*SD->Syt*SD->Sti-SD->Stih*SD->St*SD->Sytih+SD->Syt*Sqr(SD->Stih)-SD->Sth*SD->Stih*SD->Sy+SD->St*SD->Sy*SD->Sti)/d;
  C=(-SD->Sth*SD->St*SD->Sy+Sqr(SD->St)*SD->Sytih+SD->S*SD->Sth*SD->Syt-SD->S*SD->St2*SD->Sytih-SD->Stih*SD->St*SD->Syt+SD->St2*SD->Stih*SD->Sy)/d;
  ssq=A*A*SD->S + B*B*SD->St2 + C*C*SD->Sti + SD->Syy
     +2*A*B*SD->St + 2*A*C*SD->Stih + 2*B*C*SD->Sth
     -2*A*SD->Sy - 2*B*SD->Syt - 2*C*SD->Sytih;

  switch (key) {
    case 'C': return C;
    case 'B': return B;
    case 'A': return A;
    case 'X': return SD->St/SD->S;
    case 'Y': return SD->Sy/SD->S;
    case 'N': return (double)SD->n;
    case 'W': return SD->S;
    case 'y': return sqrt(ssq/(SD->n-3));
    default: return 0; }

  return 0;
}

double SDRes(char *name,char *what)
{
  for (SD=SD_head; SD; SD=SD->next)
    if (!strcmp(SD->name,name)) {
      if (*what=='?') return 1;
      if (tolower(*what)=='d') return SD_Res(tolower(what[1]));
      return SD_Res(toupper(*what)); }
  
  if (*what!='?') prt("SD: %s unknown",name);

  return 0;
}

void SDPrint(char what)
{
  for (SD=SD_head; SD; SD=SD->next) {
    prt("\nREGRESSION A+Bt+C/sqrt(t) \"%s\" %ld points (sum weight=%g)\n\
<X>=%g  <Y>=%g  dY=%g\n\
Y = %g + %g t + %g/sqrt(t)\n",
	SD->name,SD->n,SD->S,
	SD_Res('X'),SD_Res('Y'),SD_Res('y'),
	SD_Res('A'),
	SD_Res('B'),
	SD_Res('C')); }
}

void SDFree(void)
{
  while (SD_head) {
    SD=SD_head; SD_head=SD_head->next; free(SD); }
}
