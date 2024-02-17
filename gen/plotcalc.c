#include "ground.h" /* #define CALC 3 required */
#include "alloc.h"
#include "draw.h"
#include "plotcalc.h"
#include "mydraw.h"

#ifdef DOS
#define LINELEN 1024
#else
#define LINELEN 16384
#endif

double *count;
struct _Idlist_s *id;

int drawfile(char *fn,char *colx,char *coly,int maxcol,int color,int style)
{
FILE *f;
char line[LINELEN];
int i,OK;
double xy[2];

if (!fn || fn[0]==0) f=stdin;
else f=fopen(fn,"rt");

if (!f) {
  fprintf(stderr,"\"%s\" does not exist\n",fn);
  return -1; }

mysetcolor(color);
if (style>1) mysetlinestyle(0,0,2*style+1);
else if (style>0) mysetlinestyle(style,0,style);
else pointsize=(-1-style)*0.798;

myup();
*count=0;

while (fgets(line,LINELEN,f)) {
  char *c,*e;
  int n;

  for (c=line; *c!=0 && *c<=' '; c++);
  if (c[0]=='#' || c[0]=='!' || c[0]==0) {
    myup();
    if (c[0]==0) *count=0; }
  else {
    /* scan line: _Id.head->val = 1st number, etc. */
    e=c; n=0; id=_Id.head->next; /* note: _Id.head=count */
    do {
      /* patch to consider , as separator */
      while (*e==',') e++;
      c=e;
      id->val=strtod(c,&e);
      id=id->next;
      n++;
      } while (e>c && id);

    OK=1;
    loop (i,0,2) {
      double z;
      char *expr=i?coly:colx;

      e=expr; /* err */
      if (maxcol<n) z=Calc(expr,&e); else OK=0;
      if (e==expr) OK=0;
      if (OK) xy[i]=z; }
    *count+=1.0;
    if (OK) {
      if (style>=0) mydraw(xy[0],xy[1]);
      else mysmartdot(xy[0],xy[1]); }
    else
      myup(); } }

fclose(f);
return 0;
}

int minmaxfile(char *fn,int cont,char *colx,char *coly,int maxcol,
               double *X0,double *X1,double *Y0,double *Y1)
{
FILE *f;
char line[LINELEN];
int i,OK;
double xy[2];

*count=0;

if (cont==0) {
  *X0=*Y0=3e33;
  *X1=*Y1=-3e33; }

if (fn==NULL || (f=fopen(fn,"rt"))==NULL) {
  fprintf(stderr,"\"%s\" does not exist\n",fn);
  return -1; }

while (fgets(line,LINELEN,f)) {
  char *c,*e;
  int n;

  for (c=line; *c!=0 && *c<=' '; c++);
  if (c[0]=='#' || c[0]=='!' || c[0]==0) {
    if (c[0]==0) *count=0; }
  else {
    /* scan line: _Id.head->val = 1st number, etc. */
    e=c; n=0; id=_Id.head->next; /* note: _Id.head=count */
    do {
      /* patch to consider , as separator */
      while (*e==',') e++;
      c=e;
      id->val=strtod(c,&e);
      id=id->next;
      n++;
      } while (e>c && id);

    OK=1;
    loop (i,0,2) {
      double z;
      char *expr=i?coly:colx;

      e=expr; /* err */
      if (maxcol<n) z=Calc(expr,&e); else OK=0;
      if (e==expr) OK=0;
      if (OK) xy[i]=z; }
    *count+=1.0;
    if (OK) {
      Min(*X0,xy[0]) Max(*X1,xy[0])
      Min(*Y0,xy[1]) Max(*Y1,xy[1]) } } }

fclose(f);
return 0;
}
