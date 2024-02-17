#include "ground.h"
#include "solve2.h"
#include "ms.h"
/***
  solving {f(x,y)==0, g(x,y)==0}
  method=1: method of secant planes
  method=2: nested secant methods
  method=-2: nested secant methods, reversed g and r
  method=0: try 1, if fails then 2, if fails then -2
  f() and g() are called always so that f(x,y) is called before g(x,y)
    with the same parameters; thus, f() may prepare some data for g()
***/

int Solve2(func_2 f,func_2 g,
  REAL *x,REAL minx,REAL maxx,REAL epsx,
  REAL *y,REAL miny,REAL maxy,REAL epsy, int maxit,int method)
{
int m=method;
int toprt=maxit>0;

/* no detailed print if maxit<0 */
maxit=abs(maxit);

again:

if (m==-2) {
  func_2 *aux;
  aux=f; f=g; g=aux; }

switch (abs(m)) {

case 0:
case 1: { /* SECANT PLANES */

REAL x1=*x-3*epsx,y1=*y-3*epsy,x2=*x+8*epsx,y2=*y,x3=*x,y3=*y+8*epsy;
REAL f1,g1,f2=f(x2,y2),g2=g(x2,y2),f3=f(x3,y3),g3=g(x3,y3);
REAL dx2,dx3,dy2,dy3,df2,df3,dg2,dg3,ag,af,bg,bf,a,dx,dy,vp;
int done=0,it=0,isx1=0;
REAL epsq=(Sqr(epsx)+Sqr(epsy));

do {

  if (++it>maxit || x1<minx || x1>maxx || y1<miny || y1>maxy)
    if (m==0) { m=2; goto MS2; } else return -1;

  dx2=x2-x1; dy2=y2-y1;
  dx3=x3-x1; dy3=y3-y1;
  vp=dx3*dy2-dx2*dy3;
  if (fabs(vp)>(fabs(dx3*dy2)+fabs(dx2*dy3))*epsq) {
    if (!isx1) { f1=f(x1,y1); g1=g(x1,y1); isx1++; }
    df2=f2-f1; dg2=g2-g1;
    df3=f3-f1; dg3=g3-g1;
    af=dy3*df2-dy2*df3; bf=dx3*df2-dx2*df3;
    ag=dy3*dg2-dy2*dg3; bg=dx3*dg2-dx2*dg3;
    a=ag*bf-af*bg;
    if (fabs(a)>(fabs(ag*bf)+fabs(af*bg))*epsq) {
      a=vp/a;
      f3=f2; f2=f1; g3=g2; g2=g1;
      x3=x2; x2=x1; y3=y2; y2=y1;
      x1+=dx=a*(bf*g1-bg*f1);
      y1-=dy=a*(ag*f1-af*g1);
      isx1=0;
      if (toprt) prt("%20.16g %20.16g dx=%g dy=%g",Val(x1),Val(y1),Val(dx),Val(dy));
      if (fabs(dx/epsx)+fabs(dy/epsy)<1) done++;
      else done=0;
      continue; } }
  {
    REAL dx2e=dx2/epsx, dy2e=dy2/epsy;
    REAL dx3e=dx3/epsx, dy3e=dy3/epsy;
    if (Sqr(dx2e)+Sqr(dy2e)>Sqr(dx3e)+Sqr(dy3e)) { 
      if ((it+(it/5)+(it/13))&1) dx2e=-dx2e; else dy2e=-dy2e;
      x3=(x1+x2+x3)/3+epsx*dy2e*(1+it%3);
      y3=(y1+y2+y3)/3+epsy*dx2e*(1+it%3);
      f3=f(x3,y3); g3=g(x3,y3); }
    else {
      if ((it+(it/7)+(it/11))&1) dx3e=-dx3e; else dy3e=-dy3e;
      x2=(x1+x2+x3)/3+epsx*dy3e*(1+it%3);
      y2=(y1+y2+y3)/3+epsy*dx3e*(1+it%3);
      f2=f(x2,y2); g2=g(x2,y2); }
    prt("solve2: %i perturbed",it);
    done=0; }
  } while (done<1); /* NOTE: (done<2) means one more iteration */
*x=x1; *y=y1;
return 0;
}

case 2:
MS2:
{ /* NESTED SECANT METHODS */

REAL aux,X=*x,Y=*y;

MS_BEGIN(Y,epsy)
  if (MS_it>maxit)
    if (method==0 && m!=-2) { m=-2; goto again; } else return -1;

  MS_BEGIN(X,epsx)
    if (MS_it>maxit || X<minx || X>maxx || Y<miny || Y>maxy)
      if (method==0 && m!=-2) { m=-2; goto again; } else return -2;

    if (m<0) aux=g(X,Y); /* must be called before f (g:f swapped for m<0) */
    MS_f=f(X,Y);
    if (toprt) prt("%20.16g %20.16g f=%g g=%g",Val(X),Val(Y),Val(aux),Val(MS_f));
    MS_END(X,1)

  if (m>0) aux=f(X,Y); /* must be called before g(X,Y) to prepare global data */
  MS_f=g(X,Y);
  if (toprt) prt("%20.16g %20.16g f=%g g=%g",Val(X),Val(Y),Val(aux),Val(MS_f));
  MS_END(Y,1)
*x=X; *y=Y;
return 0;
}
default: ERROR(("bad method"))
}
return -3;
}

#if 0
const char *Solve2Name(int method)
/* returns constant string with method name */
{
switch (method) {
  case 0:  return "auto";
  case 1:  return "secant planes";
  case 2:  return "nested secants y:f=0 where x:g=0";
  case -2: return "nested secants y:g=0 where x:f=0"; }
return "none";
}
#endif

struct keys_s Solve2Key[] = {
  {"auto",0},
  {"SP",1},{"secant planes",1},
  {"nest1",2},{"nested secants y:f=0 where x:g=0",2},
  {"nest2",-2},{"nested secants y:g=0 where x:f=0",-2},
  {NULL,0}};
