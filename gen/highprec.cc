/*
  8/14 hyperbolic and inverse functions added

  3/09: Val reconsidered
  
  3/05: atan,asin,acos added, tan=sin/cos (less efficient) added

  6/04: sin,cos added

  bug fixed 10/97: compare of very big vs. small numbers was wrong for SEX==1.
  this error propagated to many other functions !


  added : #ifdef ABORTONERROR abort() is called if num exception

  gcc 2.96 and higher: requires -fno-strict-aliasing with optimization -O2
*/

#include "ground.h"
#include "highprec.h"

#define INTBITS 32 // # of bits in unsigned int: only INTBITS==32 expected

#define MSB (1U<<(INTBITS-1))
#define BASE (2e0*MSB)
#ifndef MAXEXP
#define MAXEXP 16384 // range is 256^-MAXEXP..256^MAXEXP
                     // smaller numbers call exception ufl and become 0
                     // bigger numbers call exception ofl
#endif

#if 4*MAXEXP < HIGHPREC
#error too low exponent limit
#endif                     

#define DOUBLEBITS 52 // accuracy of double mantissa, used by setIT() (was: 53)

// unsigned64 must be 8 bytes long integer (INTBITS==32 assumed)
#if SEX==1
// this is good for many machines...
typedef long long unsigned unsigned64;
#else
// this is SGI version; on other machines, long long unsigned may be good enough
typedef __uint64_t unsigned64;
#endif

static void err(const char *msg) /************************************** err */
{
  static int nexc,mod=1;
  nexc++;
  if (nexc%mod==0) fprintf(stderr,"\ahighprec exception: %s\n",msg);
  if (nexc>10*mod) {
    mod*=10;
    nexc=0;
    fprintf(stderr,
      "WARNING: only every %d-th highprec exception will be reported\n",
       mod); }
#ifdef ABORTONERROR
  abort();
#endif
}

#if SEX
inline void transvestit(unsigned *u)
// changes endianess of unsigned (INTBITS==32 assumed)
{
  char aux,*charu=(char*)u;
  aux=charu[0]; charu[0]=charu[3]; charu[3]=aux;
  aux=charu[1]; charu[1]=charu[2]; charu[2]=aux;
}
#endif

void Shift(unsigned u[],int s,int behind=0) /************************* Shift */
// big endian shift by +s bytes to the right or -s bytes to the left,
// padded by zeros if necessary
// small endian (SEX==1): make endian change before and after
{
  behind+=PREC*sizeof(unsigned);

#if SEX
  // swap bytes for small endians
  int i,n=(behind+3)/4;
  loop (i,0,n) transvestit(u+i);
#endif

  if (s<0) {
    memmove(u,(char*)u-s,behind+s);
    memset((char*)u+behind+s,0,-s); }
  else if (s>0) {
    memmove((char*)u+s,u,behind-s);
#if SEX
    memset(u,u[0]&128?255:0,s);
#else
    memset(u,u[0]&MSB?255:0,s);
#endif
    }

#if SEX
  // and swap back...
  loop (i,0,n) transvestit(u+i);
#endif
}

int Int(const Real &r) /************************************************ Int */
{
  int n;
  if (r.e>4) err("Int(Real) ofl");
  if (r.e<=0) return 0;

  n=(int)r.m[0];
  return n>>(8*(4-r.e));
}

void RPut(Real r) /**************************************************** RPut */
{
  int i,I,decdig,sg,decexp=0;

#if 0
  // for debugging purposes

  printf("%20.15g",Val(r));
  printf(" %4d",r.e);
  loop (i,0,PREC) printf(" %08x",r.m[i]);
  printf("\n");
#else
  sg=sign(r);
  if (sg<0) r=-r;

  if (r==0) { prt("   0.0"); return; }

  if (r<0.001) do { r*=10; decexp--; } while (r<1);
  else if (r>=1e6) do { r/=10; decexp++; } while (r>=10);

  I=Int(r);
  if (sg>0 || I) prt_("%6d.",sg*I);
  else prt_("%6s.","-0");

  decdig=(int)(PREC*INTBITS/3.3219280948873622)-1;
  loop (i,0,decdig) {
    r-=I;
    r*=10; I=Int(r);
    prt_("%01d",I);
    if (i%10==9) prtc(' '); }

  if (decexp) prt_("e%d",decexp);
  _n
#endif
}

Real::Real(double r) /***************************************** Real(double) */
{ 
  int i,s,borrow;

  e=0;
  if (r==0) {
    loop (i,0,PREC) m[i]=0;
    return; }

  if ( (s=(r<0)) ) r=-r;

  while (r>=0.5) {
    r/=256; e++;
    if (e>MAXEXP) { err("Real(double) NaN/Inf"); break; } }

  while (r<1./512.) {
    r*=256; e--;
    if (e<-MAXEXP) { err("Real(double) NaN"); break; } }

  loop (i,0,PREC) {
    r*=BASE;
    r=r-(m[i]=(unsigned4)r); }

  if (s) {
    borrow=0;
    for (i=PREC-1; i>=0; i--)
      if (borrow) m[i]=~m[i];
      else if (m[i]) { borrow++; m[i]=~m[i]+1; } }

}

Real::Real(int r)  /********************************************** Real(int) */
{ 
  *this=(double)r;
}

Real::Real(unsigned r)  /************************************ Real(unsigned) */
{ 
  *this=(double)r;
  /*.....if (r&MSB) err("Real(unsigned) ofl");*/
}

Real::Real(const char* s) /************************************ Real(char *) */
{
  const char *c=s;
  int sg=1,dp=0;
  Real dec=1;

  while (strchr(" \t\n",*c)) c++;

  if (*c=='+') c++;
  if (*c=='-') { c++; sg=-1; }

  if (!strchr(".0123456789",*c)) err("Real(char*) format");

  *this=0;

  for (; *c; c++) {
    if (strchr("0123456789",*c))
      if (dp) { dec/=10; *this += dec*((*c-'0')*sg); }
      else *this=*this*10+(*c-'0')*sg;
    else if (*c=='.')
      if (dp) err("Real(char*) two .");
      else dp++;
    else if (strchr("eE",*c)) {
      dp=atoi(c+1);
      if (dp>0) while (dp--) *this*=10;
      else if (dp<0) while (dp++) *this/=10;
      break; }
    else err("Real(char*) format"); }
}

Real::Real(int ex,unsigned ma[]) /****************** Real(exponent,mantissa) */
// this is private
{ 
  int i;
  e=ex;
  loop (i,0,PREC) m[i]=ma[i];
}

double Val(const Real &r) /************************************ Val(Real) */
{
  double d=0;
  int i,s=sign(r)<0;
  Real rr;

  if (s) rr=-r; else rr=r;
  for (i=PREC-1; i>=0; i--) d=(d+rr.m[i])*(1e0/BASE);
  if (s) d=-d;

  return d*powi(2,rr.e*8);
}

double Val(const double &r) /******************************************* Val */
{
  return r;
}

double Val(const int &r) /********************************************** Val */
{
  return (double)r;
}

double Val(const unsigned &r) /***************************************** Val */

{
  return (double)r;
}

void Normalize(int borrow,int &e,unsigned m[],int behind=0) /***** Normalize */
{
  int i;
  unsigned char *c=(unsigned char*)m;

  if (borrow) {
    Shift(m,1); e++;
#if SEX
    c[3]
#else
      c[0]
#endif
      =(m[0]&MSB)?0:255; }

  if (e<-MAXEXP) {
    err("ufl - replaced by zero");
    loop (i,0,PREC) m[i]=0;
    e=0; return; }

  loop (i,0,PREC) if (m[i]) goto nozero;
    e=0; return;
 nozero:

  while (m[0]==0 && !(m[1]&MSB)|| (m[0]==~0U) && (m[1]&MSB)) {
    Shift(m,-4,behind);
    e-=4; }

#if SEX
  while (c[3]==0 && !(c[2]&128)|| (c[3]==255) && (c[2]&128))
#else
  while (c[0]==0 && !(c[1]&128)|| (c[0]==255) && (c[1]&128))
#endif
    {
      Shift(m,-1,behind);
      e--;
    }

  if (e>MAXEXP) err("ofl");
}

int iszero(Real a) /************************************************* iszero */
{
  int i;

  loop (i,0,PREC) if (a.m[0]) return 0;
  return 1;
}

int compare(Real a,Real b) /**************************************** compare */
// private trichotomic compare: a>b gives +1, a==b gives 0, a<b gives -1
{
  int sa=sign(a),sb=sign(b);

  if (sa<sb) return -1;
  if (sa>sb) return 1;

  if (sa==0) return -sb;
  if (sb==0) return  sa;

/*.....int de=b.e-a.e;*/
/*.....*/
/*.....if (de>PREC*4) return -sb;*/
/*.....if (de<-PREC*4) return sa;*/

  // now a,b are nonzero and of the same sign 

  if (a.e<b.e) return -sa;
  if (a.e>b.e) return sa;

#if SEX
  // swap bytes for small endians
  int i;
  transvestit((unsigned*)&a.e);
  transvestit((unsigned*)&b.e);
  loop (i,0,PREC) {
    transvestit(a.m+i);
    transvestit(b.m+i); }
#endif

// WARNING: sex-dependent, this is the big-endian version

  return memcmp(&a,&b,sizeof(Real));
}

int operator==(Real a,Real b) /********************************** operator== */
{
  return compare(a,b)==0;
}

int operator!=(Real a,Real b) /********************************** operator!= */
{
  return compare(a,b)!=0;
}

int operator>(Real a,Real b) /************************************ operator> */
{
  return compare(a,b)>0;
}

int operator<(Real a,Real b) /************************************ operator< */
{
  return compare(a,b)<0;
}

int operator<=(Real a,Real b) /********************************** operator<= */
{
  return compare(a,b)<=0;
}

int operator>=(Real a,Real b) /********************************** operator>= */
{
  return compare(a,b)>=0;
}

Real operator+(Real a,Real b) /*********************************** operator+ */
{
  Real ret;
  int i,borrow=0,de;

  de=b.e-a.e;

  if (iszero(a)) return b;
  if (iszero(b)) return a;

  if (de>PREC*4) return b;
  if (de<-PREC*4) return a;

  ret.e=a.e;
  if (de>0) { Shift(a.m,de); ret.e=b.e; }
  if (de<0) Shift(b.m,-de);

  for (i=PREC-1; i>=0; i--) {
    ret.m[i]=a.m[i]+b.m[i]+borrow;
    borrow = (a.m[i]&MSB)^(b.m[i]&MSB) ? !(ret.m[i]&MSB) : !!(a.m[i]&MSB); }

  Normalize((a.m[0]&MSB)==(b.m[0]&MSB) && (a.m[0]&MSB)!=(ret.m[0]&MSB),ret.e,ret.m);

  return ret;
}


void Real::operator+=(Real b) /********************************** operator+= */
{
  *this = *this + b;
}


Real operator-(Real a,Real b) /*********************************** operator- */
{
  Real ret;
  int i,borrow=0,de;

  de=b.e-a.e;

  if (iszero(a)) return -b;
  if (iszero(b)) return a;

  if (de>PREC*4) return -b;
  if (de<-PREC*4) return a;

  ret.e=a.e;
  if (de>0) { Shift(a.m,de); ret.e=b.e; }
  if (de<0) Shift(b.m,-de);

  for (i=PREC-1; i>=0; i--) {
    ret.m[i]=a.m[i]-b.m[i]-borrow;
    borrow = (a.m[i]&MSB)^(b.m[i]&MSB) ? !(a.m[i]&MSB) : !!(ret.m[i]&MSB); }
 
  Normalize((a.m[0]&MSB)!=(b.m[0]&MSB) && (a.m[0]&MSB)!=(ret.m[0]&MSB),ret.e,ret.m);

  return ret;
}

void Real::operator-=(Real b) /********************************** operator-= */
{
  *this = *this-b;
}


Real operator-(Real b) /********************************** operator- (unary) */
{
  Real ret;
  int i,borrow=0;

  if (b.m[0]!=MSB) goto OK;
  loop (i,1,PREC) if (b.m[i]) goto OK;

    /* special case: -0.5*256^e must be shifted */
    Shift(b.m,1); b.e++;

 OK:

  ret.e=b.e;

  for (i=PREC-1; i>=0; i--)
    if (borrow) ret.m[i]=~b.m[i];
    else if (b.m[i]) { borrow++; ret.m[i]=~b.m[i]+1; }

  return ret;
}

int sign(Real r) /***************************************************** sign */
{
  int i;

  loop (i,0,PREC) if (r.m[i]) return r.m[0]&MSB ? -1:1;
  return 0;
}

Real abs(Real r) /****************************************************** abs */
{
  if (sign(r)<0) return -r; else return r;
}

Real fabs(Real r) /**************************************************** fabs */
{
  if (sign(r)<0) return -r; else return r;
}

Real operator*(Real a, Real b) /********************************** operator* */
{
  int i,j,m,phase;
  unsigned _mpl[PREC*2+1],*mpl=_mpl+2;
  unsigned b0,borrow;
  int neg=0;

  if (sign(a)<0) { a=-a; neg=!neg; }
  if (sign(b)<0) { b=-b; neg=!neg; }

  memset(_mpl,0,sizeof(_mpl));

  /* this is machine-dependent */
  unsigned64 ii;
  unsigned *ih=((unsigned*)&ii)+SEX,*il=((unsigned*)&ii)+(1-SEX);

  loop (i,0,PREC) {
    m=PREC-i+(i>0);
    /*.....  m=PREC-1;*/

    loop (phase,0,2) {
      borrow=0;
      for (j=m-phase; j>=0; j-=2) {
	ii=(unsigned64)a.m[i]*(unsigned64)b.m[j]; // machine-dependent!
	b0=mpl[i+j]&MSB;
	mpl[i+j]+=*il+borrow;
	borrow = b0^(*il&MSB) ? !(mpl[i+j]&MSB) : !!b0;
	b0=mpl[i+j-1]&MSB;
	mpl[i+j-1]+=*ih+borrow;
	borrow = b0^(*ih&MSB) ? !(mpl[i+j-1]&MSB) : !!b0; }

      j+=i;
      if (borrow) for (;;) {
	if (j<-2) err("op*: internal");
	borrow=mpl[j]&MSB;
	mpl[j]++;
	if (!borrow) break;
	if (mpl[j]&MSB) break; 
	j--; }
      } }

  i=a.e+b.e+4;
  Normalize(0,i,_mpl,8);

  if (neg) return -Real(i,_mpl);
  else return Real(i,_mpl);
}

static int IT_oneover=-1,IT_sqrt=-1;

/*
#define SETIT
  If SETIT is #defined, the best IT_sqrt and IT_oneover are determined,
  otherwise approximate formulas are used; an accuracy check of 1/r and 
  sqrt(r) is still performed.
*/

void setIT(void)
/* sets # of iterations for 1/r and sqrt(r): these functions use double seed */
{

  // sqrt: linear iterations, # of significant bin digits = (IT_sqrt+1)*DOUBLEBITS
  IT_sqrt=(PREC*INTBITS-8)/DOUBLEBITS;

  // 1/r: quadratic, # of significant bin digits = 2^IT_sqrt*DOUBLEBITS
  IT_oneover=1+(int)(log((double)(PREC*INTBITS-9)/DOUBLEBITS)/log(2.0));

  // internal accuracy check (note PREC incl. in eps)
  REAL eps=PREC*powi((REAL)0.5,PREC*INTBITS-8);
  REAL x;

#ifdef SETIT
  printf("\nformulas: IT_oneover=%d  IT_sqrt=%d\n",IT_oneover, IT_sqrt);
  IT_oneover=0; IT_sqrt=0;
 again:
#endif

  for (x=0.7; x<65000; x*=9.3) {
    REAL e=fabs((1/x)*x-1);

    if (e>eps) {
#ifdef SETIT
      IT_oneover++; goto again;
#else
      fprintf(stderr,"\nhighprec internal error: low accuracy (%g) of 1/%g\n\
use lower -O or increase IT_oneover=%d and recompile!\n",Val(e),Val(x),IT_oneover);
      exit(-1);
#endif
    }

    if (fabs(sqr(sqrt(x))-x)>eps*x) {
#ifdef SETIT
      IT_sqrt++; goto again;
#else
      fprintf(stderr,"\nhighprec internal error: low accuracy of sqrt(%g)\n\
increase IT_sqrt=%d and recompile!\n",Val(x),IT_sqrt);
      exit(-1);
#endif
      }
    }

#ifdef SETIT
  printf("\ncalculated: IT_oneover=%d  IT_sqrt=%d\n",IT_oneover, IT_sqrt);
#endif
}

Real oneover(Real r) /********************************************** oneover */
/*
  1/r for operator/
  Newton method starting from double estimate
*/
{
  if (r==0) { err("1/0"); return r; }

  int i=r.e;
  r.e=0;
  Real x(1./Val(r));
  r.e=i;
  x.e-=i;

  if (IT_oneover<0) setIT();

  loop (i,0,IT_oneover) x += x*(1-r*x);

  return x;
}

Real sqrt(Real r) /**************************************************** sqrt */
/*
  Uses LINEAR iterations avoiding Real division, using double initial estimate.
  Becomes less efficient than Newton for PREC>10 (approx.).
*/
{
  if (r<0) { err("sqrt neg arg"); return r; }
  else if (r==0) return r;

  double d=sqrt(Val(r));
  Real x(d),f(0.5/d);
  int i;

  if (IT_sqrt<0) setIT();
  
  loop (i,0,IT_sqrt) x += (r-x*x)*f;

  return x;
}

Real operator/(Real a, Real b) /********************************** operator/ */
{
  Real ret;
  ret=a*oneover(b);
  return ret;
}

void divbyint(unsigned *num,unsigned den,unsigned rprec=PREC) /*** divbyint */
// no normalization!
{
  unsigned int i;

  if (den==0) { err("div by 0"); return; }

  /* this is machine-dependent */
  unsigned64 ii=0,dden=den;
  unsigned *ih=((unsigned*)&ii)+SEX,*il=((unsigned*)&ii)+(1-SEX);

  if (num[0]&MSB) *il=den-1;
  
  loop (i,0,rprec) {
    *ih=*il; *il=num[i];
    ii=ii-dden*(num[i]=ii/dden); }
}

void mplbyint(unsigned *num,unsigned den,unsigned rprec=PREC) /*** mplbyint */
// no normalization!
{
  int i;

  /* this is machine-dependent */
  unsigned64 ii=0;
  unsigned *ih=((unsigned*)&ii)+SEX,*il=((unsigned*)&ii)+(1-SEX);

  if (num[0]&MSB) *il=den-1;

  for (i=rprec-1; i>=0; i--) {
    ii=(unsigned64)num[i]*(unsigned64)den+(unsigned64)(*ih);
    num[i]=*il; }
}

Real operator/(Real a, int b) /*************************** operator Real/int */
{
  unsigned mm[PREC+1];
  int ee;

  mm[PREC]=0;
  
  if (b>0) {
    copy(mm,a.m,PREC*sizeof(unsigned));
    ee=a.e;
    divbyint(mm,(unsigned)b,PREC+1); }
  else if (b<0) {
    a=-a;
    ee=a.e;
    copy(mm,a.m,PREC*sizeof(unsigned));
    divbyint(mm,(unsigned)(-b),PREC+1); }
  else err("div by (int)0");
  
  Normalize(0,ee,mm,4);
  
  return Real(ee,mm);
}

Real Real::operator*=(Real b) /********************************** operator*= */
{
  *this=*this*b;
  return *this;
}

void Real::operator*=(int b) /******************************** operator*=int */
{
  unsigned mm[PREC+1];
  int ee=this->e+4;
  int neg;

  if ( (neg=b<0) ) b=-b;

  copy(mm+1,this->m,PREC*sizeof(unsigned));
  mm[0]=mm[1]&MSB?~0U:0;

  mplbyint(mm,b,PREC+1);

  Normalize(0,ee,mm,4);

#if 0
  // GCC does not accept this:
  *this=Real(ee,mm);
  if (neg) *this=-*this;
#else // equivalent
  Real aux=Real(ee,mm);
  if (neg) *this=-aux; else *this=aux;
#endif
}

Real operator*(Real a, int i) /*************************** operator Real*int */
{
  a*=i;
  return a;
}

Real operator*(int i,Real a) /**************************** operator int*real */
{
  a*=i;
  return a;
}

void Real::operator/=(int i) /******************************* operator /=int */
{ 
  *this=*this/i; 
}

static Real highlog2;
static int key=2;

static double gxxbug(void)
{
  return pow(0.5,PREC*INTBITS-8);
}

Real log(Real x) /****************************************************** log */
{
  //Real eps=(Real)pow(0.5,PREC*INTBITS-8);
  Real eps=(Real)gxxbug();
  int expt; // `uninitialized' warning with not smart enough compilers
  Real sum;

  if (key==2) { key=1; highlog2=log((Real)2); key=0; }

  if (!key) {

    expt=8*x.e; x.e=0; // now 1/512<x<1/2

    if (x<=0) { err("log(<=0)"); return sum; }

    // to be optimized later...
    if (x<0.0441941738241592) { expt-=4; x*=16; }
    while (x<0.7071067811865476) { expt--; x*=2; }

    x=sqrt(x); /* to increase a bit the rate of convergence */ }

  Real y=(x-1)/(1+x),yy=y*y,y4=y*4,q=y4;
  int i=1;

  do {
    q *= yy; i += 2; sum += q/i;
    q *= yy; i += 2; sum += q/i;
    q *= yy; i += 2; sum += q/i;
  } while (fabs(q)>eps);

  sum += y4;

  /*.....printf("log: %d terms\n",i/2);*/
  
  if (key) return sum/2;
  else return sum + expt*highlog2;
}

Real exp(Real x) /****************************************************** exp */
{
  double d=Val(x);
  int expt=0,s=sign(x);

  if (key) { Real dummy=log((Real)1); /* to calculate highlog2 */ }

  expt=(int)(d*1.4426950408889634+0.5*s);

  Real y=x-highlog2*expt,sum,q=y;
  //Real eps=(Real)pow(0.5,PREC*INTBITS-8);
  Real eps=(Real)gxxbug();

  int i=2;

  do {
    q = q*y/i; i++; sum += q;
    q = q*y/i; i++; sum += q;
    q = q*y/i; i++; sum += q;
  } while (fabs(q)>eps);

  sum += y;
  sum += 1;
 
  if (expt) {
    i=expt>0?expt/8:(expt-7)/8;
    sum.e += i;
    expt -= i*8;
    sum *= (Real)(1<<expt); }

  return sum;
}

Real tanh(Real x) /**************************************************** tanh */
{
  if (fabs(x)<0.3465736) {
    return sinh(x)/cosh(x); }
  else {
    Real y=exp(2*x);

    return (y-1)/(y+1); }
}

Real sinh(Real x) /**************************************************** sinh */
{
  if (fabs(x)<0.3465736) {
    int i=2;
    Real q=x,y=x*x,sum=0;
    Real eps=(Real)gxxbug();

    do {
      q = q*y/(i*(i+1)); i+=2; sum += q;
    } while (fabs(q)>eps);
    sum+=x;

    return sum; }

  else {
    Real y=exp(x);

    return (y-1/y)/2; }
}

Real cosh(Real x) /**************************************************** tanh */
{
  Real y=exp(x);

  return (y+1/y)/2;
}

Real asinh(Real x) /************************************************* asinh */
{
  return log(x+sqrt(x*x+1)); // BUG: precision loss close to x=0
}

Real acosh(Real x) /************************************************* acosh */
{
  return log(x+sqrt(x*x-1));
}

Real atanh(Real x) /************************************************* atanh */
{
  return log((1+x)/(1-x))/2; // BUG: precision loss close to x=0, cf. atan()
}

Real pow(Real a,Real b) /*********************************************** pow */
{
  if (b==1) return a;
  else if (a==0)
    if (b>0) return (REAL)0;
    else { err("0^nonpositive"); return (REAL)0; }
  else return exp(log(a)*b);
}

#include "powi.c"

Real sqr(Real x) { return x*x; }

Real cub(Real x) { return x*x*x; }

Real pow4(Real x) { Real y=x*x; return y*y; }

Real pow5(Real x) { Real y=x*x; return y*y*x; }

Real pow6(Real x) { Real y=x*x*x; return y*y; }

#if 0 // old version

const Real PI="3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975666";

#else

// good until HIGHPREC=25 (see pi.cc)
// WARNING: some stupid compilers cannot accept longer PiMantissa
// so that this #if obfuscation
static unsigned PiMantissa[PREC]={52707178U,2290459400U
#if PREC>2
  ,3541244298U
#if PREC>3
  ,771977331U
#if PREC>4
  ,1151600952U
#if PREC>5
  ,573153073U
#if PREC>6
  ,3490197242U
#if PREC>7
  ,2565623404U
#if PREC>8
  ,2303010849U
#if PREC>9
  ,3862482963U
#if PREC>10
  ,2008962150U
#if PREC>11
  ,3476351244U
#if PREC>12
  ,1824566313U
#if PREC>13
  ,3083435088U
#if PREC>14
  ,3711927509U
#if PREC>15
  ,3048556297U
#if PREC>16
  ,395450069U
#if PREC>17
  ,3649665531U
#if PREC>18
  ,466694411U
#if PREC>19
  ,2795036597U
#if PREC>20
  ,2888826226U
#if PREC>21
  ,3687848671U
#if PREC>22
  ,3082346927U
#if PREC>23
  ,3983156862U
#if PREC>24
  ,2528803984U
#if PREC>25 && !defined(CALCPI)
#error PI not accurate enough, remove this error message, run first pi.cc for PRECISION+1 and paste the result here!
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
};

const Real PI(1,PiMantissa);
#endif

static int imove(int i)
{
  if (i>0) return (i+1)/2;
  else return i/2;
}

Real sin(Real x) /****************************************************** sin */
{
  double d=Val(x);

  if (fabs(d)<0.78539816339745) {
    int i=2;
    Real q=x,y=x*x,sum=0;
    Real eps=(Real)gxxbug();

    do {
      q = q*y/(-i*(i+1)); i+=2; sum += q;
    } while (fabs(q)>eps);
    sum+=x;

    return sum; }

  else {
    d=d/(3.14159265358979323846/4);
    int i=(int)d;

    if (d<0) i--;
    switch (i&7) {
      case 0: return sin(x-imove(i)*(PI/2));
      case 1: return cos(x-imove(i)*(PI/2));
      case 2: return cos(x-imove(i)*(PI/2));
      case 3: return -sin(x-imove(i)*(PI/2));
      case 4: return -sin(x-imove(i)*(PI/2));
      case 5: return -cos(x-imove(i)*(PI/2));
      case 6: return -cos(x-imove(i)*(PI/2));
      case 7: return sin(x-imove(i)*(PI/2));
      default: err("internal"); return 0; } }
}

Real cos(Real x) /****************************************************** sin */
{
  double d=Val(x);

  if (fabs(d)<0.78539816339745) {
    int i=3;
    Real y=x*x,qq=-y/2,q=qq,sum=0;
    Real eps=(Real)gxxbug();

    do {
      q = q*y/(-i*(i+1)); i+=2; sum += q;
      } while (fabs(q)>eps);
    sum+=qq;
    sum+=1;

    return sum; }

  else {
    d=d/(3.14159265358979323846/4);
    int i=(int)d;

    if (d<0) i--;

    switch (i&7) {
      case 0: return cos(x-imove(i)*(PI/2));
      case 1: return -sin(x-imove(i)*(PI/2));
      case 2: return -sin(x-imove(i)*(PI/2));
      case 3: return -cos(x-imove(i)*(PI/2));
      case 4: return -cos(x-imove(i)*(PI/2));
      case 5: return sin(x-imove(i)*(PI/2));
      case 6: return sin(x-imove(i)*(PI/2));
      case 7: return cos(x-imove(i)*(PI/2));
      default: err("internal"); return 0; } }
}

Real tan(Real x) /****************************************************** tan */
{
  return sin(x)/cos(x); // cheap!
}

Real atan(Real x) /**************************************************** atan */
{
  if (x>1) return (PI/2)-atan(1/x);
  if (x<-1) return -(PI/2)-atan(1/x);
  if (fabs(x)>0.4142135623731) return 2*atan(x/(sqrt(1+x*x)+1));

  int i=1;
  Real y=x*x,q=x,sum=0;
  Real eps=(Real)gxxbug();

  do {
    q = -q*y; i+=2; sum += q/i;
    } while (fabs(q)>eps);
  sum+=x;

  return sum;
}

Real asin(Real x) /**************************************************** asin */
{
  if (x==-1) return -PI/2;
  if (x==1) return PI/2;
  return 2*atan(x/(1+sqrt(1-x*x)));
}

Real acos(Real x) /**************************************************** asin */
{
  if (x==-1) return -PI;
  if (x==1) return 0;
  return 2*atan(sqrt((1-x)/(1+x)));
}

