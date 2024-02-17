/*
  NEW: see calc.c simplified
  NEW: merged with alloc.h, see ground.c
*/

#ifdef __cplusplus
extern "C" {
/*
  It is assumed that a conversion function/macro Val(X) returning double
  exists for each class X that is to be used in macros put, put2, put3,
  and similarly for Int (not used here).
*/
#else /*# __cplusplus */
#  define Val(X) (double)(X)
#  define Int(X) (int)(X)
#endif /*#!__cplusplus */

#ifdef SDS
#  error "SDS removed, use #include \"sds.h\" instead"
#endif

/*
  WARNING: because of some agressive #defines in ground.h, the system
  #include must be first
*/
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include <stdarg.h>

#include <string.h>
#include <math.h>   /* on stellar: #include <fastmath.h> */
#ifndef FP_NAN
#  define FP_NAN (1./(exp(0)-1)*(cos(0)-1))
#endif /*# FP_NAN */
#ifndef isfinite
#  define isfinite(X) (X>-1.797693e+308 && X<1.797693e+308)
#endif /*# isfinite */

extern FILE *in,*out;

#include "loop.h"
#include "int4.h"

#define PI_DOUBLE 3.14159265358979323846
#ifndef PI
#  define PI PI_DOUBLE
#endif /*# PI */

#ifndef Max
/* A=max(A,B) */
#  define Max(A,B) { if ((B)>(A)) (A)=(B); }
/* A=min(A,B) */
#  define Min(A,B) { if ((B)<(A)) (A)=(B); }
#endif /*# Max */

#ifndef max
/* see below for fmax and fmin */
#  define max(A,B) ((A)<(B)?(B):(A))
#  define min(A,B) ((A)<(B)?(A):(B))
#endif /*# max */

#ifndef copy
#  define copy memcpy
#endif /*# copy */

#define Sqr(X) ((X)*(X))
#define Cub(X) ((X)*(X)*(X))
#define Pow4(X) Sqr(Sqr(X))
#define Pow5(X) ((X)*Sqr(Sqr(X)))
#define Pow6(X) Cub(Sqr(X))

// cubrt() removed 7/2022

#define REAL double
#define FN(X) X
#include "func.h"
#undef REAL
#undef FN

/* to prevent optimizing a double in 10-real register */
double casttodouble(double x);
  
/*
  PRECISION==0 (float, not tested) and PRECISION==2 (long double) C versions.
  The get data module works in double.  Thus, there are two versions (double
  and given by PRECISION) of all these functions.
  If "prec.h" is used, double versions are #defined identical to
  float/long double; also REAL is #defined in prec.h and not here.
  PRECISION>2 (emulated) versions are in C++ and use overloading, so that no
  such naming conventions are necessary.
*/
#if PRECISION==2 || PRECISION==0
#  undef REAL
#  undef FN
#  define CAT(X,Y) X##Y
#  if PRECISION==2
#    define REAL long double
#    define FN(X) CAT(X,l)
#  else  /*# PRECISION==2 */
#    define REAL float
#    define FN(X) CAT(X,f)
#  endif /*#!PRECISION==2 */
#  include "func.h"
#  undef FN
#  undef CAT
#  undef REAL
#endif /*# PRECISION==2 || PRECISION==0 */

#define Sign(X) ((X)<0?-1:(X)>0?1:0)

/****** prt and prtc ******/

int prtc(char c);
int prt(const char *format, ...);
int prt_(const char *format, ...);
int prts(const char *c);
int prts_(const char *c);

/****** scroll ******/

void setscroll(int lines,int columns,int clear,int closegraph);
int initscroll(int buflen);
void scroll(void);

/****** legacy macro Error() code ******/
void extError(const char *msg, const char *f, int l);
#define Error(X) extError(X,__FILE__,__LINE__)

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     new more sophisticated version of error/warning handling (not PAR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

examples of usage (note double parentheses!):

if (length_of_rope<0) DISASTER(("rope length is negative"))
if (length_of_rope<1) ERROR(("too short rope: %f",length_of_rope))
if (length_of_rope>5) {
  WARNING(("too long rope: %f (cut to 5) ",length_of_rope))
  length_of_rope=5; }
*/

extern struct myErrorInfo_s {
  int Line;       /* to store value of __LINE__ macro */
  const char *File;     /* to point to __FILE__ macro string */
  int Level;      /* set by WARNING/ERROR/DISASTER : */
                       /* value status    batch mode  interactive mode */
                       /* 0     WARNING   continue    asks for cont. */
                       /* 1     ERROR     stop        asks for cont. */
                       /* 2     DISASTER  stop        stop */
} myErrorInfo;

void myError(const char *format, ...);

#define WARNING(_X) { myErrorInfo.Level=0; myErrorInfo.Line=__LINE__; myErrorInfo.File=__FILE__; myError _X; }
#define ERROR(_X) { myErrorInfo.Level=1; myErrorInfo.Line=__LINE__; myErrorInfo.File=__FILE__; myError _X; }
#define DISASTER(_X) { myErrorInfo.Level=2; myErrorInfo.Line=__LINE__; myErrorInfo.File=__FILE__; myError _X; }


/****** PUT DATA  ******/

extern const char *_putformat;
#define put(_X) prt(_putformat,#_X,Val(_X));
#define put_(_X) prt_(_putformat,#_X,Val(_X));
#define putkey(_X,_Y) { _PrtKey_(#_X,Int(_X),_Y); _n }
#define putkey_(_X,_Y) _PrtKey_(#_X,Int(_X),_Y);
#define put2(_X,_Y) { put_(_X) put(_Y) }
#define put3(_X,_Y,_Z) { put_(_X) put_(_Y) put(_Z) }
#define _n prtc('\n');
#define putv(_X) prt("%11s=[%13.6g %13.6g %13.6g ]", \
  #_X,Val(_X[0]),Val(_X[1]),Val(_X[2]));


/****** GET DATA ******/

struct keys_s {
  char *key;
  int val;
};

#ifndef _GetBufLen
#  define _GetBufLen 256
#endif /*# _GetBufLen */

typedef struct _Idlist_s {
  struct _Idlist_s *next;
  double val; /* value of the variable */
  int used;   /* set to 1 if read (used), needed for plot */
  char id[1]; /* [var len], identifier name */
} _Idlist;

extern struct _Id_s {
  char fmt[16];         /* output format for nr */
  char buf[_GetBufLen]; /*                AB[-123]XY = 2*(3+4) next... */
  double degrad;/* 1=radians, pi/180=degrees, ... */
  int list;     /* switch for creating the list of identifiers */
  int ignoreid; /* ignore unknown ID at rhs -- use zero value */
  _Idlist *head;/* list of identifiers */
  double nr;  /* numeric value            ^ ^     ^            ^       */
  char *b;    /* running pointer to buf   ^ ^     ^            *       */
  char *id;   /* pointer to identifier    * ^     ^                    */
  char *br;   /* pointer to '[' (if any)    *     ^                    */
  char *suff; /* pointer to suffix (after ']')    *                    */
  char *str;  /* string given by ID="..." */
  int indx;   /* value of subscript */
  int prt;    /* id=result of expr printed; prt=2: command `??' */
  int asg;    /* assign val to register (#); op for += _= ... */
  int echo;   /* toggle by ?=:
                 0:no echo, 1:copy input, 2:echo (print results)
                 -1: ignore some errors */
  int end;    /* set if ';' */
  int empty;  /* set if empty line or ";" only */
  int used;   /* set if match */
  int col;    /* packing identifiers to a line for the ?? function */
} _Id;


int _GetId(void);
int _GetKey(char *id,double val,struct keys_s *key);
int _GetCmp(char *id,double val);
int _GetVecCmp(char *id,char *suff,int maxindex,double val);
void _CheckData(void);
int _GetEnd(void);

#define getdata { *(_Id.b=_Id.buf)=0; _Id.list=1; do { if (_GetId()) {

#define get(_X) if (_GetCmp(#_X,Val(_X))) _X=_Id.nr;
#define getkey(_X,_Y) if (_GetKey(#_X,Val(_X),_Y)) _X=_Id.nr;
#define getmark(_X,_M) if (_GetCmp(#_X,Val(_X))) _M++,_X=_Id.nr;
#define getvec(_X,_S,_MI) if (_GetVecCmp(#_X,#_S,_MI,Val(_X[_Id.indx]_S))) _X[_Id.indx]_S=_Id.nr;
#define getvecmark(_X,_S,_MI,_M) if (_GetVecCmp(#_X,#_S,_MI,Val(_X[_Id.indx]_S))) _M++,_X[_Id.indx]_S=_Id.nr;

#define checkdata _CheckData();
#define enddata  } } while (!_GetEnd());  }

int prtkey_(int val,struct keys_s *key,int how);
void _PrtKey_(char *id,int val,struct keys_s *key);

/****** headers of tables ******/

void putline(char c,int l);
void header(char *h);
void underline(char *t);

void graph(double x,int linelen);
/*
Pseudo-graphics with a char/4 resolution.
For 0<=x<=1, a multi-character is printed at the position scaled by linelen.
*/

#ifdef __cplusplus
}
#endif /*# __cplusplus */

#ifdef MYSQRT
/*
  this is optional sqrt function replacement
  implementation requirements:
    sizeof(unsigned)=4
    sizeof(long unsigned)=8,
    sizeof(double)=8, IEEE format
  significant improvement detected on DEC Alpha processors
  no improvement found on PC (gcc/emx) nor R4400 SGI
  options:
    #define MYSQRT 1 : sqrt is a function
    #define MYSQRT 2 : sqrt is a macro (a bit faster)
*/
#  ifdef PRECISION
#    if PRECISION!=1
#      error MYSQRT requires PRECISION=1
#    endif /*# PRECISION!=1 */
#  endif /*# PRECISION */

#  if MYSQRT==1
/* function */
double mySqrt(double);
#    define sqrt mySqrt

#  else /*# MYSQRT==1 */
/* macro (inline), see ground.h */
extern unsigned sqrt_tab[64];
extern long unsigned sqrt_u;
extern double sqrt_t, sqrt_x;
extern union sqrt_un { long unsigned u; double d; } sqrt_ud;

#    define sqrt(_X)  (sqrt_ud.d = _X,\
  sqrt_u = 0x5fe80000UL-(sqrt_ud.u>>33),\
  sqrt_ud.u = sqrt_u-(long unsigned)sqrt_tab[(sqrt_u>>14) & 63UL]<<32,\
  sqrt_x = sqrt_ud.d,\
  sqrt_x *= 3.0-_X*sqrt_x*sqrt_x, \
  sqrt_x *= 12.0-_X*sqrt_x*sqrt_x, \
  sqrt_x *= 0.0625, \
  sqrt_t = _X*sqrt_x,\
  sqrt_t+sqrt_x*0.5*(_X-sqrt_t*sqrt_t))
#  endif /*#!MYSQRT==1 */
#endif /*# MYSQRT */

#include "fortran.h"

/*** alloc *** (c) J.Kolafa 1991, updated 1995, 2008
  merged with ground 04/2016 ***

  o Sophisticated alloc/free routines
  o Array range checking
  o Some string stuff (result of sprintf as pointer to a static string)

  WARNING: because of some agressive #defines, alloc.h MUST be
  #included after all system #include's !!!

  The actual functionality depends on following #define's:
    CHECKHEAP (checking level)
    FREENULL  (1: free(NULL) is skipped [new 1/2001, default #ifndef FREENULL]
               0: [approx. old version]: free(NULL) is error)
    SHM       shared memory - REMOVED

  The following macros are available for CHECKHEAP>=0:

  alloc(pointer,size)       allocates  size  of heap memory to  pointer

  alloczero(pointer,size)   allocates  size  of heap memory to  pointer
                            and fills the memory by zeros

  free(pointer)             deallocates memory pointed to by  pointer
                            and assigns NULL to it
                            #if FREENULL : free(NULL) is skipped
                            #if !FREENULL : free(NULL) is error

  In addition, the following macros are available for CHECKHEAP<0:

  ralloc(pointer,size)      as  alloc  and  alloczero  and all all alloc's are
  ralloczero(pointer,size)  recorded and can be freed at once by  release

  release(pointer)          works like in Pascal - all dynamic variables
                            allocated since  pointer  has been allocated
                            are freed.

  NOTE: Variables allocated by  alloc  must not be freed by  release  and
  those allocated by  ralloc  must not be freed by  free !

  The checking level depends on the absolute value of CHECKHEAP :
    0: no check - system  malloc, cmalloc,  and  free  used; DANGEROUS!
    1: heap overflow and bad calls to  free  are checked; safe for use, but
       no diagnostics is printed
    2: as 1 with extended error messages, possibility of array range check
       and alloc/free debug tracing.
       In addition,  copyarray()  contains element size check
       Use the following global variables:
         AllocRange  Number of bytes allocated before and after any memory
                     This area is checked for overwrites whenever any
                     alloc/free routine is called, as well as when
                     checkranges()  is called
                     The value of AllocRange must be a multiple of 16
                     WARNING: memory consuming and slow
         AllocTrace  If set, all alloc's and free's are traced
       Function to check memory overflow (in addition to alloc/free):
         checkranges(0) quiet
         checkranges(1) verbose (recommended)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "int4.h"

/* TYPEOF = typeof in C++
   TYPEOF = <nothing> in C
This obfuscation is needed because C smoothly casts void* into <any>*
but C++ does not and other ways to overcome this lead to many warnings
in the fundamentalistic gcc/g++ compiler,
*/

#ifdef __cplusplus
#  define TYPEOF(X) (typeof(X))
extern "C" {
#else /*# __cplusplus */
#  define TYPEOF(X) /* no cast needed */
#endif /*#!__cplusplus */

#ifndef FREENULL
#  define FREENULL 1
#endif /*# FREENULL */

#ifndef CHECKHEAP
#  define CHECKHEAP -1
#endif /*# CHECKHEAP */

#if CHECKHEAP<-2 || CHECKHEAP>2
#  error bad CHECKHEAP
#endif /*# CHECKHEAP<-2 || CHECKHEAP>2 */

extern int4 AllocSizeLim;

#include "loop.h"

#define CACHELINE 64  /* AMD 64, Opteron, Intel ... */
                      /* valid values are 32,64,128, applies to PARALLEL */

/* SHM removed - see shmgroundinclude.h, sys4par.h removed */
  
/***** alloc/free *****/

#if CHECKHEAP == 0 /* ==================================================== 0 */

/***** absolutely no check! *****/

#  define alloc(X,Y) X=TYPEOF(X) malloc(Y)
#  define alloczero(X,Y) X=TYPEOF(X) calloc(1,Y)

void *myfree(void *f);
#  define free(X) X=TYPEOF(X) myfree(X)
#  define checkranges(VERBOSE) /* nothing */

#elif CHECKHEAP == 1 || CHECKHEAP == -1 /* ========================== +1,-1 */

/***** heap overflow and bad calls to free are checked *****/

void *mymalloc(int4 size,int zero);
#  define alloc(X,Y) X=TYPEOF(X) mymalloc(Y,0)
#  define alloczero(X,Y) X=TYPEOF(X) mymalloc(Y,1)

void *myfree(void *f);
#  define free(X) X=TYPEOF(X) myfree(X)
#  define checkranges(VERBOSE) /* nothing */

#elif CHECKHEAP == 2 || CHECKHEAP == -2 /* ========================= +2,-2 */

/*
  bad calls to alloc and free are diagnosed with extended error messages
  range checks available:
    AllocRange bytes allocated before and after and checked
    AllocRange must be a multiple of 16
    (memory consuming: also a some control structures allocated!)
*/

extern int4 AllocRange, AllocTrace;

void *mymalloc(int4 size,int zero,char *c,char *fn,int line);
#  define alloc(X,Y) X=TYPEOF(X) mymalloc(Y,0,#X,__FILE__,__LINE__)
#  define alloczero(X,Y) X=TYPEOF(X) mymalloc(Y,1,#X,__FILE__,__LINE__)

char *mystrdup(const char *s,char *c,char *fn,int line);
#  define strdup(X) mystrdup(X,#X,__FILE__,__LINE__)

void mycheckranges(int verbose,char *fn,int line);
#  define checkranges(VERBOSE) mycheckranges(VERBOSE,__FILE__,__LINE__)

void *myfree(void *f,char *c,char *fn,int line);
#  define free(X) X=TYPEOF(X) myfree(X,#X,__FILE__,__LINE__)

#endif /*#!CHECKHEAP == 0!CHECKHEAP == 1 || CHECKHEAP == -1!CHECKHEAP == 2 || CHECKHEAP == -2 */


/***** ralloc/release *****/

#if CHECKHEAP == -1

void *myralloc(int4 size,int zero);
#  define ralloc(X,Y) X=TYPEOF(X) myralloc(Y,0)
#  define ralloczero(X,Y) X=TYPEOF(X) myralloc(Y,1)

void *myrelease(const void *p);
#  define release(X) X=TYPEOF(X) myrelease(X)

#elif CHECKHEAP == -2

void *myralloc(int4 size,int zero,char *c,char *fn,int line);
#  define ralloc(X,Y) X=TYPEOF(X) myralloc(Y,0,#X,__FILE__,__LINE__)
#  define ralloczero(X,Y) X=TYPEOF(X) myralloc(Y,1,#X,__FILE__,__LINE__)

void *myrelease(void *p,char *c,char *fn,int line);
#  define release(X) X=TYPEOF(X) myrelease(X,#X,__FILE__,__LINE__)

#endif /*#!CHECKHEAP == -1!CHECKHEAP == -2 */

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  string support, see mystring.c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

void stringinit(int dummy,int n);
/*
  dummy (not used - legacy)
  n=max number of strings available
  (if stringinit is not called, stringinit(128,3) is assumed)
*/

char *string(const char *format,...);
/*
  returns pointer to a static string containing the result of
  sprintf(<static string>,format,data)
  Example:
    FILE *f=fopen(string("config.%d",number));
*/
char *vstring(const char *format,va_list args);

char *strend(const char *str); /* returns pointer to the terminating 0 */
char *strlast(const char *str); /* returns pointer to the last char, or NULL */
char *int2sumbin(int i); /* human-readable binary: int2sumbin(5) -> 4+1 */

#ifdef __cplusplus
}
#endif /*# __cplusplus */

/* added 03/2010: optimize arrays */
#define alloc2Darray(X,N,M) do { int _I;   \
  alloc((X),(N)*sizeof((X)[0]));           \
  alloc((X)[0],(N)*(M)*sizeof((X)[0][0])); \
  for (_I=1; _I<(N); _I++) (X)[_I]=(X)[0]+_I*(M); } while (0);
#define ralloc2Darray(X,N,M) do { int _I;   \
  ralloc((X),(N)*sizeof((X)[0]));           \
  ralloc((X)[0],(N)*(M)*sizeof((X)[0][0])); \
  for (_I=1; _I<(N); _I++) (X)[_I]=(X)[0]+_I*(M); } while (0);
#define alloc2Darrayzero(X,N,M) do { int _I;   \
  alloc((X),(N)*sizeof((X)[0]));           \
  alloczero((X)[0],(N)*(M)*sizeof((X)[0][0])); \
  for (_I=1; _I<(N); _I++) (X)[_I]=(X)[0]+_I*(M); } while (0);
#define ralloc2Darrayzero(X,N,M) do { int _I;   \
  ralloc((X),(N)*sizeof((X)[0]));           \
  ralloczero((X)[0],(N)*(M)*sizeof((X)[0][0])); \
  for (_I=1; _I<(N); _I++) (X)[_I]=(X)[0]+_I*(M); } while (0);
#define free2Darray(X) do { free(X[0]); free(X); } while (0);

/* added 7/2000 */
#define allocarray(X,N) alloc((X),(N)*sizeof((X)[0]))
#define rallocarray(X,N) ralloc((X),(N)*sizeof((X)[0]))
#define allocarrayzero(X,N) alloczero((X),(N)*sizeof((X)[0]))
#define rallocarrayzero(X,N) ralloczero((X),(N)*sizeof((X)[0]))

#define allocone(X) alloc((X),sizeof((X)[0]))
#define rallocone(X) ralloc((X),sizeof((X)[0]))
#define alloconezero(X) alloczero((X),sizeof((X)[0]))
#define ralloconezero(X) ralloczero((X),sizeof((X)[0]))

/* added 1/2001 */
#if CHECKHEAP == 2 || CHECKHEAP == -2
#  define copyarray(X,Y,Z) do { \
  if (sizeof((Y)[0]) != sizeof((X)[0])) ERROR(("copyarray: element sizes disagree")) \
  copy((X),(Y),(Z)*sizeof((Y)[0])); } while (0)
#else /*# CHECKHEAP == 2 || CHECKHEAP == -2 */
#  define copyarray(X,Y,Z) copy((X),(Y),(Z)*sizeof((Y)[0]))
#endif /*#!CHECKHEAP == 2 || CHECKHEAP == -2 */

/* added 12/2003 */
#define arrayzero(X,Y) memset((X),0,(Y)*sizeof((X)[0]))
/* WARNING: this is not ANSI C portable, although I do not know about
   architecture where this trick would fail */

#include "rndgen.h"
#  define CALC 0
#include "calc.h"

/* Replacement of the old good simple deprecated gets(), see mygets.c */
extern int getsbufsize;
char *mygets(char *s);
#define gets mygets
