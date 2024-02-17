/*
  Arbitrary precision arithmetic -- (c) J.Kolafa 96

  Properties:
    ideal relative precision=1/2^(32*PRECISION-8)
    range 1e-38456..1e+38456 (for MAXEXP=16384)

  Compilation:
    highprec.cc and all modules #including highprec.h must be compiled
      by C++ with -DPRECISION=#, where #>2
    deprecated: -DHIGHPREC=# where #=2 is allowed here
      in other modules, -DPRECISION=2 means long double version
    highprec.h MUST be #included after C-style modules ground.h and alloc.h
      because it overrides some of the ground.h macros like Sqr by functions
    to check # of iterations for 1/r and sqrt(r), see setIT() and #define
      SETIT in highprec.cc

  Usage:
    high precision real numbers should be declared as e.g.
      REAL x,y[10];
    and used in the same way as double, e.g.:
      x=y[0]+sqrt(y[1]+0.5);
    available:
      - all arithmetic operators incl. conversions from double and int
      - all comparisons
      - functions sqrt,log,exp,pow,powi

  Restrictions:
    functions dVal and Int must be used for REAL-->double and REAL-->int casts:
    (macro Val(X) = (double)(X) should also exist)
      REAL res;
      int i=Int(res);
      printf("res=%g int(res)=\n",Val(res),i);
    no C++ style iostream i/o available (cannot use cout<<res;)
    trigonometric functions are not implemnented
    sqrt is tailored for PRECISION<10 and becomes inefficient for PRECISION>>10
    initializing static structures containg REAL members is not allowed
      (stupid C++ restriction), just static REAL r=INITIAL_VALUE is OK

  Implementation requirements:
    - int and unsigned int must be 32 bit long
    - 64 bit unsigned int must be available: see the definition of
      unsigned64 in highprec.cc!
    - big endian (SEX==0) version is the basis
    - small endian (SEX==1) solved by changing endianess when necessary
    - number of iterations needed for 1/r and sqrt(r) is determined from
      double accuracy 2^-53 (double result is used as seed): if not true,
      check function setIT()
    - gcc -O1 required (cannot be optimized, ?)
*/

#if defined(__cplusplus) && !defined(HIGHPRECINCLUDED)

#  define HIGHPRECINCLUDED

#  define REAL Real

#  define SEX 1
// 0=big endian, as SGI; type `__uint64_t' assumed to be 64 bit unsigned
// 1=small endian, as x86; type `long long unsigned' assumed to be 64 bit unsigned

#  ifndef HIGHPREC
#    ifdef PRECISION
#      define HIGHPREC PRECISION
#    else /*# PRECISION */
#      define HIGHPREC 1
#    endif /*#!PRECISION */
#  endif /*# HIGHPREC */

#  if HIGHPREC==1
#    define PREC 5 // default
#  else /*# HIGHPREC==1 */
#    define PREC HIGHPREC
#  endif /*#!HIGHPREC==1 */

#  if PREC<2
#    error HIGHPREC must be at least 2
#  endif /*# PREC<2 */

typedef unsigned int unsigned4;

/* undefine macro sign() */
#  undef sign

class Real {
//private:
public:
/* 
   NOTE: because of problems with "friend Val", which stopped working 
   and had to be moved outside class Real, all data were made public...
   This is ugly, but I really don't know why the compiler no longer accepts
   this construct...
*/   
  int e;              // exponent in powers of 256 (i.e. 256^e)
  unsigned m[PREC];   // mantissa in bin complement, dec.point assumed left

  friend int compare(Real a,Real b);
  friend int iszero(Real a);
  friend Real oneover(Real r);

//public:
  Real(int,unsigned[]); // exponent, mantissa
  Real() { int i; for (i=0; i<PREC; i++) m[i]=0; e=0; }
  Real(double);
  Real(int);
  Real(unsigned);
  Real(const char*);
  
  int RealExponent() { return e; }
  unsigned RealMantissa(int i) { return m[i]; }

  friend int Int(const Real&);

  friend Real operator+(Real, Real);
  friend Real operator-(Real);
  friend Real operator-(Real, Real);
  friend Real operator*(Real, Real);
  friend Real operator*(int, Real);
  friend Real operator*(Real, int);
  friend inline Real operator*(double d, Real r) { return (Real)d*r; }
  friend inline Real operator*(Real r, double d) { return r*(Real)d; }
  friend Real operator/(Real, Real);
  friend Real operator/(Real, int);
  friend inline Real operator/(Real r, double d) { return r/(Real)d; }

  void operator+=(Real);
  void operator-=(Real);
  Real operator*=(Real);
/*.....  void operator*=(unsigned);*/
  void operator*=(int);
  void operator*=(double d) { *this=*this*(Real)d; }
  void operator/=(Real b) { *this=*this/b; }
/*.....  void operator/=(unsigned);*/
  void operator/=(int);
  void operator/=(double d) { *this=*this/(Real)d; }

  friend int operator==(Real, Real);
  friend int operator!=(Real, Real);
  friend int operator>=(Real, Real);
  friend int operator>(Real, Real);
  friend int operator<=(Real, Real);
  friend int operator<(Real, Real);

  friend Real abs(Real);
  friend Real fabs(Real);
  friend int sign(Real);
  friend Real sqrt(Real);
  friend Real log(Real);
  friend Real exp(Real);

  friend Real cos(Real);
  friend Real sin(Real);
  friend Real tan(Real);
  friend Real acos(Real);
  friend Real asin(Real);
  friend Real atan(Real);

  friend Real tanh(Real);
  friend Real sinh(Real);
  friend Real cosh(Real);
  friend Real atanh(Real);
  friend Real asinh(Real);
  friend Real acosh(Real);

  friend Real pow(Real,Real);
  friend Real powi(Real,int);
  friend Real sqr(Real);
  friend Real cub(Real);
  friend Real pow4(Real);
  friend Real pow5(Real);
  friend Real pow6(Real);
  friend inline Real fmax(Real x,Real y) { if (x<y) return y; else return x; }
  friend inline Real fmin(Real x,Real y) { if (x<y) return x; else return y; }

  friend void RPut(Real);
        
};

double Val(const Real&);
double Val(const double&); // these are required by getdata
double Val(const int&);
double Val(const unsigned&);

#  define Put(R) { printf("%s = ",#R); RPut(R); }

#  ifdef Sqr
#    undef Sqr
#    define Sqr sqr
#  endif /*# Sqr */

#  ifdef Cub
#    undef Cub
#    define Cub cub
#  endif /*# Cub */

#  ifdef Pow4
#    undef Pow4
#    define Pow4 pow4
#  endif /*# Pow4 */

#  ifdef Pow5
#    undef Pow5
#    define Pow5 pow5
#  endif /*# Pow5 */

#  ifdef Pow6
#    undef Pow6
#    define Pow6 pow6
#  endif /*# Pow6 */

/* legacy cubrt eventually removed (7/2022), use cbrt = cube root by C99, POSIX.1-2001  */

#  ifdef cbrt
#    undef cbrt
#    define cbrt(X) pow((X),(REAL)(1)/3)
#  endif /*# cbrt */

#  ifdef Sign
#    undef Sign
#    define Sign sign
#  endif /*# Sign */

#  ifdef PI
#    undef PI
#  endif /*# PI */

extern const Real PI;

#endif /*# defined(__cplusplus) && !defined(HIGHPRECINCLUDED) */

