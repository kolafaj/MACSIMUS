/********** Calculator ********** (c) J.Kolafa 1991--NOW **********

* INTEGRATED to ground.c - see gen/ground.h, gen/ground.c
* USAGE without ground: see c/ev.C, c/rk.c, c/cutbin.c
* SIMPLIFIED:
  - the richest version is the default
  - the CALC switch has only two values (see below)
  - MACSIMUS rndgen (Ziff register-shift) used

#define CALC:
  0: only numbers
  1: units of measurements supported (cf. evu.c)

Synopsis:
  double Calc(char *buf,char **endscan)
  struct unit_s Calc(char *buf,char **endscan)

Parsing:
  Parses  buf , the end of expression is recognized by any character less
  than ' ', or a character that would cause a syntax error (provided that
  the previous part of the expression is correct).

  The result is returned.
  *endscan = a char after the end of the expresion parsed

Syntax:
  Numbers are in the usual format, also FORTRAN-like supported (1d1=1e1)
  Binary operators: 
    + - * / ` ^ $
    (` is modulo, ^ is power, $ is xor)
  Unary operators and functions: 
    + - \ sin cos tan asin acos atan exp 
    not(bitwise) `(logical not)
    ln log lb sqrt abs int fac
    (\ is sqrt, fac or is factorial [fac 3=1*2*3, fac 2.5=0.5*1.5*2.5],
    ln=natural log, log=decadic log)
  Parentheses: ( )
  PI or pi: 3.14159265358979323846
  Time format: [DD::]hh[:mm[:ss]]
    (the result is on hours, e.g., 2::1=49, 11:45=11.75, 0:0:36=0.01)

Precedence:
  numbers and characters (incl. : or :: in time format)
  ()
  functions and unary \ (\=sqrt)
  ^ ** (power, associative from left as * or /)
  + - ` (unary, `=logical not)
  * / `  (`=modulo, % removed)
  + - (binary)
  < <= > >= <> == (relation operators)
  & $ ($=xor)
  | (|=or)

Units (CALC&1):
  units are given after number in [], e.g.: 3[m3.mol-1]
  "more common" units are search first:
     Ma = Mega annum = 1e6 years
     Pa = Pascal (not Peta annum), because this is looked for first

Examples:
  2e3       = 2000
  'a'       = 97 (ascii value)
  2^-1/2    = (2^(-1))/2
  -2^2      = -(2^2)
  3^3^3     = (3^3)^3
  cos pi^2  = (cos (pi))^2
  cos -pi^2 = cos (-(pi^2))
  double Calc(char *buf,char **endscan)

Bugs:
  3^3^3  is not  3^(3^3)  but  (3^3)^3 (matter of taste, I guess)
  cos pi^2  is something else than  cos +pi^2
    I do not see any consistent way how to interpret such cases: use
    parentheses if you are in doubt

History:
  2025: modulo changed back to `, % removed (NB: unary ` = logical not)
        note that @ in plot and tabproc is line number
  2023: ` changed to @ (modulo), % still valid (not in ev)
  2021: sign() changed to sgn(), H() added
  2018: cbrt added; nul, null removed (replacement: (x==0))
  2017: Ctrl-D to exit, bug fixed for bad unit of a character
  2016: units version merged, CALC switch removed (as if CALC=4)
  2015: units finished (fork from calc.c)
  2014: code cleaned, rnd/rand switch, units added, DOS support removed
  2008: CALC=4 needed for some special functions (erfc, gamma, hyp.)
        new version CALC=4 is the same as old CALC=3
        (and CALC=4 becomes MACSIMUS default)
  2007: _ allowed in identifiers
  2006: binary operators ** (=^)  == >= <= <> added, floor() added
  2005: [day+]time format [DD::]hh:mm:ss (result in h)
        alternative modulo operator : changed to ` (% still OK)
        sign() added
        null,nul added (null(0)=1,..; zero() interacts with tabproc x,y,z)
  2004: operators < and > added, precedence between bitwise and +-
        (1<2 returns 1, 1>2 returns 0)
  2002: hyperbolic functions added (sh=sinh, ash=asinh etc.)
        logarithms updated (ln,log,lb)
        bitwise functions added ( & | $=xor ), not(X), bit(X) (# of bits)
        (all in CALC>2)
  2001: TABS allowed where ' '
  1997: : is modulo (the same as %)
**********************************************************************/

#ifndef CALC
#  error "CALC #undefined, use 0 for numbers only, 1 for numbers+units"
#endif /*# CALC */

#define OPSTACKLEN 128 /* operand stack */
#define NRSTACKLEN 64 /* number, or pending calculations, stack */

#include "loop.h"

static int precedence(char op) /********************************* precedence */
{
 switch (op&127) {
   case 0: return 0; /* evaluation operator */
   case '=': return 1;
   case ')': return 2;
   case '(': return 3;
   case '|': return 5;
     /* or        xor */
   case '&': case '$': return 6;
     /*                     ==        >=        <=        <> */
   case '<': case '>': case 037: case 036: case 035: case 034: return 7;
     /*                     not */  
   case '+': case '-': return op&128 ? 12 : 8; /* distinguish binary/unary */
   case '`': return op&128 ? 10 : 8; /* distinguish binary/unary */
   case '*': case '/': return 10;
   case '^': return 14;
   default: return 127; /* unary operator */ }
} /* precedence */

#if CALC&1
#  include "calcu.c"
#else /*# CALC&1 */
CALC_t retzero;
#endif /*#!CALC&1 */

CALC_t Calc(char *buf,char **endscan) /******************************** Calc */
{
  unsigned char *b,*end;
  unsigned char *op,Op[OPSTACKLEN]; /* operation stack, bit 7 set means unary operator */
#if !(CALC&1)
  double *nr,Nr[NRSTACKLEN]; /* number stack */
#endif /*# !(CALC&1) */

  b=(unsigned char*)buf;

  op=Op; *op='=';
  nr=Nr-1;

NRSCAN: /*** number, '(' or unary operator expected ***/

  if (op>Op+OPSTACKLEN-2) goto reterror;

  while (*b==' ' || *b=='\t') b++;

  if (*b==';' || *b<' ') goto reterror;

  /* # = register (old value) */
  if (*b=='#') {
    nr++; if (nr>=Nr+NRSTACKLEN) goto reterror;
    *nr=_Id.nr; b++;
    goto OPSCAN; }

  /* identifier */
  if (isalpha(*b) || *b=='_') {
    unsigned char *x=b;
    char id[64],*idp=id;
    _Idlist *l;
    char *funclist[]={"sin","cos","tan",
		      "asin","acos","atan",
		      "sinh","cosh","tanh","sh","ch","th",
		      "asinh","acosh","atanh","ash","ach","ath",
		      "erf","erfc","gamma","lgamma",
		      "not","bit", /* bitwise not, sum of bits */
		      "ln","log","lb",
		      "exp", "sqrt", "cbrt",
                      "fac","fact","rnd",
		      "pi","PI",
                      "val","unit",
		      "abs","sgn","step","int","float","floor","frac",NULL};
    /* NOTE: the Heaviside function is step(), not H(), because plot and
       tabproc interpret uppercase letters as column numbers */
    char **fptr;

    /* abc[12] allowed as identifier, abc[+2] or abc[-1] not! */
    while (isalnum(*x) || (*x && strchr("_.[]",*x))) *idp++=*x++;
    *idp=0;

    /* try to find the identifier in the list */
    for (l=_Id.head; l; l=l->next)
      if (!strcmp(l->id,id)) {
	nr++; if (nr>=Nr+NRSTACKLEN) goto reterror;
        *nr=l->val; b=x;
        l->used=1; /* needed by plot */
        goto OPSCAN; }

    /* identifier not found => will be a function, CHECK now */
    for (fptr=funclist; *fptr; fptr++) if (!strcmp(*fptr,id)) goto fnfound;

    if (_Id.ignoreid) {
      /* identifier not found (not a function), :=0 */
      nr++; if (nr>=Nr+NRSTACKLEN) goto reterror;
#if CALC&1
      nr->val=0; opunits(nr,nr,'0');
#else /*# CALC&1 */
      *nr=0;
#endif /*#!CALC&1 */
      b=x;
      goto OPSCAN; }

    goto reterror;

    fnfound:; }

  /* PI (in the list of valid keywords only pi and PI are allowed) */

  if (toupper(*b)=='P' && toupper(b[1])=='I'
   && !isdigit(b[2]) && !isalpha(b[2]) && b[2]!='_') {
    nr++; if (nr>=Nr+NRSTACKLEN) goto reterror;
#if CALC&1
    nr->val=PI; opunits(nr,nr,'0');
#else /*# CALC&1 */
    *nr=PI;
#endif /*#!CALC&1 */
    b+=2;
    goto OPSCAN; }

  /* parse a number */
  if (isdigit(*b) || *b=='.') {
    nr++; if (nr>=Nr+NRSTACKLEN) goto reterror;
#if 1
    /* numbers in the form 077=0o77 0xff 0b11
       turning this off will interpret 077=77 decadic, 0xff 0b11 are error */
    if (*b=='0' && strchr("0123456789xbo",b[1])) {
      int base=8;

      end=b+1;
#  if CALC&1
      nr->val=0; opunits(nr,nr,'0');
#  else /*# CALC&1 */
      *nr=0;
#  endif /*#!CALC&1 */
      switch (*end) {
        case 'x': base=16; goto plus;
        case 'b': base=2;
        case 'o':
plus:             end++; }

      for (;;end++) {
        int i=*end-'0';
        if (base==16) {
          if (i>='A'-'0' && i<='F'-'0') i-='A'-'0'-10;
          else if (i>='a'-'0' && i<='f'-'0') i-='a'-'0'-10; }
        if (i<0 || i>=base) break;
#  if CALC&1
        nr->val=nr->val*base+i;
#  else /*# CALC&1 */
        *nr=*nr*base+i;
#  endif /*#!CALC&1 */
      }
    }
    else
#endif /*# 1  */
      {
//	double x=strtod(b,&end);
//	double x=strtod((char*)b,(char**)&end);
/* fstrtod: accepts also FORTRAN exponent with d,D; e.g., 3.13D0 */
	double x=fstrtod((char*)b,(char**)&end);
	double xx=0;
	double q=1;

	/* time format hh:mm:ss and day+time format DD::hh:mm:ss
	   the result is always in hours */
	if (*end==':' && b!=end && end[1]==':') {
	  b=end+2;
	  x=x*24+strtod((char*)b,(char**)&end); }

	while (*end==':' && b!=end) {
	  b=end+1;
	  q*=60;
	  xx+=strtod((char*)b,(char**)&end)/q; }
#if CALC&1
	nr->val=x+xx;
#else /*# CALC&1 */
	*nr=x+xx;
#endif /*#!CALC&1 */
      }
      //      *nr=strtod((char*)b,(char**)&end);

    if (b==end) goto reterror;
    b=end; /* not in !CALC&1 */

#if CALC&1
    opunits(nr,nr,'0');
    while (*b && *b<=' ') b++;
    if (*b=='[') parseunit(nr,(char*)b+1,(char**)&end,1);
    //    fprintf(stderr,"%g*%g m^%g g^%g s^%g K^%g mol^%g A^%g\n",nr->val,nr->pow[0],nr->pow[1],nr->pow[2],nr->pow[3],nr->pow[4],nr->pow[5]);
    b=end;
#endif /*# CALC&1 */

    goto OPSCAN; }

  /* parse a character (= special form of a number, unit not accepted) */
  if (*b=='\'') {
    nr++; if (nr>=Nr+NRSTACKLEN) goto reterror;
    b++;
#if CALC&1
    if (*b=='\'') nr->val=0; else nr->val=*b++;
    opunits(nr,nr,'0');
#else /*# CALC&1 */
    if (*b=='\'') *nr=0; else *nr=*b++;
#endif /*#!CALC&1 */
    if (*b++!='\'') goto reterror;
    goto OPSCAN; }

  /* left parenthesis */
  if (*b == '(') {
    *(++op) = *b++; goto NRSCAN; }

  /* simple parsing of a function
     op-code is the 1st letter or upcase, 2nd letter for inv.func.,
     extended to hyp.func.
     (this is cumbersome, but I'm too lazy to rewrite it...) */
  if (isalpha(*b)) {
    ++op;
    *op = *b=='a' ? toupper(b[1]) : *b ; /* sin,cos,tan + inv.; also abs,exp */
    if (*b=='s' && b[1]=='q') *op='\\'; /* sqrt */
    if (*b=='c' && b[1]=='b') *op='u'; /* cbrt */
    if (*b=='s' && b[1]=='g') *op='I'; /* sgn */
    if (*b=='s' && b[1]=='t') *op='H'; /* Heaviside function step */
    /* logarithms: */
    if (*b=='l') switch (b[1]) {
      case 'g': *op='G'; break; /* lgamma */
      case 'o': *op='d'; /* log=log_10, formerly dlog */
      case 'n': break;   /* ln */
      case 'b': *op='w'; /* lb=log_2 */ }
    if (!memcmp(b,"float",5)) *op='A';
    if (!memcmp(b,"floor",5)) *op='F';
    if (!memcmp(b,"frac",4)) *op='Q';
    if (!memcmp(b,"erf",3)) *op=b[3]=='c'?'J':'j'; /* erf[c] */
    /* hyperbolic functions: (ANY function name with 'h' in) */
    if (!memcmp(b,"sh",2) || !memcmp(b,"sinh",4)) *op='x';
    if (!memcmp(b,"ch",2) || !memcmp(b,"cosh",4)) *op='y';
    if (!memcmp(b,"th",2) || !memcmp(b,"tanh",4)) *op='z';
    if (!memcmp(b,"ash",3) || !memcmp(b,"asinh",5)) *op='X';
    if (!memcmp(b,"ach",3) || !memcmp(b,"acosh",5)) *op='Y';
    if (!memcmp(b,"ath",3) || !memcmp(b,"atanh",5)) *op='Z';
    if (!memcmp(b,"val",3)) *op='v';
    if (!memcmp(b,"unit",4)) *op='V';
/*.....    fprintf(stderr,"<%s>%c\n",b,*op);*/
    *op |= 128;
    while (isalpha(*b)||*b=='_') b++;
    goto NRSCAN; }

  /* unary one-char operators like + - \ (\=sqrt) */
  *(++op) = 128 | *b++;

  goto NRSCAN;

OPSCAN: /*** binary operator or ')' or end expected ***/

  while (*b==' ' || *b=='\t') b++;
  *(++op) = *b++;

  /* multicharacter operators: ** ==  >=  <=  <>
                               ^  31  30  29  28  (dec)
                                  037 036 035 034 (oct) */
  if (*op=='*' && *b=='*') *op='^',b++;
  if (*op=='=' && *b=='=') *op=037,b++;
  if (*op=='>' && *b=='=') *op=036,b++;
  if (*op=='<' && *b=='=') *op=035,b++;
  if (*op=='<' && *b=='>') *op=034,b++;

  /* check operators */
  if (!strchr(")+-*/`^|$&<>\037\036\035\034",*op)) *op=0;

EVALSTACK: /*** stacked operations are tried to evaluate ***/

  if (op[-1]=='=' && *op==0) {
    *endscan=(char*)b-1;
    if (nr==Nr) return *nr;
    goto reterror; }

  if (*op==')' && op[-1]=='(') { op-=2; goto OPSCAN; }

  if (precedence(*op) > precedence(op[-1])) goto NRSCAN;

  op--;

#if CALC&1

#  define OP(X,C) opunits(nr,X,C)
#  define NR0 nr->val
#  define NR1 nr[1].val
#  include "calc1.c"

#else /*# CALC&1 */

#  define OP(X,C) /**/
#  define NR0 *nr
#  define NR1 nr[1]
#  include "calc1.c"

#endif /*#!CALC&1 */

  *op=op[1];

  goto EVALSTACK;

 reterror: *endscan=buf;

  return retzero;
}
