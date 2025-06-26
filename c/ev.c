/* makeev.sh

  !!! DO NOT EDIT ev.c !!!
  * edit ev.C and run C2c ev.C
  * or run makeev.sh which contains:
    C2c ev.C
    cc -o ev -O2 -Wall -DCALC=0 ev.c -I../gen -lm -lncurses
    cc -o evu -O2 -Wall -DCALC=1 ev.c -I../gen -lm -lncurses

History:
06/2024:
  character '!' now recognized (! was treated as comment even in '!')
  ` as 'logical not' was added (bitwise not is function not())
03/2024:
  m(FORMULA) added
12/2023:
  in calc.c, ` changed to @ (alt for modulo)
07/2023:
  #include <unistd.h> added (chdir())
  bug crash after "to unknown_unit" fixed
06/2022:
  added: commands cd and exe (call shell and return the result in #)
06/2020:
  empty unit [] = previous unit
03/2020:
  small fixes, hlp updated
06/2019:
  verbosity changed, now --v0 => nothing printed (useful with --r)
  --- added
02/2019:
  small fix: wrong preferred unit (command to) is (after message) removed
09/2017:
  help (evhlp.c), small fixes
04/2016:
  ev.c and evu.c unified by switch CALC, bugs fixed
07/2015:
  units added (fork from ev.c)
08/2014:
  --r added
01/2012:
  clipboard added (command clip)
11/2010:
  trailing spaces silently removed on input
4/2010:
  environment variable EV will replace 1st char to 2nd char
    (e.g., modulo is `, to use : for it, use export EV=":\`")
    NB: now modulo=@
06/2009:
  environment SCROLL removed - COLUMNS used instead
11/2007:
  support for a=b=c and a+=b
4/2007:
  code cleaned
  versions unified, DOS no longer supported
  options changed
  .evdata automatically read
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include "mygets.c"

#if CALC&1 
	 #define EV "evu"
#else 
	 #define EV "ev"
#endif

volatile int sig=0;
#include <signal.h>

int (*myprintf)(const char *format,...);
char *(*evgetline)(char *line);

#define Sqr(X) ((X)*(X))
#ifndef PI
#define PI 3.14159265358979323846
#endif

#define STDERR
#include "qmin.c"
#include "ms.h"

enum mode_e { ARGS, DUMB, PROMPT, FILTER, TEX } mode;

#include "prompt.c"

int cprintf(const char *format,...) /******************************* cprintf */
{
  int ret;
  char L[PROMPT_N*5],*tok;
  va_list args;

  va_start(args,format);
  ret=vsprintf(L,format,args);
  va_end(args);
  tok=strtok(L,"\n");
  while (tok) {
    _pr.point=_pr.pos=0;
    _pr.back=0; prompt_showline(tok); _pr.back=1;
    tok=strtok(NULL,"\n");
    addch('\n'); }
  refresh();

  return ret;
}

char *strlast(const char *str) /************************************ strlast */
/* returns pointer to the last char, or NULL if str==NULL or str="" */
{
  char *c;

  if (!str || !*str) return NULL;
  for (c=(char*)str; *c; c++);

  return c-1;
}

char *comment(char *c) /******************************************** comment */
/*
   find first ! which is comment (i.e., not character '!')
   and put 0 there or to the first consecutive space before
*/
{
  char *spc=NULL;

  if (!c) return 0; /* should not happen */
  /* line starting by ! */
  if (*c=='!') { *c=0; return c; }
  if (*c==' ') spc=c;

  while (*++c) {
    if (*c=='!' && !(c[1]=='\'' && c[-1]=='\'') ) {
      /* this ! is comment */
      if (spc) *spc=0; else *c=0; /* erase from ! or preceding space */
      return c; }
    if (*c!=' ') spc=NULL;
    else if (!spc) spc=c; }

  return NULL;
}

int texprintf(const char *format,...) /*************************** texprintf */
{
  int ret,res=0,nc=-2,inmac=0;
  char L[PROMPT_N*5],*c;
  va_list args;

  va_start(args,format);
  ret=vsprintf(L,format,args);
  va_end(args);
  for (c=L; *c; c++) {
    if (!inmac) if ((nc>44 && *c==' ') || (nc>47 && !isalnum(*c)) || nc>=50) {
      printf("\\\\\\hspace*{1em}");
      nc=0; }
    nc++;
    if (*c=='^') printf("$\\uparrow$");
    else if (*c=='_') printf("\\_");
    else if (*c=='^') printf("\\chr94 ");
    else if (*c==' ') printf("~");
    else if (*c=='{') printf("{}[");
    else if (*c=='}') printf("]");
    else if (*c=='[') {
      if (res==3) { printf("} "); res=1; }
      printf("{[\\Unit{\\rm "); inmac=1; }
    else if (*c==']') {
      printf("}]}"); inmac=0; }
    else if (*c==176) printf("$^\\circ$"); // iso-8859-2 ° 
    else if (*c=='\\' && (c>=L && c[-1]!='\\') && c[1]!='\\') printf("$\\surd$");
    else if (*c=='@') { res=1; printf("\\tt\\mbox{}\\hfill$"); } // ?
    else if (res && *c=='e' && !inmac) { res=3; printf("\\times10^{"); }
    else if (*c=='%') printf("\\%%");
    else if (*c=='.' && inmac) printf(",");
    else if (*c=='~') printf("$\\sim$");
    else printf("%c",*c); }
  if (res) {
    if (res==3) printf("}");
    printf("$"); }
  printf("\\\\\n");

  return ret;
}

int dummyprintf(const char *format,...) /*********************** dummyprintf */
{
  return 0;
}

void asksig(int signr) /********************************************* asksig */
{
  if (sig++) { /* ^C typed twice */
    if (mode==PROMPT) endwin();
    exit(-1); }
  (*signal)(SIGINT,asksig);
  myprintf("\nINTERRUPT(%d) ACCEPTED (type ^C again to exit immediately)\n",
          signr);
}

static int err,inputmod;

#if CALC&1
#  define VAL .val
#  define OP3(A,B,C) opunits(A,B,C)
static int checkto=0;
#else
#  define VAL /**/
#  define OP3(A,B,C) /**/
#endif

/* CALC_t defined here: */
#include "calc.h"

static CALC_t res;
#define NFMT 126
static char fmt[NFMT+2];

typedef struct _Idlist_s {
  struct _Idlist_s *next;
  int isconst; /* constant (cannot be changed); pi PI solved directly */
  CALC_t val;  /* value of the variable */
  int used;    /* set to 1 if read (used), needed for plot */
  char id[1];  /* [var len], identifier name */
} _Idlist; /* cf. ground.h */

/* list of variables and calc interface */
struct {
  struct _Idlist_s *head;
  double degrad;
  int ignoreid; /* 0: unknown ID = error, 1: unknown ID = silent 0 */
  CALC_t nr;
} _Id = { NULL,1 };

static int isconst;
char assignop='=';

void assign(char *b) /*********************************************** assign */
/*
  Assignment (operators =,:=) to global register res from variable named b,
  or modification (+=,-=,*=,/=,@=) of global register res by variable b
*/
{
  struct _Idlist_s *is;

  if (!strcmp(b,"pi") || !strcmp(b,"PI")) {
    myprintf("constant %s cannot be changed",b);
    return; }

  for (is=_Id.head; is; is=is->next)
    if (!strcmp(is->id,b)) {
      if (is->isconst) {
        myprintf("constant %s cannot be changed (~%s to remove)",b,b);
        return; }
      goto val; }
  is=malloc(sizeof(struct _Idlist_s)+strlen(b));
  if (!is) { myprintf("NO HEAP"); return; }
  is->next=_Id.head;
  _Id.head=is;
  strcpy(is->id,b);
  is->isconst=isconst;
 val:
  switch (assignop) {
    case '+': OP3(&res,&is->val,assignop); res VAL=is->val VAL+=res VAL; break;
    case '-': OP3(&res,&is->val,assignop); res VAL=is->val VAL-=res VAL; break;
    case '*': OP3(&res,&is->val,assignop); res VAL=is->val VAL*=res VAL; break;
    case '/': OP3(&res,&is->val,assignop); res VAL=is->val VAL/=res VAL; break;
    case '`': OP3(&res,&is->val,assignop); res VAL=is->val VAL-(int)(is->val VAL/res VAL)*res VAL; break;
    case ':':
    case '=': is->val=res; break; }

  assignop='=';
}

void deassign(char *b,char *msg) /********************************* deassign */
/* remove variable b from the list */
{
  struct _Idlist_s *is,*last=NULL;

  for (is=_Id.head; is; is=is->next) {
    if (!strcmp(is->id,b)) {
      if (last) last->next=is->next;
      else _Id.head=is->next;
      myprintf(" (%s)\n",msg); // re-added because of tevu.sh
      free(is); return; }
    last=is; }
  myprintf("variable %s not found\n",b);
}

int prompt=1;

#include "fortran.h"
#include "rndgeni.c"
#include "calc.c"

#include "powi.c"
#include "fortran.c"

static enum { DOUBLE,INT,UNSIGNED } cast;
static char *toclip; /* getenv("TOCLIP") */

static
void clip(char *s) /* ------------------------------------------------- clip */
/* copy to clipboard */
{
  if (toclip) {
    char ll[1024],*c;

    for (c=s; c[1]; c++);
    if (*c=='\n') *c=0;

    if (toclip[0]=='-')
      sprintf(ll,"echo \'%s\' | %s",s,toclip+1);
    else
      sprintf(ll,"%s \'%s\'",toclip,s);

    if (system(ll))
      myprintf("xclip: system call failed\n"); }
}

/* verbosity level, sum of flags:
   1 = print result
   2 = print f=(last value) after solve
   4 = print expanded macros
   8 = print initial constants (.evdata/.evudata)
   16= print comments (FILTER mode only)
   before 2/2019: 1: solve:f=err 2:macros 4:init.const 8:comments (FILTER only)
*/
static int verbose=1+2+16;

static
void printres(void) /* -------------------------------------------- printres */
/* print result (global res) */
{
  char s[1024];
#if CALC&1
  char *u=unit(&res);
  double q=res VAL/unitfactor;
  int donotprint=0;

  if (checkto) {
    char *c1,*c2;
  
    if (lastunit) {
      for (c1=lastunit,c2=u; *c1 && *c2; c1++,c2++) {
        if (strchr("[]",*c2)) c2++;
        if (*c1!=*c2) { donotprint++; goto fin; } }

    fin:; } }
      //  fprintf(stderr,"•%s•%s•\n",lastunit,u);

  switch (cast) {
    case DOUBLE: sprintf(s,fmt,q,u); break;
    case INT: sprintf(s,fmt,(int)q,u); break;
    case UNSIGNED: sprintf(s,fmt,(unsigned)q,u); break; }
#else
  switch (cast) {
    case DOUBLE: sprintf(s,fmt,res); break;
    case INT: sprintf(s,fmt,(int)res); break;
    case UNSIGNED: sprintf(s,fmt,(unsigned)res); break; }
#endif

#if CALC&1
  if (donotprint)
    /* do not print the value if another unit set by 'to' */
    myprintf(" unit set\n");
  else
    /* re-print value if 'to' applies to the unit of the last result */
#endif
    if (verbose&1) myprintf(s);
  clip(s);

  _Id.nr=res;
}

static
void printres2(char *fmt2,double var2) /* ------------------------ printres2 */
/* print result (global res and extra info, as error) */
{
  char s[1024],f[128];
  char *x;
#if CALC&1
  char *u=unit(&res);
  double q=res VAL/unitfactor;

  sprintf(f,"%s%s",fmt,fmt2);
  if ( (x=strchr(f,'\n')) ) *x=' '; /* why this? */
  strcat(f,"\n");
  switch (cast) {
    case DOUBLE: sprintf(s,f,q,var2,u); break;
    case INT: sprintf(s,f,(int)q,var2,u); break;
    case UNSIGNED: sprintf(s,f,(unsigned)q,var2,u); break; }

  myprintf(s);

  /* clip: only main result, not info */
  switch (cast) {
    case DOUBLE: sprintf(s,fmt,q,u); break;
    case INT: sprintf(s,fmt,(int)q,u); break;
    case UNSIGNED: sprintf(s,fmt,(unsigned)q,u); break; }
#else
  sprintf(f,"%s%s",fmt,fmt2);
  if ( (x=strchr(f,'\n')) ) *x=' '; /* why this? */
  strcat(f,"\n");
  switch (cast) {
    case DOUBLE: sprintf(s,f,res,var2); break;
    case INT: sprintf(s,f,(int)res,var2); break;
    case UNSIGNED: sprintf(s,f,(unsigned)res,var2); break; }

  myprintf(s);

  /* clip: only main result, not info */
  switch (cast) {
    case DOUBLE: sprintf(s,fmt,res); break;
    case INT: sprintf(s,fmt,(int)res); break;
    case UNSIGNED: sprintf(s,fmt,(unsigned)res); break; }
#endif
  clip(s);

  _Id.nr=res;
}

static int nop;
static CALC_t parm[4]; /* parm value ([4] is by 1 more) */
static int isparm[4]; /* if parm assigned */
static char *var,*eqn;

enum op_e { NONE,
            SOLVE,MIN,MAX,PLOT,
            INTEG,DERIV,DIFF,SUM,PROD,ITER,REPEAT,CONTFRAC, /* see below */
            DEF,EXPAND, /* see DIRTY HACK */
            UNDEF,
#if CALC&1
	TO,
#endif
            HELP /* just to mark that help has been printed */
};

static
struct mac_s {
  struct mac_s* next;
  char *X;
  char *E;
} *mhead=NULL;

static
void defmacro(char *X,char *E) /*********************************** defmacro */
/* define macro X as E (no expansion performed now) */
{
  struct mac_s *m;

  for (m=mhead; m; m=m->next) if (!strcmp(m->X,X)) goto repl;
    m=malloc(sizeof(*m));
    m->E=NULL;
    m->X=NULL;
    m->next=mhead;
    mhead=m;

 repl:
  if (m->E) free(m->E);
  if (m->X) free(m->X);
  m->X=strdup(X);
  m->E=strdup(E);
}

static
void undefmacro(char *X) /* ------------------------------------ undefmacro */
/* undefine macro X */
{
  struct mac_s **mp,*m;

  for (mp=&mhead,m=*mp; *mp; mp=&((*mp)->next),m=*mp) if (!strcmp(m->X,X)) {
    if (m->E) free(m->E);
    if (m->X) free(m->X);
    *mp=m->next;
    free(m);
    return; }
  myprintf("%s not defined",X);
}

static
void expmacro(char *expr)  /* ------------------------------------ expmacro */
/* expand expression expr using macros */
{
  char *c,*b;
  int nrepl=0;
  struct mac_s *m;
  int l,chg;

  if (!mhead) return;
  do {
    chg=0;
    for (c=expr; *c; ) {
      if (isalpha(*c) || *c=='_') {
        for (b=c; (*c && isalnum(*c)) || *c=='_'; c++);
        /* now ID starts at b and ends before c */
        for (m=mhead; m; m=m->next) {
          if (c-b==strlen(m->X) && !memcmp(b,m->X,strlen(m->X))) {
            l=strlen(expr)+strlen(m->E)-strlen(m->X)+1;
            if (l>=PROMPT_N) {
              myprintf("too long expression after macro expansion");
              return; }
            if (nrepl++>10000) {
              myprintf("more than 10000 macro expansions (recursion?)");
              return; }
            memmove(b+strlen(m->E),b+strlen(m->X),strlen(b)-strlen(m->X)+1);
            memcpy(b,m->E,strlen(m->E));
            chg=1; /* replacement has been done */
            c=expr; /* rescan string from start; NB: no c++ in 'for (c=expr..)' */
            goto rescan; }
        } /* macro list loop */
      }
      c++;
    rescan:; } /* for (c=expr; *c; ) */
  } while (chg);

  if (nrepl && verbose&4) myprintf("{%s}",expr);

  return;
}

static
void expand(char *X,char *E) /* ------------------------------------- expand */
/* define macro X as E (with all expansions performed now) */
{
  struct mac_s *m;
  char line[PROMPT_N];

  strcpy(line,E);
  expmacro(line);

  for (m=mhead; m; m=m->next) if (!strcmp(m->X,X)) goto repl;
    m=malloc(sizeof(*m));
    m->E=NULL;
    m->X=NULL;
    m->next=mhead;
    mhead=m;

 repl:
  if (m->E) free(m->E);
  if (m->X) free(m->X);
  m->X=strdup(X);
  m->E=strdup(line);
}

static
int oper(char *expr) /* ----------------------------------------------- oper */
/*
  checks expr for an operation (command), as solve,plot, etc.
  these operators are in the form OPER VAR=FROM[,TO,BY|N] EXPR
  exceptions:
    contfrac i=5,7 2++1//i ==> 2+1/(5+1/(6+1/(7)))
*/
{
  char *sp=strchr(expr,' '),*tok,*c,*eq;
  static struct op_s {
    char *name; /* keyword */
    int minparm; /* min # of parameters */
    int maxparm; /* max # of parameters */
    int dimless; /* active if UNITS: dimensioness parameter index, 0 if none
                    for check; also checked: unit(parm[0]==unit(parm[1]) */
    char *help; /* text to print if w/o parms */
    /* notes:
       * if less than maxparm parameters are give, undefined
         parameters adopt "reasonable" defaults, not previuos values
       * if minparm=0 then form "solve x x^2-2" is allowed, where x
         has been assigned before
    */
    } op[]= {
    /* index must correspond to enum op_e above */
    {"",0,0,0,""},
    {"solve",   0,3,0,"X[=INIT] EXPR (secants), X=FROM,TO[,BY] EXPR (search)"},
    {"min",     0,3,0,"X[=INIT] EXPR, X=FROM,TO[,BY] EXPR"},
    {"max",     0,3,0,"X[=INIT] EXPR, X=FROM,TO[,BY] EXPR"},
    {"plot",    2,3,2,"X=FROM,TO[,N] EXPR [N<0: 1/2-shift]"},
    {"integ",   2,3,2,"X=FROM,TO[,N] EXPR"},
    {"deriv",   0,2,0,"X[=VALUE[,D]] EXPR"},
    {"diff",    0,2,0,"X[=VALUE[,D] EXPR"}, /* = deriv */
    {"sum",     2,3,1,"N=FROM,TO[,BY] EXPR"},
    {"prod",    2,3,1,"N=FROM,TO[,BY] EXPR"},
    {"iter",    1,2,0,"X=FROM[,EPS] EXPR (repeat X=EXPR until |X-Xlast|<EPS)"},
    {"repeat",  2,2,1,"X=FROM,COUNT EXPR (COUNT* repeat X=EXPR)"},
    {"contfrac",2,3,1,"N=FROM,TO[,BY] [BEFORE++]NUM//DENOM"},
    {"def",     0,0,0,"X ANY_STRING, def X=A -> def X (A)"}, /* see DIRTY HACK */
    {"expand",  0,0,0,"= def with all macros expanded now"}, /* see DIRTY HACK */
    {"undef",   0,0,0,"X"},
#if CALC&1
	{"to",      0,0,0,"UNIT (preferred output, remove: to ~UNIT / to ~)"},
#endif
    {NULL,0,0,0,""} };
  int i,np;

  if (sp) *sp=0;

  for (nop=1; op[nop].name; nop++)
    if (!strcmp(op[nop].name,expr)) {
      if (!sp) {
        myprintf("%s %s\n",expr,op[nop].help);
        nop=HELP;
        return 1; }
      goto found; }

  if (!sp) return 0;
  *sp=' ';
  return 0;

 found:
  var=tok=strtok(sp+1," \t"); /* x=1,2,3 */

  if (!var) {
    myprintf("variable[=...] expected\n");
    return 0; }
  eq=strchr(tok,'=');
  if (eq) *eq++=0;
  res=Calc(var,&c);

  /* initialization (needed?) */
#if CALC&1
  isparm[0]=0;
  for (i=1; i<4; i++) {
    int j;
    isparm[i]=0;
    parm[i] VAL=0; loop (j,0,CALC_NUNITS) parm[i].pow[j]=0; }
#else
  for (i=0; i<4; i++) { isparm[i]=0; parm[i]=0; }
#endif
  parm[0]=res;

  np=0;

  if (eq) {
    /* = found, scan ,-separated parameters */
    for (np=0; np<4; np++) {
      char *comma=strchr(eq,',');

      if (comma) *comma=0;
      if (np>=op[nop].maxparm) {
        myprintf("superfluous parameter %d",np+1);
        return 0; }
      if (strlen(eq)==0) {
        myprintf("expression expected after ,");
        return 0; }
      parm[np]=Calc(eq,&c);
      err=eq==c;
#if CALC&1
      if (np) {
        if (np==op[nop].dimless) OP3(&parm[np],&parm[np],'i');
        else OP3(&parm[0],&parm[np],'p'); }
      if (unitserr) return 0;
#endif
      if (err) {
        myprintf("parameter %d error",np+1);
        return 0; }
      isparm[np]=1;
      if (!comma) break;
      eq=comma+1; } }
  else {
    /* no = */
    if (op[nop].minparm>0) {
      myprintf("=from,to[,by] expected\n");
      return 0; } }
  np++;

  if (np<op[nop].minparm) {
    myprintf("at least %d parameters expected but %d given",
             op[nop].minparm,np);
    return 0; }

  if (var==c) {
    /* undefined variable */
    if (np) {
      res=parm[0]; assign(var); }
    else {
      myprintf("undefined variable %s",var);
      return 0; } }

  eqn=strtok(NULL,"");
  if (!eqn) {
    if (nop!=UNDEF) {
      myprintf("expression expected\n");
      return 0; } }

  res=parm[0]; /* copy units ... should check units of parm[1], etc. */

  return nop;
}

static
double setprec(double x) /* ---------------------------------------- setprec */
/* set available precision, to be used by solve,integ, etc. */
{
  double d;
  char *c;
  double y,yy=0; /* to avoid warning only */
  int i;

  for (d=1e-33*(fabs(x)+1); d<1e99; d*=2) {
    int n=0;

    for (i=-5; i<=5; i+=2) {
      res VAL=x+d*i;
      assign(var);
      err=eqn==c;
#if CALC&1
	if (unitserr) return 0;
#endif
      if (err) { myprintf("bad equation\n"); return 0; }
      y=Calc(eqn,&c) VAL;
      if (i!=-5) n+=(yy!=y);
      yy=y; }
    if (n>2) return d; }

  return 0;
}

static
double integ(double a,double b,double d) /* -------------------------- integ */
/*
  numerical integration (of eqn), heuristic interval change, 4rd order Gauss
  WARNING - setprec bad close to x^2-1 for x approx 0 and similar
*/
{
  int st=0,sg=1;
  double s1=0,s2=0;
  char *c;


/*.....  d=pow(d,.2);*/
  d=pow(d,.3);
  if (a>b) sg=-1,d=-d;
  while (!st) {
    if (sg*(a+d-b)>=0) st++,d=b-a;
    res VAL=a+0.21132486540519*d;
    assign(var);
    s1+=Calc(eqn,&c) VAL*d;
    err=eqn==c;
#if CALC&1
	if (unitserr) return 0;
#endif
    if (err) { myprintf("integ: bad expression\n"); return 0; }
    if (sig) { myprintf("integ: interrupted\n"); return 0; }

    res VAL=a+0.78867513459481*d;
    assign(var);
    s2+=Calc(eqn,&c) VAL*d;

    a+=d;
    d*=1.003; }

  return (s1+s2)/2;
}

static
double Gauss8(double a,double b,int n) /* --------------------------- Gauss8 */
{
  double h,sum,w1,w2,h1,h2,x,s1,s2;
  int i;
  char *c;
  const double
    q1=0.430568155797026287612,
    q2=0.169990521792428132401,
    w=0.173927422568726928687;

  if (n<=0) return 0;

  sum=0;
  h=(b-a)/n;
  w1=h*w; w2=h/2-w1;
  h1=h*q1; h2=h*q2;
  loop (i,0,n) {
    x=h*(i+0.5)+a;

    res VAL=x-h1;
    assign(var);
    s1=Calc(eqn,&c) VAL;
    err=eqn==c;
#if CALC&1
	if (unitserr) return 0;
#endif
    if (err) { myprintf("integ: bad expression\n"); return 0; }
    if (sig) { myprintf("integ: interrupted\n"); return 0; }

    res VAL=x+h1;
    assign(var);
    s1+=Calc(eqn,&c) VAL;

    res VAL=x-h2;
    assign(var);
    s2=Calc(eqn,&c) VAL;
    err=eqn==c;

    res VAL=x+h2;
    assign(var);
    s2+=Calc(eqn,&c) VAL;

    sum+=s1*w1+s2*w2; }

  return sum;
}

static
int Solve(double eps) /* --------------------------------------------- Solve */
/*
  solve equation eqn=0
  eps<0 => print error messages
*/
{
  int prt=eps<0;
  char *c;

  eps=fabs(eps)*32; /*? - essentially 3e-8*/
  if (eps==0) eps=-1e-12;

  MS_BEGIN(res VAL,eps)
    if (MS_it>1000) {
      if (prt) myprintf("solve: too many iterations\n");
      return 1; }
    assign(var);
    MS_f=Calc(eqn,&c) VAL;
    err=eqn==c;
#if CALC&1
	if (unitserr) return 0;
#endif
    if (err) {
      if (prt) myprintf("solve: bad equation\n");
      return -1; }
    if (sig) {
      if (prt) myprintf("solve: interrupted\n");
      return -1; }
    MS_END(res VAL,1)

  assign(var);

  return 0;
}

static double MinMaxSg=1,MinMaxF;
static
int MinMax(double eps) /* ------------------------------------------- MinMax */
/*
  find minimum or maximum of eqn
  eps<0 => print error messages
*/
{
  int prt=eps<0;
  char *c;

  MinMaxF=0;

  qmin_prt=prt;
  eps=sqrt(fabs(eps));
  if (eps==0) eps=1e-5;

  QMIN_BEGIN(res VAL,eps,5)
    if (QMIN_it>1000) {
      if (prt) myprintf("min/max: too many iterations\n");
      return 1; }
    assign(var);
    MinMaxF=QMIN_f=MinMaxSg*Calc(eqn,&c) VAL;
    err=eqn==c;
#if CALC&1
	if (unitserr) return -1;
#endif
    if (err) {
      if (prt) myprintf("min/max: bad equation\n");
      return -1; }
    if (sig) {
      if (prt) myprintf("min/max: interrupted\n");
      return -1; }
    QMIN_END(res VAL)

  assign(var);
  MinMaxF=MinMaxSg*Calc(eqn,&c) VAL;

  return 0;
}


double froot_ext;

int Checkres(double eps) /***************************************** Checkres */
/* check numerical solution of Solve */
{
  char *c;

  assign(var);
  froot_ext=Calc(eqn,&c) VAL;

  return fabs(froot_ext)<fabs(300*eps); /* 3e-7 of typical y-value */
}

double elemmass(char *c,char **ret) /****************************** elemmass */
/* for functions M(), m() */
{
  char elem[4];
  struct mend_s {
    char *name;
    double mass; } *m,mend[]={
      {"H",1.00794},
      {"D",2.014101778},
      {"T",3.016049268},
      {"He",4.002602},
      {"Li",6.941},
      {"Be",9.012182},
      {"B",10.811},
      {"C",12.0107},
      {"N",14.0067},
      {"O",15.9994},
      {"F",18.9984032},
      {"Ne",20.1797},
      {"Na",22.989770},
      {"Mg",24.3050},
      {"Al",26.981538},
      {"Si",28.0855},
      {"P",30.973761},
      {"S",32.065},
      {"Cl",35.453},
      {"Ar",39.948},
      {"K",39.0983},
      {"Ca",40.078},
      {"Sc",44.955910},
      {"Ti",47.867},
      {"V",50.9415},
      {"Cr",51.9961},
      {"Mn",54.938049},
      {"Fe",55.845},
      {"Co",58.933200},
      {"Ni",58.6934},
      {"Cu",63.546},
      {"Zn",65.39},
      {"Ga",69.723},
      {"Ge",72.64},
      {"As",74.92160},
      {"Se",78.96},
      {"Br",79.904},
      {"Kr",83.80},
      {"Rb",85.4678},
      {"Sr",87.62},
      {"Y",88.90585},
      {"Zr",91.224},
      {"Nb",92.90638},
      {"Mo",95.94},
      {"Tc",98},
      {"Ru",101.07},
      {"Rh",102.90550},
      {"Pd",106.42},
      {"Ag",107.8682},
      {"Cd",112.411},
      {"In",114.818},
      {"Sn",118.710},
      {"Sb",121.760},
      {"Te",127.60},
      {"I",126.90447},
      {"Xe",131.293},
      {"Cs",132.90545},
      {"Ba",137.327},
      {"La",138.9055},
      {"Ce",140.116},
      {"Pr",140.90765},
      {"Nd",144.24},
      {"Pm",145},
      {"Sm",150.36},
      {"Eu",151.964},
      {"Gd",157.25},
      {"Tb",158.92534},
      {"Dy",162.50},
      {"Ho",164.93032},
      {"Er",167.259},
      {"Tm",168.93421},
      {"Yb",173.04},
      {"Lu",174.967},
      {"Hf",178.49},
      {"Ta",180.9479},
      {"W",183.84},
      {"Re",186.207},
      {"Os",190.23},
      {"Ir",192.217},
      {"Pt",195.078},
      {"Au",196.96655},
      {"Hg",200.59},
      {"Tl",204.3833},
      {"Pb",207.2},
      {"Bi",208.98038},
      {"Po",209},
      {"At",210},
      {"Rn",222},
      {"Fr",223},
      {"Ra",226},
      {"Ac",227},
      {"Th",232.0381},
      {"Pa",231.03588},
      {"U",238.02891},
      {"Np",237.05},
      {"Pu",244},
      {"Am",243},
      {"Cm",247},
      {"Bk",247},
      {"Cf",251},
      {"Es",252},
      {"Fm",257},
      {"Md",258},
      {"No",259},
      {"Lr",262},
      {"Rf",267},
      {"Db",262},
      {"Sg",269},
      {"Bh",264},
      {"Hs",269},
      {"Mt",278},
      {"Ds",281},
      {"Rg",282},
      {"Cn",285},
      {"Nh",286},
      {"Fl",289},
      {"Mc",289},
      {"Lv",293},
      {"Ts",294},
      {"Og",294},
      {NULL,0} };

  elem[0]=*c;
  elem[1]=elem[2]=elem[3]=0; /* used to support UUn etc. */
  *ret=c+1;
  if (islower(c[1])) {
    elem[1]=c[1];
    *ret=c+2;
    if (islower(c[2])) {
      elem[2]=c[2];
      *ret=c+3; } }

  for (m=mend; m->name; m++) if (!strcmp(m->name,elem)) return m->mass;

  myprintf("%s is unknown element");

  return 0;
}

void molarmass(char key,char *c) /******************************** molarmass */
/*
  In place replacement of M(FORMULA) or m(FORMULA) by MM<key>[spaces],
  where variable MM<key> (key=A..Z) contains the (molar) mass.
  char *c is a pointer to M( or m( found
  Accepted expressions:
    M(C2H6)
    M(CuSO4 H10O5)
    M(Fe0.5)
  Wrong:
    CuSO4.5H20   - accepted, but interpreted as Cu S O4.5 H2 O
    CuSO4 (H2O)5 - parentheses not supported
  molar=1: mass in g/mol (w/o units), or any compatible unit with units
  molar=0: mass in kg, or any mass unit with units
*/
{
  char *e=strchr(c,')');
  char *b=c;
  int molar=*c=='M'; /* m() or M() */
  double m,x,s=0;

  c+=2;

#if CALC&1
  OP3(&res,&res,'0');
  res.pow[1]=1;
  res.pow[4]=-molar; /* kg/mol or kg */
#endif

  if (!e) {
    myprintf("missing ) after M( or m(");
    return; }
  else if (e==c) {
    myprintf("M() or m() without compound");
    return; }

  if (key>'Z') {
    myprintf("max 26 functions M() or m() allowed in 1 line");
    return; }

  for (;;) {
    while (*c==' ') c++;
    m=elemmass(c,&c);
    if (*c==0)  {
      myprintf("M() or m(): syntax or internal error");
      return; }
    x=strtod(c,&e); /* real number, stoichiometric coeff. after element name */
    if (c==e) x=1;
    c=e;
    if (*c==0) {
      myprintf("M() or m(): syntax or internal error");
      return; }
    s+=m*x; /* molar mass summed up */

    if (*c==')') {
      /* matching ) found: s=molar mass in g/mol */
      for (e=b; e<=c; e++) *e=' ';
      b[0]=b[1]='M'; /* MM<key> used for both */
      b[2]=key;
      b[3]=0; /* temporary: b=MM<key> */
      x=res VAL; /* remember */
      if (molar) {
#if CALC&1
	res VAL=s*1e-3; /* g->kg conversion */
#else
	res VAL=s;      /* no conversion (result in g/mol) */
#endif
      } else {
#if CALC&1
	res VAL=s/6.02214076e+26; /* g->kg conversion */
#else
	res VAL=s/6.02214076e+26; /* g->kg conversion */
#endif
      }
      assign(b); /* assign res to variable b=MM<key> */
      res VAL=x; /* return back - is it needed? */
      b[3]=' '; /* restore space filling after M<key> */
      return;
    }
  }
}

void calculate(char *expr) /************************************** calculate */
/* main function to calculate an expression or line */
{
  char *c;
  char *b,*e,key,extra;

  sig=0;

  if (strchr(expr,'%')) {
    /* format */
    if (strlen(expr)>NFMT) {
      myprintf("TOO LONG FORMAT (for modulo use @)\n");
      return; }
    cast=DOUBLE;
    if (!strcmp(expr,"%")) {
#if CALC&1
	strcpy(fmt," %.8g %s\n");
#else
	strcpy(fmt," %.8g\n");
#endif
      printres();
      return; }
    if (expr[1]=='?') {
      myprintf("current format =%s",fmt);
      return; }
    if (strlen(expr)<NFMT) for (c="eEgGfiduxXoc"; *c; c++) if (strchr(expr,*c)) {
      strcpy(fmt,mode==TEX?"@":" ");
      strcat(fmt,expr);
#if CALC&1
	strcat(fmt," %s");
#endif
      if (mode) strcat(fmt,"\n");
      if (strchr("idc",*c)) cast=INT;
      if (strchr("uxXo",*c)) cast=UNSIGNED;
      printres();
      return; }
    myprintf("BAD FORMAT (for modulo use @)\n");
    return; }

 {
   /* prevent macro expansion within commands def,undef,expand */
   char x[8],*c;

   x[7]=0;
   memcpy(x,expr,7);
   for (c=x; *c; c++) if (!isalpha(*c)) *c=0;
   /* ? - def,undef,expand should not appear here... */
   if (strcmp(x,"def") && strcmp(x,"undef") && strcmp(x,"expand"))
     expmacro(expr);
 }

  /* molar mass, e.g. M(H2O) */
  key='A'; /* A,B,..*/
  for (c=expr; *c; c++) if (toupper(*c)=='M' && c[1]=='(') molarmass(key++,c);

  /* looking for id=expr */
  c=expr;

  while (*c==' ') c++;
  b=c;
  if (isalpha(*c) || *c=='_') while (isalnum(*c) || *c=='_') c++;
  e=c;
  while (*c==' ') c++;

  /* := (the same as =) += -= *= /= @= */
  extra='=';
  if (*c && strchr("+-*/:@",*c) && c[1]=='=') {
    extra=*c; *c='='; c[1]=' '; }
  if (extra==':') extra='=';

  key=*c;

  if ( (key=='=' && c[1]!='=') || key=='~')
    if (b==e) {
      /* =id */
      /* NOTE: +=id should add res to id, but is caught as error earlier */
      c++;
      while (*c==' ') c++;
      b=c;
      if (isalpha(*c) || *c=='_') while (isalnum(*c) || *c=='_') c++;
      e=c;
      *e=0;
      if (b==e) {
        struct _Idlist_s *is,*isnext;

        /* remove all id's */
        if (key=='~') {
          for (is=_Id.head; is; is=isnext) {
            isnext=is->next;
            free(is); }
          _Id.head=NULL; }

        /* print all id's - reverse order */
        if (key=='=') {
          int n,i;
          struct mac_s *m;
#if CALC&1
	struct user_s *u;
#endif

          for (n=0,is=_Id.head; is; is=is->next) n++;
          while (n) {
            n--;
            for (i=0,is=_Id.head; is; is=is->next) if (i++==n) break;
#if CALC&1
	myprintf(" %s = %.15g %s\n",is->id,is->val VAL,unit(&is->val));
#else
	myprintf(" %s = %.15g\n",is->id,is->val);
#endif
          }
#if CALC&1     
	    myprintf(" # = %.15g %s\n",res VAL,unit(&res));
#else          
	    myprintf(" # = %.15g\n",res);
#endif
          looplist (m,mhead) myprintf("def %s %s\n",m->X,m->E);
#if CALC&1
	looplist (u,userhead) myprintf("to [%s]\n",u->unit);
#endif
        } }
      else if (key=='~') deassign(b,"removed");
      else assign(b); /* assignop=extra (see above) */ }
    else {
      /* id=expr */
      calculate(c+1); /* called recursively: may handle cases as a=b=expr */
      if (!err) {
        int x=*e;
        *e=0;
        assignop=extra;
        assign(b);
        if (extra!='=') printres(); /* result of += -= *= /= */
        *e=x; } }
  else if (oper(expr)) switch (nop) {
    case SOLVE:
      if (isparm[1]) {
        int nr=0,n=500,i,j,flag;
#define NR 20
        double root[NR],froot[NR];
        double eps=1e-9*fabs(parm[1] VAL-parm[0] VAL);
        double from=parm[0] VAL,to=parm[1] VAL;
        double yrange=0;

        if (parm[1] VAL<parm[0] VAL) from=parm[1] VAL,to=parm[0] VAL;

        if (!isparm[2]) parm[2] VAL=0;
        if (parm[2] VAL) n=(int)(fabs((parm[1] VAL-parm[0] VAL)/parm[2] VAL)+.5);
        if (n<2) n=2;
        if (n>100000) {
          myprintf("solve: too fine grid\n");
          return; }

        for (i=0; i<=n; i++) {
          char *c;

          res VAL=parm[0] VAL+(double)i/n*(parm[1] VAL-parm[0] VAL);
#if CALC&1
	memcpy(res.pow,parm[0].pow,sizeof(res.pow));
#endif

          if (sig) return;
          assign(var);
          //       fprintf(stderr,"%g  %g %g %g<<<\n",res VAL,res.pow[0],res.pow[1],res.pow[2]);
          res=Calc(eqn,&c);
          err=eqn==c;
#if CALC&1
	if (unitserr) return;
#endif
          if (err) {
            myprintf("solve: bad equation\n");
            return; }
          yrange+=fabs(res VAL); }

        yrange/=n;
        if (isnan(yrange)) yrange=9e99; /* ad hoc! */

        for (i=0; i<=n; i++) {
          res VAL=parm[0] VAL+(double)i/n*(parm[1] VAL-parm[0] VAL);
#if CALC&1
	memcpy(res.pow,parm[0].pow,sizeof(res.pow));
#endif
          if (!Solve(eps) && res VAL>=from-eps && res VAL<to+eps && Checkres(eps*yrange)) {
            flag=1;
            for (j=0; j<nr; j++) if (fabs(root[j]-res VAL)<eps) {
              flag=0;
              break; }
            if (flag) {
              if (nr>=NR) {
                myprintf("solve: too many roots (>%d)\n",NR);
                return; }
              root[nr]=res VAL;
              froot[nr]=froot_ext;
              nr++; } } }

        /* sort roots */
        do {
          flag=0;
          for (i=1; i<nr; i++)
            if ( (root[i-1]-root[i])*(parm[0] VAL-parm[1] VAL)<0 ) {
              double aux;
              aux=root[i-1]; root[i-1]=root[i]; root[i]=aux;
              aux=froot[i-1]; froot[i-1]=froot[i]; froot[i]=aux;
              flag++; }
          } while (flag);

        for (i=0; i<nr; i++) {
          res VAL=root[i];
          assign(var);
          if (verbose&2) printres2(" f=%.2g",froot[i]);
          else printres(); } }
      else {
        if (Solve(-1e-9*fabs(res VAL=parm[0] VAL))<0) return;
        printres(); }
      break;

    case MAX:
      MinMaxSg=-1;
      goto MINMAX;
    case MIN:
      MinMaxSg=1;

    MINMAX:
      if (isparm[1]) {
        /* in interval */
        int n=1000,i,imin=0; /* to avoid warning only */
        double eps,X,absmin=9e99;

        if (parm[0] VAL>parm[1] VAL) X=parm[0] VAL,parm[0] VAL=parm[1] VAL,parm[1] VAL=X;
        eps=1e-8*(parm[1] VAL-parm[0] VAL);

        if (!isparm[2]) parm[2] VAL=0;
        if (parm[2] VAL) n=(int)(fabs((parm[1] VAL-parm[0] VAL)/parm[2] VAL)+.5);
        if (n<2) n=3;
        if (n>100000) {
          myprintf("min/max: too fine grid\n");
          return; }

        for (i=0; i<=n; i++) {
          //          res=parm[0]+(double)i/n*(parm[1]-parm[0]);
          X=2.*i/n-1;
          res VAL=parm[0] VAL+(parm[1] VAL-parm[0] VAL)*(X*(1-fabs(X)/2)+0.5);
#if CALC&1
	memcpy(res.pow,parm[0].pow,sizeof(res.pow));
#endif
          if (sig) return;
          assign(var);
#if CALC&1
	res=Calc(eqn,&c); res VAL*=MinMaxSg;
#else
	res=MinMaxSg*Calc(eqn,&c);
#endif
          if (res VAL<absmin) {
            imin=i;
            absmin=res VAL; } }

        //        res=parm[0]+(double)imin/n*(parm[1]-parm[0]);
        X=2.*imin/n-1;
        res VAL=parm[0] VAL+(parm[1] VAL-parm[0] VAL)*(X*(1-fabs(X)/2)+0.5);
#if CALC&1
	memcpy(res.pow,parm[0].pow,sizeof(res.pow));
#endif
        X=res VAL;
        if (MinMax(eps) || MinMaxF>absmin || res VAL<parm[0] VAL || res VAL>parm[1] VAL) {
          /* not minimized or worse or not in range */
          res VAL=X; MinMaxF=absmin; }

        assign(var);
        printres2(MinMaxSg>0?" min=%.12g":" max=%.12g",MinMaxSg*MinMaxF); }
      else {
        res=parm[0];
        if (MinMax(1e-8*fabs(res VAL))<0) return;
        printres2(MinMaxSg>0?" min=%.12g":" max=%.12g",MinMaxSg*MinMaxF); }
      break;

    case SUM: case PROD: {
      double s=nop==PROD; /* units solved later */
#if CALC&1
	int i,npow=0;
#endif
      CALC_t x;

      if (!isparm[2]) parm[2] VAL=0;
      if (parm[2] VAL==0) parm[2] VAL=parm[0] VAL<=parm[1] VAL?1:-1;
      if (parm[0] VAL>parm[1] VAL) x VAL=parm[1] VAL,parm[1] VAL=parm[0] VAL,parm[0] VAL=x VAL,parm[2] VAL*=-1;
      for (res VAL=parm[0] VAL; res VAL<=parm[1] VAL; res VAL+=parm[2] VAL) {
#if CALC&1
	memcpy(res.pow,parm[0].pow,sizeof(res.pow));
#endif
        assign(var);
        x=Calc(eqn,&c);
        err=eqn==c;
#if CALC&1
	if (unitserr) return;
#endif
        if (err) { myprintf("sum/prod: bad expression\n"); return; }
        if (sig) { myprintf("sum/prod: interrupted\n"); return; }
        if (nop==PROD) s*=x VAL; else s+=x VAL;
#if CALC&1
	npow++;
#endif
      }
#if CALC&1
	if (nop==PROD) loop (i,0,CALC_NUNITS) x.pow[i]*=npow;
#endif
      res=x;
      res VAL=s;
      assign(var);
    enough:
      printres();
      break; }

    case REPEAT: case ITER: {
      int i,n=(int)parm[1] VAL,iserr=0;
      double last=parm[0] VAL;
      double err;

      res VAL=parm[0] VAL;
      if (nop==ITER) n=100000;

      for (i=0; i<n; i++) {
        assign(var);
        res=Calc(eqn,&c);
        err=eqn==c;
#if CALC&1
	if (unitserr) return;
#endif
        if (err) { myprintf("iter/repeat: bad expression\n"); return; }
        if (sig) { myprintf("iter/repeat: interrupted\n"); return; }
        if (nop==ITER) {
          iserr++;
          if (fabs(err=res VAL-last)<=parm[1] VAL) break;
          last=res VAL; } }
      assign(var);
      if (iserr) myprintf("err=%g  ",err);
      printres();
      break; }

    case CONTFRAC: {
      /* format of command: contfrac VAR=FROM,TO[,BY] [BEFORE++]NUM//DENOM */
      double x=0,sg=1;
      double p=0,pp=1,q=1,qq=0,a,b,z;
      char *plus,*slash;

      if (!isparm[2]) parm[2] VAL=0;
      if (parm[2] VAL==0) parm[2] VAL=parm[0] VAL<=parm[1] VAL?1:-1;
      if (parm[2] VAL<0) sg=-1;
      plus=strstr(eqn,"++");
      if (plus) {
        *plus=0;
        x=Calc(eqn,&c) VAL;
        err=eqn==c;
#if CALC&1
	if (unitserr) return;
#endif
        if (err) { myprintf("contfrac: bad expression before ++\n"); return; }
        eqn=plus+2; }

      slash=strstr(eqn,"//");
      if (!slash) {
        myprintf("contfrac: missing //\nSynopsis: contfrac N=FROM,TO[,BY] [BEFORE++]NUM//DENOM\n");
        return; }
      *slash=0;
      slash+=2;

      for (res VAL=parm[0] VAL; res VAL*sg<=parm[1] VAL*sg; res VAL+=parm[2] VAL) {
        assign(var);
        b=Calc(eqn,&c) VAL;
        err=eqn==c;
#if CALC&1
	if (unitserr) return;
#endif
        if (err) { myprintf("contfrac: bad numerator\n"); return; }
        a=Calc(slash,&c) VAL;
        err=eqn==c;
        if (err) { myprintf("contfrac: bad denominator\n"); return; }
        if (sig) { myprintf("contfrac: interrupted\n"); return; }
        z=q*a+qq*b; qq=q; q=z;
        z=p*a+pp*b; pp=p; p=z; }

      res VAL=parm[0] VAL; assign(var);
      res VAL=x+p/q;
      goto enough; }

    case INTEG: {
      double s;
#if CALC&1
	int i;
#endif

      if (isparm[2])
        s=Gauss8(parm[0] VAL,parm[1] VAL,parm[2] VAL);
      else {
        double x,d;

        /* WARNING bad algorithm for d !!!!*/
        d=setprec(parm[1] VAL);
        x=setprec(parm[0] VAL);
        if (x<d) d=x; /* ??? */
        if (d==0) d=1e-5; /* const */

        if (parm[0] VAL*parm[1] VAL<0) {
          x=setprec(0);
          if (x<d) d=x;
          d*=.1;
          s=integ(0,parm[1] VAL,d)-integ(0,parm[0] VAL,d); }
        else if (fabs(parm[0] VAL)>fabs(parm[1] VAL))
          s=-integ(parm[1] VAL,parm[0] VAL,d);
        else
          s=integ(parm[0] VAL,parm[1] VAL,d);
      }

#if CALC&1
	res=Calc(eqn,&e); loop (i,0,CALC_NUNITS) res.pow[i]+=parm[0].pow[i];
#endif
      res VAL=s; }

      goto enough;

    case DIFF:
    case DERIV: {
      double d,y0[5],*y=y0+2;
      int i;

      if (isparm[1]) d=parm[1] VAL;
      else d=pow(setprec(parm[0] VAL),0.2);

      for (i=-2; i<=2; i++) if (i) {
        res VAL=parm[0] VAL+d*i;
        assign(var);
        y[i]=Calc(eqn,&c) VAL;
        err=eqn==c;
#if CALC&1
	if (unitserr) return;
#endif
        if (err) {
          myprintf("deriv: bad expression\n");
          return; } }
      res VAL=parm[0] VAL; assign(var);
#if CALC&1
	res=Calc(eqn,&e); loop (i,0,CALC_NUNITS) res.pow[i]-=parm[0].pow[i];
#endif
      res VAL=(8*(y[1]-y[-1])-(y[2]-y[-2]))/(12*d);

      goto enough; }

    case PLOT: {
      FILE *plot=fopen("/tmp/ev.tmp","wt");
      double y,unitfx=1,unitfy=1,half=0;
      int n=999; /* to avoid 0/0 in 'plot x=-10*pi,10*pi sin(x)/x' */
      int i,nto;

      if (!isparm[2]) parm[2] VAL=0;
      if (parm[2] VAL) n=(int)(parm[2] VAL);

      if (n<0) {
        if (n>-2) n=-2;
        n=-n; nto=n; half=0.5; }
      else {
        if (n<1) n=1;
        nto=n+1; half=0; }

      if (abs(n)>1000001) {
        myprintf(EV ": plot: too many points\n");
        return; }

      loop (i,0,nto) {
        res VAL=parm[0] VAL+(i+half)*(parm[1] VAL-parm[0] VAL)/n;
        assign(var);
#if CALC&1
        {
          struct unit_s RES=Calc(eqn,&c);
          y=RES VAL;
          if (i==0) {
            unitfx=unitf(&parm[1]);
            unitfy=unitf(&RES); }
        }
#else
        y=Calc(eqn,&c) VAL;
#endif
        err=eqn==c;
#if CALC&1
	if (unitserr) return;
#endif
        if (sig) {
          myprintf("plot: interrupted\n");
          return; }
        if (err) {
          myprintf("plot: bad expression\n");
          return; }
        fprintf(plot,"%.16g\t%.16g\n",res VAL/unitfx,y/unitfy);
      }
      res VAL=parm[0] VAL; assign(var);
      fclose(plot);
      if (system("plot /tmp/ev.tmp &"))
        myprintf("plot: system call failed\n");
      break; }

    case DEF:
      defmacro(var,eqn);
      deassign(var,"defined");
      //      if (mode==FILTER) myprintf(" (defined)\n");
      break;
    case UNDEF:
      undefmacro(var);
      deassign(var,"undefined");
      //      if (mode==FILTER) myprintf(" (undefined)\n");
      break;
    case EXPAND:
      expand(var,eqn);
      deassign(var,"expanded");
      //      if (mode==FILTER) myprintf(" (expanded)\n");
      break;
      //    HELP:
      //      break;
    }
  else {
    res=Calc(expr,&c);
    err=c==expr;
    if (!err) {
      while (*c==' ') c++;
      if (*c) myprintf("ERROR at %s \n",c); }
    //    if (mode) myprintf("%s =",expr);
    if (err)
      myprintf(" ERROR\n");
    else
      printres(); }
}

#include "evhlp.c"

int main(int narg,char **arg) /**************************************** main */
{
  char line0[PROMPT_N+1];
  char *line=line0+1;
  char *c,*semcol,*Xexpr="+-*/^@";
  char *evdata="/";
  FILE *f=NULL;
  int i,iarg,printed;
  int (*auxprintf)(const char *format,...);
  int retres=0;
  char *arg__c[]={"","--c"};

  if (narg<2) {

    prtsfill(Pintro);

    for (;;) {
      char s[8];

      fputs(Phelp,stdout);
      if (!fgets(s,8,stdin)) s[0]='Q';

      switch (s[0]) {
        case 'i': prtsfill(Pintro); break;
        case 'o': prtsfill(Poptions); break;
        case 'n': prtsfill(Pnumvar); break;
        case 'e': prtsfill(Penv);
          printf("\nCompile time option (see gen/prompt.c): max line length = %d\n",PROMPT_N);
          break;
        case 'x': prtsfill(Pexpr); break;
        case 'c': prtsfill(Pcommands); break;
        case 'f': prtsfill(Pfunc); break;
        case 's': prtsfill(Psolve); break;
        case 'a': prtsfill(Palg); break;
        case 'm': prtsfill(Pmac); break;
        case 'p': prtsfill(Pplot); break;
        case 't': prtsfill(Ptoggles); break;
        case 'u': prtsfill(Punits); break;
        case '.': prtsfill(Ppu); break;
        case 'd': prtsfill(Pdata); break;
        case 'r': prtsfill(Pret); break;
        case 10: arg=arg__c; narg=2; goto start__c;

        default: return 0; } } }


 start__c:;

  myprintf=printf;
  if (getenv("FMT")) {
    strcpy(fmt,getenv("FMT"));
    strcat(fmt,"\n"); }
  if (getenv("EVFMT")) {
    strcpy(fmt,getenv("EVFMT"));
    strcat(fmt,"\n"); }
  else
#if CALC&1
	strcpy(fmt,"%.8g %s\n");
#else
	strcpy(fmt,"%.8g\n");
#endif

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-' && arg[iarg][1]=='-') switch(arg[iarg][2]) {
      case 'c': {
        mode=PROMPT;
        evgetline=mygetline;
        myprintf=cprintf;
#if CALC&1
	strcpy(fmt," %.8g %s\n");
#else
	strcpy(fmt," %.8g\n");
#endif
        _pr.n=atoi(arg[iarg]+3);
        if (_pr.n==0) {
          if (getenv("COLUMNS")) _pr.n=atoi(getenv("COLUMNS")); }
        if (_pr.n<8) _pr.n=80;
        _pr.n--;
        for (i=0; i<PROMPT_H; i++) strcpy(_pr.hist[i],"?>");

        initscr(); cbreak();
        noecho(); nl();
        scrollok(stdscr,TRUE);
        intrflush(stdscr, FALSE); }
        break;
      case 'd':
        mode=DUMB;
#if CALC&1
	strcpy(fmt," %.8g %s\n");
#else
	strcpy(fmt," %.8g\n");
#endif
        myprintf=printf;
        evgetline=gets;
        goto VERBOSE;
        break;
      case 'f':
        strcpy(fmt,arg[iarg]+3);
#if CALC&1
	strcat(fmt," %s\n");
#else
	strcat(fmt,"\n");
#endif
        break;
      case 'F':
        strcpy(fmt,arg[iarg]+3);
#if CALC&1
	strcat(fmt," %s");
#endif
        break;
      case 'i':
        evdata=arg[iarg]+3;
        break;
      case 'r':
        retres=1;
        goto VERBOSE;
      case 't':
        mode=TEX;
        prompt=0;
#if CALC&1
	strcpy(fmt,"@%.8g %s\n");
#else
	strcpy(fmt,"@%.8g\n");
#endif
        myprintf=texprintf;
        evgetline=gets;
        break;
      case '-':
      case '+':
        /* -,+ not allowed as expression continuation */
        Xexpr="*/^@";
      case 'v':
      VERBOSE:
        verbose=atoi(arg[iarg]+3);
        if (arg[iarg][3]==0) verbose=1;
        break;
      case 0:
        mode=FILTER;
        prompt=0;
        myprintf=printf;
        evgetline=gets;
        break;
      default:
        if (isdigit(arg[iarg][2])) {
          mode=FILTER;
          prompt=0;
          myprintf=printf;
          evgetline=gets;
          verbose=atoi(arg[iarg]+2);
          break; }
        else {
          fprintf(stderr,"ev: wrong option (enclose expression starting with -- in ())\n");
          exit (0); }
      }
    else break;

  signal(SIGINT,asksig);

  /* reading the initial constants */
  if (!strcmp(evdata,"/")) {
    char *h=getenv("HOME");
    if (h) {
      evdata=malloc(strlen(h)+10); /* bug fixed by J. Janek 2024-08-27 */
      strcpy(evdata,h);
      strcat(evdata,"/." EV "data"); } }

  if (strlen(evdata)) f=fopen(evdata,"rt");

  auxprintf=myprintf;
  if (!(verbose&8)) myprintf=dummyprintf;

  isconst=1;
  if (f)
    while (fgets(line0,PROMPT_N+1,f)) {
      char *u=strlast(line0);

      /* removing trailing \n */
      if (u && (unsigned char)*u<=' ') *u=0;

      /* removing trailing comment (w. preceding spaces) */
      comment(line0);

      myprintf("%s\n",line0);
      if (strlen(line0)>2) {
#if CALC&1
	if (!memcmp(line0,"to ",3)) towards(line0+3); else
#endif
        calculate(line0); } }
  isconst=0;

  myprintf=auxprintf;

  strcpy(line0,"#");

  if (mode==TEX)
    printf("\
\\documentclass[twocolumn]{article}\n\
\\usepackage[utf8]{inputenc}\\usepackage[czech]{babel}\\usepackage[IL2]{fontenc}\n\
\\usepackage{/home/jiri/tex/evu}\n\
\\usepackage{color}\n\
\\parskip 0pt\\parindent 0pt\n\
\\topmargin -2cm\\oddsidemargin -1.6cm\n\
\\textwidth 19.12cm\\textheight 26cm\n\
\\columnwidth 9cm\\columnsep 0.62cm\n\
\\begin{document}\n\
\\pdfpagewidth=210mm \\pdfpageheight=297mm\n\
\\boldmath\n\
\\tt\n");

  for (;;) {
    if (prompt) {
      strcpy(_pr.prompt,"==");
      if (_Id.degrad==1) _pr.prompt[0]=toclip?'R':'r';
      else _pr.prompt[0]=toclip?'D':'d';
      _pr.prompt[1]=">,#;"[inputmod]; }
    else
      strcpy(_pr.prompt,"");

    /* ugly: for smart terminal, _pr.prompt used in mygetline (see prompt.c) */
    if (mode==DUMB) myprintf("%s",_pr.prompt);

    printed=0; /* echo only 1 comment/line: see --v8 and FILTER */

    if (iarg<narg) strcpy(line,arg[iarg++]);
    else if (!evgetline || !evgetline(line)) break;
    {
      char *c,*e=NULL;
      for (c=line; *c; c++) if ((unsigned char)*c<=' ') e=c; else e=NULL;
      if (e) *e=0;

      if (getenv("EV"))
        for (e=getenv("EV"); *e; e+=2)
          for (c=line; *c; c++) if (*c==e[0]) *c=e[1];
    }

    /* "." removed, use ctrl-d */
    if (!strcmp(line,"exit") || !strcmp(line,"quit")) break;

    if (mode==TEX) {
      if (!memcmp(line,"! #.",4)) {
        static int no=0;
        sprintf(line+1,"%2d",++no);
        line[3]='.'; }

      if (strlen(line)==0) strcpy(line,"!");
      texprintf("%s\n",line); }

    c=comment(line);
    if (c)
      if (mode==FILTER && (verbose&16) && !printed) {
        printed++;
        putchar('!');
        puts(c+1); }

    /* , -> . */
    if (inputmod&1) {
      for (c=line; *c; c++)
        if (*c==',') {
          if (c[1]==',') memmove(c+1,c+2,strlen(c+1));
          else *c='.'; } }

    /* remove 1 space */
    if (inputmod&2) {
      for (c=line; *c; c++)
        if (*c==' ' && c[1]!=' ') memmove(c,c+1,strlen(c)); }

    do { /* loop over ; separating statements */
      semcol=strchr(line,';');
      if (semcol) *semcol=0;

      if (line[0]) {
        if (line[0]=='?') myprintf("for help, run " EV " without parameters");
        else if (!strcmp(line,"rad")) _Id.degrad=1;
        else if (!strcmp(line,"deg")) _Id.degrad=PI/180;
        else if (!strcmp(line,"comma")) inputmod^=1;
        else if (!strcmp(line,"space")) inputmod^=2;
        else if (!strcmp(line,"clip")) {
          if (toclip) toclip=NULL;
          else toclip=getenv("TOCLIP"); }
        else if (!memcmp(line,"sh ",3)) {
          res VAL=system(line+3);
          printres(); }
        else if (!memcmp(line,"cd ",3)) {
          res VAL=chdir(line+3);
          printres(); }
        else if (!memcmp(line,"exe ",4)) {
          char *exe=malloc(strlen(line)+30),*c;
          FILE *evexe;

          strcpy(exe,line+4);

          if ( (c=strchr(exe,'#')) ) {
            sprintf(c,"%.16g",res VAL);
            strcat(c,strchr(line+4,'#')+1); }
          res VAL=0;
          strcat(exe," > .evexe~");
          if (system(exe)) fprintf(stderr,"exe failed\n");
          if ( (evexe=fopen(".evexe~","rt")) ) {
            if (!fscanf(evexe,"%lf",&res VAL))
              fprintf(stderr,"bad output\n");
            fclose(evexe); }
          free(exe);
          printres(); }
        else if (!memcmp(line,"write ",6)) prompt_write(line+6);
        else if (!memcmp(line,"verbose ",8)) verbose=atoi(line+8);
        else if (!strcmp(line,"write")) prompt_write("history." EV);
#if CALC&1
        else if (!memcmp(line,"to ",3)) {
          towards(line+3);
          checkto=1;
          /* last result printed even if 'to' applied to another unit */
          printres();
          checkto=0; }
#endif
        else if (!strcmp(line,"prompt")) prompt^=1;
        else
          if (strchr(Xexpr,line[0])) {
            myprintf("# %c expr =",line[0]);
            calculate(line0); }
          else {
            /* DIRTY HACK for equivalent syntax def A=x meaning def A (x)
             and expand A=x meaning expand A (x)
             BUG: cannot use def a=b in .ev[u]data */

            int dirty=0;

            if (!memcmp(line,"def ",4)) dirty=4;
            if (!memcmp(line,"expand ",7)) dirty=7;

            if (dirty) {
              char *c=line+dirty;

              while (*c==' ') c++;
              for (; *c; c++) {
                if (*c=='=') {
                  memmove(c+1,c,strlen(c)+1);
                  c[0]=' ';
                  c[1]='(';
                  strcat(c,")");
                  break; }
                if (*c==' ') break; } }

            calculate(line); }
        }
      if (semcol) memmove(line,semcol+1,strlen(semcol+1)+1);

#if CALC&1
      if (to_remove) {
        towards(to_remove);
        free(to_remove);
        to_remove=NULL; }
#endif

    } while (semcol); }

  if (mode==PROMPT) endwin();
  if (mode==TEX) printf("\\end{document}\n");

  if (retres) return res VAL==0;
  else return -err;
}
