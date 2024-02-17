/*
  07/2021: mild polishing
  04/2016: CALC simplified, DOS removed, merged with alloc.c and rndgen.c
           see rndgeni.c for the random number generator
  02/2010: PAR removed, EOF on in = just message (no error)
  04/2006: error output to out suppressed if in=NULL
  07/2005: stdarg incompatibility fixed in myError
  10/2001: minimum buffer for initscroll is 1024
  06/2000: bug fixed in commands $L $C
  11/1999: initscroll(negative) causes ERROR not to ask but stop

  This module depends on the following #defines:
  -DSCR   : scrolling of previous output available (default=no scrolling);
            also $iFILE on input requires -DSCR

  further #define variables with changable defaults are:
    STRLEN = 1 line for scrolling input (incl. $iFILE)
    GETBUFLEN = 1 line for the getdata input
  environment variables:
    LINES COLUMNS
*/

#define GETBUFLEN 1024
#define STRLEN 256

/* WARNING: alloc() etc. is #defined and free() re#defined in ground.h */
#include "ground.h"
/* here in ground.c, the original free() must apply: */

#define REAL double
#define FN(X) X
#include "func.c"
#undef REAL
#undef FN

/* to prevent optimizing a double in 10-real register */
double casttodouble(double x)
{
  double y;
  memcpy(&y,&x,sizeof(y));
  return y;
}

#if PRECISION==2
/* see ground.h */
#  define REAL long double
#  define CAT(X,Y) X##Y
#  define FN(X) CAT(X,l)
#  include "func.c"
#  undef FN
#  undef CAT
#  undef REAL
#endif /*# PRECISION==2 */

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                  in and out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

formatted i/o functions (redirectable by scroll)
~~~~~~~~~~~~~~~~~~~~~~~
prt(format,...) = fprintf(out,format,...) ended by line feed
prt_(format,...) = fprintf(out,format,...);
prts_(char *s) = fputs(s,out);  (without LF)
prts(char *s) = fputs(s,out); + LF
prtc(c) = fputc(c,out)
***/

FILE *in,*out; /*** must be initialized when program is started
                - static initializing by stdin or stdout is not portable !!! */
static FILE *masterin; /* to solve $iFILE scroll statement */

char *_S_buf=NULL; /* buffer for scroll */

#ifdef SCR
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% scroll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

(c) J.Kolafa 10/1991 (OBSOLETE, but not harmful)

This unit defines prt, prt_, and prtc so that (since initscroll(int buflen)
has been called) all data printed to out are stored in a buffer.
Any call to scroll() stops the program execution and allows looking to the
previous output up to the buffer capacity.

See the help below for available commands!

NOTES/CAVEATS/BUGS:
  assumes screen 24 lines long; this can be changed by command L or $LINES
  assumes screen 80 chars wide; this can be changed by command C or $COLUMNS
  lines longer than 80 chars are truncated and no LF is sent to display
  all control characters < Esc are interpreted as LF (tabs included!)
  does not work if data contain null characters
  may crash in extreme conditions (particularly if data of single fprintf
    call are longer than half the buffer size)
  if there is not enough memory, initscroll returns -1 and scroll simply
    does not scroll
  <stdarg.h> and ANSI standard required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

static struct {
  int L;      /* lines/screen (incl. prompt line) */
  int C;      /* no of columns */
  int ibuf;   /* _S_buf pointer */
  int buflen; /* length of _S_buf */
  int clear;  /* DOS legacy? */
  int stoponerr; /* ERROR(()), WARNING(()) don't, always stop */
#  ifdef CLOSEGRAPH
  int closegraph; /* whether to call closegraph() after quit */
#  endif /*# CLOSEGRAPH */
} S={25,80,0,0,1,0};

#  define NONE 0
#  define ANSI 1
#  define LINEFEEDS 2
#  define SPACES 3

static void _S_getscreen(int *size,char *c) /****************** _S_getscreen */
{
  if (c) *size=atoi(++c);
}

void setscroll(int lines,int columns,int clear,int closegraph) /** setscroll */
{
  if (lines) S.L=lines;
  if (columns) S.C=columns;
  S.clear=clear;
#  ifdef CLOSEGRAPH
  S.closegraph=closegraph;
#  endif /*# CLOSEGRAPH */
}

int initscroll(int buflen) /************************************* initscroll */
{
  char *LINES=getenv("LINES");
  char *COLUMNS=getenv("COLUMNS");

  if (in==NULL) in=stdin;
  if (out==NULL) out=stdout;

  if (buflen<=1023) {
    S.stoponerr=buflen;
    return 0; }

  if (LINES) S.L=atoi(LINES);
  if (COLUMNS) S.C=atoi(COLUMNS);

  buflen -= buflen & 1; /* must be even !!! */
  if (_S_buf) free(_S_buf);
  alloczero(_S_buf,buflen);
  /* if ( (_S_buf=(char*)malloc(buflen))==NULL) return -1; */
  /*  memset(_S_buf,0,buflen); */
  S.buflen=buflen;

  _S_buf[0]='\n';
  S.ibuf=1; /* to point after \n */
  _Id.echo=(in==stdin && out==stdout)+1;

  return 0;
}

static void _S_swapbuf(void) /*********************************** _S_swapbuf */
{
  int i,bh=S.buflen/2;
  char *b,c;

  if (S.ibuf>bh) {
    b=_S_buf+bh;
    loop (i,0,bh) { c=_S_buf[i]; _S_buf[i]=b[i]; b[i]=c; }
    S.ibuf-=bh; }
}

int prtc(char c) /***************************************************** prtc */
{
  int i;

  i=putc(c,out);
  if (out==stdout) if (_S_buf) {
    _S_buf[S.ibuf++]=c; _S_buf[S.ibuf]=0;
    _S_swapbuf(); }

  return i;
}

int prts_(const char *c) /********************************************* prts */
{
  int i;

  i=fputs(c,out);
  if (out==stdout) if (_S_buf) {
    while (*c) _S_buf[S.ibuf++]=*c++;
    _S_buf[S.ibuf]=0;
    _S_swapbuf(); }

  return i;
}

static int vprt_(const char *format,va_list args) /******************* vprt_ */
{
  int i;
  char *s;

  if (out==stdout && _S_buf) {
    i=vsprintf(s=_S_buf+S.ibuf,format,args);
    if (fputs(s,out)==EOF) i=EOF;
    else { i=0; while (*s) s++,i++; /* to return # of chars */ }
    while (_S_buf[S.ibuf]) S.ibuf++;
    _S_swapbuf(); }
  else
    i=vfprintf(out,format,args);

  return i;
}

void scroll(void) /************************************************** scroll */
{
  char x[STRLEN], *xnl;
#  define B(I) _S_buf[I+S.buflen*(I<0)]

  if (_S_buf) {
    static char s[STRLEN];
    static int l;
    FILE *f=NULL;
    int pos=S.ibuf-1,i,nl=0,key=0;

    if (S.buflen) if (_S_buf[pos]>=27) prtc('\n');
    for (;;) {
    again:
      fputs(key ? key>0 ? "bot$" : "top$" : "$",stdout);
      if (!fgets(x,STRLEN,stdin)) return;
      if ( (xnl=strchr(x,'\n')) ) *xnl=0; /* remove possible LF */

      switch (*x) {
        case 'i': masterin=in;
          if ( !(in=fopen(x+1,"rt")) ) {
            if (_Id.echo>=0)
              ERROR(("getdata: $i%s: no include file",x+1))
              in=masterin; }
        case '.':
        case 'q': return;
        case  0 : nl=1; break;
        case 'u': nl= -(S.L-2); break;
        case 'd': case ' ': nl=S.L-2; break;
        case '/': xnl=x+1; nl=strlen(xnl);
                  if (nl) { strcpy(s,xnl); l=nl; }
                  i=pos;
                  /* screen-line up to start search from 2nd line */
                  nl= -(S.L-2);
                  while (nl<0) {
                    while (i--,B(i)>=27);
                    if (!B(i)) { i++; break; }
                    nl++; }
                  for ( ; B(i); i++)
                    if (!memcmp(s,&B(i),l)) { pos=i; nl=S.L-1; break; }
                  break;
        case 'w': if (x[1]<=' ') { nl=0; break; }
                  f=fopen(x+1,"wt");
        case 't': case 'g': nl= -32000; break;
        case 'L': _S_getscreen(&S.L,x); Min(S.L,200) Max(S.L,4) break;
        case 'C': _S_getscreen(&S.C,x); Min(S.C,200) Max(S.C,10) break;
        case '!': _S_getscreen(&S.clear,x); break;
        case 'e': fputc('\n',out);
                  if (!strcmp(x,"exit"))
#  ifdef CLOSEGRAPH
                  if (S.closegraph) closegraph();
#  endif /*# CLOSEGRAPH */
                  exit(1);
        case 'b': case 'G': nl=32000; break;
        case '?': case 'h': fprintf(stderr,"\n\
%i lines x %i cols  clrscr=%d\n\n\
$t $g    top\n\
$b $G    bottom\n\
$<Enter> line down\n\
$-       line up\n\
$u       page up\n\
$d $<Spc>page down\n\
$#       by # lines\n\
$L#      set lines\n\
$C#      set cols\n\
$!#      set clrscr mode\n\
$wFILE   write\n\
$iFILE   include\n\
$/STRING find\n\
$/       last find\n\
$q $.    quit scroll\n\
$exit    exit prog\n\n\
",S.L,S.C,S.clear);
                  goto again;
        case '-': if (x[1]==0) { nl=-1; break; }
        default: nl=atoi(x); }

      key=0;
      /* by nl lines forward */
      while (nl>0) {
        while (pos++,B(pos)>=27);
        if (!B(pos)) { pos--; nl=0; key++; break; }
        nl--; }

      /* subtract screenlen to move cursor 1 screen up */
      nl-=(S.L-1);

      /* by nl lines backward */
      while (nl<0) {
        while (pos--,B(pos)>=27);
        if (!B(pos)) { pos++; key--; break; }
        nl++; }

      switch (S.clear) {
        case NONE:
          /* DOS legacy */
          break;
        case ANSI:
          printf("\n\33[2J\n");
          break;
        case SPACES:
          loop (i,0,S.C*S.L) putc(' ',stdout);
          putc('\n',stdout);
          break;
        case LINEFEEDS:
          loop (i,0,S.L) putc('\n',stdout); }

      if (f) { /* to file */
        i=pos;
        while (B(i)) {
          putc(B(i),f);
          i++; }
        fclose(f); f=NULL; }

      nl=S.L-1;

      while (nl>0) {
        i=0;
        while (pos++,B(pos)>=27) if (i++<S.C) putc(B(pos),stdout);
        if (i<S.C) putc('\n',stdout);
        if (!B(pos)) { pos--; break; }
        nl--; }
      } /* for (;;) */
    }
  else {
    /* batch mode: only $i solved here */
    if (!fgets(x,STRLEN,in)) ERROR(("getdata: cannot read in"))
    if ( (xnl=strchr(x,'\n')) ) *xnl=0; /* remove possible LF */
    if (*x=='i') {
      masterin=in;
      if ( !(in=fopen(x+1,"rt")) ) {
        if (_Id.echo>=0) ERROR(("getdata: $i%s: no include file",x+1))
        in=masterin; } }
    }
#  undef B
#  undef STRLEN
} /*scroll*/

#else /*# SCR */
/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% no scroll %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/
int initscroll(int buflen) /************************************* initscroll */
{
  if (in==NULL) in=stdin;
  if (out==NULL) out=stdout;
  if (buflen) _S_buf=(char*)1; /* to mark only */
  _Id.echo=(in==stdin && out==stdout)+1;

  return 0;
}

int prtc(char c) /***************************************************** prtc */
{
  return putc(c,out);
}

int prts_(const char *c) /******************************************** prts_ */
{
  return fputs(c,out);
}

static int vprt_(const char *format,va_list args) /******************* vprt_ */
{
  return vfprintf(out,format,args);
}

void scroll(void) /************************************************** scroll */
{
  char x[8];
  if (_S_buf) fgets(x,8,in);
}
#endif /*#!SCR */

/*** common for SCR and !SCR ***/

int prts(const char *c) /********************************************** prts */
{
  int i=prts_(c);

  if (i>=0) return prtc('\n');
  else return i;
}

static int vprt(const char *format,va_list args) /********************* vprt */
{
  int i,j;

  i=vprt_(format,args);
  j=prtc('\n');

  return j<0 ? j : i;
}

int prt_(const char *format,...) /************************************* prt_ */
{
  int ret;
  va_list args;

  va_start(args,format);
  ret=vprt_(format,args);
  va_end(args);

  return ret;
}

int prt(const char *format,...) /*************************************** prt */
{
  int ret;
  va_list args;

  va_start(args,format);
  ret=vprt(format,args);
  va_end(args);

  return ret;
}


/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    Error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/* for macros ERROR DISASTER WARNING */

struct myErrorInfo_s myErrorInfo;

void myError(const char *format, ...) /***************************** myError */
{
  va_list args;
  const char *errfmt="\n*** %s in %s:%d\n*** ";
  int ierr=20;

  if (myErrorInfo.Level<0 || myErrorInfo.Level>2) {
    fprintf(stderr,"\n?! TOTAL BREAKDOWN: unknown errorlevel");
    myErrorInfo.Level=2; }

  if (out!=stdout
#ifdef SCR
                  || S.stoponerr
#endif /*# SCR */
                                  ) {
    fprintf(stderr,errfmt,
            "WARNING\0ERROR\0  DISASTER"+myErrorInfo.Level*8,
            myErrorInfo.File,myErrorInfo.Line);
    va_start(args,format);
    vfprintf(stderr,format,args);
    va_end(args);
    fputc('\n',stderr);
  }

  prt_(errfmt,
       "WARNING\0ERROR\0  DISASTER"+myErrorInfo.Level*8,
       myErrorInfo.File,myErrorInfo.Line);
  va_start(args,format);
  vprt(format,args);
  va_end(args);
  _n

  if (
#ifdef SCR
      !S.stoponerr &&
#endif /*# SCR */
                      in==stdin && out==stdout) {
    char s[32];

  prompt:
#ifdef SCR
    fprintf(stderr,"*** select:  (a)bort  (e)xit%s%s\n",
      myErrorInfo.Level==0?"  (c)ontinue":myErrorInfo.Level==1?"  (i)gnore":"",
      S.buflen>1?"  (s)croll":"");
#else  /*# SCR */
    fprintf(stderr,"*** select:  (a)bort  (e)xit  %s\n",
      myErrorInfo.Level==0?"  (c)ontinue":myErrorInfo.Level==1?"  (i)gnore":"");
#endif  /*#!SCR */
    if (fgets(s,32,in)==NULL) { prt("*** EOF ***"); s[0]='e'; }
    if (!ierr--) s[0]='e';
    switch (s[0]) {
      case 'a': fclose(out);
#ifdef CLOSEGRAPH
                closegraph();
#endif /*# CLOSEGRAPH */
                abort();
      case 'e': fclose(out);
#ifdef CLOSEGRAPH
                closegraph();
#endif /*# CLOSEGRAPH */
                exit(myErrorInfo.Level);
      case 'c': if (myErrorInfo.Level) goto prompt;
                prt("*** continuing ***"); return;
      case 'i': if (myErrorInfo.Level==2) goto prompt;
                prt("*** ignoring error ***"); return;
#ifdef SCR
      case 's': case '$': scroll();
#endif /*# SCR */
      default : goto prompt; }
    }
  else if (myErrorInfo.Level){
    fclose(out);
    exit(myErrorInfo.Level);  }
}

/* legacy version, kept for compatibility reasons (and PAR): */

void extError(const char *msg, const char *f, int l) /************* extError */
{
  static char e[]="*** %s:%i ERROR %s ***\n";

  fprintf(stderr,e,f,l,msg);

  if (in) { /* added 4/2006 */
    prtc('\n');
    prt_(e,f,l,msg); }
  if (in==stdin && out==stdout) {
    char s[32];
  prompt:
#ifdef SCR
    fprintf(stderr,"select:  (a)bort  (e)xit  (r)esume  (s)croll\n");
#else /*# SCR */
    fprintf(stderr,"select:  (a)bort  (e)xit  (r)esume\n");
#endif /*#!SCR */
    if (fgets(s,32,in)==NULL) {
      prt("*** EOF ***");
      s[0]='e'; }
    switch (*s) {
      case 'a': fclose(out);
#ifdef CLOSEGRAPH
                closegraph();
#endif /*# CLOSEGRAPH */
                abort();
      case 'e': fclose(out);
#ifdef CLOSEGRAPH
                closegraph();
#endif /*# CLOSEGRAPH */
                exit(msg[0]+(msg[0]==0));
      case 'r': prt("*** execution resumed ***");
                return;
#ifdef SCR
      case 's': case '$': scroll();
#endif /*# SCR */
      default : goto prompt; } }
  else {
    fclose(out);
#ifdef CLOSEGRAPH
    closegraph();
#endif /*# CLOSEGRAPH */
    exit(msg[0]+(msg[0]==0)); }
}

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           C-version of the PL/I GET DATA and PUT DATA statements            %
%                            (c) J.Kolafa 1991                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
* input stream should contain a list of assignments of the form
    identifier=expression
  or
    identifier[decimal_const]suffix=expression
  ended by a semicolon (suffix may be empty)
* identifiers must begin with 'A'..'Z', 'a'..'z' ( not '_' )
* input is read from file `in'
* output is written to file `out'
*   `out' and `in' must be defined by the user or by initscroll()
* use in the source must start by `getdata' and end by `enddata'
* macro get(S) reads scalar S ( statement in input stream is S=expression )
* getvec(V,S,maxindex) reads one item from a vector of structures
  ( statement in input stream is V[index]S=expression )
  ( use getvec(V,,maxindex) for a simple vector, V[index]=expression )
* checkdata checks for misspelled identifiers, to be used just before enddata
* statement `?' single on line prints a brief help
* `?=' toggles echo mode (auto print of expressions results)
* `??' prints a list of all identifiers
* #ifdef SCR then using $ in the input stream enables scrolling
  $ must be followed by a scroll command (single letter or number)
* any text following '!' until end of line is discarded (treated as a comment)
* line beginning by * is copied and ignored (comment, or WARNING passed)

* example of use

  #include "ground.h"

  double dd=7.7;
  double z[3]={1.1,2.2,3.3};
  double abc[2];

  int main()
  {
    int i=11, J3=33;
    double xx;
    int xxset=0;
    struct { int i; double d; } s[10];
    struct keys_s initkeys[]={
      {"cont",0},  {"continue",0},
      {"append",1},
      {"start",2},
      {"random",3},
      {NULL,0}};

    initscroll(0); // selects in=stdin, out=stdout, echo on, no scrolling

    getdata
      get(i) get(J3) get(abc[0]) get(abc[1])
      getvec(z,,3) getvec(s,.i,10) getvec(s,.d,10)
      getmark(dd,ddset)
    checkdata enddata
    if (!ddset) ERROR(("dd has not appeared in input data"))

    put(dd)
    put2(i,sqrt((double)J3))
    put3(z[0]*2,z[1]*3,z[2]*4)
    put2(abc[0],abc[1])
    put2(s[3].i,s[3].d)

    return 0;
  }

* example of input data

  dd=acos(-1)^2 J3=123.4^2
  i=2
  i = #*3*4*5*6
  z[1] = pi/3
  s[3].i=123 s[3].d=(1+1/(2+3))^.5
  abc[1]=1;

* will generate the following output

i=2
  i = 2
i = #*3*4*5*6
  i = 720
z[1] = pi/3
  z[1] = 1.04719755
s[3].i=123 s[3].d=(1+1/(2+3))^.5
  s[3].i = 123
  s[3].d = 1.09544512
abc[1]=1;
  abc[1] = 1

         dd=7.7
          i=720           sqrt((double)J3)=5.74456
     z[0]*2=2.2                z[1]*3=3.14159            z[2]*4=13.2
     abc[0]=0                  abc[1]=1
     s[3].i=123                s[3].d=1.09545

* NOTES/CAVEATS/BUGS:
- Indexed variables cannot be used in expressions; however, the following
  workaround may be used:
    L[1]
    x=# ! now x=L[1]
- One assignment statement in the input stream must be within one line of
  maximum length GETBUFLEN
- In the case of more dimensions the varying index must be in the first []
- Error messages appear after the whole line is read and (if _Id.echo==1)
  copied to out - the position of the error is not marked
- An expression is evaluated in double, the result is undefined if it cannot
  be cast into the type requested
- If the variable is int, result printed as float
- pi, PI, and standard function names cannot be used as identifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "calc.c"

const char *_putformat="%11s=%-13.6g ";

/* this structure must be global - used in getdata macros */
struct _Id_s _Id = { " = %-16.9g", "", 1.0, 0};

static void prtecho(void) /***************************************** prtecho */
{
  prt("<echo is %s>",_Id.echo?_Id.echo==1?"copy input":"on":"off");
}

int _GetId(void) /*************************************************** _GetId */
/*
  gets one statement from `in'
  reads a whole line (until `;') if necessary
  scans lot of statements
  returns 1 if get() statements are to be scanned, otherwise 0
*/
{
  int c,sg,excl=0,hack;
  char *x,*end;
  static char zero=0;

  double lhs=_Id.nr;

  _Id.str=NULL;

  /* 1st pass just to scan all scalar variables to create a list with values */
  if (_Id.list) return 1;

  /* += -= /= *= ^= %= pending */
  _Id.used=0;
  if (_Id.asg>1) goto evalexpr;

  /* reads one line (until `;') if necessary */
  if (*_Id.b==0) {
    x=_Id.buf;
    hack=0;
    do {
      c=getc(in);
      hack++;
      if (hack==1 && c=='*')
        /* line starting with * is copied and skipped 
           (hack to handle *** WARNINGs, *** ERRORs in ble-files) */
        for (;;) {
          prtc(c);
          c=getc(in);
          if (c<0) break;
          if (c=='\n') { prtc(c); goto eoln; } }
      if (c<0) {
        /* EOF found: swap back to master `in' if this is EOF on file
           included by $iFILE (this is handled by scroll() */
        fclose(in);
        in=masterin; masterin=NULL;
        if (in==NULL) {
          fprintf(stderr,"EOF on in - exit\n");
          if (out) prt("EOF on in - exit");
          exit(0); }
        c='\n'; }
#ifdef SCR
      while (c=='$') {
        scroll();
        c=getc(in); }
#endif /*# SCR */
      if (_Id.echo==1) if (c>0) prtc(c);
      if (c<' ') { if (c=='\n') goto eoln; c=' '; }
      if (c=='!') excl=1; /* comment: skip to eoln */
      if (!excl) *x++=c;
      if (x-_Id.buf>=GETBUFLEN-1) {
        ERROR(("getdata: too long input line\n\
*** 4 bytes in HEX: %08x",*(unsigned int*)_Id.b))
        *(_Id.b=_Id.buf)=0;
        return 0; }
    } while (c!=';' || excl); /* `;' after `!' ignored - not end getdata set */
    eoln:
    *x++=0;
    *x=0;
    _Id.b=_Id.buf; }

/* remove trailing spaces */
  while (*_Id.b==' ') _Id.b++;

  _Id.prt=_Id.asg=0;
  _Id.id=_Id.b;
  _Id.suff= &zero;
  _Id.br=NULL;

#include "question.c"

  /* empty line or `;' */
  if ((_Id.empty = (*_Id.id==';' || *_Id.id==0))) {
    noscan:
    _Id.id= &zero;
    _Id.used=1;    /* suppress check `identifier found' */
    goto ret; }

  if (_Id.echo==2) _Id.prt=1; /* ID name to be printed */


  /*** identifier scan, incl. # and indexed id's ***/

  if (isalpha(*_Id.b))
    while (*_Id.b && (isalnum(*_Id.b) || strchr("_.[",*_Id.b))) {
      if (!_Id.br && *_Id.b=='[') {
        /* scanning [INDEX] (the 1st []) */
        _Id.br=_Id.b++; sg=1; _Id.indx=0;
        if (*_Id.b=='-') { sg= -1; _Id.b++; }
        while (isdigit(*_Id.b)) { _Id.indx=_Id.indx*10+sg*(*_Id.b++-'0'); }
        _Id.suff=_Id.b;
        if (*_Id.b!=']') { 
          ERROR(("getdata: missing ] in \"line\":\n\
\"%s\"",_Id.buf))
          goto next; }
        _Id.suff=_Id.b+1; }
      next: _Id.b++; }

  while (*_Id.b==' ') *_Id.b++=0;

  /* postponing operators += -= *= ... */
  if (strchr("+-*/%^|&$",*_Id.b))
    if (*(_Id.b+1)=='=') {
      _Id.asg=*_Id.b; *_Id.b=0; _Id.b+=2;
      goto ret; }

  if (_Id.b==_Id.id || *_Id.b!='=') {
    /* anything else than `=' after ID : ID value or expression printed */
    _Id.prt=1;  /* to print ID name */
    _Id.used=0; /* to check if ID is OK */
    _Id.asg=1;  /* to assign to _Id.nr (= # in expr), CALC>=2 only */
    if (_Id.b==_Id.id
     || (*_Id.b && strchr("+-*/^%:()",*_Id.b))
     || !strcmp("pi",_Id.id) || !strcmp("PI",_Id.id)  ) {
      /* expression to evaluate and print */
      /* if (_Id.br) Error("indexed variable in expression"); */
      _Id.used=2; /* suppress check + mark not id=expr */
      _Id.b=_Id.id;
      goto evalexpr; }
    goto ret; }

  /* there was = (or += +- ...) after ID */
  *_Id.b++=0; /* mark string end (replace `=' etc. by zero) */

  /*** expression evaluation ***/
 evalexpr:

  while (*_Id.b==' ' || *_Id.b=='\t') _Id.b++;

  if (*_Id.b=='\"') {
    /* NEW: support of ID="string" -- _Id.str exported instead of _Id.nr */
    _Id.str=end=++_Id.b;
    _Id.nr=_Id.str[0]; /* ascii of 1st char */
    _Id.used=0; /* force error */
    while (*end) {
      if (end-_Id.b>64) 
        ERROR(("getdata: missing right \" or string too long near \"strings\":\n\
\"%s\"\n\
\"%s\"",_Id.buf,_Id.b))
      if (*end=='\"') { *end++=0; break; }
      end++; } }
  else {
    _Id.str=NULL;
    _Id.nr=Calc(_Id.b,&end);

    /* += -= /= *= solved */
    switch (_Id.asg) {
      case '+': _Id.nr+=lhs; break;
      case '*': _Id.nr*=lhs; break;
      case '^': _Id.nr=pow(lhs,_Id.nr); break;
      case '-': _Id.nr=lhs-_Id.nr; break;
      case '/': _Id.nr=lhs/_Id.nr; break;
      case '%': _Id.nr=lhs-(int)(lhs/_Id.nr)*_Id.nr; break; }

    _Id.asg=0; }

  if (end==_Id.b) {
    if (_Id.echo>=0)
      ERROR(("getdata: bad number or expression at \"id\":\n\
\"%s\"\n\
*** offending \"string\":\n\
\"%s\"\n\
*** 4 bytes in HEX: %08x %s",
             _Id.buf,
             _Id.b,*(unsigned int*)_Id.b,
             (unsigned char)_Id.b[0]&128?" (ILLEGAL NON-ASCII CHARACTER - check coding)":""))
    *(_Id.b=_Id.buf)=0;
    return 0; }

  _Id.b=end;

 ret:
  _Id.end= *_Id.b==';';

  return 1;
}

static struct keys_s *passkey;
static int errorkey;

static int findkey(void) /****************************************** findkey */
{
  struct keys_s *key=passkey;

  errorkey=0;
  while (key->key) {
    if (!strcmp(_Id.str,key->key)) return key->val;
    key++; }
  errorkey++;

  return 0;
}

int _GetCmp(char *id,double val) /********************************** _GetCmp */
{
  _Idlist *l;

  if (_Id.list) {
    l=(_Idlist *)malloc(sizeof(_Idlist)+strlen(id));
    alloc(l,sizeof(_Idlist)+strlen(id));
    /* =(_Idlist *)malloc(sizeof(_Idlist)+strlen(id)); */
    if (!l) {
      fprintf(stderr,"no heap for _Id.list\n");
#ifdef CLOSEGRAPH
      closegraph();
#endif /*# CLOSEGRAPH */
      exit(-2); }
    l->next=_Id.head; _Id.head=l;
    strcpy(l->id,id);
    l->val=val;
    return 0; }

  if (_Id.prt>=2) {
    int x=strlen(id)+1;
    if ((_Id.col+=x)>79) { _Id.col=x; prtc('\n'); }
    prt_("%s ",id);
    return 0; }

  if (_Id.br!=NULL) *_Id.br='[';
  if (strcmp(id,_Id.id)) return 0;

  if (_Id.str) {
    if (passkey) {
      _Id.nr=findkey();
      if (errorkey) ERROR(("getdata: unknown string value in expression:\n\
%s=\"%s\"",id,_Id.str)) }
    else
      ERROR(("getdata: near \"id\":\n\
\"%s\"\n\
*** variable %s accepts only numbers, not strings",_Id.str,id)) }

  if (_Id.asg) _Id.nr=val;
  _Id.used=1;
  if (_Id.prt) prt_("  %s",id);

  /* try to find the identifier in the list to update value */
  for (l=_Id.head; l; l=l->next)
    if (!strcmp(l->id,id)) l->val=_Id.nr;

  return 1;
}

int _GetVecCmp(char *id,char *xsuff,int maxindex,double val) /*** _GetVecCmp */
{
  /* this patch is needed for compilers which make string " " and not ""
     from an empty macro argument (e.g. Cray) */
  char suff[32];
  strcpy(suff,xsuff);
  if (*suff==' ') *suff=0; /* replace suffix " " by "" */

  /* nothing instead of list creating - cannot use indexed variables in expr. */
  if (_Id.list) return 0;

  if (_Id.prt>=2) {
    int x=strlen(id)+strlen(suff)+8;
    if (maxindex>999) x+=2;
    if ((_Id.col+=x)>79) { _Id.col=x; prtc('\n'); }
    prt_("%s[0..%d]%s ",id,maxindex-1,suff); return 0; }

  if (_Id.br==NULL) return 0;
  *_Id.br=0;
  if (strcmp(id,_Id.id) || strcmp(suff,_Id.suff)) return 0;
  if (_Id.asg) _Id.nr=val;
  _Id.used=1;
  if (_Id.prt) prt_("  %s[%i]%s",id,_Id.indx,suff);
  if (_Id.indx>=maxindex || _Id.indx<0) {
    ERROR(("getdata: index in %s[%d] out of range",id,_Id.indx))
    return 0; }

  return 1;
}

void _CheckData(void) /****************************************** _CheckData */
{
  if (!_Id.used)
    if (!_Id.list) {
      if (_Id.echo>=0) {
        if (_Id.b==NULL || _Id.b[0]==0)
          ERROR(("getdata: bad identifier %sat \"%s\"",
                 strpbrk(_Id.buf,"[]")?"or unexpected [] ":"",
                 _Id.buf))
        else
          ERROR(("getdata: bad identifier or unexpected [] at \"%s\"\n\
*** or syntax error at \"%s\" (HEX=%08x%s)",
                 _Id.buf,
                 _Id.b,*(unsigned int*)_Id.b,
                 (unsigned char)_Id.b[0]&128?": NON-ASCII CHARACTER - check coding":"")) }
        *(_Id.b=_Id.buf)=0; }
}

int _GetEnd(void) /************************************************* _GetEnd */
{
  if (_Id.list) {
    /* was just after creating the list of identifiers */
    _Id.list=0;
    return 0; }

  if (_Id.prt) {
    if (_Id.prt>=2) prtc('\n');
    else prt(_Id.fmt,_Id.nr); }

  if (_Id.end) {
    _Idlist *Q,*QQ;

    if (_Id.echo==2) prtc('\n');

    /* cleaning the list */

    QQ=_Id.head;
    while (QQ) {
      Q=QQ; QQ=QQ->next; free(Q); }
    _Id.head=NULL; }

  return _Id.end;
}

int _GetKey(char *id,double val,struct keys_s *key) /*************** _GetKey */
{
  int ret;

  passkey=key;
  ret=_GetCmp(id,val);
  passkey=NULL;

  return ret;
}

int prtkey_(int val,struct keys_s *key,int how) /******************* prtkey_ */
/*
  prints string ="NAME" according to val and key
  in case of several NAMEs with the same val:
    how=0:  prints the first occurrence in key
    how=1:  prints the last occurrence in key
    how=2:  prints all occurrences in key
  returns the number of characters printed
*/
{
  char *found=NULL;
  int i=0;

  while (key->key) {
    if (val==key->val) {
      if (how==1) found=key->key;
      else i+=prt_("=\"%s\"",key->key);
      if (how==0) break; }
    key++; }

  if (found) i+=prt_("=\"%s\"",found);

  return i;
}

void _PrtKey_(char *id,int val,struct keys_s *key) /*************** _PrtKey_ */
/*
  for putkey(id): id=VALUE="NAME" (1st occurrence) is printed
  filled to 26 chars
*/
{
  int i=prt_("%11s=%d",id,val);

  prtkey_(val,key,0);

  while (i<26) prtc(' '),i++;
}

/****** headers of tables ******/

void putline(char c,int l) /**************************************** putline */
{
  int i;
  /*** prints a line of l characters c on out ***/
  /* new: negative l means no LF */
  if (l) loop (i,0,abs(l)) prtc(c);
  if (l>0) _n
}

void header(char *h) /*********************************************** header */
/***
  to make headers of tables. Usage:
    header("AA BB ");
    ...
    header("");
  will produce
  ======
  AA BB
  ------
  ...
  ======
  on out
  if the 1st character of the arg of header is special, it is used in front of
  hlines
  new 8/99: h may contain \n
*/
{
  static int len=0;
  static char c1,c2;
  char *c;

  if (h && h[0]) {
    if (strchr("!@#$%^&*+[]:<>()./\\",h[0])) c1=c2=h[0];
    else { c1='='; c2='-'; }

    len=0;
    for (c=h; *c; c++) if (*c<' ') len=0; else len++;
    _n
    putline(c1,len);
    prt_("%s\n%c",h,c2);
    putline(c2,len-1); }
  else
    putline(c1,len);
}

void underline(char *t) /***************************************** underline */
/*
 produces

title
^^^^^

in BLEND, prepended by !
*/
{
#ifdef BLEND
  prt("\n! %s",t);
  prt_("! "); putline('^',strlen(t));
#else  /*# BLEND */
  prt("\n%s",t);
  putline('^',strlen(t));
#endif  /*#!BLEND */
}

#define MAXGRAPHLEN 125

void graph(double x,int linelen) /************************************ graph */
/*
  Pseudo-graphics with a char/4 resolution.
  For 0<=x<=1, a multi-character is printed at the position scaled by linelen.
*/
{
  int c,i;
  static char a[]="X X>><<X";
  static char l[MAXGRAPHLEN+3];

  Min(linelen,MAXGRAPHLEN)

/*.....char *l; alloc(l,linelen+3);*/

  memset(l,' ',linelen+2);
  l[1]=l[linelen]='|';
  if (linelen&1) l[linelen/2+1]='|';
  c=(int)((linelen-1)*4*x+4.5);
  i=(c+256)%4; c=(c-i)/4;
  if (c>=0 && c<=linelen) memcpy(l+c,a+2*i,2);
  l[linelen+1]=0;
  prt("%s",l+1);
/*.....free(l);*/
}

#ifdef MYSQRT
/* OBSOLETE: optional sqrt function replacement, see ground.h */

unsigned sqrt_tab[64] = {
   5376U,  12024U,  19815U,  27394U,  34750U,  41877U,  48762U,  55398U,
  61770U,  67867U,  73677U,  79186U,  84377U,  89240U,  93748U,  97887U,
 101635U, 104961U, 107845U, 110218U, 112068U, 113409U, 114654U, 115341U,
 115422U, 114907U, 113267U, 110876U, 107701U, 103741U,  98918U,  93152U,
  92222U,  96728U, 100941U, 104850U, 108425U, 111685U, 114529U, 117129U,
 119149U, 120699U, 122335U, 123565U, 124351U, 124648U, 124500U, 123853U,
 122666U, 120373U, 117526U, 114220U, 110158U, 105438U,  99982U,  93742U,
  86655U,  78666U,  69713U,  59729U,  48641U,  36365U,  22820U,   7901U };

#  if MYSQRT==1
/* function */
double mySqrt(double _X) /******************************************* mySqrt */
{
  long unsigned sqrt_u;
  double sqrt_t, sqrt_x;
  union { long unsigned u; double d; } sqrt_ud;

  sqrt_ud.d = _X;
  sqrt_u = 0x5fe80000UL-(sqrt_ud.u>>33);
  sqrt_ud.u = sqrt_u-(long unsigned)sqrt_tab[(sqrt_u>>14) & 63UL]<<32;
  sqrt_x = sqrt_ud.d;
  sqrt_x *= 3.0-_X*sqrt_x*sqrt_x;
  sqrt_x *= 12.0-_X*sqrt_x*sqrt_x;
  sqrt_x *= 0.0625;
  sqrt_t = _X*sqrt_x;

  return sqrt_t+sqrt_x*0.5*(_X-sqrt_t*sqrt_t);
} /* end mySqrt */

#  else /*# MYSQRT==1 */
/* macro (inline), see ground.h */
long unsigned sqrt_u;
double sqrt_t, sqrt_x;
union sqrt_un sqrt_ud;
#  endif /*#!MYSQRT==1 */
#endif /*# MYSQRT */

#include "fortran.c"
#include "rndseed.c"
#include "rndgeni.c"

/* original alloc.c merged 4/2016 ========================================== */
/* see ground.h for help */

/* my alloc/free needed: */
#include "mystring.c"
#include "mygets.c"

/* changes 11/2009:
   FNLEN increased (from 28 to 40)
   output pointer hex only
   5/2012: small bug fixed - no trailing 0 in id
*/

/* WARNING: original system free must be valid here, do not use my alloc/free */
#undef free

int4 AllocSizeLim=0x40000000; /* 1 GiB 6/2012; NB: to be 64 bit, not int4 */

#if CHECKHEAP > 2 || CHECKHEAP < -2
#  error CHECKHEAP out of range
#endif /*# CHECKHEAP > 2 || CHECKHEAP < -2 */

/* SHM removed - see shmgroundinclude.c, sys4par.h removed */

/***** alloc/free *****/

#if CHECKHEAP == 0 /* ==================================================== 0 */

/***** no check *****/

void *myfree(void *f) /********************************************** myfree */
{
#  if FREENULL
  if (f)
#  endif /*# FREENULL */
    free(f);
  return NULL;
} /* myfree */

#elif CHECKHEAP == 1 || CHECKHEAP == -1 /* ========================== +1,-1 */

/***** heap overflow and bad calls to free are checked *****/

void *mymalloc(int4 size,int zero) /******************************* mymalloc */
{
  void *a;

  if (size>AllocSizeLim) {
    Error("AllocSizeLim exceeded");
    return NULL; }
  if ((a=(void*)malloc(size)) == NULL) {
    Error("malloc failed");
    return NULL; }
  if (zero) memset(a,0,size);

  return a;
} /* mymalloc */

void *myfree(void *f) /********************************************** myfree */
{
  if (f) { free(f); f = NULL; }
#  if !FREENULL
  else Error("free(NULL)");
#  endif /*# !FREENULL */
  return f;
} /* myfree */


#elif CHECKHEAP == 2 || CHECKHEAP == -2 /* ========================== +2,-2 */

/* IDLEN-1 chars of id's recorded */
#  define IDLEN 12
#  define IDFMT "%12s"
#  define FNLEN 40

int4 AllocRange, AllocTrace;

static unsigned char AllocRangeData[16]={
  85,122,49,97,165,202,35,27,109,181,59,65,206,108,162,170};

struct AllocRange_s {
  unsigned char check[16];
  void *p;
  struct AllocRange_s *next;
  int4 size;
  int line;
  char fn[FNLEN];
  char id[1]; /* variable length */
} *Range0;

static int AllocCount=0;

static struct AllocRange_s *staticR;

static void RError(const char *format,...) /************************* RError */
{
  va_list args;
  char enough[256];

  sprintf(enough,
    "RANGE CHECK for `%s':\n*** allocated with %ld B in %s:%d\n*** ",
    staticR->id,(long)staticR->size,staticR->fn,staticR->line);

  va_start(args,format);
  vsprintf(enough+strlen(enough),format,args);
  myError(enough);
  va_end(args);
}

static void checkone(struct AllocRange_s *R) /********************* checkone */
{
  if (!AllocRange)
    return;
  else {
    unsigned char *head,*tail;
    int4 size;
    unsigned int i;

    staticR=R;

    if (AllocRange<0 || AllocRange&15) {
      RError("AllocRange smashed");
      AllocRange=16; }

    head=(unsigned char*)R->p-(AllocRange+16);
    size=*(int4*)(head+8);
    if (size<1 || size>AllocSizeLim) {
      RError("size overwritten");
      return; }
    if (*(int4*)(head+12)!=AllocRange) {
      RError("write %ld B before",AllocRange+8L);
      *(int4*)(head+12)=AllocRange; }
    head+=16;
    loop (i,0,AllocRange) if (head[i]!=AllocRangeData[i&15]) {
      RError("write %ld B before",(long)AllocRange-i);
      head[i]=AllocRangeData[i&15]; }
    tail=head+size+AllocRange;
    loop (i,0,AllocRange) if (tail[i]!=AllocRangeData[i&15]) {
      RError("write %ld B after",(long)i+1L);
      tail[i]=AllocRangeData[i&15]; } }
}

void mycheckranges(int verbose,char *fn,int line) /*********** mycheckranges */
{
  struct AllocRange_s *R;
  int n=0;

  myErrorInfo.Level=1; myErrorInfo.Line=line; myErrorInfo.File=fn;

  if (verbose)
    prt("ALLOC and RANGE check called from %s:%d %s",
        fn,line,Range0?"":"- empty!");

  for (R=Range0; R; R=R->next) {
    staticR=R;
    if (verbose)
      prt("%3d " IDFMT " %08p %lu B by %s:%d",
         ++n,R->id,R->p,R->size,R->fn,R->line);
    if (memcmp(R->check,AllocRangeData,16))
      RError("bad control data");
    checkone(R); }
}

void *mymalloc(int4 size,int zero,char *c,char *fn,int line) /***** mymalloc */
/* note: zero&2 marks call from ralloc */
{
  void *a;

  myErrorInfo.Level=1; myErrorInfo.Line=line; myErrorInfo.File=fn;

  if (size<1 || size>AllocSizeLim)
    myError("alloc(%s,%ld): size out of range",c,(long)size);

  if (AllocRange) {
    unsigned char *head,*tail;
    unsigned int i;
    struct AllocRange_s *R;

    /* 12/2004: change AllocRange smoothly to a multiple of 16: */
    AllocRange=(AllocRange+15)/16*16;

    mycheckranges(0,fn,line);

    head=(unsigned char*)malloc(size+2*AllocRange+16);
    if (!head)
      myError("alloc(%s,%ld): no heap",c,(long)size);
    *(int4*)(head+8)=size;
    *(int4*)(head+12)=AllocRange;
    head+=16;
    tail=head+size+AllocRange;
    loop (i,0,AllocRange) head[i]=tail[i]=AllocRangeData[i&15];
    a=(void*)(head+AllocRange);
    R=malloc(sizeof(struct AllocRange_s)+strlen(c));
    if (!R)
      myError("alloc(%s,%ld): no heap (check)",c,(long)size);
    else {
      R->next=Range0;
      R->p=a;
      R->size=size;
      copy(R->check,AllocRangeData,16);
      strcpy(R->id,c);
      /* copy only last FNLEN-1 chars of fn if longer than FNLEN */
      strcpy(R->fn,fn+max((int)strlen(fn)-(FNLEN-1),0));
      R->line=line;
      Range0=R; }
    }
  else
    if ((a=(void*)malloc(size)) == NULL)
       myError("alloc(%s,%ld): no heap",c,(long)size);

  if (zero&1) memset(a,0,size);

  if (AllocTrace)
    prt("%3d:" IDFMT " %08p %calloc %luB%s in %s:%d",
         ++AllocCount,
         c,
         a,
         zero&2?'r':' ',
         (long)size,
         zero?":=0":"",
         fn,line);
  return a;
} /* mymalloc */

char *mystrdup(const char *s,char *c,char *fn,int line) /********** mystrdup */
{
  char *a=mymalloc(strlen(s)+1,0,c,fn,line);
  strcpy(a,s);
  return a;
}

void *myfree(void *f,char *c,char *fn,int line) /******************** myfree */
{
  myErrorInfo.Level=1; myErrorInfo.Line=line; myErrorInfo.File=fn;

  if (f == NULL) {
#  if FREENULL
    prt("free(%s==NULL) in %s:%d",c,fn,line);
#  else /*# FREENULL */
    myError("free(%s==NULL)",c);
#  endif /*#!FREENULL */
    return NULL; }

  if (AllocRange) {
    struct AllocRange_s *R,**C=NULL;

#  if 0
    /* check of the memory being freed first */
    staticR=malloc(sizeof(struct AllocRange_s)+strlen(c));
    if (!staticR) myError("no heap for control data");
    strcpy(staticR->id,c);
    staticR->p=*f;
    /* statement that has allocated it is not known */
    strcpy(staticR->fn,"?");
    staticR->line=-1;
    checkone(staticR);
    free(staticR);
#  endif /*# 0 */

    mycheckranges(0,fn,line);

    for (R=Range0; R; R=R->next)
      if (f==R->p) {
        if (C) *C=R->next;
        else Range0=R->next;
        if (AllocTrace)
	  prt("%3i:" IDFMT " %08p as %s in %s:%d freed",
	      AllocCount--,c,f,
	      R->id,R->fn,R->line);
        free(R);
        break; }
      else
        C=&R->next;
    free((char*)f-AllocRange-16); }
  else
    free(f);

  if (AllocTrace & !AllocRange)
    prt("%3i:" IDFMT " %08p freed",
         AllocCount--,c,f);

  return NULL;
}

#endif /*#!CHECKHEAP == 0!CHECKHEAP == 1 || CHECKHEAP == -1!CHECKHEAP == 2 || CHECKHEAP == -2 */

/***** ralloc/release *****/

#if CHECKHEAP < 0

#  ifndef RHEAPTABLEN
#    define RHEAPTABLEN 32
#  endif /*# RHEAPTABLEN */

typedef struct {
  void *p;
#  if CHECKHEAP == -2
  char id[IDLEN];
#  endif /*# CHECKHEAP == -2 */
} rheap_item;

typedef struct rheap_s {
  rheap_item i[RHEAPTABLEN];
  int pos;
  struct rheap_s *next;
} rheap_t;

rheap_t *rheap;

#  if CHECKHEAP == -1

void *myralloc(int4 size,int zero) /******************************* myralloc */
{
  rheap_t *h;
  void *a;

  if (rheap==NULL || rheap->pos>=RHEAPTABLEN) {
    h=rheap;
    alloc(rheap,sizeof(rheap_t));
    rheap->next=h;
    rheap->pos=0; }

  if ((a=(void*)malloc(size)) == NULL) Error("ralloc");
  if (zero) memset(a,0,size);

  rheap->i[rheap->pos++].p=a;
  
  return a;
}

void *myrelease(const void *p) /********************************** myrelease */
{
  void *x;
  rheap_t *h;

  while (rheap!=NULL) {
    x=rheap->i[--rheap->pos].p;
    rheap->i[rheap->pos].p=myfree(rheap->i[rheap->pos].p);
    if (rheap->pos <= 0) {
      h=rheap; rheap=rheap->next;
      h=TYPEOF(h) myfree(h); }
    if (x==p) return NULL; }

  Error("release");
  return NULL;
}

#  else /* CHECKHEAP = -2 */ /*# CHECKHEAP == -1 */

void *myralloc(int4 size,int zero,char *c,char *fn,int line) /***** myralloc */
{
  rheap_t *h;
  void *a;

  myErrorInfo.Level=1; myErrorInfo.Line=line; myErrorInfo.File=fn;

  if (rheap==NULL || rheap->pos>=RHEAPTABLEN) {
    h=rheap;
    alloc(rheap,sizeof(rheap_t));
    rheap->next=h;
    rheap->pos=0; }

  a=mymalloc(size,zero|2,c,fn,line);

  if (AllocTrace) {
    copy(rheap->i[rheap->pos].id,c,IDLEN);
    rheap->i[rheap->pos].id[IDLEN-1]=0; }

  rheap->i[rheap->pos++].p=a;

  return a;
}

void *myrelease(void *p,char *c,char *fn,int line) /************** myrelease */
{
  void *x;
  rheap_t *h;

  myErrorInfo.Level=1; myErrorInfo.Line=line; myErrorInfo.File=fn;

  while (rheap!=NULL) {
    x=rheap->i[--rheap->pos].p;

    if (AllocTrace)
      rheap->i[rheap->pos].p=myfree(rheap->i[rheap->pos].p,rheap->i[rheap->pos].id,fn,line);
    else
      rheap->i[rheap->pos].p=myfree(rheap->i[rheap->pos].p,c,fn,line);

    if (rheap->pos <= 0) {
      h=rheap; rheap=rheap->next;
      h=myfree(h,"rheap",fn,line); }
    if (x==p) return NULL; }

  myError("release %s",c);
  return NULL;
}

#  endif /*#!CHECKHEAP == -1 */
#endif /*# CHECKHEAP < 0 */
