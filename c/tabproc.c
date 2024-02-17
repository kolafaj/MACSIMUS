/* \make tabproc
*/

#include "ground.h"
// double nan(const char *tagp); /* header missing in some older systems */

// #define N 27 /* max # of columns; see also plotcalc.c */

#include "parm.h"

#define N 1000 /* max # of columns (excl. n=0); see also plotcalc.c */
double *count;
#include "tabinc.c"

struct {
  char *expr; /* expression with A, B, ... instead of #1, #2, ... */
  char *fmt;  /* output format */
  int intfmt; /* 1 if the format accepts int, not double */
  int maxcol; /* max column in expression */
} tab[N+1];

int main(int narg, char **arg)
{
  char s[10240];
  char *c,*a,*e,*p,*sum;
  int i,n,iarg,ierr=0;
  int base=1;
  struct _Idlist_s *id;
  int isNaN=0;
  double NaN;
  char *FMT=getenv("FMT");

  if (!FMT) FMT="%.8g"; /* default format */
  
  initscroll(0);

  if (narg<2) {
    fprintf(stderr,"\
Spreadsheet filter.  (c) J.Kolafa 1991-2022. Call by:\n\
  tabproc EXPR1[:FMT1] [EXPR2[:FMT2]] ... < INPUT > OUTPUT\n\
where\n\
  EXPR is expression, e.g.: \"A/2+sin(B)\"  \"#1/2+sin(#2)\"  \"pi^z\"\n\
  may include sum(FROM,TO), where FROM,TO are column numbers: \"sum(1,9)/9\"\n\
  1st column=A=#1=c1=x, .., 26th column=Z=#26=c26, 555th column=#555=c555,...\n\
  line number (from 0, reset by blank line)=@=#0=n=c0\n\
  FMT is int or double format; if missing, the previous applies\n\
  Examples: \"7.4f\"  \"%%7.4f\" \"err=%%.2e\"  \"file-%%c.dat\"\n\
Lines beginning with \"#\" or \"!\" are copied as comments\n\
Max %d columns (A..Z,c27..c%d) allowed\n\
Environment variables:\n\
  FMT  : default format for output; if not given, %%.8g applies\n\
  RNDSEED : seed for random numbers (default=0=time)\n\
  NOLF : (Any value) suppresses line feeds on output\n\
  BASE : Assumes input of integers with given BASE (see strtol; 1=float)\n\
  NAN  : If undefined (default), any malformed number on input produces NaN\n\
         If defined, garbage at end of number is ignored (123GARBAGE -> 123)\n\
         and garbage gives the value of NAN (with NAN=7, z1 -> 7)\n\
  a b .. j : will be available as variables a b .. j (as in plot)\n\
See also:\n\
  ev evu (for expression syntax, list of functions, etc.)\n\
  tab mergetab mergeg field prettab shifttab transtab data2tab\n\
  derivtab filttab latextab maketab selcol difxmin\n"
          ,N,N);
  exit(0); }

  NaN=nan("0x8000000000000");
  if (getenv("NAN")) isNaN++,NaN=strtod(getenv("NAN"),&e);

  if (getenv("BASE")) base=atoi(getenv("BASE"));
  if (base<0 || base>36) Error("bad environment BASE");

  loop (i,0,PARM) {
    char ID[8]="`",*ge;
    
    ID[0]=i+'a';
    if ( (ge=getenv(ID)) ) {
      
      alloczero(id,sizeof(struct _Idlist_s)+1);
      id->next=_Id.head;
      _Id.head=id;
      id->id[0]='a'+i;
      id->val=atof(ge); } }

  columnidlist();

  /* parsing args: extract expressions and formats */
  loop (iarg,1,narg) {
    i=iarg-1;
    if (i>N) ERROR(("too many output columns"))
    if (i) tab[i].fmt=tab[i-1].fmt;
    else tab[i].fmt=FMT;
    /* expression sum(FROM,TO) */
    if ( (sum=strstr(arg[iarg],"sum")) ) {
      char *lp=strchr(sum,'(');
      char *rp=strchr(sum,')');
      char *cm=strchr(sum,',');
      int from,to,j;
      
      if (!lp) Error("tabproc: sum: missing ( after sum");
      if (!cm) Error("tabproc: sum: missing , after sum(");
      if (!rp) Error("tabproc: sum: missing ) after sum(,");
      from=atoi(lp+1);
      to=atoi(cm+1);
      if (from>N) Error("tabproc: sum: too many columns");
      if (to>N) Error("tabproc: sum: too many columns");
      if (to<1) Error("tabproc: sum: wrong to column");
      alloc(tab[i].expr,strlen(arg[iarg])+5*to);
      strcpy(tab[i].expr,arg[iarg]);
      tab[i].expr[sum-arg[iarg]]=0;
      strcat(tab[i].expr,"(");
      loop (j,from,to) strcat(tab[i].expr,string("c%d+",j));
      strcat(tab[i].expr,string("c%d",to));
      strcat(tab[i].expr,rp); }
    else {    
      alloc(tab[i].expr,strlen(arg[iarg])+1);
      strcpy(tab[i].expr,arg[iarg]); }

    //    fprintf(stderr,"%d: %s\n",i,tab[i].expr);
    
    if ( (a=to_colon(tab[i].expr)) ) {
      a++; /* bug fixed 01/2022 */
      alloc(tab[i].fmt,strlen(a)+1); 
      if (strchr(a,'%'))
        strcpy(tab[i].fmt,a);
      else {
        strcpy(tab[i].fmt+1,a);
        tab[i].fmt[0]='%'; } }
    
    tab[i].maxcol=max_column;

    p=strchr(tab[i].fmt,'%');
    if (!p) Error("missing % in format");
    do {
      p++;
      if (!*p) Error("bad format"); 
    } while (!isalpha(*p));

    if (*p=='l') Error("l (long) specifier in format not supported");
    if (*p=='h') Error("h (short) specifier in format not supported");
    if (strchr("cxXodiu",*p)) tab[i].intfmt=1;
    else if (strchr("gGeEf",*p)) tab[i].intfmt=0;
    else Error("unknown or unsupported format"); }

  getsbufsize=10240;
  
  while (gets(s)) {
    for (c=s; *c!=0 && *c<=' '; c++);

    if (c[0]=='#' || c[0]=='!' || c[0]==0) {
      if (c[0]==0) *count=0;
      prts(s); }

    else {
      /* scan line: _Id.head->val = 1st number, etc. */
      e=c; n=0; id=_Id.head->next; /* note: _Id.head=count */
      do {
        /* skip white (needed for malformed number) */
        while (*e && *e<=' ') e++;
        c=e;
        if (base==1) id->val=fstrtod(c,&e);
        else id->val=strtol(c,&e,base);

        if (c==e) {
          id->val=NaN; /* was not valid number */
          while (*e>' ') e++; }
        else if (*e>' ') {
          /* garbage at end of number: skip to end of field */
          while (*e>' ') e++; 
          if (!isNaN) id->val=NaN; /* 123GARBAGE is NaN */ }

        id=id->next;
        n++;
      } while (e>c && id);

      loop (i,0,narg-1) {
        double y=0;
        
        if (i!=0) prtc(' ');
        e=tab[i].expr; /* err */

//    fprintf(stderr,"tab[%d].maxcol=%d n=%d expr=%s\n",i,tab[i].maxcol,n,tab[i].expr);
        if (tab[i].maxcol<n) y=Calc(tab[i].expr,&e);
        if (e==tab[i].expr) {
          prt_("???");
          ierr=-1; }
        else
          if (tab[i].intfmt) prt_(tab[i].fmt,(int)y); 
          else prt_(tab[i].fmt,y); }
      if (!getenv("NOLF")) _n
      *count+=1.0; } }

  return ierr;
}
