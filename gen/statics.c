/***************** statics.c ******* (c) J. Kolafa ****************

      autocorrelation analysis of serially correlated data

  01/2017 max # of points increased from 2^31-1 to 2^32-1
          BUG: it is assumed without check that int/unsigned is 4 byte long
  02/2011 error estimate (automatic) improved
  01/2010 dt added, StaAdd changed
  12/2009 key @ added (.tcf), code cleaned
  02/2007 bug (stuck if very large number) fixed
  08/2006 function StaStdErr added
  08/2005 to produce files of time covariances
  10/1996 to produce files of time correlation functions
  11/1995 for better portability, old version=statics0.c
  02/1992 2^n subaverage (block) analysis added
     1991 first C-version statisti.c
     1990 PL/1 version statics.pli, Pascal version statics.pas

******************************************************************/

#define WILD '*'
/*
  StaMean, StaN, StaVar accept WILD at the end of the name as a wildcard *
*/

#include "ground.h"
#include "varfile.h"
#include "statics.h"

/*** StaType and StaPtr ***/
typedef struct { double FirstVal,LastVal,SumCov; } StaType;
typedef struct QStruct *StaPtr;
struct QStruct {
  StaPtr    next;     /* -> to the next set of measurements */
  unsigned4 no;       /* # of measurements: saved from this point */
  int4      lag;      /* for original data */
  int4      blklag;   /* for block-averaged data */
  int4      noblk;    /* # of blocks */
  char      name[StaNameLen];
  StaType   s[1/*var len*/];
};

/***
  format of data in StaPtr->s[] :

  subscript_______________________contents_________________
  0                               FirstVal=min x, LastVal=max x, SumCov=SUM x^3
  1 .. 1+lag                      statistics upto lag of original data
  2+lag .. 2+lag+blklag           statistics of averages by 2 pieces of data
  3+lag+blklag .. 3+lag+2*blklag  statistics of averages by 4 pieces of data
  etc.                            etc.

  one block of statistics (of length of lag+1 or blklag+1) :

  rel.index___0_____________1______________2_________..._______lag_______
  FirstVal  SUM y          y[1]           y[2]               y[lag]
  LastVal   DT/Input*      y[n]           y[n-1]          y[n-lag+1]
  SumCov    SUM y^2   SUM y[i]*y[i-1]  SUM y[i]*y[i-2]  SUM y[i]*y[i-lag]
  * 1st block only: single value passed to StaAddBlk, then replaced by DT
    other blocks: values summed here by 2, then passed to StaAddBlk

  where y is x or a subaverage of 2,4,... x's  and n the number of y's recorded

  see statics.h for the functions
***/

/*** StaInfo ***/
static struct StaInfo_s {
  double   dt;     /* time interval of 1 datum */
  unsigned lag;    /* lag for original data */
  unsigned blklag; /* lag for data averaged in blocks of lengths of */
  unsigned noblk;                                 /* 2,4,...2^noblk */
  StaPtr   head;   /* first item in the list */
  StaPtr   tail;   /* last item in the list (see StaFind) */
  char     tname[StaNameLen];
                   /* truncated name (see StaFind) */
} StaInfo = { 1,32, 2, 12, NULL, NULL, "<ndef>" };


/*** local functions ***/

static unsigned StaDblSize(int lag,int noblk,int blklag) /* # of doubles */
{
  return (lag+2+noblk*(blklag+1))*3;
}

static unsigned StaSize(unsigned qdsize) /* size in bytes */
{
  return sizeof(struct QStruct) + sizeof(double)*(qdsize-3);
}

int StaError; /* set to 1 if name not found */

static StaPtr StaFind(const char *name,int msglev)
{
  StaPtr Q,QQ=NULL;
  char *wild=strchr(name,WILD);
  int l=wild?wild-name:0;

  StaError=0;

  /* truncate a too long name */
  if (strlen(name)>=StaNameLen) {
    copy(StaInfo.tname,name,StaNameLen);
    StaInfo.tname[StaNameLen-1]=0; }
  else
    strcpy(StaInfo.tname,name);

  looplist (Q,StaInfo.head) {
    if ( (l && memcmp(Q->name,StaInfo.tname,l)==0)
      || (l==0 && strcmp(Q->name,StaInfo.tname)==0) ) {
      if (QQ) prt("StaFind WARNING: %s multiple match",name);
      QQ=Q; }
    StaInfo.tail=Q; }

  if (!QQ) {
    StaError++;
    if (msglev) prt("StaFind WARNING: %s not found",StaInfo.tname); }

  return QQ;
}

static double StaAddBlk(StaType *q, unsigned4 n, int lag)
{
  int i;
  double x=q->LastVal; /* input value = block subaverage */

  q->FirstVal += x; /* sum x */
  if (n<=lag) q[n].FirstVal=x;
  loopto (i,0,lag) q[i].SumCov += x*q[i].LastVal;
  for (i=lag; i>0; i--) q[i].LastVal = q[i-1].LastVal;
  q->LastVal=0;

  return x;
}

static char fmt[]=" %3u%7.4f %6.3f";

static int errmode;

static int lstd;
static double av,dt;
static struct {
  double err[3]; /* blocked w/o covariance, blocked with c_1, blocked upto c_2 */
  int n;
} *laststd;

static enum {TCF,COV,NONE,ONE,STA,PRECSTA} statics; /* do not change order! */
  /*  ONE=one line of summary statistics to out
      STA=standard verbose table to out
      PRECSTA=more precise standard verbose table to out
      TCF=time correlation function, to out or file
      COV=covariances, to out or file */

static double StaPrtBlk(StaType *q, unsigned n, unsigned lag, unsigned blk)
{
  double sf=q->FirstVal;
  double sl=q->FirstVal;
  double t,ti,cor=1,lastcor=1,vari,varn,mx=0;
  unsigned i,m=n/3;

  if (n<2) return 0;
  blk /= n;
  if (m>lag) m=lag;
  av=q->FirstVal/n;
  vari=q->SumCov/n - av*av;
  if (vari<0) vari=0.0;
  varn=vari/((n-1)*Sqr((double)blk));
  if (vari==0) return 0;
  if (statics==COV) cor=vari; /* otherwise 1 */
  t=cor;

  for (i=0;;) {
    ti=i;
    if (dt>0) ti*=dt;

    switch (statics) {
      case COV:
        if (dt<1e-5) prt("%.8g %.8g %.8g",ti,cor,varn*t);
        else prt("%9.6f %.8g %.8g",ti,cor,varn*t);
        break;
      case TCF:
        if (dt<1e-5) prt("%.8g %9.6f %9.6f",ti,cor,t);
        else prt("%9.6f %9.6f %9.6f",ti,cor,t);
        break;
      default: { /* error estimates */
        double e;
        int pr=statics>=STA;

        if (pr) {
          fmt[0]=' ';
          if (i%3==0) {
            _n
            if (blk!=1) fmt[0]='*'; }
          fmt[13]='2'+(t<100)+(statics==PRECSTA)*2;
          fmt[5]='7'+(statics==PRECSTA)*2;
          prt_(fmt,i,cor,t); }

        /* the default: uncorrelated stdev multiplied by sqrt(1+2*sum c_i) */
        if (t<0) e=0;
        else e=sqrt(varn*t); /* sqrt(var*(1+2*tau_i) */
        
        if (errmode) {
          if (i==1)
            /* block-block: 1st order subseries process assumed */
            e=sqrt(varn*(t+4*Cub(cor)/((1+2*cor)*(1-cor))));
          if (i>1)
            /* block-...-block: 1st order process assumed */                                                 
            e=sqrt(varn*(t+2*Sqr(cor)/(lastcor-cor)));
        }

        Max(mx,e)

        if (laststd) {
          laststd[lstd].n=n;
          if (i<3) laststd[lstd].err[i]=e;
        }

        if (pr) {
          if (statics==PRECSTA)
            prt_("%13.5e",e);
          else {
            /* e-format is not compatible => it is created here from f and i */
            int p=0;

            if (e) {
              while (e<0.995) { e*=10; p--; if (p<-999) break; }
              while (e>9.95) { e/=10; p++; if (p>999) break; } }
            if (abs(p)>999) prt_(" ???????");
            else prt_("%4.1fe%c%i%i",e,p>=0?'+':'-',abs(p)/10,abs(p)%10); } } } }

    if (++i>m) break;

    sf=sf-q[i].FirstVal;
    sl=sl-q[i].LastVal;
    lastcor=cor;
    cor=(q[i].SumCov-sf*sl/(n-i))/(n-i);
    if (statics!=COV) cor/=vari;
    t+=cor+cor; }

  return mx;
}

static
double StaPrintErr(const char *name,const char *HOW) /********** StaPrintErr */
/*
  Prints (if pr) statistics of one quantity.
  See statics.h for HOW
  String HOW determines the mode:
    "/"    nothing printed
    NULL   only one line of mean+stderr is printed to out (statics=ONE)
    ""     compact table with c_t, tau, StDev printed to out (STA)
    "+"    as above with more decimal digits
    "-"    t:c_t table printed to stdout (TCF)
    "-@"   t:c_t table printed to file name.tcf (TCF)
    OTHER: t:c_t table printed to file HOW.name.tcf (TCF)
           (special characters in `name' are ignored - legacy mode)
    "-t"   in front of HOW = as above (TCF)
    "-T"   in front of HOW = as above + special chars kept except '/' -> '_'
    "-c"   in front of HOW = as -t except covariances are printed (COV)
    "-C"   in front of HOW = as above + special chars kept except '/' -> '_'

    REMOVED: "*" in front of HOW - is ERROR now!
*/
{
  StaPtr Q=StaFind(name,1);
  double av,vari /*,varn,sdev0*/ ;
  int nostd=0;
  int CHAR; /* 0: change / -> _ only ; 1: remove all problematic chars */
  double mx=0;
  unsigned i,prec;
  unsigned4 shift;
  StaType *q;
  FILE *oldout=NULL;
  static int pass=1;

  if (!Q) return 0;

  /*** set output format ***/
  if (!HOW)
    /* only 1 line of statistics (mean, stderr), low-precision */
    statics=ONE;
  else if (!strcmp(HOW,"/"))
    /* no output */
    statics=NONE;
  else if (!strcmp(HOW,""))
    /* "standard" low-prec output (this used to be default) */
    statics=STA;
  else if (!strcmp(HOW,"+"))
    /* "standard" high-prec output */
    statics=PRECSTA;
  else {
    if (HOW[0]=='*')
      Error("\"*\" is no longer a valid key in HOW: use \"-c\" or \"-C\"  instead");

    statics=TCF; CHAR=0; /* default */
    if (!memcmp(HOW,"-t",2)) { statics=TCF, CHAR=0; HOW+=2; }
    if (!memcmp(HOW,"-c",2)) { statics=COV, CHAR=0; HOW+=2; }
    if (!memcmp(HOW,"-T",2)) { statics=TCF, CHAR=1; HOW+=2; }
    if (!memcmp(HOW,"-C",2)) { statics=COV, CHAR=1; HOW+=2; }

    if (!strcmp(HOW,"-"))
      /* covariances/correlations to out */;
    else {
      /* covariances/correlations to file */
      char *fn,*c1;
      const char *c2;

      alloc(fn,strlen(HOW)+strlen(name)+6);
      if (strcmp(HOW,"-@")) { /* used to be "@" */
        strcpy(fn,HOW);
        if (*fn) strcat(fn,"."); }
      else
        fn[0]=0;
      for (c2=name; *c2==' '; c2++);

      for (c1=fn+strlen(fn); *c2; c2++) {
        if (CHAR) {
          /* keep only "good" characters */
          if (isalnum(*c2) || strchr("_-+~!@#.",*c2)) *c1++=*c2; }
        else {
          /* only / changed to _ */
          if (*c2=='/') *c1++='_';
          else *c1++=*c2; } }
      strcpy(c1,statics==COV?".cov":".tcf");
      oldout=out;
      out=fopen(fn,"wt");
      if (!out) { ERROR(("cannot write to %s",fn)) out=oldout; }
      free(fn); } }

  /*** 1st block of (not blocked) data ***/
  q=Q->s; q++;
  dt=q->LastVal;

  if (statics!=ONE && statics!=NONE)
    prt("%s%-24s  No=%lu  dt=%g",
         statics>=STA?"":"# ",
           Q->name,
                 (long unsigned)Q->no,
                               dt);

  /*** statistics of original data (no subaverages) ***/
  av=q->FirstVal/Q->no;
  vari=q->SumCov/Q->no - av*av;
  if (vari<0) vari=0.0;

  prec=9+(statics!=STA)*6;

  if (statics!=ONE && statics!=NONE)
    prt_("%s  Mean=%.*g",statics>=STA?"":"# ", prec,av);

  if (Q->no<2) {
    if (statics!=ONE && statics!=NONE) _n
    nostd++;
    goto ret0; }

  if (statics!=ONE && statics!=NONE) {
    prt("  Var=%.*g  skewness=%.*g",
               prec,vari,
                               prec, (Q->s->SumCov/Q->no
                                    - 3*(q->SumCov/Q->no)*(q->FirstVal/Q->no)
                                    + 2*Cub(q->FirstVal/Q->no))/(vari*sqrt(vari)));

    prt("%s  range=[%.*g,%.*g]=%.*g",
         statics>=STA?"":"# ",
                    prec,Q->s->FirstVal,
                         prec,Q->s->LastVal,
                               prec,Q->s->LastVal-Q->s->FirstVal); }

  if (vari==0 || Q->no<3) {
    nostd++;
    goto ret0; }

  switch (statics) {
    case COV: prt("#___t__Cov[t]__Var+2*SUM(Cov)"); break;
    case TCF: prt("#___t____c[t]__1+2tau"); break;
    case NONE:
    case ONE: break;
    case STA: loop (i,0,3) prt_("  _t__c[t]__1+2tau__StDev_"); break;
    case PRECSTA: loop (i,0,3) prt_("  _t____c[t]___1+2tau_______StDev_"); }

  mx=StaPrtBlk(q,shift=Q->no,Q->lag,Q->no);
  /* no blocking => max of Var*(1+2*tau) encountered (often pessimistic) */

  if (statics<=COV) goto ret;

  q += Q->lag+1;
  /*** statistics of data in 2,4,8,... blocks ***/

  if (Q->noblk) allocarrayzero(laststd,Q->noblk);

  loopto (i,1,Q->noblk) {
    lstd=i-1;
    shift >>= 1;
    StaPrtBlk(q,shift,Q->blklag,Q->no);
    q += Q->blklag+1; }

  if (Q->noblk>2 && Q->no>31) {
    int from=log(Q->no)-1.5; // heuristic (consider to derive from StaInfo.tau)
    int to=log(Q->no)/log(2)-3.4; // min ~10 blocks
    int j;
    double mxsum=0,wsum=0,w;

    if (from>to) to=from;
    if (from>=Q->noblk) from=to=Q->noblk-1;
    if (to>=Q->noblk) to=Q->noblk-1;

    switch (errmode) {

      case 0:
        loopto (j,from,to)
          if (laststd[j].n) {

            /* uncorrelated */
            if (from==to) w=1;
            else w=(double)(j+1-from)/(to-from+1);
            w*=1-1/sqrt(laststd[j].n);

            if (statics>=STA)
              prt_("\nD %2d n=%-3d  stderr(0)=%-10.4g w=%5.3f  ",
                   j+1,laststd[j].n,laststd[j].err[0],w);
            wsum+=w; mxsum+=w*Sqr(laststd[j].err[0]);

            /* c_1 */
            w=1.5-2/sqrt(laststd[j].n);
            if (w<0) w=.1;
            if (statics>=STA) prt(" stderr(1)=%-10.4g w=%5.3f",laststd[j].err[1],w);
            wsum+=w; mxsum+=w*Sqr(laststd[j].err[1]); }

        if (mxsum) {
          mx=sqrt(mxsum/wsum); }
        break;

      case 1:
        mx=laststd[to].err[1];
        break;

      case 2:
        mx=laststd[to].err[2];
        break;

      default:
        mx=0;
    }

  } /* enough data */

  if (Q->noblk) free(laststd);
  laststd=NULL;

 ret0:
  if (statics<=COV) goto ret;
  if (statics!=ONE && statics!=NONE) _n

  if (statics!=NONE) {
    if (pass) {
      pass=0;
      prt("#________average__std.err___rel.err______no_=name="); }

    if (nostd)
      prt(statics==PRECSTA
          ? "# %.15g  -      -     %d =%s="
          : "#%15.9g     -          -   %7d =%s=",
          av,Q->no,name);
    else
      prt(statics==PRECSTA
          ? "# %.15g %.6g %.6g %d =%s="
          : "#%15.9g %8.3g %9.3g %7d =%s=",
              av,    mx,   mx/av,Q->no,name);
    if (statics!=ONE) _n }

 ret:
  if (oldout) { fclose(out); out=oldout; }

  return mx;
} /* StaPrintErr */


/*** user functions ***/

void StaSet(double DT,int lag,int blklag,int noblk) /**************** StaSet */
{
  if (DT>0) StaInfo.dt=DT; /* remembers the old value */

  StaInfo.lag=lag;
  StaInfo.blklag=blklag;
  StaInfo.noblk=noblk;
}

void StaAdd(const char *name, double x) /**************************** StaAdd */
/*** next value x of quantity tx is processed ***/
{
  StaPtr Q=StaFind(name,0);
  unsigned i;
  unsigned4 shift;
  StaType *q;

  if (!Q) {
    /*** new data set appended to the list ***/
    int size=StaDblSize(StaInfo.lag,StaInfo.noblk,StaInfo.blklag);

    alloc(Q,StaSize(size));
    Q->next=NULL;
    Q->no=0;
    Q->lag=StaInfo.lag;
    Q->blklag=StaInfo.blklag;
    Q->noblk=StaInfo.noblk;
    copy(Q->name,StaInfo.tname,StaNameLen);
    if (StaInfo.head==NULL) StaInfo.head=Q;
    else StaInfo.tail->next=Q;
    loop (i,2,size) ((double*)Q->s)[i]=0.0;
    Q->s[1].LastVal=StaInfo.dt; /* 1st block */
    Q->s[0].FirstVal=Q->s[0].LastVal=x; }

  /*** x added to the data set ***/

  Q->no++;
  if (Q->no==0) ERROR(("StaAdd(\"%s\",%g):\n\
*** implementation limitation: number of measurements overflow",name,x))

  /* header */
  q=Q->s;
  if (x<q->FirstVal) q->FirstVal=x; /* MIN(x) */
  if (x>q->LastVal) q->LastVal=x;   /* MAX(x) */
  q->SumCov += x*x*x;               /* SUM(x^3) */

  /*** statistics of original data (no subaverages) ***/
  q++;

  q->LastVal=x; /* here: input only */
  StaAddBlk(q,shift=Q->no,Q->lag);
  q->LastVal=StaInfo.dt; /* to be stored */

  q += Q->lag+1;
  /*** statistics of data in 2,4,8,... blocks ***/

  loopto (i,1,Q->noblk) {
    x = q->LastVal += x;
    if (shift & 1) break;
    shift >>= 1;
    StaAddBlk(q,shift,Q->blklag);
    q += Q->blklag+1; }

} /* StaAdd */


void StaFree(void) /*********************************************** StaFree */
/* clears all measurements (not necessary in the first call) */
{
  StaPtr Q;
  while (StaInfo.head) {
    Q=StaInfo.head;
    StaInfo.head=StaInfo.head->next;
    free(Q); }
} /* StaFree */


double StaPrint(const char *name,const char *how) /**************** StaPrint */
{
  return StaPrintErr(name,how);
}

void StaPrintAll(const char *HOW) /***************************** StaPrintAll */
/* see StaPrint for `HOW' */
{
  StaPtr Q=StaInfo.head;

  while (Q) {
    StaPrintErr(Q->name,HOW);
    Q=Q->next; }
} /* StaPrintAll */


static void StaXSave(char *fn)
/* saves all measurements to file named fn */
{
  StaPtr Q=StaInfo.head;

  VarOpen(fn,"wb");
  while (Q) {
    VarPut(&Q->no,
           StaSize( StaDblSize(Q->lag,Q->noblk,Q->blklag))
           -(int)((char*)&Q->no-(char*)Q) );
    Q=Q->next; }
} /* StaXSave */


void StaSave(char *fn) /******************************************** StaSave */
{
  StaXSave(fn);
  VarClose();
}


void StaKeySave(char *fn,int4 key) /******************************** StaSave */
{
  StaXSave(fn);
  KeyClose(key);
}

void StaLoad(char *fn) /******************************************** StaLoad */
/*
  loads all measurements from file named fn until EOF
*/
{
  StaPtr Q=StaInfo.head;
  unsigned oldsize,size;

  VarOpen(fn,"rb");
  while (!VarFile.eof) {

    /* find list end */
    looplist (Q,StaInfo.head) StaInfo.tail=Q;

    /* normal read - no conversion */
    alloc(Q,oldsize=VarFile.size+(int)((char*)&Q->no-(char*)Q));
    VarRead(&Q->no,VarFile.size);

    size=StaSize(StaDblSize(Q->lag,Q->noblk,Q->blklag));
    if (size != oldsize) {
      ERROR(("StaLoad: size read=%u != expected = %u (try STAOLDFMT?)",
             oldsize,size))
      return; }

    /* incl. to list end */
    if (StaInfo.head==NULL) StaInfo.head=Q;
    else StaInfo.tail->next=Q;
    Q->next=NULL; }

  VarClose();
} /* StaLoad */


unsigned4 StaN(const char *name) /****************************************** StaN */
/* # of recorded measurements of quantity named name */
{
  StaPtr Q=StaFind(name,1);
  if (Q) return Q->no; else return 0;
} /* StaN */


double StaMean(const char *name) /********************************** StaMean */
/* mean value of quantity named name */
{
  StaPtr Q=StaFind(name,1);

  if (Q) return Q->s[1].FirstVal/Q->no; else return 0.0;
} /* StaMean */


double StaVar(const char *name) /************************************ StaVar */
/* variance of quantity named name */
{
  StaPtr Q=StaFind(name,1);

  if (Q==NULL) return 0.0;
  else {
    double av=Q->s[1].FirstVal/Q->no;
    return Q->s[1].SumCov/Q->no - av*av; }
} /* StaVar */


double StaStdErr(const char *name) /****************************** StaStdErr */
/*
  estimated standard error: the performance depends on StaSet() settings!
  (the same as StaPrint last line)
*/
{
  return StaPrintErr(name,"/");
} /* StaStdErr */

double StaLag(const char *name,char key) /************************** StaLag */
/*
  returns the lag according to the key:
  'l': lag
  'k': blklag
  'n': noblk
  't': time = lag*DT
  'd': DT
  other: -9e9;
*/
{
  StaPtr Q=StaFind(name,1);
  switch (key) {
    case 'l': return Q->lag; break;
    case 'k': return Q->blklag; break;
    case 'n': return Q->noblk; break;
    case 'd': return Q->s[1].LastVal; break;
    case 't': return Q->lag*Q->s[1].LastVal; break;
    default: return -9e9; }
} /* StaLag */

void StaErrMethod(int m) /************************************* StaErrMethod */
{
  errmode=m;
} /* StaErrMethod */
