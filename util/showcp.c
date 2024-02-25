/* make showcp
  2024/02: V2.0f: bug of V2.0e (selection of columns) fixed
  2024/01: V2.0e: option -z extended, -e removed, lin.regr. of Etot removed
  2023/11: V2.0d: option -- added
  2022/07: V2.0c: option -z extended
  2022/05: V2.0c: option -z added
  2022/03: V2.0b: linear regression uses t=DT*n, not n
  2019/11: V2.0a: new CPmarks, reverse endian removed, time support (cf. option -i)
  2019: SHOWCPGEOMETRY added
  2018: Lx*Ly*Lz and 1/rho added
  2015: Ekflow added, -a-# (add column Etot-Etot0) removed

  option -r changed, now = running (cumulative) convergence profile
  REMOVED reverse endian = environment variable REVERSEENDIAN

  printing MACSIMUS convergence profiles, as produced by cook* etc.

  use `make' to compile showcp - see makefile

  IF the first record in the .cp file starts with CPmark as the 1st
  number, then the 2nd number contains NCP either in the 1st or 4th byte
  (to recognize endianess), the 4th or 1st beeing zero.  The remaining
  (NCP-2) 4-byte fields in the record contain max 4 letter info; the 1st
  two info's are assumed to be "Etot" and "Tkin".
  The minimum NCP is 2.

  Very old format: 5 columns without a header
*/

#define VERSION "2.0f"
#include "options.h"
#define x badoption
int optionlist[32] =
/* ` a b c d   e f g h i j k  l   m  n  o p q r s t  u v w x y z { | } ~   */
  {0,0,0,0,500,0,1,0,0,x,x,1,29,127,15,-1,0,x,0,0,0,-1,0,x,0,x,0,x,0,x,x,x};
/* @ A B C D E F G H I J K  L   M  N  O P Q R S T  U V W X Y Z [ \ ] ^ _ */
#undef x
#define INTFMT "%d"

#include "ground.h"
#include <time.h>
#include <unistd.h> /* unsigned sleep(unsigned); */
#include "cpmark.h"
#include "pakcp.h"

#include "varfile.h"
#include "statics.h"

// #include "linregr.h"

#include "showcphlp.c"

struct opto_s {
  struct opto_s *next;
  int n;
  int col[1]; /* [n] */
} *head=NULL;

void addopto(char *opt) /******************************************* addopto */
{
  char *c;
  int n=1;
  struct opto_s *o;

  for (c=opt; *c; c++) n+=*c==',';

  alloc(o,sizeof(*o)+(n-1)*sizeof(o->col[0]));

  o->col[0]=atoi(opt)-1;
  if (o->col[0]<0) ERROR(("-o: syntax or nonpositive arg"))
                     n=1;
  for (c=opt; *c; c++) if (*c==',') {
    o->col[n]=atoi(c+1)-1;
    if (o->col[n]<0) ERROR(("-o: syntax or nonpositive arg"))
                       n++; }
  o->next=head;
  o->n=n;
  head=o;
}

FILE *CP;
int NCP=5;
char *simname,*plotname,*dot;
int CPL=79;   /* chars per line -1 */
int LINES=22; /* lines per page for -# options -2 */

/* pseudograph: */
#define STEPS "_-^"
#define CHARLEV 3 /* 2 or 3 */

char defCPkey[5]=  {'E',   'T',    'U',   'i',     'P'};
char deftit[5][8]= {"Etot","Tkin", "Epot","Eintra","P   "};
char *CPkey;
char (*tit)[8];
int argname;

static double sumT,sumTT;
static double *sumr,*sumrT,*drdT;
double averT;
double t=0,DT=0,lastDT=0; /* cook t,DT reconstructed */

int block,from,to,ncp,n;

int check=0, file=0;

float *r,*minr,*maxr,*firstr,*lastr;
int *zero; /* if set, 1st record is subtracted */
double *sum;
int    i,c;
char   *l; /* [CPL+1] */

static int rewound=1, packed;
static int verbose;
static char *zlist;

#define MAXNCP option('m')
char *cols;
int ncols;
unsigned char **gr; /* pseudo-graph */

int eread(void) /***************************************************** eread */
/* read one record of CP */
{
  int i,j;
  static int next=0;

  if (rewound) { /* newly opened or rewound */
    char key='@';

    next=0;
    r=getcpheader(CP,&NCP,verbose);
    if (!r) Error("showcp: read CP/CPZ");

    if (!tit) {
      allocarrayzero(tit,NCP);
      alloc(CPkey,NCP);
      loop (i,0,2) {
        strcpy(tit[i],deftit[i]);
        CPkey[i]=defCPkey[i]; }
      loop (i,2,NCP) {
        memcpy(&tit[i],&r[i],4);
        CPkey[i]=tit[i][0];
        loop (j,0,i) if (CPkey[j]==CPkey[i]) {
          if (i<9) CPkey[i]='0'+i+1;
          else CPkey[i]=0; } }

      loop (i,0,NCP) if (!CPkey[i]) {
      advancekey: key++;
        loop (j,0,NCP) if (CPkey[j]==key) goto advancekey;
        CPkey[i]=key; }

      allocarrayzero(zero,NCP);
      allocarrayzero(minr,NCP);
      allocarrayzero(maxr,NCP);
      allocarrayzero(lastr,NCP);
      allocarrayzero(firstr,NCP);
      allocarrayzero(sum,NCP);
      allocarrayzero(sumr,NCP);
      allocarrayzero(sumrT,NCP);
      allocarrayzero(drdT,NCP); }

    rewound=0;

    return 1; }

  if (!packed) {
    int ret=fread(r,sizeof(float),NCP,CP)==NCP;
    return ret; }
  else { /* packed */
   again:
    if (!next) {
      if (!getpakcp(CP)) return 0;
      else next=1; }
    if (nextcprec(r))
      return 1;
    else {
      endgetpakcp(CP);
      next=0;
      goto again; } }
}

static void initgr(int n) /****************************************** initgr */
{
  int i,j;

  if (option('u')) {
    alloczero(gr[n],(CPL+1)*LINES*3);
    loop (j,0,CPL)
      loop (i,0,LINES) {
        gr[n][(i*(CPL+1)+j)*3+0]=226;
        gr[n][(i*(CPL+1)+j)*3+1]=160;
        gr[n][(i*(CPL+1)+j)*3+2]=128; } }
  else {
    alloczero(gr[n],(CPL+1)*LINES);
    loop (j,0,CPL)
      loop (i,0,LINES) gr[n][i*(CPL+1)+j]=' '; }
}

static void addbraille(unsigned char *grn,int x,int y) /********* addbraille */
{
  unsigned char *G=grn+((y/4)*(CPL+1)+(x/2))*3;
  int x1=x&1;
  int y3=y&3;

  if (G-grn>=(CPL+1)*LINES*3 || G-grn<0) ERROR(("addbraille %d",G-grn))

  if (y3==0) G[1]|=1<<x1;
  else G[2]|=1<<((3-y3)+x1*3);
}

int Pipe=0;
void PIPE(void) /****************************************************** PIPE */
{
  if (Pipe) prtc('#');
}

int main(int narg, char **arg) /*************************************** main */
{
  FILE *CPA=NULL, **xfile;
  int pp=0,Block,icp,iscpa=0,warn=0,usetimemarks=1;
  int passC=0,simplepipe=0; /* option -- */
  double t0,tfrom=-9.9e99,sumt;
  char *plotcommand; /* must accept -pPARENT_ID */
  double Ekflow=0,sumEkflow=0,rho0=0,P0=0;
  char *geometry=getenv("SHOWCPGEOMETRY");

#if CHECKHEAP==2
  AllocRange=96; AllocTrace=1;
#endif

  plotcommand=getenv("PLOTCOMMAND");
  if (!plotcommand) plotcommand="plot";
  if (!geometry) geometry="900x343";

  initscroll(0);

  /*** option analysis ***/
  loop (i,1,narg)
    if (arg[i][0]=='-') {
      if (arg[i][1]==0)
        Pipe=1;
      else if (isdigit(arg[i][1])) {
        int n=atoi(arg[i]+1)-1;

        if (n<0 || n>=MAXNCP) Error("showcp: option -# (column) out of range");
        if (!cols) allocarrayzero(cols,MAXNCP);
        cols[n]++;
        Max(pp,n) }
      else if (arg[i][1]=='-')
        simplepipe++;
      else {
        if (arg[i][1]=='C') passC++;
        if (arg[i][1]=='o' && arg[i][2]!='-' && arg[i][2]!='0')
          addopto(arg[i]+2);
        getoption(arg[i])
        if (arg[i][1]=='z') zlist=arg[i]+2;
        if (arg[i][1]=='h') DT=atof(arg[i]+2); } }
    else if (!simname) {
      char *c;

      alloc(simname,strlen(arg[i])+StaNameLen+8);
      *simname++='-'; /* for .cov */
      *simname++="cC"[passC]; /* for .cov */
      strcpy(simname,arg[i]);
      for (c=simname; *c; c++) if (*c=='.') dot=c;

      if (dot) {
        int unpacked=!strcmp(dot,".cp");

        packed=!strcmp(dot,".cpz");
        if (packed||unpacked) {
          CP=fopen(simname, "rb");
          if (!CP) ERROR(("%s not found",simname)); }
        else {
          PIPE();
          prt("%s: unexpected extension, trying to append .cp or .cpz",
              simname); } }

      if (!CP) {
        dot=c;
        strcpy(dot,".cp");
        PIPE();
        prt("no extension, trying %s",simname);
        CP=fopen(simname, "rb");
        if (CP) {
          FILE *tst;
          strcpy(dot,".cpz");
          tst=fopen(simname, "rb");
          strcpy(dot,".cp");
          if (tst) {
            WARNING(("both %s and .cpz files present: %s will be used",
                     simname,simname))
              fclose(tst); } }
        else {
          strcpy(dot,".cpz");
          CP=fopen(simname, "rb");
          packed=1;
          PIPE();
          prt("not found, trying %s",simname);
          if (!CP) {
            ERROR(("no %s nor .cp file",simname))
            exit(0); } } } }
    else /* !simname */
      if (!argname) argname=i;

  /* MAXNCP=option('m') */
  allocarrayzero(xfile,MAXNCP);
  if (!cols) allocarrayzero(cols,MAXNCP);
  allocarrayzero(gr,MAXNCP);

  if (option('p')) option('a')=option('p');
  if (option('a')<0) option('a')=1;
  if (option('a')>1) option('a')=2;
  if (DT<0) DT=-DT,usetimemarks=0;
  //  option('z')--;

  if (option('o')<0) addopto("2,6,7");

#if 0
  if (!option('@')) {
    Pipe++;
    PIPE();
    Pipe--;
    prt(" showcp: view and analyze MACSIMUS convergence profiles"); }
#endif

  if (!simname) {

    prtsfill(Pintro);

    for (;;) {
      char s[8];

      prts_(Phelp);
      if (!fgets(s,8,stdin)) s[0]=0;
      switch (s[0]) {
        case 'a': prtsfill(Pactions); break;
        case 'c': prtsfill(Pcolumn); break;
        case 'e': prtsfill(Penv); break;
        case 'g': prtsfill(Pgui); break;
        case 'f': prtsfill(Pfiles); break;
        case 'i': prtsfill(Pintro); break;
        case 'o': prtsfill(Poptions); break;
        case 'r': prtsfill(Prange); break;
        case 't': prtsfill(Ptime); break;
        default: return 0; } }

  exit(1); }

  if (getenv("LINES")) LINES=atoi(getenv("LINES"))-2;
  if (getenv("COLUMNS")) CPL=atoi(getenv("COLUMNS"))-1;

  option('k')&=3;
  if (option('k')==3) LINES--;
  allocarrayzero(l,CPL+1);

  if (option('u')<0) {
    char *e=getenv("LOCALE");

    if (e) option('u')=!!strcmp("C",e);
    e=getenv("LC_ALL");
    if (e) option('u')=!!strcmp("C",e);
    if (option('u')<0) option('u')=1; }

  if (!CP) Error("showcp: CP/CPZ file not found");

  if (option('s')==1) option('s')=80;
  if (option('s')) initscroll(option('s')*1024);
  from=option('f');
  to=option('t');
  block=option('b');

  if (!eread()) Error("showcp: empty CP/CPZ file or bad header");
  rewind(CP);

  if (option('z')) {
    unsigned i;

    if (zlist[0]) {
      i=atoi(zlist);

      while (i) {
        if (i<=NCP) zero[i-1]++;
        zlist=strchr(zlist,',');
        if (zlist++) i=atoi(zlist);
        else i=0; } }
    else
      loop (i,0,NCP) zero[i]++; }

  if (argname)
    loop (i,argname,narg) if (arg[i][0]!='-') {
      int n;

      loop (n,0,NCP) if (!memcmp(arg[i],tit[n],min(strlen(arg[i]),4))) {
        cols[n]++;
        Max(pp,n)
        goto found; }
      ERROR(("field \"%s\" not found",arg[i]))
   found: ; }

  loop (i,0,NCP) if (cols[i]) ncols++;

  /* no action and any column selected forces pseudograph */
  if (!(option('p') || option('e') || option('T') || option('a') || option('c')))
    if (ncols) option('g')=1;

  /* no column selected -> all columns */

  if (!ncols) loop (i,0,NCP) cols[i]++;
  ncols=NCP;

  if (option('g')) loop (i,0,NCP) if (!gr[i] && cols[i]) initgr(i);

  PIPE();
  verbose=!option('@') && !simplepipe;
  if (verbose)
    prt("%s from=%d to=%d (0=eof) block=%d (0=auto)",simname,from,to,block);
  /* note: because of verbose=1, next output will be ncp= from eread() */

  PIPE();

  rewound=1;

  ncp=n=0;

  if (simplepipe) {
    while (eread()==1) {
      if (r[0]>CPmark) {
        if (++ncp>=from) {
          if (to) if (ncp>to) break;
          if (zlist) loop (i,0,NCP) if (zero[i]) { lastr[i]=r[i]; zero[i]=0; }
          loop (i,0,NCP) if (cols[i]) printf("%11.8g ",r[i]-lastr[i]);
          putchar('\n'); } } }
    return 0; }

  /* determining min/max and non-blocked statistics */
  while (eread()==1) {
    if (r[0]>CPmark) {
      if (++ncp>=from) {
        double V=1,irho=0;
        int key=0;
        static int notStaSet=1;

        if (option('l')) {
          if (notStaSet) {
            StaSet(DT,option('l'),2,option('n'));
            notStaSet=0; } }

        if (to) if (ncp>to) break;
        // LRAdd("Etot",1,n*DT,r[0]);

        if (option('l'))
          loop (i,0,NCP) if (cols[i]) {
            StaAdd(tit[i],r[i]);
            if (!strcmp(tit[i],"rho ")) irho=1/r[i];
            if (!strcmp(tit[i],"Lx  ")) V*=r[i],key|=1;
            if (!strcmp(tit[i],"Ly  ")) V*=r[i],key|=2;
            if (!strcmp(tit[i],"Lz  ")) V*=r[i],key|=4; }

        if (irho) StaAdd("1/rho",irho);
        if (key==7) StaAdd("LxLyLz",V);

        sumT+=r[1]; sumTT+=Sqr(r[1]);
        loop (i,2,NCP) { sumr[i]+=r[i]; sumrT[i]+=r[i]*r[1]; }

        if (n++)
          loop (i,0,NCP) {
            Min(minr[i],r[i])
            Max(maxr[i],r[i])
            lastr[i]=r[i]; }
        else {
          copy(firstr,r,sizeof(float)*NCP);
          copy(minr,r,sizeof(float)*NCP);
          copy(maxr,r,sizeof(float)*NCP); } } }
    else if (r[0]==CPmark) {
      /* legacy time in ASCII format;
         expr. after && indicates ASCII, not integer number<=0xffff */
      if (!option('@') && ((char*)(r+1))[0]*((char*)(r+1))[3]) {
        unsigned char t[25],*c;

        copy(t,r+1,24); t[24]=0;
        for (c=t; *c; c++) if (*c<' ' || *c>127) {
            strcpy((char*)t,"-");
            break; }
        PIPE(); prt("# %7d: %s",n,t); } }
    else if (r[0]==CPmarkT) {
      /* real time in time_t format */
      if (!option('@')) {
        time_t tt;
        char *ct;

        if (NCP<=2) {
          int4 u;
          copy(&u,r+1,4);
          tt=u; }
        else
          copy(&tt,r+1,8);

        ct=ctime(&tt);
        PIPE(); prt_("# %7d: %s",n,ct); } }
    else if (r[0]==CPmarkU) {
      /* simulation time t in ps */
      if (usetimemarks) {
        if (NCP<=2) {
          float u;
          copy(&u,r+1,4);
          t=u; }
        else
          copy(&t,r+1,8);
        prt("#         t=%.6f ps",t); } }
    else if (r[0]==CPmarkV) {
      /* timestep */
      if (usetimemarks) {
        lastDT=DT;
        if (NCP<=2) {
          float u;
          copy(&u,r+1,4);
          DT=u; }
        else
          copy(&DT,r+1,8);
        if (lastDT && DT!=lastDT)
          prt("WARNING: DT changed %.15g -> %.15g",lastDT,DT), warn++;
        if (!option('@'))
          prt("#         dt=h*noint=%.6f ps",DT); } }
    else if (r[0]==CPmarkW) {
      /* simulation time t and timestep */
      if (usetimemarks) {
        if (NCP>=5) {
          copy(&t,r+1,8);
          lastDT=DT;
          copy(&DT,r+3,8);
          if (lastDT && DT!=lastDT)
            prt("WARNING: DT changed %.15g -> %.15g",lastDT,DT), warn++;
          if (!option('@'))
            prt("#         t=%.6f ps  DT=h*noint=%.6f ps",t,DT); } } }
    else
      Error("showcp: unknown CPmark type");
  } /* eread() */

  if (option('@')) {
    prt("%d",ncp);
    return 0; }

  prt("\n# %d data records read\n",ncp);

  if (DT==0) {
    fprintf(stderr,"\
=========================== WARNING =============================\n\
Time cycle (DT=h*noint) information is not present in the cp-file\n\
and it has not been specified using option -h. DT=1 will be used.\n\
=================================================================\n");
    DT=1; }

  *dot=0; /* no extension */

  if (option('l')) {
    if (!option('c'))
      /* one line of statistics */
      StaPrintAll(NULL);
    else {
      if (option('c')&1) StaPrintAll(simname); /* TCF */
      if (option('c')&2) StaPrintAll(simname-2); /* COV */
      if (option('c')&4) StaPrintAll(option('c')&32?"+":""); }
    StaFree(); }


  rewind(CP);
  rewound=1;

  if (block==0) {
    if (option('|')) block = (n+99)/100; /* cca 100 lines */
    else if (n<=1000) block=1;
    else if (n<=10000) block=10;
    else if (n<=100000) block=100;
    else if (n<=1000000) block=1000;
    else if (n<=10000000) block=10000;
    else if (n<=100000000) block=100000;
    else block=1000000;
    prt("### block autoset: block=%d (%d %sdata points)",block,n/block,n==1?"":"blocked "); }

  if (option('g')) {
    if (option('u')) block = (n+CPL*2-1)/(CPL*2);
    else block = (n+CPL-1)/CPL; }

  plotname=simname;

  if (option('a')) {
    strcpy(dot,option('r')?".cpr":".cpa");
    CPA=Pipe?stdout:fopen(simname,"wt");
    if (!CPA) {
      CPA=fopen(plotname=option('r')?"/tmp/showcp.cpr":"/tmp/showcp.cpa","wt");
      if (!CPA) ERROR(("cannot write to %s nor /tmp/showcp.cpa/.cpr",simname))
      else prt("\
WARNING: cannot create %s, using /tmp/showcp.cpa/.cpr\n\
(do not run more instances of showcp)",simname); }
    iscpa++;
    fprintf(CPA,"# block="INTFMT"\n#__",block);
    if (option('a')>1) fprintf(CPA,"    1:t/ps ");
    loop (i,0,NCP) fprintf(CPA,"%d:%-10s",i+1+(option('a')>1),tit[i]);
    if (option('v')) fprintf(CPA,"   %d:%-8s\n",NCP+1,"Ekflow");
    else fprintf(CPA,"\n"); }

  PIPE();
  put(averT=sumT/n)
  _n PIPE();
  prt("colm_ky_name__________min____________max________max-min_");
  loop (i,0,NCP) {
    PIPE();
    prt_("%2d %c  %c %4s %13.6g  %13.6g  %13.6g",
         i+1,i+'A',CPkey[i],
         tit[i],minr[i],maxr[i],maxr[i]-minr[i]);
    if (zero[i]) prt_("   1st item :=0");
    _n }

  if (block>1) {
    PIPE();
    prt("averaged in blocks by %i",block); }

  if (option('|')) putline('=',CPL);

  ncp=n=0;
  memset(sum,0,sizeof(sum[0])*NCP);
  sumt=0;
  verbose=0;

  if (option('x')) {
    loop (i,0,NCP) if (cols[i]) {
      char *c=dot,*t=tit[i];

      *c++='.';
      while (*t>' ') *c++=*t++;
      strcpy(c,".xy");
      xfile[i]=fopen(simname,"wt");
      if (!xfile[i])
        ERROR(("cannot create %s",simname))
      else {
        PIPE();
        prt("writing %s",simname); } } }

  /* T-dependence, blocked statistics, graph */
  t0=0;
  icp=0; /* for sim. time, may be restarted */
  while (eread()==1) {
    if (r[0]>CPmark) {
      /* valid data record */
      icp++;
      t=t0+icp*DT;

      if (++ncp>=from) {

        if (option('z'))
          loop (i,0,NCP) {
            /* column shifted to start from 0 */
            if (zero[i]) r[i]-=firstr[i]; }

        if (tfrom<-9e99) tfrom=t;
        if (to) if (ncp>to) break;

        if (option('l') && option('c')&16)
          loop (i,2,NCP) if (cols[i]) {
            char st[64];

            strcpy(st,"linT-"); strcat(st,tit[i]);

            StaAdd(tit[i],r[i]-drdT[i]*(r[1]-averT)); }

        loop (i,0,NCP) sum[i] += r[i];
        sumt+=t;

        if (option('v')) {
          if (rho0>0) Ekflow-=(r[4]-P0)*(0.5/r[3]+0.5/rho0);
          rho0=r[3]; P0=r[4];
          sumEkflow+=Ekflow*(option('v')/1e6); }

        if ((++n)%block==0) {

          Block=option('r')?n:block;

          memset(l,' ',CPL); l[CPL]=0; l[0]='|'; l[CPL-1]='|';

          loop (i,0,NCP) if (cols[i]) {
            if (block!=1 && option('c')&8 && option('l') && !option('r')) {
              StaAdd(string("blocked-%s",tit[i]),sum[i]/block); }

            c=floor( (sum[i]/Block-minr[i])/(maxr[i]-minr[i]+1e-30)*CPL );
            if (c<0) c=0;
            if (c>=CPL) c=CPL-1;
            l[c] = l[c]==' ' || l[c]=='|'  ? CPkey[i] : 'X';
            if (gr[i]) {
              if (option('u')) {
                c=(sum[i]/Block-minr[i])/(maxr[i]-minr[i]+1e-30)*(LINES*4);
                Max(c,0) Min(c,LINES*4-1)
                addbraille(gr[i],(n-1)/block,c); }
              else {
                c=(sum[i]/Block-minr[i])/(maxr[i]-minr[i]+1e-30)*(LINES*CHARLEV);
                Max(c,0) Min(c,LINES*CHARLEV-1)
                  gr[i][(CPL+1)*(c/CHARLEV)+(n-1)/block]=STEPS[c%CHARLEV]; } } }

          if (CPA) {
            if (option('a')>1) fprintf(CPA,"%13.10g ",sumt/block);
            loop (i,0,NCP) fprintf(CPA,"%11.8g ",sum[i]/Block);
            if (option('v')) fprintf(CPA,"%11.8g\n",sumEkflow/Block);
            else fprintf(CPA,"\n"); }

          loop (i,0,NCP) if (xfile[i]) {
            int ier;

            if (DT!=0) ier=fprintf(xfile[i],"%.8g %.8g\n",ncp*DT,sum[i]/Block);
            else ier=fprintf(xfile[i],INTFMT" %.8g\n",ncp,sum[i]/Block);
            if (ier<=0) ERROR(("cannot write to xy file - disk full?")) }

          if (option('|')) { PIPE(); prts(l); }

          if (!option('r')) loop (i,0,NCP) sum[i]=0;
          sumt=0;
          sumEkflow=0; } } }
    else if (usetimemarks) {
      if (r[0]==CPmarkU) {
        /* simulation time t in ps */
        if (NCP<=2) {
          float u;
          copy(&u,r+1,4);
          t=u; }
        else
          copy(&t,r+1,8);
        t0=t; icp=0; }
      else if (r[0]==CPmarkV) {
        /* timestep */
        if (NCP<=2) {
          float u;
          copy(&u,r+1,4);
          DT=u; }
        else
          copy(&DT,r+1,8); }
      else if (r[0]==CPmarkW) {
        /* simulation time t and timestep */
        if (NCP>=5) {
          copy(&t,r+1,8);
          t0=t; icp=0;
          copy(&DT,r+3,8); } } }
  }

  fclose(CP);
  if (CPA) fclose(CPA);
  loop (i,0,NCP) if (xfile[i]) fclose(xfile[i]);

  if (option('|')) putline('=',CPL);

  if (option('c')&24) {
    _n PIPE();
    StaPrintAll(option('c')&4?"":NULL); }

  // LRPrint(' ');

  if (iscpa) {
    PIPE();
    prt("column numbers:");
    if (option('a')>1) prt_(" 1=t");
    loop (i,0,NCP) prt_(" %d=%s",i+1+(option('a')>1),tit[i]);
    if (option('v')) prt_(" %d=%s",NCP+1,"Ekflow");
    _n }

  if (option('g')) putline('=',CPL);
  loop (i,0,NCP) if (gr[i]) {
    if (maxr[i]>minr[i]) {
      if (option('u')) for (n=LINES-1; n>=0; n--) prts((char*)(gr[i]+n*(CPL+1)*3));
      else for (n=LINES-1; n>=0; n--) prts((char*)(gr[i]+n*(CPL+1))); }
    else
      prt("<<<EMPTY RANGE>>>");
    if (option('k')&1) putline('=',
      CPL-prt_("=== %d:%s === [%g,%g] === dif=%g === by "INTFMT" ",
                   i+1,tit[i],minr[i],maxr[i],maxr[i]-minr[i],block));
    if (option('k')&2) putline('=',
      CPL-prt_("=== %d:%s === 1st->last=%g->%g === dif=%g by "INTFMT" ",
                   i+1,tit[i],firstr[i],lastr[i],lastr[i]-firstr[i],block)); }

  if (option('s')) {
    prt("use $COMMAND ($? for help) or ; to end");
    getdata enddata }

  if (option('p')) {
    struct opto_s *o;
    int pid=getpid();
    char *fn=strdup(string("/tmp/plot.%d",pid));
    int nplots=0;

    if (option('z'))
      loop (i,0,NCP) if (zero[i]) {
        minr[i]-=firstr[i];
        maxr[i]-=firstr[i]; }

    loop (i,0,NCP) if (cols[i]) {
      looplist (o,head) {
        int j;
        loop (j,0,o->n) if (o->col[j]==i) goto colfound; }
      /* not plotted for empty range unless specified in -o */
      if (maxr[i]>minr[i]) addopto(string("%d",i+1));
    colfound:; }

    loop (i,0,NCP) if (cols[i]) looplist (o,head) if (o->col[0]==i) nplots++;
    if (nplots>1) remove(fn);

    loop (i,0,NCP) if (cols[i]) looplist (o,head) if (o->col[0]==i) {
      char s[1024];
      int j;

      if (nplots>1)
      sprintf(s,"%s -p%d \"[:][%.8g:%.8g]\" %s",
              plotcommand,
              pid,
              minr[i],maxr[i],
              plotname);
      else
        sprintf(s,"%s \"[:][%.8g:%.8g]\" %s",
                plotcommand,
                minr[i],maxr[i],
                plotname);

      loop (j,0,o->n)
        if (option('a')>1)
          sprintf(s+strlen(s),"\":1:%d\" ", o->col[j]+2);
        else
          sprintf(s+strlen(s),"\":%.15g+%.15g*((n+.5)*%d-.5):%d\" ", tfrom,DT,block,o->col[j]+1);

      {
        static char genv[64],nenv[256];
        useconds_t sl=option('d')*1000;

        sprintf(genv,"PLOTGEOMETRY=%s+%d+%d",geometry,8*(i+1),30*(i+1));

        putenv(genv);
        sprintf(nenv,"PLOTNAME=%s %d ",simname,i+1);
        *dot='.';
        loop (j,0,o->n) if (o->col[j]<NCP) {
          strcat(nenv," ");
          strcat(nenv,tit[o->col[j]]); }
        putenv(nenv);
        strcat(s,"&");
        //        fprintf(stderr,"%s\n",s);
        if (system(s)) fprintf(stderr,"%s failed",s);
        usleep(sl); }
      } }

  if (warn) WARNING(("%d changes of time cycle (DT=h*noint) have been detected or wrong -h used.\n\
*** This is OK for showing graphs with options -a2|p2.\n\
*** For correct statistics, timestep-uniform ranges must be selected\n\
*** by options -f|-t.",warn))

  fclose(out);

  return 0;
}
