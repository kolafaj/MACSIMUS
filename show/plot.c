/* make plot

  number color point       line
  1      W / K full circ   ----
  2      Y   B open circ   -- --
  3      C   R full sq     --.--.
  4      M   G open sq     ......
  5      G   M full diam   - - - -
  6      R   C open diam   -..-..
  7      B   Y full ^      -...-...
  8      brown open ^      .. .. ..
  9      cyan  full v      ... ...
  10     mag   open v      -- - -- -
  11     gray  +           . . . . .
  12     green x           . ... . ...
  13     red   *           . .. . ..
  14     blue  odot        -- - -- -
  15  =1, etc

  12/2023 (2.0s) ! missing column 2 filled by the prev value in plot file (w/o :)

  02/2023 buggy hot key U fixed

  10/2022 -n changed to NAME (ps.def -> NAME.def, plot.eps -> NAME.eps)
    active variables added to ps output

  5/2020 (2.0q) bit cleaned, hotkeys [ ] 0..9 fixed, hotkey U added

  11/2016 interface changed to xdraw.c, LAZYX11 removed

  7/2009 interface to fit.c extended, some hot keys changes

  01/2007 interface to fit.c - hot key t,T,ctrl-t
    must be called as (example):
      plot 2fit.dat:A:B:o:C :"a+b*A+c*A^2"
    file 2fit.dat must have x in 1st column, y in 2nd column, and
      optional std in 3rd column
    max PARM parms, far from fool-proof
    uses external routine exp2pow to convert ^ to pow or multiplication
    hot key '=' can serve to change parameters sent to minimize
                (see FIT4PLOT/fit.c and gen/minimize.h)

  06/2004 double leftclick emits coordinates and if env TOCLIP is set
    then puts it to the clipboard
  TOCLIP is the command to paste the text (appended 1 arg) to the clipboard
  e.g.:
    "dcop --user jiri klipper klipper setClipboardContents" (for KDE)
    "xclip.sh", where xclip.sh = shell script "echo $* | xclip"
  The result will be called as by system()

  04/2004 environment variables a b c d e f (see PARM) now provide default
          values for parameters
  02/2004 mouse interface changed from drawplus to mousekey

  07/2000 update: plot and plotc unified
  for the old version ("LEVEL=1") see plot0.c
  10/2000: Any redraw-causing event that comes during drawing is ignored.
           This may cause that the graph is not redrawn when a window
           is minimized or eclipsed during long drawing, but
           suppresses unnecessary redrawing on systems that
           occasionally generate several Expose events on start,
           maximize or so.
           drawx11: environment variable LAZYX11 causes waiting
           LAZYX11 seconds after creating a window (otherwise more
           Expose events can come later causing redraws)
environment variable PLOTINIT can contain the inital string to execute
e.g., to make PostScript and quit, use (csh/tcsh)
  setenv PLOTINIT "#q"
or (sh/bash)
  export PLOTINIT="#q"
*/

#define VERSION "2.0t"

/* allows changing param a by hotkeys a A etc. */
#include "parm.h"
/* maximum 10 (params=a..j): code has to be changed at several places */

#ifndef PARM
#  error "PARM not #defined (see gen/parm.h)"
#endif /*# PARM */

/* PostScript definitions to be placed in file "ps.def"
   reading see psdef.c
   all lengths are in pt if not stated otherwise
   values in [] are optional or, if inside <>, defaults
   case insensitive
  . EOF
  ! <comment>
  # <comment>
  b B/W mode: replaces color lines/points by different style black
  s <font size [14]> <encoding> [<default font [Helvetica]>]
  m <mode: Landscape Portrait Encapsulated [L]>
    !!! Portrait swaps xsize <--> ysize if used after command w
  x <left margin [100]> <right margin [5]> <xaxis text []>
  y <bottom margin [100]> <top margin [5]> <yaxis text []>
  w <x window size [576=8*72]> <y window size [432=6*72]>
    (negative = in cm)
    !! command m p (used after w) swaps xsize <--> ysize
  t <line thickness [1]> <frame thickness [0.5]>
  f <approx. # of ticks in x [5]> <# of ticks in y [5]> <tick length [5]>
    negative: abs. tick span
    this option may repeat, but labeling applies to the first f
  c R G B color
  # R G B[,<dash>,<linewidth>[,<pointsize>,<point>]
    (set PS color and style, #=1,2,...;
    <dash>=line pause [line pause ...]
    <point>={circle,fullcircle...}
    example: 1 1 0 0,5 3 1 3,0.5,2,opensquare )
    BUG: numbering is wrong (cf. line[16] below), b recommended instead
  r <rotation angle for l and y (y needs EXTPSSTRING)
  l x y STRING : string to write; can use _Index ^Exponent \Symbol
  l x y $#-STRING as above, prepend by a line of style # (=1,2,3..)
  l x y $#.STRING as above, prepend by a symbol of style # (=1,2,3..)
  l x y $#+STRING as above, both keys (line+symbol)
  l x y $#=STRING as above, errbar
        WARNING: in scripts (incl. cat > ps.def <<EOF), protect $ by \$
  l x +[dy] STRING : (next y by dy below)
*/

#include "ground.h"
#include <float.h>
#include "xdraw.h"
#include "ploterr.h"
#include "mydraw.h"
#include "tabinc.c"

#if 0
// attempt to prevent "fatal IO error 11.." - does not work
#  include <signal.h>

void MYHUP(int i)
{
  if (display) {
    XFreeGC(display,gc);
    XCloseDisplay(display); }
  exit(1);
}

//  signal(SIGHUP,MYHUP);
#endif /*# 0 */

#define fint(R) (REAL)(int)(R)

#define LXWIDTH 6 /* # of digits in the y-axis numbers */
int lbl=1; /* label and  grid style */
int maxxn,maxyn;

char *psshow(char *text) /******************************************* psshow */
/* see also psstring.c */
{
  static char *line;
  static int nline=254;
  int i,j;

 again:
  if (!line) allocarray(line,nline+2);

  for (i=0,j=0; text[i];) {
    if (j>=nline) {
      free(line);
      nline+=16;
      goto again; }
    if (strchr("()%\\",text[i])) line[j++]='\\';
    line[j++]=text[i++]; }
  line[j]=0;

  return line;
}

void makegrid(int xaxis,REAL x0,REAL x1,REAL y0,REAL y1) /********* makegrid */
/*
  NEW version for any font, new interface xdraw.c
  old BGI version (still valid for NSK): gen/makegrid.c
  lbl: 0=no grid
       1=light grid
       2=denser and brighter grid
  xaxis: 1=x-axis
         0=y-axis (x and y are swapped: x0,x1=the axis to be labeled)
*/
{
  if (lbl) {
    REAL dd=1,x1x0=x1-x0;
    int ifmt=0, j=0;
    static char fmt[16]; /* 8 enough, but the compiler shouts loud */
    int l,lfrom,lto,at;
    REAL ladd,decticks,minspan;

    if (x1x0<=0) return;

    if (xaxis) {
      at=getmaxy()-1;
      decticks=(2+lbl*4)*sqrt(getmaxy()/480.);
      minspan=16*(6-lbl*2)*sqrt(getmaxy()/480.);
      settextjustify(1,0); }
    else {
      at=1;
      decticks=(2+lbl*4)*sqrt(getmaxx()/640.);
      minspan=16*(6-lbl*2)*sqrt(getmaxy()/480.);
      settextjustify(0,1); }

    /* density of the grid */
    while (x1x0 > 10*decticks*dd) { dd*=10; ifmt--; }
    while (x1x0 <    decticks*dd) { dd/=10; ifmt++; }

    if (ifmt>16) ifmt=16;
    sprintf(fmt,"%%.%df",ifmt);
    if (ifmt<0 || (fabs(x0)<1e-5 && fabs(x1)<1e-5) ) strcpy(fmt,"%g");

    /* this is strange, but it seems to work */
    ladd=fint(x0/(dd*10))*10;
    lfrom=Int(x0/dd-1e-4-ladd);
    lto=Int(x1/dd+1e-4-ladd);

    // if (xaxis) put2(dd*scaling.x,getmaxx()/(dd*scaling.x)) else put2(dd*scaling.y,getmaxx()/(dd*scaling.y))

    loopto (l,lfrom,lto) {
      REAL xi=(l+ladd)*dd;

      if (xi<x0-dd*0.001 || xi>x1+dd*0.001) continue;

      setcolor(DARKGRAY);
      setlinestyle(2,lbl==1?0x0601:0x0301,1);

      if (xaxis) {
        //        put2(dd*scaling.x,minspan)
        if (dd*scaling.x>minspan) lline(xi+dd/2,y0, xi+dd/2,y1); }
      else {
        //        put2(dd*scaling.y,minspan)
        if (dd*scaling.y>minspan) lline(y0,xi+dd/2, y1,xi+dd/2); }

      setcolor(lbl==1?DARKGRAY:LIGHTGRAY);
      setlinestyle(2,0x0301,1);

      if (l%5==0) {
        setcolor(lbl==1?LIGHTGRAY:WHITE);
        if (l%10==0) setlinestyle(2,0x0202,1); }

      if (l%10==0 || (l%5==0 && lto-lfrom<33) || lto-lfrom<11) {
	char nr[24];
	int ii;

	sprintf(nr,fmt,Val(xi));
	if (xaxis) {
	  ii=SX(xi);
          line(ii,at-xfont.height-3,ii,at-xfont.height+3);
          if (ii-outtextwidth(nr)/2-xfont.width>j)
            j=ii+outtextxy(ii,at,nr)/2; }
	else
	  outtextxy(at,SY(xi),nr); }

      if (xi>x0-0.0001*x1x0) {
	if (xaxis) lline(xi,y0,xi,y1);
	else lline(y0,xi,y1,xi); } }

    setlinestyle(0,0,1);
    settextjustify(0,2); }
}

static int color[14]={
  WHITE,YELLOW,LIGHTCYAN,LIGHTMAGENTA,LIGHTGREEN,LIGHTRED,LIGHTBLUE,
  BROWN,CYAN,MAGENTA,LIGHTGRAY,GREEN,RED,BLUE};
int column,font=2;

#include "plotbit.c"
#define palcolor(COL) COL
#include "xdrawhelp.c"
#include "plothelp.c"
#include "plothlp.c"

#include "plotdef.c"

#if PARM
struct parm_s {
  double *parm;
  double step; /* changed by hotkeys [/] [*] and use by a A b B .. */
  double min,max; /* for the slider, initially min=*parm-step, max=*parm+step */
  int on;   /* changed or present in environment */
  int used; /* used in formulas */
} parm[PARM];
int lastparm=-1,maxparm=-1; /* NB: maximum maxparm is PARM-1 */
static char name[8]="aname";
double sliderpos; /* premistit do xdraw ! */
void minmaxset(int iparm)
{
  if (iparm>=0 && iparm<PARM) {
    parm[iparm].min=*parm[iparm].parm-parm[iparm].step;
    parm[iparm].max=*parm[iparm].parm+parm[iparm].step;
    parm[iparm].on=1; }
}

/* see FIT4PLOT/fit.c for info on parms */
char *FIT4PLOT; /* where the code for fitting is located */
double D=0,eps=0,err=-2;
int method=-2,maxit=10000,nerr=0;
#endif /*# PARM */

int percnr=0; /* counter if integer format in file name;
                 -0x7fffffff: if % does not mean format */
char *perc=NULL; /* !NULL if % found in file name */
char *numberedfile(char *template)
{
  if (perc)
    return string(template,percnr);
  else
    return template;
}

#include "plotmenu.c"

void roundrg(double *x,double *X,int enlarge) /********************* roundrg */
{
  double dx=*X-*x,q=1;

  if (dx<=0) return;
  while (dx>50) dx/=10,q*=10;
  while (dx<5) dx*=10,q/=10;

  *x=((int)(*x/q+0.5001*(Sign(*x)-enlarge)))*q;
  *X=((int)(*X/q+0.5001*(Sign(*X)+enlarge)))*q;
}

int smartcol(char **pT) /******************************************* smartcol */
/*
  changes column numbers 0,1,...,26 into @,A,...Z
  new: also 27 into c27 (then, allocates a new copy)
  01 changed into A, 1. or +1 etc. not changed (=expression)
  returns whether T started with 0 (turns off hotkeys 1..9)
  NB: @ is changed into n in to_colon() in tabinc.c
  BUG: for columns >26, a string is reallocated
       => and undo from hotkeys 1..0 does not work!
*/
{
  char *T=*pT;
  int i=atoi(T),r=T[0]=='0';
  char *c;

  for (c=T; *c; c++) if (!isdigit(*c)) return r; /* no digit => no change */
  if (i<=26) T[0]=i+'@',T[1]=0;
  else if (r) T[0]='c';
  else {
    alloc(T,strlen(T)+2);
    strcpy(T+1,*pT);
    T[0]='c';
    *pT=T; }

  return r;
}

char *coltok(char *s) /********************************************** coltok */
{
  static char *S;
  char *o;

  if (s) S=s;
  if (!S) return NULL;

  o=S;
  S=strchr(S,':');
  if (S) *S++=0;

  return o;
}

int main(int narg,char **arg) /**************************************** main */
{
  int iarg,I,I0=1,Iopt,advI=1,col,aspect=0;
  int fixcoly=0; /* see struct f_s */
  int pass=0,verbose=0,prtenv=0,batch=0,recompile=1;
  int difcol;
  static char colA[]="A",colB[]="B";
  char *colx=colA,*coly=colB;
  char *coldy1=NULL,*coldy2=NULL; /* error bars, either +dy-dy or +-dx+-dy */
  static struct rgfn_s {
    struct rgfn_s *next;
    char name[L_tmpnam];
  } *rgfn0,*rgfn;
  double x;
#define RGUNDO 16
  static struct rg_s
    rginit={DBL_MAX,-DBL_MAX,DBL_MAX,-DBL_MAX},
    rg    ={DBL_MAX,-DBL_MAX,DBL_MAX,-DBL_MAX},
    rgset ={DBL_MAX,-DBL_MAX,DBL_MAX,-DBL_MAX},
    rgundo[RGUNDO]={{DBL_MAX,-DBL_MAX,DBL_MAX,-DBL_MAX}};
  int lblx=0,lbly=0,lblshow=0,psty=0,lsty=1,fix=0,i;
  int colorsty=0,errsty=1,newerrsty;
  struct f_s {
    struct f_s *next;
    char *colx;    /* expression for the x-column */
    char *coly;    /* expression for the y-column */
    char *coldy1;  /* expression for the -dy1 error bar */
    char *coldy2;  /* expression for the +dy2 error bar (if none, dy2=dy1) */
    int errbar;    /* if error bars */
    int maxcol;    /* max. column (to be read in by drawfile or drawerrbar) */
    int fixcoly;   /* column with leading 0 => not changed for hotkeys 1..0 */
    int col;       /* color */
    int lsty;      /* linestyle */
    int psty;      /* pointstyle */
    int fix;       /* fixed style (SPACE inactive) */
    char *colon[5];/* pos of : in fn (or NULL) ??? OUT OF ORDER BECAUSE OF strdup */
    char *colydup; /* original value of column y */
    char fn[3];    /* filename [variable length]: strlen(fn) added */
  } *head=NULL,*f,*last=NULL;
  char *tok1,*tok2,*tok3,*toks;
  char *tok4,*tok5;
  char fn[ARGLEN], A[ARGLEN];
  int key,lastkey=-99999;
  static FILE *resp;
  int isstdin=0,iserrbar=1;

#if CHECKHEAP==-2 || CHECKHEAP==2
  AllocRange=1024;
#endif

  loop (i,1,RGUNDO) rgundo[i]=rgundo[0];

  initscroll(0);

  //  signal(SIGTERM,MYHUP);

  /*** print help if no argument ***/
  if (narg<2) {
    prtsfill(Pintro);

    for (;;) {
      char s[8];

      prts_(Phelp);
      if (!fgets(s,8,stdin)) return 0;
      switch (s[0]) {
        case 'i': prtsfill(Pintro); break;
        case 'g': prtsfill(Pgui); break;
        case 'o': prtsfill(Poptions); break;
        case 'r': prtsfill(Prange); break;
        case 'f': prtsfill(Pfiles); break;
        case 'c': prtsfill(Pcolumns); break;
        case 'l': prtsfill(Pcolors); break;
        case 's': prtsfill(Pstyle); break;
        case 'e': prtsfill(Penv); break;
        case 'm': prtsfill(Pmore); break;
        case 't': prtsfill(Pfit); break;
        case 'p': prtsfill(Pprint); break;
        case 'x': prtsfill(Pexpr); break;
        case 'G': prtsfill(Pgeo); break;
        default: return 0; } }
  } /* help */

  xwindowhints.geometry=getenv("PLOTGEOMETRY");
  xwindowhints.winname=getenv("PLOTNAME");
  xwindowhints.icowidth=plotbit_width;
  xwindowhints.icoheight=plotbit_height;
  xwindowhints.icobits=(char*)plotbit_bits;
  if ((tok1=getenv("PLOTINIT"))) addkbdstring((unsigned char*)tok1);

/* making a list of id's for Calc() */
#if PARM
  FIT4PLOT=getenv("FIT4PLOT");
  if (!FIT4PLOT) {
    char *USER=getenv("USER");
    if (!USER) USER=getenv("USERNAME");
    if (!USER) USER="plot";
    FIT4PLOT=dupstr(string("/home/%s/macsimus/c/fit4plot",USER)); }

  loop (i,0,PARM) {
    char *ge;
    static char cstep[8]="?step";

    alloc(id,sizeof(struct _Idlist_s)+1);
    id->next=_Id.head; _Id.head=id;
    parm[i].parm=&id->val;
    parm[i].step=1;
    parm[i].on=0;
    parm[i].used=0;
    id->id[0]='a'+i;
    id->val=0;
    id->used=0;
    if ( (ge=getenv(id->id)) ) {
      id->val=atof(ge);
      parm[i].on=1; }
    cstep[0]='a'+i;
    if ( (ge=getenv(cstep)) ) {
      parm[i].step=atof(ge);
      parm[i].on=1; }
  }
#endif /*# PARM */

  columnidlist();

  /* parse -options (must be first), skip empty arg */
  loop (iarg,1,narg) if (arg[iarg][0]) {
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'I': addkbdstring((unsigned char*)arg[iarg]+2); break;
        case 'a': aspect++; break;
        case 'b': batch++; break;
        case 'c': colorsty=256; break;
        case 'd': plotopt.slowdraw=atof(arg[iarg]+2)*1e6+0.5; break;
        case 'e': prtenv++; break;
        case 'h': percnr=-0x7fffffff;
                  if (arg[iarg][2]) percnr=atoi(arg[iarg]+2);
                  break;
        case 'g': xwindowhints.geometry=arg[iarg]+2; break;
        case 'k': plotopt.key=arg[iarg]+2; break;
        case 'n':
          if (arg[iarg][2]) {
            psdef=dupstr(string("%s.def",arg[iarg]+2));
            plotps=dupstr(string("%s.ps",arg[iarg]+2));
            ploteps=dupstr(string("%s.eps",arg[iarg]+2)); }
          else
            Error("empty option -n\n*** (for color change on blank line, use -c)");
          break;
        case 'N':
          if (arg[iarg][2]) xwindowhints.winname=arg[iarg]+2;
          break;
        case 'p': parent=atoi(arg[iarg]+2); break;
        case 't': recompile=0;
                  fprintf(stderr,"plot: WARNING: fitting function not recompiled because of option -t\n");
                  break;
        case 'v': verbose++; break;
        case 'z': rndinit(7,atoi(arg[iarg]+2)); break;
        case 0:
        case ':': goto exopt; /* -: is pipe */
        default: fprintf(stderr,"plot: option %s unknown\n",arg[iarg]); }
    else break; }

 exopt:
  Iopt=iarg;
  I0=iarg;

  if (!xwindowhints.winname) {  
    int l=6;

    loop (i,1,narg) if (arg[i][0]!='-') l+=1+strlen(arg[i]);
    alloc(xwindowhints.winname,l);
    strcpy(xwindowhints.winname,"plot");
    loop (i,1,narg) if (arg[i][0]!='-') {
      strcat(xwindowhints.winname," ");
      strcat(xwindowhints.winname,arg[i]); } }

  if (strlen(xwindowhints.winname)<16)
    xwindowhints.iconame=xwindowhints.winname;
  else {
    static char iconame[16];
    memcpy(iconame,xwindowhints.winname,15);
    xwindowhints.iconame=iconame; }

  if (parent>=0) {
    char *fn=string("/tmp/plot.%d",parent);
    FILE *f=fopen(fn,"at");
    if (f) { fprintf(f,"%d\n",getpid()); fclose(f); } }

  for (I=iarg=Iopt-1;;) {

    /* read arg from a response file */
    if (resp) do {
      if (!fgets(A,ARGLEN,resp)) {
        fclose(resp);
        resp=NULL;
        goto nextarg; }
      } while (A[0]=='!');
    else {
     nextarg:
      if (++iarg>=narg) break;
      if (arg[iarg][0]=='@' && arg[iarg][1]!=':') {
        /* @FILE is a response file */
        resp=fopen(arg[iarg]+1,"rt");
        if (!resp) fprintf(stderr,"plot: no file %s\n", arg[iarg]+1);
        continue; }
      else if (!arg[iarg][0]) /* ignore empty argument */ ;
      else {
        /* string starts with filename, [###], or @: (dummy filename) */
        strcpy(A,arg[iarg]); } }

    /* remove trailing whites  */
    while ( *(toks=strend1(A))<=' ' && toks>=A) *toks=0;

    I+=advI;

    if (iarg==Iopt && A[0]=='[' && *strend1(A)==']') {
      /* (optional) initial range: must be 1st arg and must be of form [...] */
      I0=iarg+1;
      toks=A+1;
      if (*toks!=':') rgset.x=atof(toks);
      if ( (toks=strchr(toks,':')) ) {
        toks++;
        if (!strchr("[],",*toks)) rgset.X=atof(toks);
        toks=strchr(tok1=toks,'[');
        if (!toks) toks=strchr(tok1,',');
        if (toks) {
          toks++;
          if (*toks!=':') rgset.y=atof(toks);
          if ( (toks=strchr(toks,':')) ) {
            toks++;
            if (*toks!=']' && *toks) rgset.Y=atof(toks); } } } }

    else {
      /* arg starts by a filename, @:, :, or a file of consecutive numbers [###] */
      if (A[0]=='[') {
        /* [###] (# of intervals and optional further parameters) */
        static int n=1000;
        static double r=DBL_MAX,R=-DBL_MAX;
        FILE *rgf;

        toks=A+1;
        if ( !(tok1=strchr(toks,']')) ) {
          fprintf(stderr,"plot: no matching ]\n");
          goto shut; }
        if (!strchr(":],",*toks)) n=atoi(toks);
        if ( (toks=strchr(toks,':')) && toks<tok1) {
          toks++;
          if (!strchr(":],",*toks)) r=atof(toks);
          if ( (toks=strchr(toks,':')) && toks<tok1) {
            toks++;
            if (!strchr(":],",*toks)) R=atof(toks); } }

        if (r>1e308) r=rgset.x;
        if (R<-1e308) R=rgset.X;

        if (r>1e308) r=rg.x;
        if (R<-1e308) R=rg.X;

        if (r>1e308 || R<-1e308 || n<=0 || n>1000000) {
          fprintf(stderr,"plot: no RANGE or bad # of points (max 1000000)\n");
          goto shut; }

        alloc(rgfn,sizeof(*rgfn));
        if (!tmpnam(rgfn->name)) Error("plot: cannot create temporary file of [...]");
        rgfn->next=rgfn0;
        rgfn0=rgfn;
        rgf=fopen(rgfn->name,"wt");
        loopto (i,0,n) fprintf(rgf,"%.20g\n", i*(R-r)/n+r);
        fclose(rgf);
        strcpy(fn,rgfn->name);
        tok1++;
        if (tok1[0]!=':') strcat(fn,":");
        strcat(fn,tok1); }
      else if (!last && A[0]==':') {
        /* first plot argument starts with : -> @: */
        fn[0]='@';
        strcpy(fn+1,A); }
      else
        /* A starts with : FILENAME or @ (dummy filename) */
        strcpy(fn,A);

      /* not counted for colors */
      if (fn[0]=='@') I--;

      i=sizeof(struct f_s)+strlen(fn);
      if (fn[0]==':') {
        /* :(columns) construct */
        i+=strlen(last->fn); /* only filename */
        alloczero(f,i);
        strcpy(f->fn,last->fn);
        strcat(f->fn,fn); }
      else {
        alloczero(f,i);
        strcpy(f->fn,fn); }

      if (f->fn[0]=='-' && (f->fn[1]==0 || f->fn[1]==':') && !isstdin) {
        /* pipe: file '-' created */
        FILE *ff;
        int ch;

        isstdin++;
        ff=fopen("-","wt");

        while ( (ch=getchar())>=0 ) fputc(ch,ff);
        fclose(ff); }

      if (percnr!=-0x7fffffff) perc=strchr(f->fn,'%');

      /* dirty hack: change FILE:Y into FILE::Y (NB: f->fn by 1 byte longer) */
      tok1=strchr(f->fn,':');
      // max_column  should have been set
      // this is wrong for plot :0:1 file.dat:
      // if (!tok1) max_column=2; else max_column=-1;
      if (tok1 && !strchr(tok1+1,':')) memmove(tok1+1,tok1,strlen(tok1)+1);

      tok1=coltok(f->fn); /* filename or [###], thrown away here */
      tok2=tok3=tok4=tok5=NULL;
      tok1=coltok(NULL);

      max_column=0; /* ??? -1 was here */

      if (tok1) {
        f->colon[0]=tok1-1;
        tok2=coltok(NULL);
        if (tok2) {
          f->colon[1]=tok2-1;
          tok3=coltok(NULL);
          if (tok3) {
            f->colon[2]=tok3-1;
            tok4=coltok(NULL);
            if (tok4)
              f->colon[4]=tok4-1;
            tok5=coltok(NULL); } } }

      if (tok1 && *tok1) {
        smartcol(&tok1);
        to_colon(tok1); // max_column set here
        colx=tok1; }

      if (tok2 && *tok2) {
        fixcoly=smartcol(&tok2);
        to_colon(tok2); // max_column set here, may reallocate the string */
        coly=tok2;
        f->colydup=dupstr(tok2); /* memory leak + dirty patch */
      }

      if (tok3 && *tok3) {
        int II=I0;

        psty=0; lsty=-1; fix=1;
        newerrsty=0;
        for (toks=tok3; *toks; toks++) switch (*toks) {
          case 's': fix=0; break;
          case 'x': if (newerrsty==1) newerrsty=3; else newerrsty=2; break;
          case 'y': if (newerrsty==2) newerrsty=3; else newerrsty=1; break;
          case 'n': colorsty=256-colorsty; break;
          case 'c': I=II++; break;
          case 'C': advI^=1; break;
          case '.': psty=-1; break;
          case 'p': psty=-2; break;
          case 'P': psty=-3; break;
          case 'o': psty=-4; break;
          case 'O': psty=-6; break;
          case '-': lsty=1; break;
          case '=': lsty=2; break;
          case 'd': lsty=0; break;
          default:
            if (isdigit(*toks)) {
              advI=0;
              II=atoi(toks);
              if (II>9) toks++;
              I=I0-1+II++; }
            else fprintf(stderr,"plot: `%c' bad key in STYLE\n",*toks); }
        if (newerrsty) errsty=newerrsty;
        if (psty==0 && lsty==-1) lsty=1; }

      if (tok4) {
        if (*tok4) {
          smartcol(&tok4);
          to_colon(tok4); // max_column set here
          coldy1=tok4; }

        /* FILE:A:B:-: removes previous setting */
        else coldy1=NULL; }

      if (tok5) {
        if (*tok5) {
          smartcol(&tok5);
          to_colon(tok5); // max_column set here
          coldy2=tok5; }

        /* FILE:A:B:-:C: removes previous setting */
        else coldy2=NULL; }

      /* missing : means :1:2 -- solved above
         if (max_column<0) max_column=2; */

      if (f->fn[0]!='@') {
        f->colx=colx;
        f->coly=coly;
        f->errbar=0;
        /* new version: plot :A:B:-:C FILEs now shows ALL error bars */
        if (coldy1)  {
          f->errbar++;
          f->coldy1=coldy1;
          if (coldy2) f->coldy2=coldy2;
          else f->coldy2=coldy1; }
        f->fixcoly=fixcoly;
        f->psty=psty;
        f->lsty=lsty;
        f->fix=fix;
        f->col=color[(I-I0)%14];
        f->maxcol=max_column;

        if (!minmaxfile(numberedfile(f->fn),colx,coly,f->maxcol, &rg.x,&rg.X,&rg.y,&rg.Y)) {
          if (head) last->next=f,last=f;
          else head=last=f;
          f->next=NULL;
          if (rgset.x< 1e308) rg.x=rgset.x;
          if (rgset.y< 1e308) rg.y=rgset.y;
          if (rgset.X>-1e308) rg.X=rgset.X;
          if (rgset.Y>-1e308) rg.Y=rgset.Y; }
        rginit=rg; } }
  } /* iarg & @resp-file */

  looplist (f,head) f->colydup=dupstr(f->coly);
  
  if (0) looplist (f,head) {
      fprintf(stderr,"DEBUG LIST1 %s x=%s y=%s max=%d :%p:%p\n",
              f->fn,f->colx,f->coly,f->maxcol,f->colon[0],f->colon[1]);
    }

#ifdef DEBUG
  fprintf(stderr,"plot: rgx=%.9g %.9g  rgy=%.9g %.9g\n",rg.x,rg.X,rg.y,rg.Y);
#endif /*# DEBUG */

#if 1
  if (rg.x>=rg.X || rg.y>=rg.Y) {
    fprintf(stderr,"plot: rgx=%.9g %.9g  rgy=%.9g %.9g\n",rg.x,rg.X,rg.y,rg.Y);
    fprintf(stderr,"not enough points (in given ranges)\n");
    goto shut; }
#endif /*# 1 */

  psty=0;
  lsty=1;
  repeatprefix('`');

  
 notify:

  if (display) {
    /* X display already on */
    selectfont(font);
    setmenusize();
    setviewport(0,0,maxxn-1,maxyn-1,1);
    clearviewport();
    setviewport(0,0,maxxn-1,maxyn-1,0);
    redrawmenu=1; }
  else {
    /* start graphics */
    if (startgraph(batch?0:-9)<0) exit(-1);
    //    if (getmaxy()<479) font=0; else font=2;
    font=2;
    if (getmaxy()<248 || getmaxx()<400) font=0;
    selectfont(font);
    if (getenv("GUI")) {
      menu=!!strchr(getenv("GUI"),'p');
      KILL=!!strchr(getenv("GUI"),'K'); }
    setmenusize(); }

  lwindow(lbl?LXWIDTH*xfont.width+3:4,maxxn-4,4,getmaxy()-(lbl?xfont.height+3:4));
  if (menu) makemenu(1);

  for (;;) { /* main loop */

#ifdef DEBUG
    fprintf(stderr,"drawing...\n");
#endif /*# DEBUG */

    if (rg.x>=rg.X) rg.x-=1,rg.X+=1;
    if (rg.y>=rg.Y) rg.y-=1,rg.Y+=1;

    scale(rg.x,rg.X,rg.y,rg.Y,aspect);
    setviewport(0,0,maxxn-1,maxyn-1,1);
    clearviewport();
    setwritemode(0);
    setfillstyle(1,BLACK);
    makegrid(1,rg.x,rg.X,rg.y,rg.Y);
    makegrid(0,rg.y,rg.Y,rg.x,rg.X);

    if (!menu) makemenu(0);

    if (my_style.ps && lbl) myclip(); /* contains gsave ! */

#ifdef PARM
    if (getenv("PLOTCMD")) {
      int ii;
      static char *env; /* max 32 B using format %c=%.16g */

      if (!env) allocarray(env,32*PARM);
      loop (ii,0,PARM) {
        sprintf(env+32*ii,"%c=%.16g",'a'+ii,parm[ii].parm[0]);
        putenv (env+32*ii); }
      if (system(getenv("PLOTCMD")))
        fprintf(stderr,"plot: PLOTCMD failed\n"); }
#endif /*# PARM */

    looplist (f,head) if (f->fn[0]!='@') {

      if (0) fprintf(stderr,"DEBUG LIST2 %s x=%s y=%s max=%d :%p:%p\n",
                     f->fn,f->colx,f->coly,f->maxcol,f->colon[0],f->colon[1]);

      if (f->errbar && iserrbar)
        drawerrbar(numberedfile(f->fn),
                   f->colx,f->coly,f->coldy1,f->coldy2,f->maxcol,
                   f->col|colorsty,f->psty<0?f->psty:psty<0?psty:-6,errsty);

      if (f->fix) {
        if (f->lsty>=0)
          drawfile(numberedfile(f->fn),
            f->colx,f->coly,f->maxcol, f->col|colorsty,f->lsty);
        if (f->psty<0)
          drawfile(numberedfile(f->fn),
            f->colx,f->coly,f->maxcol, f->col|colorsty,f->psty); }
      else {
        if (psty<0)
          drawfile(numberedfile(f->fn),f->colx,f->coly,f->maxcol,
                   f->col|colorsty,psty);
        if (lsty>=0)
          drawfile(numberedfile(f->fn),f->colx,f->coly,f->maxcol,
                   f->col|colorsty,lsty);
        if (lsty<0 && psty>=0)
          drawfile(numberedfile(f->fn), f->colx,f->coly,f->maxcol,
                   f->col|colorsty,0); }
      setlinestyle(0,0,1); /* default if style=1 (?)*/ }

#if PARM
    /* checking which parameters have been used in formulas */
    if (pass==0) {
      char id[2]="a";
      struct _Idlist_s *l;

      loop (i,0,PARM) {
        id[0]='a'+i;
        looplist (l,_Id.head)
          if (!strcmp(l->id,id) && l->used) parm[i].used=1; }
      pass=1;
      makemenu(1); }
#endif /*# PARM */

    if (lblshow) {
      int aty=lbly-xfont.top;
#if PARM
      char as[64],*env;
      int ii;

      Max(maxparm,lastparm)

      aty-=(maxparm+1)*(xfont.height+xfont.top+2);
      loopto (ii,0,maxparm) {
        name[0]='a'+ii;
        if ( (env=getenv(name)) )
          sprintf(as,"%c: %s = %g",'a'+ii,env,*parm[ii].parm);
        else
          sprintf(as,"%c=%.12g",'a'+ii,*parm[ii].parm);
        aty=outbgtextxy(lblx,aty,font,
                        parm[ii].used?WHITE:LIGHTGRAY,
                        ii==lastparm?RED:BLACK,as);
        if (my_style.ps)
          fprintf(my_style.ps,"%.2f %.2f moveto (%s) show\n",
                  lblx/(double)maxxn*my_style.xsize,
                  (1-lbly/(double)maxyn)*my_style.ysize+(maxparm-ii)*my_style.fontsize*1.15,
                  psshow(as));


      }
#endif /*# PARM */

      for (col=0,f=head; f; col++,f=f->next) if (f->fn[0]!='@') {
        int i;
        loop (i,0,5) if (f->colon[i]) *(f->colon[i])=':';
        mysetcolor(f->col);
        aty=outbgtextxy(lblx,aty,font,f->col,BLACK,numberedfile(f->fn));
        if (my_style.ps)
          fprintf(my_style.ps,"%.2f %.2f moveto (%s) show\n",
                  lblx/(double)maxxn*my_style.xsize,
                  (1-lbly/(double)maxyn)*my_style.ysize-(col+1)*my_style.fontsize*1.15,
                  psshow(numberedfile(f->fn)));
        loop (i,0,5) if (f->colon[i]) *(f->colon[i])=0; } }

#ifdef PARM
    if (menu && lastparm>=0 && lastparm<PARM) {
      int atx=maxxn+2,aty=maxyn-3*xfont.height-2, i,len;
      char s[32],ss[48],*env,senv[2];

      name[0]='a'+lastparm;
      env=getenv(name);

      if (!env || strlen(env)>8 ) {
        senv[0]=name[0]; senv[1]=0;
        env=senv; }

      sprintf(ss,"[%.14g,%.14g]",parm[lastparm].min,parm[lastparm].max);

      setviewport(0,0,maxxn-1,maxyn-1,0);
      setfillstyle(1,BLUE);
      bar(atx,aty,getmaxx()-2,maxyn-3);

      atx+=2;
      len=strlen(ss);
      if (len>18) aty+=3;
      else aty+=xfont.height/2+1;

      setcolor(parm[lastparm].used?WHITE:LIGHTGRAY);

      i=17; do { i--;
        sprintf(s,"%s=%.*g",env,i,*parm[lastparm].parm); }
      while (strlen(s)>18);
      outtextxy(atx,aty,s);

      aty+=xfont.height;
      if (len>18) {
        i=17; do { i--;
          sprintf(s,"[%.*g,",i,parm[lastparm].min); }
        while (strlen(s)>18);
        outtextxy(atx,aty,s);
        aty+=xfont.height;
        i=17; do { i--;
          sprintf(s,"%.*g]",i,parm[lastparm].max); }
        while (strlen(s)>18);
        outtextxy(atx,aty,s); }
      else
        outtextxy(atx,aty,ss);
      setviewport(0,0,maxxn-1,maxyn-1,1); }
#endif /*# PARM */

    if (my_style.ps) {
      struct lett_s *l;
      int it;

      myup();
      if (lbl) fprintf(my_style.ps,"grestore\n");
      loop (it,0,my_style.nticks) {
        myaxis(it,1,1);
        myaxis(it,1,-1);
        myaxis(it,-1,1);
        myaxis(it,-1,-1); }

      for (l=lett; l; l=l->next) {
        fprintf(my_style.ps,"\n%%----- %s -----\n",l->s);
        if (l->r)
          fprintf(my_style.ps,
                  "gsave\n%.2f %.2f translate %.2f rotate\n0 0 moveto\n",
                  l->x,l->y,l->r);
        else
          fprintf(my_style.ps,"%.2f %.2f moveto\n",l->x,l->y);
        it=0;
        if (l->s[0]=='$') {
          it=atoi(l->s+1)-1;
          if (it<0) {
            fprintf(stderr,"plot: %s: negative $color changed to 0\n",l->s);
            it=0; } }
        mypssetxy(l->x,l->y,color[it%14]);
        if (l->cmd) fprintf(my_style.ps,"%s\n",l->cmd);
        psstring(l->s,-1);
        if (l->r) fprintf(my_style.ps,"grestore\n"); }

      fprintf(my_style.ps,"grestore\n");
      mycloseplotps();
      my_style.ps=NULL; }

#ifdef DEBUG
    fprintf(stderr,"plot: drawn\n");
#endif /*# DEBUG */

/*.....   do { key=extch(); } while (key==X_REDRAW_ME);*/

    /* all extra X_REDRAW_ME (typically Expose) coming just now are ignored */
    /* This is not the best solution :-( */
    setviewport(0,0,maxxn-1,maxyn-1,0);

    if (lastparm>=0) {
      sliderpos=(*parm[lastparm].parm-parm[lastparm].min)/(parm[lastparm].max-parm[lastparm].min);
      drawslider(0); }

#if 0
    for (;;)
      if (kbhit()) {
        if ( (key=readkey())==X_REDRAW_ME) {
          continue; }
        else
          break; }
      else {
        key=readkey();
        break; }
#endif /*# 0 */

    //    if (!menu) makemenu(0);

    if (display) XFlush(display); /* added because of problems with empty graph */

    if (lastkey==F4) key='Q'; /* causes stop */
    else key=readkey();

#ifdef DEBUG
    fprintf(stderr,"plot: switch key=%d\n",key);
#endif /*# DEBUG */

    if (key=='u' || key==('Z'&31)) {
      if (rgundo[0].x<rgundo[0].X) {
        rg=rgundo[0];
        loop (i,1,RGUNDO) rgundo[i-1]=rgundo[i]; } }
    else {
      if (rgundo[0].x!=rg.x || rgundo[0].y!=rg.y || rgundo[0].X!=rg.X || rgundo[0].Y!=rg.Y) {
        for (i=RGUNDO-1; i>0; i--) rgundo[i]=rgundo[i-1];
        rgundo[0]=rg; }

      switch (key) {
        case 'O'&31:
        case F5:
          font=(font+2)%4; // to extend to font 4 : %6
          erasebuttons();
          redrawmenu=1;
          goto notify;

        case F10:
          menu=!menu;
          erasebuttons();
          redrawmenu=1;
          goto notify;

        case F9:
        case 'U'&31:
          rg.x=rg.y=DBL_MAX;
          rg.X=rg.Y=-DBL_MAX;

          looplist (f,head) if (f->fn[0]!='@')
            minmaxfile(numberedfile(f->fn),f->colx,f->coly,f->maxcol, &rg.x,&rg.X,&rg.y,&rg.Y);
          break;

        case '?': case F1: case 'h': case 'H':
          maxxn=getmaxx()+1;
          setviewport(0,0,maxxn-1,maxyn-1,1);
          help();
          setmenusize();
          redrawmenu=1;
          goto notify;

        case 'L'&31:
        case X_REDRAW_ME: /* event forcing redraw */
#ifdef DEBUG
          fprintf(stderr,"plot: redraw requested\n");
#endif /*# DEBUG */
          goto notify;
        case X_UNKNOWN_KEY:
          break; /* some more non-alpha keys */
        case 'K':
          if (parent>=0) {
            char *fn=dupstr(string("/tmp/plot.%d",parent));
            FILE *f=fopen(fn,"rt");
            char line[16];
            while (fgets(line,16,f)) {
              if (atoi(line)!=getpid())
                if (system(string("kill %s",line))) fprintf(stderr,"plot: kill %s failed",line); }
            fclose(f);
            //            fprintf(stderr,"%s\n",fn);
            remove(fn);
            return 0; }
          else if (KILL) {
            /* is dirty, because the caller is supposed to be killed */
            if (system("killall plot")) fprintf(stderr,"plot: cannot killall plot\n"); }
          break;

        /*** MOUSE ***/
        case LEFTCLICK:
#if PARM
          if (lastparm>=0 && isinslider()) {
            drawslider(0);
            *parm[lastparm].parm=parm[lastparm].min+sliderpos*(parm[lastparm].max-parm[lastparm].min);
            break; }
          else
#endif /*# PARM */
          if (dosmouse.x<maxxn) {
            /* print x,y [dx,dy] */
            char clipit[128],*name;
            double x=iSX(dosmouse.x),y=iSY(dosmouse.y);
            static double lastx,lasty;
            static int col=LIGHTCYAN;

            printf("%14.10g %14.10g %14.10g %14.10g CLICK\n",x,y,x-lastx,y-lasty);
            lastx=x; lasty=y;
            if ( (name=getenv("TOCLIP")) ) {
              if (name[0]=='-') sprintf(clipit,"echo \'%.9g %.9g\' | %s",x,y,name+1);
              else sprintf(clipit,"%s \'%.9g %.9g\'",name,x,y);
 	      if (system(clipit)) fprintf(stderr,"plot: %s failed\n",clipit); }
            sprintf(clipit,"%.9g %.9g",x,y);
            if (col==LIGHTCYAN) col=YELLOW; else col=LIGHTCYAN;
            setcolor(col);
            line(dosmouse.x-8,dosmouse.y,dosmouse.x+8,dosmouse.y);
            line(dosmouse.x,dosmouse.y-8,dosmouse.x,dosmouse.y+8);
            if (display) {
              XFlush(display);
              usleep(200000); } }
          break;

        case MIDCLICK:
          if (dosmouse.x<maxxn) {
            /* print x,y [dx,dy] */
            char clipit[128];
            double x=iSX(dosmouse.x),y=iSY(dosmouse.y);
            static int col=LIGHTCYAN;

            my_rg=rg;
            readpsdef();
            sprintf(clipit,"%.2f %.2f",mysx(x),mysy(y));
            prt("%s moveto () show",clipit);
            //            prt("l %s ",clipit);
            if (col==LIGHTCYAN) col=YELLOW; else col=LIGHTCYAN;
            setcolor(col);
            line(dosmouse.x-8,dosmouse.y,dosmouse.x+8,dosmouse.y);
            line(dosmouse.x,dosmouse.y-8,dosmouse.x,dosmouse.y+8);
            if (display) {
              XFlush(display);
              usleep(200000); } }
          break;

        case 'Q':
        case 'Q'&31:
          goto TheEnd;

        case 'q':
        case F12:
        case ESC:
          while (kbhit()) readkey();
          i=readkey();
          if (i=='q' || i==ESC || i==F12) goto TheEnd;
          break;

        case RIGHTCLICK: /* show info */
          lbly=dosmouse.y;
          lblx=dosmouse.x;
          lblshow=1;
          break;

        case MIDDRAG: {
          double q=dosmouse.dx/scaling.x;

          rg.x-=q; rg.X-=q;
          q=dosmouse.dy/scaling.y;
          rg.y+=q; rg.Y+=q; }
          break;

        case LEFTDRAG:
#if PARM
          if (lastparm>=0 && isinslider()) {
            drawslider(0);
            *parm[lastparm].parm=parm[lastparm].min+sliderpos*(parm[lastparm].max-parm[lastparm].min);
            break; }
          else
#endif /*# PARM */
          if (dosmouse.x<maxxn) {
            setcolor(WHITE);
            setwritemode(1);
            setlinestyle(0,0,1);
            for (;;) {
              rectangle(dosmouse.x0,dosmouse.y0,dosmouse.xl,dosmouse.yl);
              rectangle(dosmouse.x0,dosmouse.y0,dosmouse.x,dosmouse.y);
              if (dosmouse.drag<2) {
                key=readkey();
                if (key!=LEFTDRAG) break; /* should never happen */ }
              else {
                if (dosmouse.x0!=dosmouse.x && dosmouse.y0!=dosmouse.y) {
                  rg.x=iSX(dosmouse.x0); rg.X=iSX(dosmouse.x);
                  if (rg.x>rg.X) x=rg.x, rg.x=rg.X, rg.X=x;

                  rg.y=iSY(dosmouse.y0); rg.Y=iSY(dosmouse.y);
                  if (rg.y>rg.Y) x=rg.y, rg.y=rg.Y, rg.Y=x; }
                break; } }
            break; }

        case RIGHTDRAG: {
          double x,y;

          x=(rg.X-rg.x)*dosmouse.dx*(2./getmaxx()); rg.X-=x; rg.x+=x;
          y=(rg.Y-rg.y)*dosmouse.dy*(2./getmaxy()); rg.Y+=y; rg.y-=y; }
          break;

        case 'U':
          looplist (f,head) strcpy(f->coly,f->colydup);
          break;

        case 'k':
          rg=rginit;
          break;

        case '.': psty=lsty=-1; break;

#ifdef PARM
        case F8: {
          //          loopto (i,0,maxparm) *parm[i].parm=0;
          loop (i,0,PARM-1) *parm[i].parm=0;
          maxparm=lastparm=-1;
          goto notify; }
        case 'N'&31: /* key duplicated */
        case 136: nerr=0; break;
        case 137: nerr=10; break;
        case 138: nerr=100; break;
        case 139: nerr=999; break;
        case 'N': if (nerr==0) nerr=200; else nerr+=12+nerr*3;
        case 'n': if (nerr>5) nerr/=2; else nerr=0;
          prt("nerr=%d (# of samples to calc. stderr of fit) will be used by hotkeys t ^T T",nerr);
 	  break;
        case 'T'&31:
        case 't':
        case 'T': {
	  FILE *sh=fopen("plotfit.sh","wt");
	  int np=0,prec=(key=='t')+2*(key=='T');
	  char *c;
          static char *precname[3]={"double","long","high"}; /* [key] */
	  static int lastnp[3]={-1,-1,-1}; /* [key] */
          static int defaultmethod[3]={0,0,2}; /* [key] */
          static int defaultmaxit[3]={10000,10000,12}; /* [key] */
          static double defaultD[3]={2e-5,1e-5,1e-8}; /* [key] */
          static double defaulteps[3]={1e-7,1e-8,1e-13}; /* [key] */

	  if (!head->next) break;

	  fprintf(sh,"#!/bin/sh\n\
# this script has been written by plot\n");
	  fprintf(sh,"echo method=%d par=1 D=%g eps=%g maxit=%d nerr=%d err=%g > fit.get\n",
                  method==-2?defaultmethod[prec]:method,
                  D==0?defaultD[prec]:D,
                  eps==0?defaulteps[prec]:eps,
                  method==-2?defaultmaxit[prec]:maxit,
                  nerr,err);

	  /* determine np */
	  for (c=head->next->coly; *c; c++)
	    if (*c>='a' && *c<='`'+PARM) {
	      /* to extend if _ used */
	      if (isalnum(c[1])) continue; /* part of longer id */
	      if (c>head->next->coly) if (isalnum(c[-1]) || c[-1]=='.') continue; /* part of exp number */
	      Max(np,*c-'`') }

 	  loop (i,0,np)
	    fprintf(sh,"echo \"A[%d]=%.15g\" >> fit.get\n",i,*parm[i].parm);
	  fprintf(sh,"echo \"; method=9;\" >> fit.get\n");
	  fprintf(sh,"cd %s\n",FIT4PLOT);
	  if (np!=lastnp[prec] && recompile) {
	    fprintf(sh,"echo \"// generated by plot\" > fitinc.c\n");
	    fprintf(sh,"echo \"#define N %d\" >> fitinc.c\n",np);

            fprintf(sh,"echo >> fitinc.c\n");
            fprintf(sh,"echo 'REAL func(REAL *P,REAL A) {' >> fitinc.c\n");
	    loop (i,0,np)
              fprintf(sh,"echo \"#define %c (P[%d])\" >> fitinc.c\n",i+'a',i);
            fprintf(sh,"echo return >> fitinc.c\n");
            fprintf(sh,"echo -n '  ' >> fitinc.c\n");
	    fprintf(sh,"exp2pow \"%s\" >> fitinc.c\n",head->next->coly);
	    fprintf(sh,"echo ';}' >> fitinc.c\n");

	    fprintf(sh,"echo >> fitinc.c\n");
            fprintf(sh,"echo 'REAL expr(REAL *P,REAL A) {' >> fitinc.c\n");
            fprintf(sh,"echo return >> fitinc.c\n");
            fprintf(sh,"echo -n '  ' >> fitinc.c\n");
            if (head->next->next)
              fprintf(sh,"exp2pow \"%s\" >> fitinc.c\n",head->next->next->coly);
            else
              fprintf(sh,"echo 0 >> fitinc.c\n");
	    fprintf(sh,"echo ';}' >> fitinc.c\n");

	    fprintf(sh,"echo >> fitinc.c\n");
            fprintf(sh,"echo 'REAL (*funcexpr)(REAL *P,REAL A)=%s;' >> fitinc.c\n\n",head->next->next?"expr":"func");

	    fprintf(sh,"echo >> fitinc.c\n");
	    loop (i,0,np)
              fprintf(sh,"echo \"#undef %c\" >> fitinc.c\n",i+'a');

	    fprintf(sh,"cp makefile.%s makefile\n\
rm fit.%s\n\
if ! make fit.%s ; then echo plotfit.sh: COMPILATION ERROR ; exit; fi\n", precname[prec],precname[prec],precname[prec]);
	    lastnp[prec]=np; }
	  fprintf(sh,"cd -\n");
	  fprintf(sh,"%s/fit.%s %s < fit.get\n",FIT4PLOT, precname[prec], head->fn);
 	  fclose(sh);
	  if (system("bash plotfit.sh")) sh=NULL;
          else sh=fopen("fit.fit","rt");
	  if (sh) {
	    char line[32];

	    loop (i,0,np) {
	      if (!fgets(line,32,sh)) Error("plot: read fit.fit");
	      sscanf(line+4,"%lf",parm[i].parm); }
	    fclose(sh); }
	  break; }
#endif /*# PARM */

        /* DUMP */
#define FNAME "plot0000"
        case F2:
        case 'P'&31: {
          /* see also 'P' of blend, note that dark-on-white differs! */
          static char fn[32]=FNAME,fn0[32]; /* 16 enough, but the compiler shouts loud.. */
          FILE *f;
          int x,y,wbg,c,ch;
          unsigned4 pix4;
          unsigned char *pix=(unsigned char *)&pix4;
          unsigned char col[16][3]={
            {0,0,0}, {0,0,127}, {0,127,0}, {0,127,127},
            {127,0,0}, {127,0,127}, {127,127,0}, {153,153,153},
            {77,77,77}, {0,0,255}, {0,255,0}, {0,255,255},
            {255,0,0}, {255,0,255}, {255,255,0}, {255,255,255} };
          XImage *ximage=XGetImage(display,win,0,0,maxxn,maxyn,AllPlanes,ZPixmap);
          unsigned long xcol;
          extern long unsigned *xcoltab;

          if (menu) {
            makemenu(2);
            redrawmenu=1; }
          else {
            fprintf(stderr,"\
plot: PrintScreen:\n\
in the graphical window, select output format:\n\
  P=PPM      p=PPM inverted (white background)\n\
  E=EPS      e=EPS inverted\n\
  O=gray PS  o=gray PS inverted\n\
  C=PS       c=PS inverted\n"); }

          ch=readkey();

          if (!strchr("pPoOcCeE",ch)) goto notify;
          wbg=islower(ch);
          ch=tolower(ch);

          strcpy(fn+8,".ppm");

          fprintf(stderr,"plot: writing %s\n",fn);
          f=fopen(fn,"wb");
          if (wbg) loop (x,0,16) loop (y,0,3) col[x][y]=~col[x][y];

          fprintf(f,"P6\n%d %d\n255\n",maxxn,maxyn);
          loop (y,0,maxyn) loop (x,0,maxxn) {
            xcol=XGetPixel(ximage,x,y);
            loop (c,0,16) if (xcoltab[c]==xcol) break;
            if (c>=16) c=0;

            copy(pix,col[c],3);
            fwrite(pix,3,1,f); }

          fclose(f);

          c=0;
          strcpy(fn0,fn);
          switch (ch) {
            case 'e':
              strcpy(fn+8,".eps");
              c=system(string("ppm2ps -e -c %s %s",fn0,fn));
              break;
            case 'c':
              strcpy(fn+8,".ps");
              c=system(string("ppm2ps -c %s %s",fn0,fn));
              break;
            case 'o':
              strcpy(fn+8,".ps");
              c=system(string("ppm2ps -2 -c %s %s",fn0,fn));
              break; }

          if (c)
            fprintf(stderr,"plot: Conversion to %s failed\nCheck whether \'ppm2ps\' is installed\n",fn);
          if (ch!='p') unlink(fn0);

          sprintf(fn+5,"%03d",atoi(fn+5)+1);
          makemenu(1);
          break; }

        case F4: /* will stop later */
        case F3:
        case '#':
          my_rg=rg;
          readpsdef();

          myopenplotps(my_style.mode=='E'?ploteps:plotps, xwindowhints.winname);
          break;

        case 'p': if (psty<-1) psty++; break;
        case 'P': psty--; break;
        case 'L': if (lsty<4) lsty++; break;
        case 'l': if (lsty>0) lsty--; break;
        case ' ':
          if (lsty==-1) lsty=1,psty=-4;
          else if (psty<0) psty=0,lsty=1;
          else lsty=-1,psty=-3;
          break;
        case '=':
#if PARM
          {
            double a=*parm[0].parm,b=*parm[1].parm,c=*parm[2].parm,d=*parm[3].parm,e=*parm[4].parm,f=*parm[5].parm,g=*parm[6].parm,h=*parm[7].parm,i=*parm[8].parm,j=*parm[9].parm;
            double astep=parm[0].step,bstep=parm[1].step,cstep=parm[2].step,dstep=parm[3].step,estep=parm[4].step,fstep=parm[5].step,gstep=parm[6].step,hstep=parm[7].step,istep=parm[8].step,jstep=parm[9].step;
            static double x;
            struct _Idlist_s *oldhead=_Id.head;
            _Id.head=NULL;
            closegraph();
            prt("Enter data in the getdata format:\n\
  a..j = the parameters (of formula, fit, etc.)\n\
  astep..jstep = sensitivity (by aA, etc.)\n\
  x = auxiliary variable\n\
  method = -2: default ([std]=t and [fast]: amoeba, [high]=t: Newton-Raphson)\n\
           -1: greedy, 0: amoeba, 1: conj.grad 2: Newton-Raphson 3: MC\n\
  eps = precision (default=0=precision-dependent)\n\
  D = step for numerical derivatives (default=0=precision-dependent)\n\
  maxit = minimization steps (-> precision-dependent for method=-2)\n\
  nerr = to calculate errors of parameters by sampling (hotkeys n N)\n\
  err = -1: uses stdev of data provided (assuming 1 if not provided)\n\
        -2: as err=-1 if stdev of data available, err=0 otherwise (default)\n\
        -3: stdev of dy_i is proportional to 3rd column,\n\
            rescaled by stdev calculated from residual sum of squares\n\
            (for data with known weights~1/dy_i^2 but unknown precision)\n\
       >=0: uses stdev calc. from resid. sum sq., multiplied by y^err:\n\
         0: the same abs.err. (of stdev above) of all data\n\
         1: the same rel.err. of all data\n\
       0.5: error is proportional to y^0.5 (suited for histogram data)\n\
  end data by ; (control returns to the plot window)");
            getdata
	      get(x)
              get(a) get(b) get(c) get(d) get(e) get(f) get(g) get(h) get(i) get(j)
              get(astep) get(bstep) get(cstep) get(dstep) get(estep) get(fstep) get(gstep) get(hstep) get(istep) get(jstep)
              get(method) get(maxit) get(D) get(eps) get(nerr) get(err)
              /* no get(par) because given by t/T */
              checkdata
            enddata
            *parm[0].parm=a,*parm[1].parm=b,*parm[2].parm=c,*parm[3].parm=d,*parm[4].parm=e,*parm[5].parm=f,*parm[6].parm=g,*parm[7].parm=h,*parm[8].parm=i,*parm[9].parm=j;
            parm[0].step=astep,parm[1].step=bstep,parm[2].step=cstep,parm[3].step=dstep,parm[4].step=estep,parm[5].step=fstep,parm[6].step=gstep,parm[7].step=hstep,parm[8].step=istep,parm[9].step=jstep;
            _Id.head=oldhead;
            startgraph(batch?0:-9);
            break;
          }
#endif /*# PARM */
        case '\\': psty=0,lsty=1; break;
        case ':': case ';': psty=0; lsty=0; break;
        case 'z': lbl=(lbl+1)%3; break;
        case 'Z': lbl=(lbl+2)%3; break;

        case F11:
        case '^': iserrbar=!iserrbar; break;

        case 'R'&31: /* ^R */
          key=0;
          goto round;

        case 'R': key=1;
      round:
          roundrg(&rg.x,&rg.X,key);
          roundrg(&rg.y,&rg.Y,key);
          //        case 'L'&31:
        case 'r': break;

        case '+': x=(rg.X-rg.x)*0.2; rg.X-=x; rg.x+=x;
                  x=(rg.Y-rg.y)*0.2; rg.Y-=x; rg.y+=x; break;
        case '-': x=(rg.X-rg.x)*0.5; rg.X+=x; rg.x-=x;
                  x=(rg.Y-rg.y)*0.5; rg.Y+=x; rg.y-=x; break;

        case 'X': x=(rg.X-rg.x)*0.45; rg.X-=x; rg.x+=x; break;
        case 'x': x=(rg.X-rg.x)*4.5; rg.X+=x; rg.x-=x; break;
        case 'x'&31:
          rg.x=fabs(rg.x); rg.X=fabs(rg.X);
          if (rg.X>rg.x) rg.x=-rg.X; else rg.X=rg.x,rg.x=-rg.x;
          break;

        case WHEELFORWARD: x=(rg.Y-rg.y)*0.1; rg.Y-=x; rg.y+=x; break;
        case WHEELBACKWARD: x=(rg.Y-rg.y)*0.125; rg.Y+=x; rg.y-=x; break;

        case 'Y': x=(rg.Y-rg.y)*0.45; rg.Y-=x; rg.y+=x; break;
        case 'y': x=(rg.Y-rg.y)*4.5; rg.Y+=x; rg.y-=x; break;
        case 'y'&31:
          rg.y=fabs(rg.y); rg.Y=fabs(rg.Y);
          if (rg.Y>rg.y) rg.y=-rg.Y; else rg.Y=rg.y,rg.y=-rg.y;
          break;

        case RIGHT:
          x=(rg.X-rg.x)*0.1; rg.X+=x; rg.x+=x; break;
        case LEFT:
          x=(rg.X-rg.x)*-0.0618; rg.X+=x; rg.x+=x; break;
        case UP:
          x=(rg.Y-rg.y)*0.1; rg.Y+=x; rg.y+=x; break;
        case DOWN:
          x=(rg.Y-rg.y)*-0.0618; rg.Y+=x; rg.y+=x; break;

        case INS:
        case 'W':
          if (lblshow) {
            lblx=(lblx+33)%(maxxn*3/4);
            lbly=(lbly+13)%(maxyn-xfont.height); }
          else {
            lblshow=1;
            lblx=lbly=10; }
          break;
        case DEL:
        case 'w':
          lblshow=0; break;

        case PGDN:
        case '{': if (perc) { percnr--; redrawmenu=1; goto notify; } break;
        case PGUP:
        case '}': if (perc) { percnr++; redrawmenu=1; goto notify; } break;

        case '[': difcol=-1; goto plus;
        case ']': difcol=1;
      plus:
          column=0;
          looplist (f,head)
            if (!f->fixcoly) {
              if (isupper(f->coly[0])) {
                f->coly[0]+=difcol;
                Max(f->coly[0],'A')
                Min(f->coly[0],'Z') }
              if (!column) column=f->coly[0]; }
          if (menu) goto notify;
          else break;

        case 'v':
          prt("x-range = [%.9g %.9g] = %.9g",rg.x,rg.X,rg.X-rg.x);
          prt("y-range = [%.9g %.9g] = %.9g",rg.y,rg.Y,rg.Y-rg.y);
#if PARM
          {
            loop (i,0,PARM)
              if (parm[i].on) prt("%c=%.16g",i+'a',*parm[i].parm);
            loop (i,0,PARM)
              if (parm[i].on) prt("%cstep=%.12g range=[%g %g]",i+'a',parm[i].step,parm[i].min,parm[i].max);
          }
#endif /*# PARM */
          break;

#if PARM
        case 'V': {
          FILE *f=fopen("plotenv.sh","wt");
          loop (i,0,PARM)
            if (parm[i].on) fprintf(f,"export %c=%.16g\n",i+'a',*parm[i].parm);
          fclose(f); }
          break;

          /* see also show/old+misc/plot-slash+star.c w/o recentering */
        case '/': if (lastparm>=0 && lastparm<PARM) parm[lastparm].step*=.5;
          minmaxset(lastparm);
          continue;

        case '*': if (lastparm>=0 && lastparm<PARM) parm[lastparm].step*=2;
          minmaxset(lastparm);
          continue;

        case '<': if (lastparm>=0) *parm[lastparm].parm-=parm[lastparm].step;
          minmaxset(lastparm);
          break;
        case '>': if (lastparm>=0) *parm[lastparm].parm+=parm[lastparm].step;
        case '|':
          minmaxset(lastparm);
          break;

        case 'o': if (lastparm>=0) *parm[lastparm].parm=0;
          minmaxset(lastparm);
          break;
#endif /*# PARM */

        default: if (isdigit(key)) {
          if (key=='0') {
            char r=readkey();

            if (isupper(r)) key=r;
            else if (strchr("xyz",r)) key=r+'A'-'x';
            else if (r=='n') key='n';
            else if (isdigit(r)) key=r*10+readkey()+('@'-11*'0');
            else key='0';

            if (key>'Z' && key!='n') key='0'; }
          else
            key+='@'-'0';

          looplist (f,head) {
            f->coly[0]=key; /* only 1 char -> max 26 columns */
            f->coly[1]=0; }
          } /* default */
#if PARM
          else if ( (key>=('A'&31) && key<('A'&31)+PARM)
                    || (key>='a' && key<'a'+PARM)
                    || (key>='A' && key<'A'+PARM)
                    || (key>=141 && key<141+PARM) ) {

            int newparm;

            if (key>=141) {
              newparm=key-141; }
            else if (key>='a') {
              newparm=key-'a';
              *parm[newparm].parm-=parm[newparm].step; }
            else if (key>='A') {
              newparm=key-'A';
              *parm[newparm].parm+=parm[newparm].step; }
            else {
              newparm=key-('A'&31);
              *parm[newparm].parm=0; }

            minmaxset(newparm);
            if (newparm!=lastparm) {
              lastparm=newparm; goto notify; }
          }
#endif /*# PARM */
        }
      lastkey=key; }
    }

 TheEnd:
  closegraph();
  colx=getenv("GUI");
  if (verbose || (colx && strchr(colx,'v')) ) {
    prt("x-range = [%.9g %.9g] = %.9g",rg.x,rg.X,rg.X-rg.x);
    prt("y-range = [%.9g %.9g] = %.9g",rg.y,rg.Y,rg.Y-rg.y);
#if PARM
    loop (i,0,PARM) if (parm[i].on) prt("%c=%.12g",i+'a',*parm[i].parm);
#endif /*# PARM */
  }

 shut:
  if (isstdin) unlink("-");
  for (rgfn=rgfn0; rgfn; rgfn=rgfn->next) remove(rgfn->name);

  if (prtenv) {
    FILE *outf=fopen("plot.env","wt");
    loop (i,0,PARM) if (parm[i].on) fprintf(outf,"export %c=%.15g\n",i+'a',*parm[i].parm);
    fclose(outf); }

  return 0;
}
