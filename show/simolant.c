/* KEY=.MORE.outtext.hh.makebutton. kostnice2 -8 -u < simolant.utf8.c > simolant.c ; make simolant
   for LANG=2 version (Czech without accents):
   KEY=.MORE.outtext.hh.makebutton. kostnice2 -8 -a < simolant.utf8.c >| simolant.c ; make simolant
   compile all versions: simolant.sh
   ====================================================================
   >>>>> DO NOT EDIT simolant.c: USE simolant.utf8.c AND CONVERT <<<<<<
   ====================================================================
 */

/* have a look at simolant.sh before editing this line: */
#define LANG 2 /* 0=English
                  1=Czech (iso-8859-2 fonts)
                  2=Czech without accents */
//#define DEBUG /* timing + hotkey codes + hotkeys i/I INTERNALDELAY */

#include "ground.h"
#include "xdraw.h"
#include "alloc.h"
#include "rndgen.h"
#include "getdata.c"

int th=5; /* line thickness, 5=ximage discs - cf. hotkey a */

void redraw();
void drawcfg();

/* if slow drawing so slow that accepting events is delayed */
int INTERNALDELAY=10;
int internaldelay=0;
int demo;

XImage *ximage,*ximage0;
extern GC gc;
extern long unsigned *xcoltab;

int maxxn,maxyn;

typedef float myreal;

myreal gravity;  /* acceleration of gravity */
myreal walldens=1.5; /* density of smoothed atoms on walls */
#define MAXN 1000     /* max # of atoms */

typedef struct {
  myreal x,y;
} vector;

#if LANG==0
#  include "simhelpen.c"
#elif LANG==1
#  include "simhelpcz.c"
#elif LANG==2
#  include "simhelpcz0.c"
#else /*#!LANG==0!LANG==1!LANG==2 */
#  error unknown LANG
#endif /*#!LANG==0!LANG==1!LANG==2 */

vector r[MAXN],v[MAXN],a[MAXN]; /* configuration, velocities, accelerations */
vector rt,ri,incirc; /* trial position, r[i], random vector */
struct {
  int ix,iy;
} scr[MAXN]; /* screen coordinates: to erase old atom positions */
int n=300; /* # of atoms */
int nn;
int ii,i,j,k,rad,maxrad,nrep;
myreal L,Lh,Lorig; /* box size, Lh=L/2, Lorig=starting value of L */
myreal T; /* temperature as NVT parameter */
myreal Tk,Tm,TT; /* Tk=kinetic T, Tm is summed during a sweep, TT is auxiliary */
myreal Tmax=2.593743;
myreal deltaU; /* MC: Utrial-U */
myreal bag; /* Creutz daemon bag */
myreal ff,vv,rr,xx,maxvv;
myreal qT=1.1;
myreal t=0;

enum bctype {BOX,SLIT,CYCLIC} bc=BOX; /* boundary conditions */

enum cfgtype {GAS,LIQUID,CRYSTAL,DEFECT,VACANCY,INTERSTITIAL,CAPILLARY,DUMP} cfg=GAS; /* initial cfg. */

enum stat_e {half,one/*,two*/,Lhalf,Auto,fixed} dstat; /* status of MC displacement */
myreal dadj; /* displacement for MC: automatically adjustable */
myreal d; /* displacement for MC: actual */
int accepted; /* # of accepted MC moves */
myreal ar; /* acceptance ratio */

myreal acc=0.3333333; /* acc ratio for auto set */
myreal tau=1; /* MD thermostat time constant */
myreal range=2.25; /* range^2 for counting neighbors (hotkey 'n') */

myreal h=0.01; /* MD timestep */
int hauto=0; /* status of MD timestep */
myreal hadj=0.01; /* MD timestep: automatically adjustable */
myreal P=0; /* pressure, as measured from the virial of force (only if pressure) */
int pressure=1;

int md; /* true: MD, false:MC */
int mauto; /* MC first, then switch to MD */
int constt; /* true:NVT, false:NVE */

myreal sx,sy; /* scaling for the screen output */
int aspx,aspy; /* pixel aspect ratio */
int col;

char ch;
char st[32];
#ifdef DEBUG
double ttime0,ttime; /* to control measure sampling rate, in seconds */
#endif /*# DEBUG */
int sweeps; /* # of sweeps / one cycle (lasting ~ sample s) */
#define MSWEEPS 32
int msweeps=MSWEEPS;
int hot=1;
int ndelay=1;

FILE *plb;

#ifdef DEBUG
#  include <time.h>

#  ifdef CHEAPTIME
/* real time in resolution by 1 s */
double mytime(void) /************************************************ mytime */
{
  time_t t0;
  time(&t0);
  return (double)t0;
}
#  else /*# CHEAPTIME */
#    include <sys/time.h>
/* real time in better resolution (10^-6 s if provided by the system) */
double mytime(void) /************************************************ mytime */
{
  struct timeval tv;
  struct timezone tz;

  gettimeofday(&tv,&tz);

  return (unsigned)tv.tv_usec/1e6+tv.tv_sec;
}
#  endif /*#!CHEAPTIME */
#endif /*# DEBUG */

/* atom-atom potential u(r)=r^8-r^4 (2D Lennard-Jones), as the FUNCTION of rr=r^2 */
myreal u(myreal rr)
{
  rr=1/Sqr(rr);
  return Sqr(rr)-rr;
}

/* site-site forces (to be multiplied by vector r), as the FUNCTION of rr=r^2 */
myreal f(myreal rr)
{
  myreal r4;

  r4=Sqr(rr);

  return (8/r4-4)/rr/r4;
}

/* attractive walls: site-wall potential */
myreal uwalla(myreal r)
{
  myreal rr;

  if (r<=0) return -1e30*r;
  else {
    rr=Sqr(r);
    return (walldens*0.375/Sqr(rr)-walldens*0.5)/(r*rr);
  }
}

/* attractive walls: site-wall force (to be multiplied by vector r) */
myreal fwalla(myreal r)
{
  myreal rr;

  rr=Sqr(Sqr(r));
  return (walldens*2.625/rr-walldens*1.5)/rr;
}

/* repulsive walls: site-wall potential */
myreal uwallr(myreal r)
{
  myreal rr;

  if (r<=0) return -1e30*r;
  else {
    rr=Sqr(r);
    return (walldens*0.375/Sqr(rr))/(r*rr); }
}

/* repulsive walls: site-wall force (to be multiplied by vector r) */
myreal fwallr(myreal r)
{
  myreal rr;

  rr=Sqr(Sqr(r));

  return (walldens*2.625/rr)/rr;
}

typedef myreal (*walltype)(myreal r);
walltype uwallx,fwallx,uwally,fwally,uwallxL,fwallxL,uwallyL,fwallyL;

/* squared vector in nearest image convention */
myreal Sqrni(myreal r)
{
  return Sqr(Lh-fabs(Lh-fabs(r)));
}

myreal histgrid=10.0;
#define histmax 50
unsigned hist[histmax+1];

void addP(double rr)
{
  int ir;

  P=P+rr*f(rr);
  ir=(int)(sqrt(rr)*histgrid);
  if (ir<=histmax) hist[ir]++;
}

/* support to print the hotkey panel */
int atx,aty,tox;

/* support to print info on variables */
void info(char *a, myreal x, int c,int y)
{
  char st[32];

  if (y>=0) {
    if (x==0) strcpy(st,"0.0");
    else if (a[1]=='g') sprintf(st,"%4.2f",x);
    else sprintf(st,x<0.01?"%4.4f":"%4.3f",x); }
  else {
    y=0;
    st[0]=0; }

  y=y*19+1;
  setcolor(c);
  selectfont(2);
  setfillstyle(1,BLACK);
  bar(atx,y,tox,y+15);
  outtextxy(atx,y,string("%s%s",a,st));
  selectfont(0);
}

/* print the current gravity info */
void prtgravity(myreal add)
{
  if (bc!=CYCLIC) {
    gravity=gravity+add;
    if (fabs(gravity)<0.001) gravity=0;
    aty=40;
    tox=getmaxx()-50;
    info(" g=",gravity,WHITE,3); }
}

/* draws walls: dotted=penetrable (cyclic)
               solid green = attractive
               solid red = repulsive */
void drawwalls(void)
{
  setlinestyle(0,0,th);
  if (bc==BOX) {
    if (uwallx==uwalla) setcolor(LIGHTGREEN);
    else setcolor(LIGHTRED); }
  else {
    setlinestyle(1,0,th); setcolor(LIGHTGRAY); }
  line(0,0,0,(int)(sy*L));
  if (bc==BOX) {
    if (uwallxL==uwalla) setcolor(LIGHTGREEN);
    else setcolor(LIGHTRED); }
  line((int)(sx*L),0,(int)(sx*L),(int)(sy*L));

  setlinestyle(0,0,th);
  if (bc!=CYCLIC) {
    if (uwally==uwalla) setcolor(LIGHTGREEN);
    else setcolor(LIGHTRED); }
  else { setlinestyle(1,0,th); setcolor(LIGHTGRAY); };
  line(0,0,(int)(sx*L),0);
  if (bc!=CYCLIC) {
    if (uwallyL==uwalla) setcolor(LIGHTGREEN);
    else setcolor(LIGHTRED); }
  line(0,(int)(sy*L),(int)(sx*L),(int)(sy*L));
}

/* support to print the help screen */
void hh(char *s)
{
  outtextxy(3,aty,s);
  aty+=xfont.height+3;
}

#if LANG==0
#  define MORE "once more to quit"
#else /*# LANG==0 */
#  define MORE "stiskni jeste jednou"
#endif /*#!LANG==0 */

/* finish, with optional ASCII dump of cfg */
void finish(void)
{
  int i;

  fprintf(stderr,"MD t=%g\n",t);

  setcolor(YELLOW);
  settextstyle(0,0,1);
  outtextxy(2,2,MORE);
  XFlush(display);

  if (ESC!=readkey()) {
    setcolor(BLACK);
    outtextxy(2,2,MORE);
    return; }

  closegraph();
  if (plb) fclose(plb);

  if (0) if (demo || yes("write configuration to file simolant.dup")) {
    out=fopen("simolant.dup","wt");
    if (!out) Error("cannot write to simolant.dup");
    put3(n,L,walldens)
    put3(T,Tmax,bc)
    put2(h,dadj)
    prt("; ! x         y");
    loop (i,0,n) prt("%9.5f %9.5f",r[i].x,r[i].y);
    fclose(out); }
  exit(0);
}

void makebuttons(void)
{
  aty=253;
#define DY (xfont.height+15)
#define DYS (xfont.height+12)

   setcolor(LIGHTGREEN);
   outtextxy(atx,getmaxy()-41,"(c) J.Kolafa 1993-2013");
   outtextxy(atx,getmaxy()-27,"www.vscht.cz/fch/");
   outtextxy(atx,getmaxy()-13,"    software/simolant");

   makebutton(getmaxx()-47,58,'s',"-",SPEEDHELP);
   makebutton(getmaxx()-19,58,'S',"+",SPEEDHELP);

#if LANG==0
   makebutton(atx+114,6,F1,"help",HELP);

   setcolor(YELLOW); outtextxy(atx,aty+3,"T:");
   makebutton(atx+113,aty,'T',  "warm",THELP);
   makebutton(atx+58,aty, 't',  "cool",THELP);
   aty+=DY;

   setcolor(YELLOW); outtextxy(atx,aty+3,"method:");
   makebutton(atx+58,aty, 'm',  "MC/MD",METHODHELP);
   makebutton(atx+113,aty,'u',   "auto",METHODHELP); aty+=DY;

   setcolor(YELLOW);
   outtextxy(atx,aty+3,"gravity:");
   makebutton(atx+88,aty, 'g',"\31",GHELP);
   makebutton(atx+111,aty,'G',"\30",GHELP);
   makebutton(atx+134,aty,'z',"0",GHELP);
   aty+=DY;

   aty+=DYS;
   setcolor(YELLOW); outtextxy(atx,aty+3,"walls:");
   makebutton(atx+86,aty-DYS,'Y', "top", WALLHELP);
   makebutton(atx+57,aty,    'x', "left",  WALLHELP);
   makebutton(atx+106,aty,   'X', "right", WALLHELP);
   makebutton(atx+76,aty+DYS,'y', "bottom", WALLHELP);;
   aty+=DYS+DY;

   setcolor(YELLOW); outtextxy(atx,aty+3,"prepare:");
   makebutton(atx+70,aty,F5,"gas",GASHELP); aty+=DYS;
   makebutton(atx+70,aty,F6,"liquid",LIQUIDHELP); aty+=DYS;
   makebutton(atx+70,aty,F7,"capillary",CAPILHELP); aty+=DYS;

   makebutton(atx+50,aty,F9,"crystal",IDHELP); aty+=DYS;

   makebutton(atx+50,aty,F10,"dislocation",EDGEHELP); aty+=DYS;
   makebutton(atx+50,aty,F11,"vacancy",VACANCYHELP); aty+=DYS;
   makebutton(atx+50,aty,F12,"interstitial",INTERHELP); aty+=DYS;

   setcolor(YELLOW); outtextxy(atx,aty+3,"show:");
   makebutton(atx+48,aty,'n',"neighbors",NBRHELP); aty+=DY;
#else /*# LANG==0 */
   makebutton(atx+107,6,F1,"navod",HELP);

   setcolor(YELLOW); outtextxy(atx,aty+3,"T:");
   makebutton(atx+106,aty,'T',   "ohrej",THELP);
   makebutton(atx+36,aty, 't',  "ochlad",THELP);
   aty+=DY;

   setcolor(YELLOW); outtextxy(atx,aty+3,"metoda:");
   makebutton(atx+58,aty, 'm',  "MC/MD",METHODHELP);
   makebutton(atx+113,aty,'u',   "auto",METHODHELP); aty+=DY;

   setcolor(YELLOW);
   outtextxy(atx,aty+3,"gravitace:");
   makebutton(atx+88,aty, 'g',"\31",GHELP);
   makebutton(atx+111,aty,'G',"\30",GHELP);
   makebutton(atx+134,aty,'z',"0",GHELP);
   aty+=DY;

   aty+=DYS;
   setcolor(YELLOW); outtextxy(atx,aty+3,"steny:");
   makebutton(atx+80,aty-DYS,'Y', "horni", WALLHELP);
   makebutton(atx+57,aty,    'x', "leva",  WALLHELP);
   makebutton(atx+106,aty,   'X', "prava", WALLHELP);
   makebutton(atx+80,aty+DYS,'y', "dolni", WALLHELP);;
   aty+=DYS+DY;

   setcolor(YELLOW); outtextxy(atx,aty+3,"priprav:");
   makebutton(atx+70,aty,F5,"plyn",GASHELP); aty+=DYS;
   makebutton(atx+70,aty,F6,"kapalinu",LIQUIDHELP); aty+=DYS;
   makebutton(atx+70,aty,F7,"kapilaru",CAPILHELP); aty+=DYS;

   makebutton(atx+40,aty,F9,"krystal",IDHELP); aty+=DYS;

   makebutton(atx+40,aty,F10,"dislokaci",EDGEHELP); aty+=DYS;
   makebutton(atx+40,aty,F11,"vakanci",VACANCYHELP); aty+=DYS;
   makebutton(atx+40,aty,F12,"intersticial",INTERHELP); aty+=DYS;

   setcolor(YELLOW); outtextxy(atx,aty+3,"ukaz:");
   makebutton(atx+48,aty,'n',"sousedy",NBRHELP); aty+=DY;
#endif /*#!LANG==0 */
}

int nbrs[10];

/* count and draw neighbors (hotkey n) */
void neighbors(int show)
{
  int i,j,ix,iy;
  int nbr;
  char ch[2];
  vector rt;
  unsigned char col[11]={BLACK,
                         DARKGRAY,
                         RED,
                         MAGENTA,
                         LIGHTBLUE,
                         LIGHTCYAN,
                         WHITE, /* 6 */
                         LIGHTRED,
                         LIGHTMAGENTA,
                         LIGHTGREEN,
                         LIGHTGRAY};

  loopto (nbr,0,7) nbrs[nbr]=0;
  if (show) {
    settextstyle(0,0,1);
    setfillstyle(1,BLACK); };

  loop (i,0,n) {
    rt=r[i];
    nbr=0;
    switch (bc) {
    case BOX:
      loop (j,0,n) if (i!=j)
        if (Sqr(rt.x-r[j].x)+Sqr(rt.y-r[j].y)<range) nbr++;
      break;
    case SLIT:
      loop (j,0,n) if (i!=j)
        if (Sqrni(rt.x-r[j].x)+Sqr(rt.y-r[j].y)<range) nbr++;
      break;
    case CYCLIC:
      loop (j,0,n) if (i!=j)
        if (Sqrni(rt.x-r[j].x)+Sqrni(rt.y-r[j].y)<range) nbr++;
      break;
    }
    if (show) {
      ix=(int)(rt.x*sx)-4;
      iy=(int)(rt.y*sy)-4;
      ch[0]='0'+nbr;
      ch[1]=0;
      if (nbr>9) ch[0]='X',nbr=10;
      if (nbr) {
        setcolor(col[nbr]);
        bar(ix-1,iy-1,ix+9,iy+10);
        outtextxy(ix+1,iy,ch); } }
    else
      if (nbr<=7) nbrs[nbr]++;
  }

  if (show) {
    XFlush(display);
    sleep(3);
    if (kbhit()) if (readkey()==ESC) finish();
    redraw(); }
}

/* prints the panel with variable info and sets min and max T and some other things */
void infopanel(int all)
{
  tox=getmaxx()-50;
  strcpy(st,constt?"NVT":"NVE");
  strcat(st,md?" MD":" MC");
  //  strcat(st,string("      %+d",ndelay));
  info(st,0,YELLOW,-1);

  tox=getmaxx()-50;
  if (T<9.51337e-4) T=9.51337e-4;
  if (T>Tmax) T=Tmax;
  if (constt) bag=T;
  Lh=L/2;

  info(" T=",T,8+7*constt,1);
  if (md) {
    if (hauto) {
      i=LIGHTRED; h=hadj; }
    else i=WHITE;
    info("dt=",h,i,4);
  }
  else {
    if (dstat==Auto) {
      i=LIGHTRED; d=dadj; }
    else i=WHITE;
    info("dr=",d,i,4); };

  setcolor(YELLOW);
  setfillstyle(1,BLACK);
  bar(getmaxx()-40,40,getmaxx(),50);
  outtextxy(getmaxx()-40,40,string("r=%+d",ndelay));

  if (all) {
    if (md) info("Tk=",Tm,LIGHTCYAN,2);
    else info("Td=",Tm,LIGHTCYAN-3*constt,2);
    if (demo) {
      setcolor(YELLOW);
      setfillstyle(1,BLUE);
      settextstyle(0,0,1);
      bar(atx,217,getmaxx()-1,245);
      neighbors(0);
#define L1 219
#define L2 232
#if LANG==0
      if ((nbrs[6]<=1) && (T>0.25)) {
        hot=1;
        outtextxy(atx,219," \20 try to cool");
        outtextxy(atx,232,"   to about T=0.15"); }
      else if ((nbrs[5]+nbrs[4]>=n/2) && (T>0.1)) {
        if (hot) {
          outtextxy(atx,219," \20 try to cool");
          outtextxy(atx,232,"   to about T=0.05"); }
        else {
          outtextxy(atx,219," \20 try to warm");
          outtextxy(atx,232,"   to about T=1.0"); } }
      else if ((nbrs[6]>n/2) && (T<0.08)) {
        hot=0;
        outtextxy(atx,219," \20 try to warm");
        outtextxy(atx,232,"   to about T=0.15"); }
      else if (Tm>1.2*T) outtextxy(atx,227,"  refrigerating...");
      else if (Tm<0.8*T) outtextxy(atx,227,"  heating...");
#else /*# LANG==0 */
      if ((nbrs[6]<=1) && (T>0.25)) {
        hot=1;
        outtextxy(atx,219,"   zkus ochladit");
        outtextxy(atx,232,"   na asi T=0.15"); }
      else if ((nbrs[5]+nbrs[4]>=n/2) && (T>0.1)) {
        if (hot) {
          outtextxy(atx,219,"   zkus ochladit");
          outtextxy(atx,232,"   na asi T=0.05"); }
        else {
          outtextxy(atx,219,"   zkus ohrat");
          outtextxy(atx,232,"   na asi T=1.0"); } }
      else if ((nbrs[6]>n/2) && (T<0.08)) {
        hot=0;
        outtextxy(atx,219,"   zkus ohrat");
        outtextxy(atx,232,"   na asi T=0.15"); }
      else if (Tm>1.2*T) outtextxy(atx,227,"   chladim...");
      else if (Tm<0.8*T) outtextxy(atx,227,"   topim...");
#endif /*#!LANG==0 */
    }
    info("ar=",ar,LIGHTCYAN-3*md,6);
  }
}

void redraw(void) /************************************** redraw everything */
{
  clearviewport();
  selectfont(0);
  setfillstyle(1,BLACK);
  drawwalls();
  setcolor(WHITE);
  aty=getmaxy()-xfont.height*13-17;
  settextstyle(0,0,1);
  makebuttons();
  prtgravity(0);
  infopanel(0);
  drawcfg();
}

void setscr(void) /***************** calculates screen coordinates of atoms */
{
  loop (i,0,n) {
    scr[i].ix=(int)(r[i].x*sx);
    scr[i].iy=(int)(r[i].y*sy); }
}

/* draw all atoms, using scr[] */
void drawcfg(void)
{
  int i;

  if (ximage) {
    int j,x,y,x0,x9,y0,y9,rr=rad*(rad+1);

    if (col==BLACK) {
      loop (i,1,maxyn-2)
        loop (j,1,maxyn-2) XPutPixel(ximage,i,j,xcoltab[BLACK]); }
    else if (col==BLUE) {
      loop (i,1,maxyn-2)
        loop (j,1,maxyn-2)
          if (XGetPixel(ximage,i,j)==xcoltab[WHITE])
            XPutPixel(ximage,i,j,xcoltab[BLUE]); }
    else {
      loop (i,1,maxyn-2)
        loop (j,1,maxyn-2)
          if (rnd()<.02) XPutPixel(ximage,i,j,xcoltab[BLACK]); }

    loop (i,0,n) {
      x0=scr[i].ix-rad,x9=scr[i].ix+rad;
      y0=scr[i].iy-rad,y9=scr[i].iy+rad;
      Max(x0,1) Min(x9,maxyn-1)
        Max(y0,1) Min(y9,maxyn-1)
        loopto (x,x0,x9) loopto(y,y0,y9)
        if (Sqr(scr[i].ix-x)+Sqr(scr[i].iy-y)<=rr)
          XPutPixel(ximage,x,y,xcoltab[WHITE]); }
    setwritemode(0);
    XPutImage(display,win,gc,ximage,2,2,2,2,maxyn-4,maxyn-4);
    XFlush(display); }
  else {
    setcolor(WHITE);
    setlinestyle(0,0,th);
    loop (i,0,n) circle(scr[i].ix,scr[i].iy,rad); }
}

/* resize the simulation box + redraw */
void resize(myreal newL)
{
  int i;
  myreal q;

  if (newL>Lorig) newL=Lorig;
  if (newL<sqrt(n/3)) newL=sqrt(n/3);
  q=newL/L;
  L=newL;
  loop (i,0,n) { r[i].x*=q; r[i].y*=q; };
  setscr();
  redraw();
}

/* read ASCII dump of cfg */
void readcfg()
{
  int i;
  FILE *f=in;
  char line[128];

  in=fopen("simolant.dup","rt");
  getdata
    get(n) get(L) get(walldens)
    get(T) get(Tmax) get(bc)
    get(h) get(dadj)
  checkdata enddata
  fgets(line,128,in);
  loop (i,0,n) {
    double x,y;
    if (!fgets(line,128,in) || 2!=sscanf(line,"%lf%lf",&x,&y)) {
      fprintf(stderr,"not enough or bad data in simolant.dup\n");
      exit(0); }
    r[i].x=x; r[i].y=y; }
  fclose(in);
  in=f;
}

void help(void) /******************************* the help screen (hotkey h) */
{
  clearviewport();
  aty=4;
  settextstyle(0,0,1);
  setcolor(LIGHTGREEN);
  selectfont(2);
#if LANG==0
  hh("SYMBOLS");
  setcolor(WHITE);
  hh("NVE MC  Creutz Monte Carlo (microcanonical ensemble, NVE)");
  hh("NVT MC  Metropolis Monte Carlo (canonical ensemble, NVT)");
  hh("NVE MD  microcanonical molecular dynamics");
  hh("NVT MD  molecular dynamics with Berendsen thermostat");
  hh("T       temperature (parameter of the NVT ensemble)");
  hh("Tk      temperature (from kin.en.)  Td    temperature (from Creutz demon energy)");
  hh("dr      max MC displacement         dt    MD timestep");
  hh("P       virial pressure             ar    acceptance ratio of MC displacements");
  aty+=8;
  setcolor(LIGHTGREEN);
  hh("HOTKEYS (functions marked by - are not available by mouse)");
  setcolor(WHITE);
  hh("m       toggle between Monte Carlo and molecular dynamics");
  hh("u       toggle auto MC/MD selection");
  hh("e      -toggle between NVE and NVT ensembles");
  hh("t T     decrease/increase temperature (NVT only), also wheel") ;
  hh("* /    -decrease/increase step for the above");
  hh("g G z   decrease/increase/zero gravity (negative=down)");
  if (bc==BOX) hh("x X y Y toggle between repulsive (red) and attractive (WHITE) walls");
  else if (bc==SLIT) hh("y Y     toggle between repulsive (red) and attractive (WHITE) walls");
  hh("d D    -MC: change displacement, red=auto; MD: toggle fixed timestep/auto (red)");
  hh("c C    -cycle between different colors and repaint of atom trajectories");
  hh("r R    -decrease/increase atom radius on screen (pot.unchanged)");
  hh("a A    -toggle thickness of circles/lines; try in case of flickering");
  hh("l L o  -decrease/increase/original box size + rescale configuration accordingly");
  hh("p      -toggle virial pressure calculation");
  hh("n       count neighbors of atoms");
  hh("cursors-move the whole configuration");
  hh("0      -freeze (pause), continue by any key   SPACE  -redraw screen");
  hh("=      -set more parameters (from the console)");
  hh("F      -toggle trajectory dump to file (only if caled with parm, HUGE FILES!)");
  hh("ESC    -return (press twice)");
  aty+=8;
  setcolor(LIGHTMAGENTA);
  hh("...press any key");
#else /*# LANG==0 */
  hh("SYMBOLY");
  setcolor(WHITE);
  hh("NVE MC  Creutzovo Monte Carlo (microkanonicky soubor, NVE)");
  hh("NVT MC  Metropolisovo Monte Carlo (kanonicky soubor, NVT)");
  hh("NVE MD  microkanonicka molekulova dynamika");
  hh("NVT MD  molekulova dynamika s Berendsenovym termostatem");
  hh("T       teplota (parametr NVT souboru)");
  hh("Tk      teplota (z kineticke enegie)   Td   teplota z energie Creutzova demona");
  hh("dr      max delka MC pohybu            dt   casovy krok v MD");
  hh("P       tlak (z virialu sil)           ar   zlomek prijatych konfiguraci v MC");
  aty+=8;
  setcolor(LIGHTGREEN);
  hh("KLAVESY (funkce ozn. - nelze ovladat mysi)");
  setcolor(WHITE);
  hh("m       prepina mezi metodami Monte Carlo a molekulova dynamika");
  hh("u       zapne/vypne automaticky vyber metody MC/MD");
  hh("e      -prepina mezi soubory NVE a NVT");
  hh("t T     snizi/zvysi teplotu (jen NVT), tez kolecko mysi");
  hh("* /    -snizi/zvysi krok pro zmenu teploty");
  hh("g G z   snizi/zvysi/vynuluje gravitaci (zaporna=dolu)");
  if (bc==BOX) hh("x X y Y prepina mezi odpudivou (cervenou) a pritazlivou (zelenou) stenou");
  else if (bc==SLIT) hh("y Y     prepina mezi odpudivou (cervenou) a pritazlivou (zelenou) stenou");
  hh("d D    -MC: zmena delky pohybu, cervena=auto; MD: casovy krok konstantni delky/auto");
  hh("c C    -vyznaceni trajektorie (modre, bile, zadne)");
  hh("r R    -zmensuje/zvetsuje velikost atomu (potencial se nemeni)");
  hh("a A    -prepina tloustku car/krouzku; zkus, pokud zobrazeni blika");
  hh("l L o  -zmensi/zvetsi/puvodni simulacni krabice, konfigurace se preskaluje");
  hh("p      -zapne/vypne vypocet tlaku");
  hh("n       vyznaci pocet sousedu atomu");
  hh("kurzory-pohyb cele konfigurace");
  hh("0      -zastavi simulaci (pauza), pokracuj cimkoliv   SPACE  -prekresli obrazovku");
  hh("=      -nastav parametry (na terminalu)");
  hh("F      -zapne/vypne vypis trajektorie do souboru (obrovske soubory!)");
  hh("Esc    -konec (stiskni 2x)");
  aty+=8;
  setcolor(LIGHTMAGENTA);
  hh("...press any key");
#endif /*#!LANG==0 */
  if (readkey()==ESC) finish();
  redraw();
}

/* set some extra parameters (hotkey `=') */
void setparms()
{
#if LANG==0
  if (md) {
    h=GetReal("integration time step dt",h);
    Tmax=GetReal("max NVT temperature",Tmax);
    tau=GetReal("time constant for MD friction thermostat",tau);
  }
  else {
    acc=GetReal("acceptance ratio for automatic MC displacement",acc);
    dadj=GetReal("MC displacement (if not auto)",dadj);
  };
  range=GetReal("range for counting neighbors",sqrt(range));
#else /*# LANG==0 */
  if (md) {
    h=GetReal("integrační krok dt",h);
    Tmax=GetReal("max NVT teplota",Tmax);
    tau=GetReal("časová konstanta termostatu",tau);
  }
  else {
    acc=GetReal("zlomek přijetí pro automatick nastavení MC pohybu",acc);
    dadj=GetReal("MC pohyb (není-li auto)",dadj);
  };
  range=GetReal("dosah definující sousedy",sqrt(range));
#endif /*#!LANG==0 */
  range=range*range;
}

/* to move the configuration (using cursors) */
void movexy(myreal dx, myreal dy)
{
  loop (i,0,n) {
    r[i].x+=dx; if (r[i].x>L) r[i].x-=L; if (r[i].x<0) r[i].x+=L;
    r[i].y+=dy; if (r[i].y>L) r[i].y-=L; if (r[i].y<0) r[i].y+=L;  }
  setscr();
  redraw();
}

void initcfg(enum cfgtype cfg)
{
  int i,ii,j,k;

  msweeps=MSWEEPS;
  if (demo) mauto=1;

  n=nn;

  Lh=L/2;
  if (cfg==DUMP) return;

  hauto=1;
  dstat=Auto;
  constt=1;
  pressure=1;

  md=0;
  gravity=0;
  uwallx=uwallr; fwallx=fwallr;
  uwallxL=uwallr; fwallxL=fwallr;
  uwally=uwallr; fwally=fwallr;
  uwallyL=uwallr; fwallyL=fwallr;

  if (cfg==GAS) {
    /* initial random configuration */
    loop (i,0,n) {
      r[i].x=0.5+rnd()*(L-1); r[i].y=0.5+rnd()*(L-1);
    }
    T=1;
    return;
  }

  if (cfg==CAPILLARY) {
    /* initial random configuration */
    loop (i,0,n) {
      r[i].x=0.5+rnd()*(L-1); r[i].y=0.5+(rnd()*0.3+0.7)*(L-1);
    }
    T=0.15;
    gravity=-0.01;
    uwallx=uwalla; fwallx=fwalla;
    uwallxL=uwalla; fwallxL=fwalla;
    uwallyL=uwalla; fwallyL=fwalla;
    return;
  }

  T=0.05;
  if (cfg==LIQUID) T=0.15;

  k=1;
  while (3*k*(k+1)+1<n) k++;
  if (3*k*(k+1)+1>=MAXN) k--;
  n=3*k*(k+1)+1;

  switch (cfg) {
    case DEFECT: n-=k+1; break;
    case VACANCY: n--; break;
    case INTERSTITIAL: n++; break;
    default:;
  }

  ii=0;
  rt.x=1.13;
  if (cfg==LIQUID) rt.x=1.35;
  rt.y=0.866*rt.x;

  for (i=k; i>=0; i--) {
    for (j=-k; j<=k-i; j++) {
      if (cfg==VACANCY && i==0 && j==1) ii--;
      r[ii].x=Lh+(i*0.5+j)*rt.x;
      r[ii].y=Lh+i*rt.y;
      if (i!=0) {
        ii++;
        r[ii].x=r[ii-1].x;
        r[ii].y=L-r[ii-1].y;
      }
      ii++;
    } }

  if (cfg==INTERSTITIAL) {
    r[ii].x=Lh+0.5;
    r[ii].y=Lh+0.5; }

  if (cfg==LIQUID) ff=1; else ff=sqrt(T)*0.3;

  loop (i,0,n)  {
    r[i].x+=rnd()*ff; r[i].x+=rnd()*ff; }

  if (cfg==DEFECT) loop (i,0,n) if (r[i].x>Lh) {
    ff=sqrt((r[i].x-Lh)/k)*0.43;
    if (r[i].y>Lh) { r[i].x+=ff; r[i].y+=ff; }
    else { r[i].x-=ff; r[i].y+=ff; };
  }
}

void reinitcfg(enum cfgtype cfg)
{
  initcfg(cfg);
  setscr();
  redraw();
  drawcfg();
  infopanel(1);
}

void randomv(void)
{
  int i;
  loop (i,0,n) {
    v[i].x=vv*(rnd()-rnd()+rnd()-rnd()+rnd()-rnd());
    v[i].y=vv*(rnd()-rnd()+rnd()-rnd()+rnd()-rnd()); }
}

/* read hotkey switch */
void hotkey(void)
{
  usleep(internaldelay);

  while (kbhit()) {
    unsigned int ch=readkey();

#ifdef DEBUG
    fprintf(stderr,"hotkey: %c=%d\n",ch,ch);
#endif /*# DEBUG */
    switch (ch) {

      case X_REDRAW_ME:
        redraw();
        break;

      /* toggle MC/MD */
      case 'm':
      case 'M':
        mauto=0;
        vv=sqrt(T*2);
        md=!md;
        if (md)
          /* switching to molecular dynamics:
             reset velocities to Maxwell-Boltzmann */
          randomv();
        else {
          /* switching to Monte Carlo: reset daemon's bag */
          bag=T;
        }
      break;

      /* temperature */
      case WHEELFORWARD:
      case 'T':
        if (constt) T=T*qT;
        break;
      case WHEELBACKWARD:
      case 't':
        if (constt) T=T/qT;
        break;

      case '*': qT*=qT; break;
      case '/': qT=sqrt(qT); break;

      /* set and draw gravity */
      case 'g':
        prtgravity(-0.01);
        break;
      case 'G':
        prtgravity(0.01);
        break;
      case 'z':
        gravity=0;
        prtgravity(0);
        break;

      /* toggle repulsive/attractive walls */
      case 'x': {
        if (uwallx==uwalla) { uwallx=uwallr; fwallx=fwallr; }
      else { uwallx=uwalla; fwallx=fwalla; }
        drawwalls(); }
        break;
      case 'X': {
        if (uwallxL==uwalla) { uwallxL=uwallr; fwallxL=fwallr; }
        else { uwallxL=uwalla; fwallxL=fwalla; }
        drawwalls(); }
        break;
      case 'Y': {
        if (uwally==uwalla) { uwally=uwallr; fwally=fwallr; }
        else { uwally=uwalla; fwally=fwalla; }
        drawwalls(); }
        break;
      case 'y': {
        if (uwallyL==uwalla) { uwallyL=uwallr; fwallyL=fwallr; }
        else { uwallyL=uwalla; fwallyL=fwalla; }
        drawwalls(); }
        break;

      /* cycle between MC displacements/MD timesteps */
      case 'd':
        if (md) hauto=!hauto;
        else if (dstat<fixed) dstat++; else dstat=half;
        break;
      case 'D':
        if (md) hauto=!hauto;
        else if (dstat>half) dstat--; else dstat=fixed;
        break;

      /* toggle const T / const E */
      case 'e': case 'E':
        constt=!constt;
        break;

      /* trace color */
      case 'c':
      case 'C':
        if (col==LIGHTMAGENTA) col=BLACK;
        else if (col==BLACK) col=BLUE;
        else if (col==BLUE) col=LIGHTMAGENTA;
        break;

      /* slow/fast */
      case 's':
        ndelay--;
        if (ndelay<-6) ndelay=-6;
        break;
      case 'S':
        ndelay++;
        if (ndelay>9) ndelay=9;
        break;

      /* atom size */
      case 'r':
        rad=rad*4/5;
        if (rad<=0) rad=1;
        redraw();
        break;
      case 'R':
        rad=maxrad;
        redraw();
        break;

      case 'A': case 'a':
        th+=2;
        if (th>6) th=1;
        if (th==5) {
          ximage=ximage0;
          setlinestyle(0,0,3); }
        else {
          ximage=NULL;
          setlinestyle(0,0,th); }
        redraw();
        break;

      case 'u':
        mauto=!mauto;
        break;

      /* box resize */
      case 'l':
        resize(L/1.03);
        break;
      case 'L':
        resize(L*1.03);
        break;
      case 'o': case 'O':
        resize(Lorig); /* original one */
        break;

      case 'L'&31:
      case ' ':
        if (ximage) {
          int i,j;
          loop (i,1,maxyn-2)
            loop (j,1,maxyn-2) XPutPixel(ximage,i,j,xcoltab[BLACK]); }
        redraw();
        break;

      /* freeze */
      case '0':
        if (readkey()==ESC) finish();
        break;

      /* toggle: calculate and show virial pressure */
      case 'p': case 'P': {
        pressure=!pressure;
        if (!pressure) {
          info(" P=",P,DARKGRAY,5); }
        break;
      }

      case 'n': case 'N':
        neighbors(1);
        break;

      case 'h': case 'H':  case '?':
        help();
        break;

      case '=':
        setparms();
        break;

      case F1:
        outbwtextxy(2,2,2,HELP);
        readkey();
        redraw();
        break;

      case F5: reinitcfg(GAS);  break;
      case F6: reinitcfg(LIQUID); break;
      case F7: reinitcfg(CAPILLARY); break;
      case F9: reinitcfg(CRYSTAL); break;

      case F10: reinitcfg(DEFECT); break;
      case F11: reinitcfg(VACANCY); break;
      case F12: reinitcfg(INTERSTITIAL); break;

      /* move the cfg */
      case DOWN: movexy(0,0.1+(bc==CYCLIC));  break;
      case UP: movexy(0,-0.1-(bc==CYCLIC)); break;
      case LEFT: movexy(-0.1-(bc>BOX),0); break;
      case RIGHT: movexy(0.1+(bc>BOX),0); break;

      case 'F': /* if (!demo) */ {
        static int nplb=0;

        /* toggle trajectory dump */
        settextstyle(0,0,1);
        setcolor(WHITE);
        setfillstyle(1,GREEN);
        bar(atx-3,getmaxy()-33,getmaxx(),getmaxy());
        if (plb) {
#if LANG==0
          outtextxy(atx,getmaxy()-21," trajectory closed");
#else /*# LANG==0 */
          outtextxy(atx,getmaxy()-30,"   trajektorie");
          outtextxy(atx,getmaxy()-15,"    uzavrena");
#endif /*#!LANG==0 */
          fclose(plb);
          plb=NULL; }
        else {
          plb=fopen(string("sim%04d.traj",nplb),"wt");
          nplb++;
#if LANG==0
          outtextxy(atx,getmaxy()-30," trajectory opened");
          outtextxy(atx,getmaxy()-15," F again to close");
#else /*# LANG==0 */
          outtextxy(atx,getmaxy()-30,"zaznam trajektorie");
          outtextxy(atx,getmaxy()-15,"stiskni F pro konec");
#endif /*#!LANG==0 */
          fprintf(plb,"#  x    y    vx    vy  N=%d L=%f\n",n,L); } }
        break;

#ifdef DEBUG
      case 'i': INTERNALDELAY/=2; break;
      case 'I': INTERNALDELAY=INTERNALDELAY*3/2+2; break;
#endif /*# DEBUG */

      case ESC: finish(); break;
      default: ;
    }
    infopanel(0);
  }
}

void MDerror(void) /********************************************************/
{
  setfillstyle(1,RED);
  selectfont(3);
  setcolor(YELLOW);
  bar(0,0,660,45);
#if LANG==0
  outtextxy(10,12,"MD integration problem, switching to NVT MC");
#else /*# LANG==0 */
  outtextxy(10,12,"Problem MD integrace, prepinam na NVT MC");
#endif /*#!LANG==0 */
  XFlush(display);
  sleep(2);
  if (mauto) msweeps=MSWEEPS;
  md=0; constt=1;
  redraw();
}

int main(int narg,char **arg)
{
  XArc *lastxcirc,*xcirc;

#if CHECKHEAP==2
  AllocRange=64000;
#endif /*# CHECKHEAP==2 */

  selectfont(0);

  bag=0;
  sweeps=1;
  pressure=1;

  initscroll(0);
  rndinit(0,0);
  demo=narg==1;

  _n
#if LANG==0
  prt("Monte Carlo and molecular dynamics simulations of");
  prt("2D Lennard-Jones-like particles in a square simulation cell");
  prt("(c) J. Kolafa 1993-2012");
  _n
  prt("Atom-atom potential: u(r)=1/r^8-1/r^4");
  prt("Attractive atom-wall potential: uwall(r)=pi*D*(3/8/r^7-1/2/r^3)");
  prt("Repulsive atom-wall potential: uwall(r)=pi*D*3/8/r^7");
  _n
  if (demo) prt("Need more options? Use a paramater:\n  simolant x");
#else /*# LANG==0 */
  prt("Simulace dvourozměrných částic v čtvercové simulační buňce");
  prt("metodami Monte Carlo a molekulová dynamika");
  prt("(c) J. Kolafa 1993-2012");
  _n
  prt("Potenciál atom-atom: u(r)=1/r^8-1/r^4");
  prt("Přitažlivost atom-stěna: uwall(r)=pi*D*(3/8/r^7-1/2/r^3)");
  prt("Odpudivost atom-stěna: uwall(r)=pi*D*3/8/r^7");
  _n
  if (demo) prt("Potřebuješ víc možností? Volej s parametrem:\n  simolant x");
#endif /*#!LANG==0 */
  _n

  mauto=1;

  if (!demo) {
#if LANG==0
    if (yes("special initial conditions")) {
      prt("select initial conditions:");
      prt("  0=gas (the default: short timestep, high Tmax)");
      prt("  1=liquid drop (longer timestep, lower Tmax)");
      prt("crystals (even longer timestep and lower Tmax):");
      prt("  2=ideal crystal  3=with edge defect  4=with vacancy  5=with interstitial");
      prt("  6=read from file simolant.dup");
      cfg=(enum cfgtype)(GetInteger("",0)); }
#else /*# LANG==0 */
    if (yes("speciální počáteční podmínky")) {
      prt("zvol počáteční podmínky:");
      prt("  0=plyn (default: krátký časový krok, vysoké Tmax)");
      prt("  1=kapka kapaliny (delší časový krok, nižší Tmax)");
      prt("krystaly: (ještě delší časový krok a nižší Tmax):");
      prt("  2=ideální krystal  3=s hranovou dislokací  4=s vakancí  5=s intersticiálem");
      prt("  6=čti soubor simolant.dup");
      cfg=(enum cfgtype)(GetInteger("",0)); }
#endif /*#!LANG==0 */
  if (cfg!=DUMP) {
#if LANG==0
    prt("select boundary conditions:");
    bc=(enum bctype)(GetInteger("  0=hard walls  1=slit pore (periodic in x)  2=periodic",0));
#else /*# LANG==0 */
    prt("zvol okrajové podmínky:");
    bc=(enum bctype)(GetInteger("  0=pevné zdi  1=štěrbina (periodická v ose x)  2=periodické",0));
#endif /*#!LANG==0 */
    }
  }

  /* some defaults */
  L=1;
  dstat=one; dadj=1; d=1;

  h=0.02;

  if (cfg!=GAS) {
    L=0.9;
    dstat=Auto; dadj=0.5; d=0.5;
    h=0.05; hadj=0.05; Tmax=0.3; n=140;
    pressure=0;
  }

  if (cfg>LIQUID) {
    L=0.8;
    dadj=0.1; hadj=0.05; d=0.1;
    h=0.1; Tmax=0.16; n=200;
  }

  if (cfg!=DUMP) {
    if (demo)
      L*=(sqrt(n*6)+1-(bc)/2.);
    else {
#if LANG==0
      n=GetInteger(string("number of atoms (max %d)",MAXN),n);
#else /*# LANG==0 */
      n=GetInteger(string("počet atomů (max %d)",MAXN),n);
#endif /*#!LANG==0 */

      if (n>MAXN || n<1) exit(0);
#if LANG==0
      L=GetReal("initial box size L",L*(sqrt(n*6)+1-(bc)/2.));
      if (bc!=CYCLIC) walldens=GetReal("density of atoms on walls D",1.5);
#else /*# LANG==0 */
      L=GetReal("počáteční velikost boxu L",L*(sqrt(n*6)+1-(bc)/2.));
      if (bc!=CYCLIC) walldens=GetReal("hustota atomů stěny D",1.5);
#endif /*#!LANG==0 */
    };
    walldens*=PI;
  }
  else
    readcfg();

  nn=n;

  hauto=1;
  dstat=Auto;
  constt=1;

  Lorig=L;

  pressure=1;

  uwallx=uwallr; fwallx=fwallr;
  uwally=uwallr; fwally=fwallr;
  uwallxL=uwallr; fwallxL=fwallr;
  uwallyL=uwallr; fwallyL=fwallr;
  gravity=0;

  col=BLACK;

  /*T=GetReal("initial temperature T",1);*/

  initcfg(cfg);

  xwindowhints.winname=xwindowhints.iconame="simolant";
  xwindowhints.geometry="840x680+80+20";
  //  xwindowhints.geometry="640x490+80+20";

  startgraph(-9);
  maxxn=getmaxx()+1;
  maxyn=getmaxy()+1;

  /* turning this off selects direct draw */
  ximage=XGetImage(display,win,0,0,maxxn,maxyn,AllPlanes,ZPixmap);
  ximage0=ximage;

  atx=getmaxy()+6;

  aspx=aspy=1;
  sy=getmaxy()/Lorig;
  sx=sy*aspy/aspx;
  maxrad=(int)(sx/2.0+0.01);
  //maxrad=(int)(sx/2.0+.5);
  rad=maxrad;

  setscr();

  md=0; constt=1;
#ifdef DEBUG
  ttime0=mytime();
#endif /*# DEBUG */
  redraw();

  for (;;) { /* main loop */

    /* (?) otherwise X11 is slow ... */
    if (sweeps>0) internaldelay=INTERNALDELAY*sqrt(n)/sweeps;
    if (ndelay<0) internaldelay=10;

    if (col>15) col=0;
    if (col<0) col=15;

    Tm=0;
    accepted=0;

    if (ndelay>0) sweeps=pow(1.618034,2+ndelay)/2.236068+.5;
    else sweeps=1;

    if (!ximage) {
      lastxcirc=startcircles(n);
      loop (i,0,n) circles(lastxcirc,scr[i].ix,scr[i].iy,rad);

      xcirc=startcircles(n); }

    loop (ii,0,sweeps) { /* <<<<<< */

      if (md) { /* ========== molecular dynamics ========== */

        if (hauto) h=hadj;

        /* calculate forces (accelerations) caused by walls and gravity */
        loop (i,0,n)
          if (bc==CYCLIC) {
            a[i].x=a[i].y=0; }
          else {
            ri=r[i];
            a[i].y=fwally(ri.y)-fwallyL(L-ri.y)-gravity;
            if (bc==BOX) a[i].x=fwallx(ri.x)-fwallxL(L-ri.x);
            else a[i].x=0; }

        /* calculate atom-atom pair forces */
        switch (bc) {
          case BOX:
            loop (i,0,n) {
              ri=r[i];
              loop (j,i+1,n) {
                rt.x=ri.x-r[j].x; rt.y=ri.y-r[j].y;
                ff=f(Sqr(rt.x)+Sqr(rt.y));
                a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
                a[j].x-=rt.x*ff; a[j].y-=rt.y*ff; }
            }
            break;
          case SLIT:
            loop (i,0,n) {
              ri=r[i];
              loop (j,i+1,n) {
                rt.x=ri.x-r[j].x;
                if (rt.x>Lh) rt.x=rt.x-L;
                if (rt.x<-Lh) rt.x=rt.x+L;
                rt.y=ri.y-r[j].y;
                ff=f(Sqr(rt.x)+Sqr(rt.y));
                a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
                a[j].x-=rt.x*ff; a[j].y-=rt.y*ff;  }
            }
            break;
          case CYCLIC:
            loop (i,0,n) {
              ri=r[i];
              loop (j,i+1,n) {
                rt.x=ri.x-r[j].x;
                if (rt.x>Lh) rt.x=rt.x-L;
                if (rt.x<-Lh) rt.x=rt.x+L;
                rt.y=ri.y-r[j].y;
                if (rt.y>Lh) rt.y=rt.y-L;
                if (rt.y<-Lh) rt.y=rt.y+L;
                ff=f(Sqr(rt.x)+Sqr(rt.y));
                a[i].x+=rt.x*ff; a[i].y+=rt.y*ff;
                a[j].x-=rt.x*ff; a[j].y-=rt.y*ff;  }
            }
        }

        /* leap-frog integrator, calculate kinetic energy, and draw atoms */
        vv=0;
        maxvv=0;

        loop (i,0,n) {
          /* crash integrator test */
          if ((fabs(a[i].x)>1000) || (fabs(a[i].y)>1000)) {
            put(h)
            MDerror();
            goto escape; };

          /* leap-frog */
          v[i].x+=h*a[i].x; v[i].y+=h*a[i].y;
          xx=Sqr(v[i].x)+Sqr(v[i].y);
          vv+=xx;  /* Ekin is calculated at time t+h/2 */
          if (xx>maxvv) maxvv=xx;

          r[i].x+=h*v[i].x;
          r[i].y+=h*v[i].y;

          /* prepare draw */
          if (ii==sweeps-1) {
            scr[i].ix=(int)(r[i].x*sx); scr[i].iy=(int)(r[i].y*sy);
            if (xcirc) circles(xcirc,scr[i].ix,scr[i].iy,rad);
          }
        } /*n*/

        /* adjust the timestep */
        maxvv=sqrt(maxvv);
        xx=maxvv*h;
        if (xx>0.15) {
          h=0.06/sqrt(maxvv);
          hadj=h;
          MDerror(); };
        if (xx>0.12) hadj=0.1/maxvv;
        else hadj=hadj*0.95+(0.09*0.05)/maxvv;
        if (hadj>0.1) hadj=0.1;

        vv/=2; /* Ekin */
        Tk=vv/n; /* kinetic temperature */

        if (constt) {
          /* friction thermostat: rescale velocities, correlation time ~ tau */
          ff=exp((T-Tk)*h/tau);
          loop (i,0,n) { v[i].x*=ff; v[i].y*=ff; } }

        Tm+=vv;

        t+=h;
      } /* md */

      else { /* ========== Monte Carlo ========== */

        /* max displacement = Lh=L/2 */
        if (dadj>Lh) dadj=Lh;

        /* determining the current displacement */
        switch (dstat) {
          case half: d=0.5; break;
          case one:  d=1; break;
          case Lhalf: d=L/2;  break;
          default: d=dadj; }

        loop (i,0,n) {
          /* trial move of atom i */

          ri=r[i]; /* old position */

          /* rnd() displacement in a 1-circle */
          do {
            incirc.x=rndcos();
            incirc.y=rndcos(); }
          while (Sqr(incirc.x)+Sqr(incirc.y)>1);

          /* new (trial) position (not too much efficient but simple) */
          rt.x=ri.x+incirc.x*d;
          /* normalize to the box,
             even for box+slit b.c. which prevents out-of-box problems */
          if (rt.x<0) rt.x+=L;
          else if (rt.x>L) rt.x-=L;

          /* ... and the same for y coord */
          rt.y=ri.y+incirc.y*d;
          if (rt.y<0) rt.y=rt.y+L;
          else if (rt.y>L) rt.y=rt.y-L;

          /* energy difference caused by the trial move: walls and gravity */
          deltaU=0;
          if (bc<=SLIT) {
            deltaU=uwally(rt.y)-uwally(ri.y)
              +uwallyL(L-rt.y)-uwallyL(L-ri.y)
              -gravity*(ri.y-rt.y);
            if (bc==BOX)
              deltaU=deltaU+
                uwallx(rt.x)-uwallx(ri.x)
                +uwallxL(L-rt.x)-uwallxL(L-ri.x);
          }

          /* : all other atoms ... not efficient, old energies calculated
             again and again though they might have been stored */
          switch (bc) {
            case BOX:
              loop (j,0,n) if (i!=j)
                deltaU+=u(Sqr(rt.x-r[j].x)+Sqr(rt.y-r[j].y))-u(Sqr(ri.x-r[j].x)+Sqr(ri.y-r[j].y));
              break;
            case SLIT:
              loop (j,0,n) if (i!=j)
                deltaU+=u(Sqrni(rt.x-r[j].x)+Sqr(rt.y-r[j].y))-u(Sqrni(ri.x-r[j].x)+Sqr(ri.y-r[j].y));
              break;
            case CYCLIC:
              loop (j,0,n) if (i!=j)
                deltaU+=u(Sqrni(rt.x-r[j].x)+Sqrni(rt.y-r[j].y))-u(Sqrni(ri.x-r[j].x)+Sqrni(ri.y-r[j].y));
          }

          if ( (constt && -deltaU/T>log(rnd()+1e-99)) /* Metropolis test */
               || (!constt && deltaU<bag) /* Creutz microcanonical test */ ) {
            /* move accepted */
            if (!constt) bag-=deltaU;
            accepted++;
            r[i]=rt;

            /* automatic adjustment of displacement size to reach
               acceptance ratio = acc; tiny violation of
               microreversibility is usually acceptable */
            if (dstat==Auto) dadj*=1+(1-acc)*1e-3; }
          else {
            /* move rejected */
            /* automatic adjustment of displacement size */
            if (dstat==Auto) dadj*=1-acc*1e-3; }

          if (ii==sweeps-1) {
            scr[i].ix=(int)(r[i].x*sx); scr[i].iy=(int)(r[i].y*sy);
            if (xcirc) circles(xcirc,scr[i].ix,scr[i].iy,rad);  }

          Tm+=bag;
        } /*n*/
      } /*MC*/
    } /* sweeps */

    if (!ximage) {
      setcolor(col);
      showcircles(lastxcirc);
      free(lastxcirc);

      setcolor(WHITE);
      showcircles(xcirc);
      free(xcirc);
      XFlush(display); }
    else
      drawcfg();

    if (ndelay<0) delay(30/pow(1.6,ndelay));

    if (plb) {
      loop (i,0,n)
        if (md) fprintf(plb,"%g %g %g %g\n",r[i].x,r[i].y,v[i].x,v[i].y);
        else fprintf(plb,"%g %g\n",r[i].x,r[i].y);
      fprintf(plb,"\n"); }

    /* fool-proof out-of-box check */
  escape:
    loop (i,0,n) {
      if ((r[i].x<=0) || (r[i].x>=L)) {
        if ((!md) || (bc==BOX)) MDerror();
        while (r[i].x<0) r[i].x+=L;
        while (r[i].x>L) r[i].x-=L; }
      if ((r[i].y<=0) || (r[i].y>=L)) {
        if ((!md) || (bc!=CYCLIC)) MDerror();
        while (r[i].y<0) r[i].y+=L;
        while (r[i].y>L) r[i].y-=L; } }

    hotkey();

    if (msweeps>0) msweeps--;

    if (mauto && !md && msweeps==0) {
      md=1;
      vv=sqrt(T*2);
      randomv(); };

    Tm/=(n*sweeps);
    ar=(double)accepted/(n*sweeps);

    /* calculate virial pressure and RDF */
    if (pressure) {
      P=0; /* P=sum r*f = -virial */
      loopto (i,0,histmax) hist[i]=0;
      loop (i,0,n) {
        ri=r[i];
        if (bc<CYCLIC) {
          P+=ri.y*fwally(ri.y)+(L-ri.y)*fwallyL(L-ri.y);
          if (bc==BOX) P+=ri.x*fwallx(ri.x)+(L-ri.x)*fwallxL(L-ri.x); }
        switch (bc) {
          case BOX:
            loop (j,i+1,n)
              addP(Sqr(ri.x-r[j].x)+Sqr(ri.y-r[j].y));
            break;
          case SLIT:
            loop (j,i+1,n)
              addP(Sqrni(ri.x-r[j].x)+Sqr(ri.y-r[j].y));
            break;
          case CYCLIC:
            loop (j,i+1,n)
              addP(Sqrni(ri.x-r[j].x)+Sqrni(ri.y-r[j].y));
        } }

      if (md) TT=Tk; else TT=T;
      P=(n*TT+P/3)/Sqr(L);
      info(" P=",P,LIGHTCYAN,5);
      vv=0;
      loopto (i,0,histmax) {
        rr=(double)hist[i]/(2*i+1);
        if (rr>vv) vv=rr; }
      vv=80.0/vv;

      setfillstyle(1,GREEN);
      bar(atx,214-80,atx+3*histmax+1,214-1);
      setcolor(WHITE);
      line(atx,214,atx,214-80);
      line(atx+35,214,atx+35,214-80);
      line(atx+70,214,atx+70,214-80);
      line(atx+105,214,atx+105,214-80);
      line(atx+140,214,atx+140,214-80);
      setcolor(RED);
      setlinestyle(0,0,3);
      loopto (i,0,histmax) {
        j=214-(int)(vv*hist[i]/(2*i+1)+0.5);
        if (j!=214) line(atx+i*3,j,atx+i*3,214); }
      setlinestyle(0,0,th); }

    infopanel(1);

#ifdef DEBUG
    {
      static double av=0;
      ttime=mytime();
      av=av*.99+.01*((ttime-ttime0)/sweeps);
      fprintf(stderr,"%7.5f %6.4f %d\n",av,(ttime-ttime0)/sweeps,INTERNALDELAY);
      ttime0=ttime;
    }
#endif /*# DEBUG */

  } /* main loop */
}
