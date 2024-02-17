/* make show
 */

#define VERSION "2.2j"

/* Notes on raytracing:

  `show' does not do ray tracing but generates data for two raytracers.

  1/ `ray', the "reasonably intelligent raytracer" by Mark
     VandeWettering, with some improvements.  In the ball or
     bar-and-ball mode, use hot key `N', followed by `o' (one) or `m'
     (merge) or `s' (series).  To render the picture, use
       ray -n FILE#### [MORE OPTIONS]
     where FILE####.nff is the file generated
     use
       ray -h
     to get help
     Example:
       ray -n test0000 -S.5 -vxv
     will calculate test0000.ppm from test0000.nff of half the
     original size (-S.5) with fast but reasonable quality
     antialiasing and show it using xv

     Advantages:
       better `ambient light' (`inside' of molecules)
       easy change of picture aspect ratio
       tailored to this type of scenes
     Disadvantages:
       no textures and special effects
       sometimes cumbersome control
       radiosity calculation not supported

  2/ PovRay, "Persistence Of Vision" raytracer.   In the ball or
     bar-and-ball mode, use hot key `V', followed by `o' (one) or `m'
     (merge) or `s' (series).  To render the picture, use
       ray +IFILE####.pov [MORE OPTIONS]
     where FILE####.pov is the file generated
     File "show.inc" defining colors and textures is needed

     Advantages:
       much richer set of options and special effects, textures,...
       widely used - a sort of standard
     Disadvantages (may be fixed in newer versions):
       (?: perhaps somebody tells me how to overcome these problems...)
       uses nonstandard left-handed coordinate system
       difficult to change picture aspect ratio
       cannot render large molecules (core dump!)
       bad ambient light (`inside' of molecules poorly rendered)
     ?: `radiosity' is supported, but according to my experience it
       does not work and just makes ugly random patches -- perhaps
       will be fixed in future versions...

     Example:
       povray res320.ini +A +Itest0001.pov +D
     will calculate test0001.<what-is-your-default> from test0001.pov
     in size 320x200 with fast but reasonable quality
     antialiasing (+A) displaying picture in progress (+D)

******************************************************************************/


// #define BETTERBALLS /* slightly more precise but slower ball algorithm */
                    /* not needed for large ZBUFDEPTH (X11) */

#define CHARGES     /* recoloring by charges possible: takes space! */
#define SHELL       /* allows option -j and hotkeys 'j','J': show waters close
                       to protein BUG: 'j','J' takes effect after next frame
                       is reloaded */
#define DEPTH       /* depth cue */
#define SLICES      /* enable z-slices */

#define ASPECT 1
#define CONVERT "-loop 0 -delay 5 -limit memory 4GiB -limit map 4GiB -limit disk 8GiB"

#ifndef XIMAGE
#  define XIMAGE 2
/* X11 mode:
   0=draw image (in ball/stick modes) directly pixel-by-pixel (slow, shaky)
   1=put image (faster, avoids flickering), ximage copied from window
     (problems with full screen, moving of window screen out of screen)
     INCOMPATIBLE WITH MENU
   2=put image (faster, avoids flickering), ximage created from scratch
     (not fully compatible, fails on newer CygWins)
*/
#endif /*# XIMAGE */

#include "ground.h"

int SCRLEN,oldSCRLEN=-1;
#include "xdraw.h"

#if XIMAGE
XImage *ximage;
extern GC gc;
extern long unsigned *xcoltab;
#endif /*# XIMAGE */

#include "options.h"
#define x badoption
int optionlist[32] =
/*` a b c d e f  g h i j k l  m n o p  q r   s t  u v w x y z  { | } ~ ? */
{-1,4,0,0,0,1,-2,7,0,0,0,2,0,-1,0,0,0,50,x,100,0,-1,1,0,0,0,0,-5,1,1,0,0};
/*@ A B C D E F  G H I J K L  M N O P  Q R   S T  U V W X Y Z  [ \ ] ^ _*/
/* NB: -r is used by option -r#g#b */
#undef x

#define VEC(R) R[0],R[1],R[2]
/* WARNING: povray uses left-handed coordinate system */
#define iVEC(R) R[0],R[1],-R[2]

/* ZBUFDEPTH<=max(unsigned SITEINT), normally equal */
#define ZBUFDEPTH 65535
#define SITEINT short

typedef struct site_s {
#ifdef CHARGES
  float charge;
#endif /*# CHARGES */
  unsigned SITEINT c;       /* internal (palette) atom color (for spheres) */
  unsigned SITEINT r;       /* radius in 0.01 AA; temporary color info */
  unsigned SITEINT mark;    /* 0=none, 1=H, 2=O (for HB), 4=click-marked */
  unsigned SITEINT cluster; /* for function "find clusters" */
  SITEINT molcol;           /* to calculate col,rad if no other info */
#ifdef SHELL
  SITEINT shell;            /* 1=water, 2=water Oxygen, 4=to draw */
#  define SHELLTEST(I) if (site[I].shell&4)
#else /*# SHELL */
#  define SHELLTEST(I) /*empty*/
#endif /*#!SHELL */
  char Crad[4];             /* new 4 bytes instead of .gol */
  char *id;
  char *type;
} site_t;

char *shellkey="HOH";
char *parameter_set="DUMMY";

typedef float fvector[3];
typedef double real;
#include "vector3d.h" /* needed by dihangle.c */
typedef int ivector[2];

int bodyframe; /* rotations WRT fixed frame/body frame */
int presenter; /* F4 = left presenter = 'i' */

double drot[3][3]={
  {1,0,0},
  {0,1,0},
  {0,0,1} };

fvector rot[3]={
  {1,0,0},
  {0,1,0},
  {0,0,1} };

fvector *cfg; /* [ns] */

/* modes BOND and BOND2:
  xy are screen coordinates, xy0 are the old ones
  (because the old bonds are erased by xor-redrawing)
  similarly bond and bond0
*/

ivector *xy,*xy0; /* [ns], screen positions of atoms (to be swapped) */
struct forclick_s {
  float z; /* z-position of the center (real units) */
  int rr; /* radii for ball modes (screen units) */
} *forclick; /* [ns] */

int nxygold,nxygold0;
ivector *xygold,*xygold0;
char boxcolor='O'; /* color of box, one of WYRGBCMO */
int nwalls; /* 0,1,2 */
int wallatomr=-300; /* vdW radius in 0.01 AA (site[].r-compatible) */
float wallz[2]={-9e9,-9e9};

site_t *site; /* [ns] */
static fvector center;
static float typicalsize;

int ns; /* number of sites */
int nc,nc0; /* number of bonds (constraints), incl. possible hydrogen bonds */
            /* nc0 refers to bond0 = erase by xor redrawing */
int ncfix;  /* number of bonds (constraints), as read from .mol */

char *donor="H*",*acceptor="O*";
float HBDIST=2.4,hbdist=0;
/* NOTE(?): a hot key to change hbdist is OK, however,
   changing hbdist=0 to hbdist>0 will break the bond allocation algorithm */
void xhbdist(void) {
  float x=hbdist;
  hbdist=HBDIST;
  HBDIST=x;
  if (hbdist && HBDIST) hbdist=0;
  fprintf(stderr,"H-bonds turned %s %g %g\n",hbdist?"on":"off",hbdist,HBDIST);
}
int ihbcolor=3; /* GREEN: index of turbocolor */
int hbcolor;
void hbonds(void);
void readGOL(void);
unsigned key;

char
  *molname, /* MOLNAME (without extension) */
  *plbname, /* playback file name or format accepting int (without extension) */
  *plbfn;   /* playback file name or format accepting int (with extension) */

static char *Fn(char *ext) /******************************************* Fn */
{
  static char fn[256];

  strcpy(fn,molname);
  strcat(fn,".");
  strcat(fn,ext);

  return fn;
}

unsigned char *scrbuf;

unsigned del=0;
int font=2;
static struct timeb { int time,millitm; } t0;
void ftime(struct timeb *x){}

char *li;
typedef int bond_t[3]; /* [2]&0xf0 = bond color (unless BOND2), permanent
                             &0x0f = bond color to show (removed if over box) */
bond_t *bond,*bond0; /* nbond,nbond0 : this and last (cf. xy0) */
                     /* warning: bond may reallocate */
int nbond,nbond0;
int bondeqbond0; /* 1 if bond=bond0 (identical array) */
int maxx,maxxn,maxy,maxyn;

int enlarge(bond_t **bond,int from,int to) /*********************** enlarge */
/* change size of *bond
 the new size is returned
 sizes are in sizeof(bond_t)
 NB: realloc may not be compatible with module alloc.h/alloc.c
*/
{
  bond_t *b;
  int n=min(from,to);

  allocarray(b,to);
  copyarray(b,*bond,n);
  free(*bond);
  *bond=b;

//  prt("bond reallocated: %d->%d",from,to);

  return to;
}

#define NOCOL 8 /* number of ball colors */
int MODE=31*NOCOL+2; /* # of palette entries, max. 256
                        background=0 (used as black, for wire modes, as set by GUI=..#HEXCOLOR)
                        orange=MODE-2
                        background=MODE-1 (for bar and ball modes, as set by -bg)
                        for byte-deep displays minus some for other applications */
int COLLEN; /* length of palette for one ball colors */
int wscale100=100;

/*
  available atom colors
  WARNING: different color code than turbocolor
  don't change the order or names of colors! - see myline2 !!!
*/
struct colortab_s {
  char *col;      unsigned char rgb[4]; char key,turbocolor; char *NFFcol;} colortab[NOCOL+1] = {
/*      col       {r   g b n.a.  }   -@   for .gol NFF:diff refl spec refr ior */
/*0*/ {"WHITE"   ,{255,255,255,0},  'W', WHITE,       "1 0 20 0 0"},
/*1*/ {"YELLOW"  ,{255,255,0,  0},  'Y', YELLOW,      "1 0 21 0 0"},
/*2*/ {"RED"     ,{255,0,  0,  0},  'R', LIGHTRED,    "1 0 22 0 0"},
/*3*/ {"GREEN"   ,{0,  255,0,  0},  'G', LIGHTGREEN,  "1 0 23 0 0"},
/*4*/ {"BLUE"    ,{0,  0,  255,0},  'B', LIGHTBLUE,   "1 0 24 0 0"},
/*5*/ {"CYAN"    ,{0,  223,223,0},  'C', LIGHTCYAN,   "1 0 25 0 0"},
/*6*/ {"MAGENTA" ,{255,0,  255,0},  'M', LIGHTMAGENTA,"1 0 26 0 0"},
/*7*/ {"ORANGE"  ,{255,159,0,  0},  'O', BROWN,       "1 0 26 0 0"},
/*8*/ {""        ,{0,0,0,0},        'X',0,NULL} /* omit atom */
};

/* color->gray conversion for postscript hardcopy */
unsigned char *psgray;
float col2gray[3]={0.3,0.59,0.11};

/* Turbo color -> palette color (offset in the palette table) */
/* AVAILABLE COLORS
   0..8 not available
    9=BLUE
   10=GREEN
   11=CYAN
   12=RED
   13=MAGENTA
   14=YELLOW
   15=WHITE;
   */
unsigned char turbo2pal[16];

struct collist_s {
  struct collist_s *next;
  int colno;
  char *pattern;
} *collisthead=NULL;

static char aminoacids[]="GLY ALA VAL LEU ILE PRO SER THR ASP GLU ASN GLN LYS ARG HIS PHE TYR TRP MET CYS";

static unsigned char (*pal)[3]; /* the palette: RGB */

char *ballbackground; /* hex number, as passed by -bg */

void setmypalette(void) /*************************************** setmypalette */
{
  int c,i,p,ci;
  unsigned char pa;
  float psg;
  int COLOFF;

  if (pal) return; /* palette already set */

  COLOFF=COLLEN/2;

  allocarrayzero(psgray,MODE);
  allocarrayzero(pal,MODE);

  /* pal[0]=background, set by xwindowhints.background */
  ci=1;
  loop (c,0,NOCOL) {
    turbo2pal[(int)colortab[c].turbocolor]=1+COLLEN*c;
    loop (i,0,COLLEN) {
      psg=0;
      loop (p,0,3) {
        pa=colortab[c].rgb[p]*(i+COLOFF)/(COLOFF+COLLEN-1);
        pal[ci][p]=pa;
        psg+=pa*col2gray[p]; }
      psg=255.0f-(255.0f-psg)*wscale100/100; /* bit more white */
      if (psg>255.0f) psg=255.0f;
      psgray[ci]=psg;
      ci++; } }

  if (xwindowhints.background) {
    /* background for wire modes, as passed from GUI, should be dark */
    unsigned x[3];
    char *c=xwindowhints.background;

    if (*c=='#') c++;
    sscanf(c,"%2x%2x%2x",x,x+1,x+2);
    psg=0;
    loop (p,0,3) {
      pa=pal[0][p]=x[p];
      psg+=pa*col2gray[p]; }
    psg=255.0f-(255.0f-psg)*wscale100/100; /* bit more white */
    if (psg>255.0f) psg=255.0f;
    psgray[0]=psg; }

  /* background for ball and bond modes, as defined by -bg */
  if (ballbackground) {
    unsigned x[3];
    char *c=ballbackground;

    if (*c=='#') c++;
    sscanf(c,"%2x%2x%2x",x,x+1,x+2);
    psg=0;
    loop (p,0,3) {
      pa=pal[MODE-1][p]=x[p];
      psg+=pa*col2gray[p]; }
    psg=255.0f-(255.0f-psg)*wscale100/100; /* bit more white */
    if (psg>255.0f) psg=255.0f;
    psgray[MODE-1]=psg; }
  else {
    pal[MODE-1][0]=0; pal[MODE-1][1]=0; pal[MODE-1][2]=1;
    /* true black 0,0,0 will become background */
    psgray[MODE-1]=255; }

  //  loop (i,0,MODE) printf("%d  %d %d %d\n",i,pal[i][0],pal[i][1],pal[i][2]);
  loop (i,0,MODE) setpal(i,pal[i][0],pal[i][1],pal[i][2],256);
}

unsigned SITEINT *zbuf;

fvector L,Lh; /* box size, if option('l')) */
int varL;

fvector shift;
char shiftinL[3];
int isshift;

FILE *dump,*file,*plb;
enum { NONE,PLB,PPM,ZBUF,NFF,ATM,POV,PS,COLORPS,EPS,BONDPS } dumpmode;
char *dumpinfo[]={"none","PLB","PPM","ZBUF","NFF","ATM","POV","PS","color PS","EPS","bond PS"};
/* NFF must be the 1st TEXT output file */
enum { ONE,SERIES,START,CONT } dumpstat;
char *dumpext[]={"???","plb","ppm","zbuf","nff","atm","pov","ps","ps","eps","ps"};
char *dumpname=NULL;
int dumpindex;
float hdr[2];
float NFFfloor;
#ifdef SLICES
int slicecolor;
unsigned SITEINT slicez;
#endif /*# SLICES */
fvector eye={0,0,-3e33}; /* eye position in AA; -3e33=auto setup
                            NB: convention fixed/changed in V2.2e */
fvector Scale; /* cfg -> screen scaling factors; now ASPECT=1 so SCALE[1]=-SCALE[0] */
int ids[3]; /* offset, pixels */
float zmin,zmax,zbufmin,zbufmax;
#ifdef DEPTH
float zmind,zmaxd,cue,icue;
#endif /*# DEPTH */
int dumpprepared;
int iplb,lastiplb,percent;
int menu=1,redrawmenu=1;
int central=1; /* central projection */

int fastmove;
int step=0,laststep=1,swing=0;
int pos;        /* frame number: 1st frame=1 (warning - this has changed!) */
int maxpos;     /* = number of frames */
int lastpos=-1; /* last frame */
double sliderpos,trace=0;
char percentfn[256];

void newplb(int advance) /******************************************* newplb */
{
  FILE *f;
  char *c,*d=NULL;

  iplb+=advance;
  sprintf(percentfn,plbfn,iplb);
  f=fopen(percentfn,"rb");
  if (!f) {
    if (advance) { iplb-=advance; return; }
    else ERROR(("open %s",percentfn)); }
  fprintf(stderr,"showing %s\n",percentfn);

  for (c=percentfn; *c; c++) if (*c=='.') d=c;
  if (d) *d=0;
  if (plb) fclose(plb);
  plb=f;
}

static float *nmf;

void readnmf(int ns) /********************************************** readnmf */
{
  if (!nmf && getenv("SHOWNMF")) {
    FILE *nmff;

    allocarray(nmf,3*ns);
    if ( (nmff=fopen(getenv("SHOWNMF"),"rt")) ) {
      /*
#normal mode frequencies
#matrix n=3420   d2U/drdr=-df/dr, dr=2*1e-05   Jacobi diag.err.=9.9303e-12
#i        [THz]          [1/cm]
   0      0.00003634      0.00121221
      */
#define NMFLINE 48
      char line[NMFLINE];
      int i=0;
      float x;

      while (fgets(line,NMFLINE,nmff)) if (line[0]!='#') {
        sscanf(line,"%d%f%f",&i,&x,&x);
        if (i>=0 && i<ns*3) nmf[i]=x; }
      fclose(nmff); } }
}

#include "showdump.c"

#ifdef SHELL
#  include "inshell.c"
#endif /*# SHELL */

#ifdef CHARGES
#  include "groupc.c"
int negC,posC;
#endif /*# CHARGES */

int bigpoint,first=1;

#define ORANGE 16 // as "turbocolor"
#define BG 17
int palcolor(int turbo) /****************************************** palcolor */
/*
  convert TURBO color (BLACK=0,..,WHITE=15) into the (closest) palette color,
  to be used in setcolor() for lines (bonds) and text
*/
{
  /* dark colors but BLUE should not occur */
  switch (turbo) {
    case BLACK: return 0;
    case BG: return MODE-1;
    case WHITE: return COLLEN;
    case LIGHTGRAY: return COLLEN*2/3;
    case DARKGRAY: return COLLEN/3;
    case YELLOW: return COLLEN*2;
      //    case BROWN: return COLLEN*3/2;
    case BROWN: return COLLEN*15/2;
    case RED: return COLLEN*5/2;
    case LIGHTRED: return COLLEN*3;
    case GREEN: return COLLEN*7/2;
    case LIGHTGREEN: return COLLEN*4;
    case BLUE: return COLLEN*9/2;
    case LIGHTBLUE: return COLLEN*5;
    case CYAN: return COLLEN*11/2;
    case LIGHTCYAN: return COLLEN*6;
    case MAGENTA: return COLLEN*13/2;
    case LIGHTMAGENTA: return COLLEN*7;
    default: /* ORANGE */ return COLLEN*8; }
}

void myline(ivector *xy,int i,int j,int col) /********************** myline */
/*
  draw bond i-j with color col (in the palette; DOS: turbo color)
  col=0: do not draw
  xy is a precalculated array of 2D screen coordinates
  bigpoint>0: i=j (point)
*/
{
  if (!col) return;
#ifdef SHELL
  if (xy[i][0]==-9999 || xy[j][1]==-9999) return;
#endif /*# SHELL */

  setcolor(col);

  if (bigpoint==0)
    line(xy[i][0],xy[i][1],xy[j][0],xy[j][1]);
  else if (bigpoint==1)
    putpixel(xy[i][0],xy[i][1]+bigpoint,col);
  else {
    line(xy[i][0]-bigpoint+1,xy[i][1]-bigpoint+1,xy[i][0]+bigpoint,xy[i][1]+bigpoint);
    line(xy[i][0]-bigpoint+1,xy[i][1]+bigpoint,xy[i][0]+bigpoint,xy[i][1]-bigpoint+1); }
}

void myline2(ivector *xy,int i,int j,int bcol) /******************* myline2 */
/*
  bcol=0: draw bond i-j divided into 2 parts using atom colors
  bcol>0: draw bond i-j with color col
  bigpoint>0: i=j (point)
  if colors of i,j are the same, myline is called
  (note that if a bond composed of 2 parts is redrawn by 1 color, it must be done in 2 parts, too)
  NB: bcol is the palette color (DOS: TURBO color)
*/
{
  int col,xx,yy;

#ifdef SHELL
  if (xy[i][0]==-9999 || xy[j][1]==-9999) return;
#endif /*# SHELL */

  if (!bcol)
    col=site[i].c+(COLLEN-1);
  else
    col=bcol;

  if (site[i].c==site[j].c)
    myline(xy,i,j,col); /* this is called also if bigpoint set */
  else {
    setcolor(col);

    xx=(xy[i][0]+xy[j][0])/2;
    yy=(xy[i][1]+xy[j][1])/2;

    line(xy[i][0],xy[i][1],xx,yy);
    if (!bcol)
      col=site[j].c+(COLLEN-1);
    setcolor(col);
    line(xx,yy,xy[j][0],xy[j][1]); }
}

/* INFRONT = what to do if there is something to show in front of the screen */
#define INFRONT return; /* do not show and return */

void showsphere(struct forclick_s *fci,
                fvector r,float rad,int c) /********************* showsphere */
/* returns info in *fci (if not NULL):
   rr for marking in no-wire modes, BUG: aspect=1 assumed
   z of atom */
{
  fvector cr,R;
  unsigned map;
  int x0,x1,y0,y1;
  int sx01=ids[0]+1;
  int i,j;
  float xx,yy,rr,d,dr,dx,cscale,epsrr,q;
  int h,col,dotcol;

  if (fci) fci->rr=0,fci->z=0;

  if (c>=NOCOL*COLLEN || r[0]>99999.0) return;

  if (rad==0) return;
  VVV(cr,=r,-center)
  R[0]=rot[0][0]*cr[0]+rot[0][1]*cr[1]+rot[0][2]*cr[2];
  R[1]=rot[1][0]*cr[0]+rot[1][1]*cr[1]+rot[1][2]*cr[2];
  R[2]=rot[2][0]*cr[0]+rot[2][1]*cr[1]+rot[2][2]*cr[2];

  if (dumpmode==PLB) {
    opendump();
    fwrite(R,sizeof(fvector),1,dump); }

  if (dumpmode==NFF) {
    dumpnffcol(c);
    calcNFFfloor(R[1]-rad);
    fprintf(dump,"s %.4f %.4f %.4f %.4f\n",VEC(R),rad); }

  if (dumpmode==ATM) {
    int ic=(float (*)[3])r-cfg;
    opendump();
    if (ic>=0 && ic<NS) /* BRRRRRRR */
      fprintf(dump,"%c %lf %lf %lf\n",site[ic].type[0],cfg[ic][0],cfg[ic][1],cfg[ic][2]); }

  if (dumpmode==POV) {
    calcNFFfloor(R[1]-rad);
    fprintf(dump,"sphere {<%f,%f,%f>, %.4f texture {%s%s}}\n",
            iVEC(R),rad,colortab[5].rgb[1]<4?"ALT":"SHOW",colortab[c/COLLEN].col); }

  R[0]-=eye[0];
  R[1]-=eye[1];

  if (central) {
    q=eye[2]-R[2];
    if (q<=0) INFRONT
    q=eye[2]/q;
    rad*=q;
    R[0]*=q; R[1]*=q; }

  rr=rad*rad;
  epsrr=rr*0.999999;

  {
    Min(zmin,R[2])
    Max(zmax,R[2]+rad)
    Min(zbufmin,R[2]-rad)
    Max(zbufmax,R[2]+rad)
#ifdef DEPTH
    Min(zmind,R[2])
    Max(zmaxd,R[2]+rad)
#endif /*# DEPTH */
  }

  cscale=
#ifdef DEPTH
    (icue+cue*(R[2]-zmind)/(zmaxd-zmind))*
#endif /*# DEPTH */
                                         (COLLEN*2-1)/rr;
  dotcol=cscale*rr/2;
#ifdef SLICES
  if (slicecolor)
    slicecolor=(int)(
#  ifdef DEPTH
               (icue+cue*((slicez-ids[2])/Scale[2]-zmind)/(zmaxd-zmind))*
#  endif /*# DEPTH */
                                                          ((COLLEN*2-1)/2) )+c;
#endif /*# SLICES */

  if (option('f')<0 && rad*Scale[1]>-1.5 /* note that Scale[1]<0 */ ) {
    /* dot instead too small a sphere */
    x0=(int)(R[0]*Scale[0])+ids[0]; if (R[0]*Scale[0]<0) x0--;
    y0=(int)(R[1]*Scale[1])+ids[1]; if (R[1]*Scale[1]<0) y0--;
    if (y0>=0 && y0<maxy && x0>=0 && x0<maxx) {
      map=y0*maxxn+x0;
      h=(int)(R[2]*Scale[2])+ids[2];
#ifdef SLICES
      if (!slicez || slicez>h)
#endif /*# SLICES */
        if (h>zbuf[map]) {
          zbuf[map]=h;
          scrbuf[map]=c+dotcol; } } }
  else {

    y0=(int)((R[1]+rad)*Scale[1])+ids[1]; if (y0<0) y0=0;
    y1=(int)((R[1]-rad)*Scale[1])+ids[1]; if (y1>maxy) y1=maxy;

    loopto (j,y0,y1) {
      yy=(j-ids[1])/Scale[1]-R[1]; yy=yy*yy;
      d=epsrr-yy;
      if (d>0) dx=sqrt(d); else dx=0;
      x0=(int)((R[0]-dx)*Scale[0])+sx01; if (x0<0) x0=0;
      x1=(int)((R[0]+dx)*Scale[0])+ids[0]; if (x1>maxx) x1=maxx;
      map=j*maxxn+x0;
#ifdef SLICES
      if (slicez) {
        int hh;

        loopto (i,x0,x1) {
          xx=(i-ids[0])/Scale[0]-R[0]; xx=xx*xx;
          dr=d-xx; if (dr<0) dr=0;
          h=(int)((sqrt(dr)+R[2])*Scale[2])+ids[2];
          if (slicez>h) {
            if (h>zbuf[map]) {
              col=(int)(cscale*dr);
              zbuf[map]=h;
              scrbuf[map] = c+col/2;
              if (col&1) if ((i+j)&1) scrbuf[map]++; } }
          else {
            hh=(int)((-sqrt(dr)+R[2])*Scale[2])+ids[2];
            if (slicez>hh) {
              zbuf[map]=slicez;
              scrbuf[map]=slicecolor?slicecolor:MODE-2; } }
          map++; } }
      else
#endif /*# SLICES */
        loopto (i,x0,x1) {
          xx=(i-ids[0])/Scale[0]-R[0]; xx=xx*xx;
          dr=d-xx; if (dr<0) dr=0;
          h=(int)((sqrt(dr)+R[2])*Scale[2])+ids[2];
          if (h>zbuf[map]
#ifdef BETTERBALLS
              /* more precise but more time consuming */
              || (h==zbuf[map] && ((sqrt(dr)+R[2])*Scale[2])+ids[2])>zbuf[map]
#endif /*# BETTERBALLS */
              ) {
            col=(int)(cscale*dr);
            zbuf[map]=h;
            scrbuf[map] = c+col/2;
            if (col&1) if ((i+j)&1) scrbuf[map]++; }
          map++; }
      }
  }

  if (fci) fci->rr=Sqr(rad*Scale[1]+0.5),fci->z=R[2];
}

void voidsphere(struct forclick_s *fci,
                fvector r,float rad,int c) /********************* voidsphere */
/* the same as showsphere without showing (just eturn fci) */
{
  fvector cr,R;
  float q;

  if (fci) fci->rr=0,fci->z=0;
  else return;

  if (c>=NOCOL*COLLEN || r[0]>99999.0) return;

  if (rad==0) return;
  VVV(cr,=r,-center)
  R[2]=rot[2][0]*cr[0]+rot[2][1]*cr[1]+rot[2][2]*cr[2];

  if (central) {
    q=eye[2]-R[2];
    if (q<=0) INFRONT
    q=eye[2]/q;
    rad*=q; }

  if (fci) fci->rr=Sqr(rad*Scale[1]+0.5),fci->z=R[2];
}

float grad=0.15; /* thickness of bars (modes 1-4) in AA */
float gradscale=-1.5; /* >0: absolute thickness of box and walls
                         <0: relative to grad */

void showbar(fvector r0,fvector r1,int c) /************************* showbar */
{
  fvector cr[2],R[2];
  unsigned map;
  int x0,x1,y0,y1;
  int k,i,j,col;
  float xx,yy,rr,DD,q,lp,lm,lR0,lR1,ff,tt,cscale;
  int h;
  float D[3],d[2];
  float dd,ss,rad=grad,scal;

  if (r0[0]>99999.0 || r1[0]>99999.0 || !c) return;

  VVV(cr[0],=r0,-center)
  VVV(cr[1],=r1,-center)

  loop (k,0,2) {
    R[k][0]=rot[0][0]*cr[k][0]+rot[0][1]*cr[k][1]+rot[0][2]*cr[k][2];
    R[k][1]=rot[1][0]*cr[k][0]+rot[1][1]*cr[k][1]+rot[1][2]*cr[k][2];
    R[k][2]=rot[2][0]*cr[k][0]+rot[2][1]*cr[k][1]+rot[2][2]*cr[k][2]; }

  if (dumpmode==NFF) {
    dumpnffcol(c);
    fputc('c',dump);
    loop (k,0,2) {
      fprintf(dump," %.4f %.4f %.4f %.4f",VEC(R[k]),rad);
      calcNFFfloor(R[k][1]-rad); }
    fputc('\n',dump); }

  if (dumpmode==POV) {
    fprintf(dump,"cylinder {");
    loop (k,0,2) {
      fprintf(dump,"<%f,%f,%f>,",iVEC(R[k]));
      calcNFFfloor(R[k][1]-rad); }
    fprintf(dump,"%f open texture{%s%s}}\n",
            rad,colortab[5].rgb[1]<4?"ALT":"SHOW",colortab[c/COLLEN].col); }

  loop (k,0,2) {
    R[k][0]-=eye[0];
    R[k][1]-=eye[1]; }

  loop (k,0,2) {
    if (central) {
      q=eye[2]-R[k][2];
      if (q<=0) INFRONT
      q=eye[2]/q;
      rad*=sqrt(q);
      R[k][0]*=q; R[k][1]*=q; }

    Min(zmin,R[k][2])
    Max(zmax,R[k][2])
    Min(zbufmin,R[k][2]-rad)
    Max(zbufmax,R[k][2]+rad)
#ifdef DEPTH
    Min(zmind,R[k][2])
    Max(zmaxd,R[k][2]+rad)
#endif /*# DEPTH */
  }

  VVV(D,=R[1],-R[0])

  loop (i,0,2) if (D[i]==0) D[i]=rad*0.01;

  DD=Sqr(D[0])+Sqr(D[1]);

  if (D[1]<0) loop (i,0,3) {
    D[i]=-D[i];
    xx=R[0][i],R[0][i]=R[1][i],R[1][i]=xx; }

  d[0]=D[1]; d[1]=-D[0]; /* rot -90 deg = clockwise */

  xx=rad/sqrt(DD);
  d[0]*=xx; d[1]*=xx;
  rr=rad*rad;
  cscale=
#ifdef DEPTH
    (icue+cue*((R[0][2]+R[1][2])/2-zmind)/(zmaxd-zmind))*
#endif /*# DEPTH */
                                                         (2*COLLEN-1)/rr;

  y1=(R[0][1]-fabs(d[1]))*Scale[1]+ids[1]; y0=(R[1][1]+fabs(d[1]))*Scale[1]+ids[1];

  if (y0<0) y0=0;
  if (y1>maxy) y1=maxy;

  loopto (j,y0,y1) {
    yy=(j-ids[1])/Scale[1]-R[0][1];
    lp=R[0][0]+d[0]+(yy-d[1])*(D[0]/D[1]);
    lm=R[0][0]-d[0]+(yy+d[1])*(D[0]/D[1]);
    lR0=R[0][0]+     yy*(d[0]/d[1]);
    lR1=R[1][0]+    (yy-D[1])*(d[0]/d[1]);
    if (D[0]>0) ff=fmax(lm,lR0),tt=fmin(lp,lR1);
    else ff=fmax(lm,lR1),tt=fmin(lp,lR0);

    x0=ff*Scale[0]+ids[0]; if (x0<0) x0=0;
    x1=tt*Scale[0]+ids[0]; if (x1>maxx) x1=maxx;
    loopto (i,x0,x1) {
      xx=(i-ids[0])/Scale[0]-R[0][0];
      scal=D[0]*xx+D[1]*yy;
      dd=rr-(xx*xx+yy*yy-Sqr(scal)/DD);
      if (dd>0) {
        ss=scal/DD;
        /* WARNING: without term +sqrt(dd), the bond is shaded + flat
           term +sqrt(dd) is incorrect, but it looks a bit better */
        h=(int)((R[0][2]*(1-ss)+R[1][2]*ss +sqrt(dd) )*Scale[2])+ids[2];
        map=j*maxxn+i;
        if (h>zbuf[map]) {
          col=(int)(cscale*dd);
          zbuf[map]=h;
          scrbuf[map] = c+col/2;
          if (col&1) if ((i+j)&1) scrbuf[map]++; } }
    }
  }
}

void showbar2(fvector r0,fvector r1,int c0,int c1) /*************** showbar2 */
{
  fvector r;
  int i;

  if (c0==c1)
    showbar(r0,r1,c1);
  else {
    loop (i,0,3) r[i]=(r0[i]+r1[i])/2;
    showbar(r0,r,c0);
    showbar(r,r1,c1); }
}

static enum mode_e { UNDEF, BAR, BAR2, DUMBELL, DUMBELL2, BALL, BOND, BOND2 } curmode=UNDEF,lastmode;
static double mouserot;
static int clickmode; /* 0: one best match marked, 1,2,3: all marked (increasing radius) */

static char *mygetline(void) /************************************ mygetline */
/***
  gets one line `li' form `file',
  skipping comments (lines beginning with !)
***/
{
  char *c;

  li[0]=0;
  do {
    c=fgets(li,256,file);
    if (!c) return c;
  } while (li[0]=='!');

  return c;
}

void setmenusize(void);
void makemenu(int which);

void setgmode(enum mode_e mode) /********************************** setgmode */
{
  static int gstarted=0;

  selectfont(font); /* global */
  setmenusize();

  if (curmode!=UNDEF)
    if ((mode<=BALL && curmode<=BALL) || (mode>=BOND && curmode>=BOND)) goto ret;

  // if (mode!=BALL) NEW: always needed, no optimization...
  {
    /* arrays needed for showing bonds */
    static int pass;

    if (!pass) {
      if (ns<=0) ERROR(("ns=%d illegal",ns))
      allocarray(xy,ns);
      allocarray(xy0,ns);
      allocarray(forclick,ns);
      pass=1; } }

  if (mode<=BALL) {
    /* 3D modes (XImage) */

    /* here is the 1st call to startgraph: MODE is max color */
    if (!gstarted) {
      startgraph(MODE);
      if (getmaxy()<479) font=0; else font=2;
      selectfont(font);
      setmenusize(); }

    redrawmenu=1;

    SCRLEN=maxxn*maxyn;
    if (SCRLEN!=oldSCRLEN) {
      //oldSCRLEN not set!
      if (zbuf) {
        free(zbuf);
        free(scrbuf); }
#if XIMAGE
      if (ximage) XDestroyImage(ximage);
#endif /*# XIMAGE */
      oldSCRLEN=SCRLEN;

      allocarray(zbuf,SCRLEN);
      alloc(scrbuf,SCRLEN);
#if XIMAGE
#  if XIMAGE==1
      /*
        read existing (just having created) window
        - probably correct code
        - slow over slow network
        - crashes if the window is shifted out of screen
        - fullscreen problems
      */
      ximage=XGetImage(display,win,0,0,getmaxx()+1,maxyn,AllPlanes,ZPixmap);
#  else /*# XIMAGE==1 */
      {
        /* create ximage from scratch - probably system-dependent
           cf. jkv.c and testx.c
        */
        unsigned *data;
        Visual *visual=NULL;
        data=malloc(SCRLEN*sizeof(data[0])); /* freed from XDestroyImage by system free */
        if (!data) ERROR(("no heap for ximage->data"))
        ximage=XCreateImage(display,visual,24,ZPixmap,0,(char*)data,maxxn,maxyn,32,0);
      }
#  endif /*#!XIMAGE==1 */
#endif /*# XIMAGE */
    }
    setmypalette();
  } /* mode<=BALL */

  if (mode>=BOND) {
    if (!gstarted) startgraph(-9);
    setmenusize();
    setcolor(palcolor(WHITE));
    setlinestyle(0,0,1);
    setfillstyle(0,0); }

  mouserot=2*PI/(maxxn+maxy);
  gstarted=1;

 ret: curmode=mode;
}

void project(ivector xyi,fvector r) /****************************** project */
{
  fvector cr;
  float q,sxq,syq;

  VVV(cr,=r,-center)

  if (central) {
    q=eye[2]-(rot[2][0]*cr[0]+rot[2][1]*cr[1]+rot[2][2]*cr[2]);
    if (q<=0) INFRONT
    q=eye[2]/q;
    sxq=Scale[0]*q; syq=Scale[1]*q; }
  else {
    sxq=Scale[0]; syq=Scale[1]; }

  xyi[0]=(int)(sxq*(rot[0][0]*cr[0]+rot[0][1]*cr[1]+rot[0][2]*cr[2]-eye[0]))+ids[0];
  xyi[1]=(int)(syq*(rot[1][0]*cr[0]+rot[1][1]*cr[1]+rot[1][2]*cr[2]-eye[1]))+ids[1];
}

void rotate(int axis,double angle) /******************************** rotate */
{
  double sa=sin(angle), ca=cos(angle);
  double o[3][3], rot2[3][3];
  int i=(axis+1)%3,j=(i+1)%3,k=(j+1)%3,a,b;

  if (bodyframe) sa=-sa;

  o[i][i]=ca; o[i][j]=-sa; o[i][k]=0;
  o[j][i]=sa; o[j][j]= ca; o[j][k]=0;
  o[k][i]= 0; o[k][j]=  0; o[k][k]=1;

  loop (a,0,3) loop (b,0,3)
    rot2[a][b]=o[a][0]*drot[0][b]+o[a][1]*drot[1][b]+o[a][2]*drot[2][b];

  loop (a,0,3) loop (b,0,3) {
    drot[a][b]=rot2[a][b];
    if (bodyframe) rot[b][a]=drot[a][b];
    else rot[a][b]=rot2[a][b]; }
}

void endian(char *a) /********************************************** endian */
{
  char x;

  x=a[0]; a[0]=a[3]; a[3]=x;
  x=a[1]; a[1]=a[2]; a[2]=x;
}

int readcfg(void) /************************************************* readcfg */
/* returns 1 if new cfg read, position=pos (1st frame=1) */
{
  int i,j,aux;

  if (pos==lastpos && iplb==lastiplb) return 0;

  fseek(plb,0,SEEK_END);
  maxpos=(ftell(plb)-option('k')*sizeof(float))/(sizeof(fvector)*(ns+varL));
  if (pos>maxpos) pos=maxpos;
  if (pos<1) pos=1;
  lastpos=pos;
  lastiplb=iplb;

  if (fseek(plb,(long)sizeof(fvector)*(ns+varL)*(pos-1)+option('k')*sizeof(float),SEEK_SET))
    ERROR(("read plb"))
  if (varL) {
    if (fread(L,sizeof(fvector),1,plb)!=1)
      ERROR(("read plb (frame=%d, # of frames=%d): cannot read L",pos,maxpos))
    Lh[0]=L[0]/2; Lh[1]=L[1]/2; Lh[2]=L[2]/2; }
  //  if (L[0]+L[1]+L[2]==0 && defaultslit) slit=incube=0;

  i=fread(cfg,sizeof(fvector),ns,plb);
  if (i!=ns)
    ERROR(("read plb (frame=%d, # of frames=%d): only %d sites read",pos,maxpos,i))

  for (j=2,aux=option('c'); j>=0; j--,aux/=10)
    if (aux%10) {
      float center=0;

      loop (i,0,ns) center+=cfg[i][j];
      center/=ns;
      loop (i,0,ns) cfg[i][j]-=center; }

  if (isshift) {
    fvector Lshift;

    loop (j,0,3)
      if (shiftinL[j]) Lshift[j]=L[j]*shift[j];
      else Lshift[j]=shift[j];

    loop (i,0,ns) loop (j,0,3) cfg[i][j]+=Lshift[j]; }

#ifdef SHELL
  inshell();
  if (option('j')<=0)
#endif /*# SHELL */

  if (option('l')) { /* periodic b.c. */
    loop (i,0,ns) loop (j,0,3) if (L[j]) {
      while (cfg[i][j]>L[j]) cfg[i][j]-=L[j];
      while (cfg[i][j]<0)    cfg[i][j]+=L[j]; }

    loop (i,0,nc) { /* bonds over periodic box removed (=replaced by 0) */
      aux=bond[i][2]&0xf0;
      bond[i][2]=aux|(aux>>4);
      loop (j,0,3) if (L[j]) if (fabs(cfg[bond[i][0]][j]-cfg[bond[i][1]][j])>Lh[j])
        bond[i][2]&=0xf0; } }

  return 1;
}

#define OPTION_C option('c')
#include "match.c"

#define fullangle 0 /* 0:dihedral potential in (-180,180); 1: (0,360) */
#include "dihangle.c"
#undef fullangle

#include "showmove.c"
#include "showmol.c"
#include "xdrawhelp.c"
#include "showhelp.c"
#include "showhlp.c"

#include "showmenu.c"

int main(int narg,char **arg) /**************************************** main */
{
  int i,j,i0,j0;
  static int bondcolor,number=1;
  double an;
  int b,ncmax;
  fvector maxr,minr;
  void *swap;
  float rad=0.01,SXY=0,cylsph=3;
  int redraw=1,reread=1,bookmark=0;
  struct collist_s *cl;
  enum mode_e newmode;
  char nr[32],*ext;
  int boxpalcolor,boxturbocolor;

#if CHECKHEAP==-2
  AllocRange=1024;
  AllocTrace=0;
#endif /*# CHECKHEAP==-2 */

  initscroll(0);
  repeatprefix('`');
  menu=getenv("GUI")==NULL || !!strchr(getenv("GUI"),'s');

  /* colors */
  if (getenv("SHOWPAL")) MODE=atoi(getenv("SHOWPAL"));
  if (MODE>255) MODE=255;
  COLLEN=(MODE-2)/NOCOL;
  if (COLLEN<2) COLLEN=2;
  MODE=NOCOL*COLLEN+2;
  markcolor=COLLEN*(NOCOL-1)+1; // start of oraNge
  boxpalcolor=markcolor;
  boxturbocolor=ORANGE;
  dosmouse.button[0]=palcolor(BLACK);
  dosmouse.button[1]=palcolor(DARKGRAY);
  dosmouse.button[2]=palcolor(LIGHTGRAY);
  dosmouse.button[3]=palcolor(WHITE);
  dosmouse.button[4]=palcolor(LIGHTCYAN);

  {
    char *name;
    /* extended ... */
    if ((name=getenv("SHOWGEOMETRY"))) xwindowhints.geometry=name;
    if ((name=getenv("SHOWINIT"))) addkbdstring((unsigned char*)name);
  }

  /*** option analysis ***/
  loop (i,1,narg)
    if (arg[i][0]=='-') {

      if (arg[i][1]=='Z') {
        char *comma=strchr(arg[i],',');

        wallz[0]=atof(arg[i]+2);
        if (comma++) {
          wallz[1]=atof(comma);
          comma=strchr(comma,',');
          if (comma++) {
            wallatomr=atof(comma)*100.001; } } }
      else if ( (arg[i][1]&31) == ('|'&31) ) {
        char *c=arg[i]+2;
        if (*c) {
          if (isalpha(*c)) boxcolor=*c++;
          gradscale=atof(c); }
        else
          option('|')=0; }

      else if (arg[i][1]=='I')
        addkbdstring((const unsigned char*)(arg[i]+2));
      else if (arg[i][1]=='S') {
        SXY=1/atof(arg[i]+2);
        typicalsize=10; }
      else if (arg[i][1]=='P') {
        char *comma=strchr(arg[i],',');
        if (comma) {
          eye[0]=atof(arg[i]+2);
          if (comma) {
            eye[1]=atof(comma+1);
            comma=strchr(comma+1,',');
            if (comma) {
              eye[2]=atof(comma+1); } } }
        else presenter++; }
      else if (arg[i][1]=='K')
        colortab[5].rgb[0]=colortab[5].rgb[1]=colortab[5].rgb[2]=128;
      else if (arg[i][1]=='J')
        shellkey=arg[i]+2;
      else if (arg[i][1]=='D')
        donor=arg[i]+2;
      else if (arg[i][1]=='A')
        acceptor=arg[i]+2;
      else if (arg[i][1]=='H') {
        hbdist=atof(arg[i]+2);
        if (hbdist==0) hbdist=HBDIST;
        else HBDIST=hbdist;
        loop (j,0,NOCOL) if (strchr(arg[i]+2,colortab[j].key)) ihbcolor=j;
        /* see def. of bond_t why *0x11 */
        hbcolor=colortab[ihbcolor].turbocolor*0x11; }
      else if (arg[i][1]=='L')
        animated.options=arg[i]+2;
      else if (arg[i][1]<='Z') {
        loop (j,0,NOCOL+1)
          if (arg[i][1]==colortab[j].key) {
            alloc(cl,sizeof(struct collist_s));
            cl->colno=j;
            cl->pattern=arg[i]+2;
            cl->next=collisthead; collisthead=cl;
            goto OK; }
        ERROR(("%s is unknown color: use one of {W,Y,R,G,B,C,M,O,X} (X=omit)",arg[i])) }
      else if (arg[i][1]=='r') {
        /* option -r%g%b% */
        char c[64];

        /*.....prt("%s",arg[i]); delay(2000);*/
        sprintf(c,"%s%s",arg[i]+2,"g59b11");
        col2gray[0]=(float)atoi(c)*0.01f;
        col2gray[1]=(float)atoi(strchr(c,'g'))*0.01f;
        col2gray[2]=(float)atoi(strchr(c,'b'))*0.01f; }
      else if (arg[i][1]=='b' && arg[i][2]=='g') {
        if (arg[i][3]) ballbackground=arg[i]+3;
        else if (++i>=narg) {
          fprintf(stderr,"WARNING: misplaced -bg option\n");
          break; }
        else
          ballbackground=arg[i]; }
      else if (arg[i][1]=='g') {
        if (arg[i][2]) xwindowhints.geometry=arg[i]+2;
        else if (++i>=narg) {
          xwindowhints.geometry="100x100";
          fprintf(stderr,"WARNING: misplaced -g option\n");
          break; }
        else
          xwindowhints.geometry=arg[i]; }
      else if (strchr("xyz",arg[i][1])) {
        isshift=arg[i][1]-'x';
        shift[isshift]=atof(arg[i]+2);
        if (strend(arg[i])[-1]=='L') shiftinL[isshift]++;
        isshift++; }
      else if (arg[i][1]=='m') {
        matchfile=strdup(arg[i]+2); }
      else
        getoption(arg[i])

#ifdef CHARGES
      if (arg[i][1]=='n') {
        negA=optA(arg[i]);
        loop (j,0,NOCOL)
          if (negA) if (negA[0]==colortab[j].key) { negC=j; goto OK; } }
      if (arg[i][1]=='p') {
        posA=optA(arg[i]);
        loop (j,0,NOCOL)
          if (posA) if (posA[0]==colortab[j].key) { posC=j; goto OK; } }
#endif /*# CHARGES */
    OK:; }
    else {
      /* name or file name */
      if (!molname) molname=arg[i];
      else if (!plbfn) plbfn=arg[i];
      else ERROR(("%s: superfluous file argument",arg[i])) }

  an=(PI/8)/option('a');
  wscale100=100+option('w')*(option('w')<0);

  /* show box and wall setup */
  loop (j,0,NOCOL) if (boxcolor==colortab[j].key) {
    boxturbocolor=colortab[j].turbocolor;
    boxpalcolor=1+COLLEN*j; }
  nwalls=(wallz[0]>-9)+(wallz[1]>-9);

#ifdef SHELL
  distq=Sqr((double)option('j'));
#endif /*# SHELL */

#ifdef DEPTH
  cue=option('q')*0.01; icue=1-cue;
#endif /*# DEPTH */

  if (option('u')==1) grad=.4f,cylsph=1.f;

  if (!molname) {
    prtsfill(Pintro);

    for (;;) {
      char s[8];

      prts_(Phelp);
      if (!fgets(s,8,stdin)) s[0]=0;
      switch (s[0]) {
        case 'i': prtsfill(Pintro); break;
        case 'g': prtsfill(Pgui); break;
        case 'f': prtsfill(Pfiles); break;
        case 'b': prtsfill(Pbasic); break;
        case 'c': prtsfill(Pcolor); break;
        case 'o': prtsfill(Pfileopt); break;
        case 'x': prtsfill(Pbox); break;
        case 's': prtsfill(Pshow); break;
        case 'p': prtsfill(Pplb); break;
        case 'G': prtsfill(Pgeo); break;
        default: return 0; } }

    exit(0); }

  /* removing extension from molname */
  j=strlen(molname);
  ext=molname+j-4;
  if (j>4 && (!strcmp(ext,".mol")||!strcmp(ext,".gol")||!strcmp(ext,".plb"))) {
    char *a;

    j-=4;
    alloczero(a,j+1);
    copy(a,molname,j);
    molname=a; }

  if (option('v')) prt("MOLNAME=%s PLBFN=%s",molname,plbfn);

  if (!plbfn) {
    plbname=molname;
    alloczero(plbfn,j+5);
    sprintf(plbfn,"%s.plb",plbname); }
  else {
    j=strlen(plbfn);
    if (j>4 && !strcmp(plbfn+j-4,".plb")) {
      char *a;

      j-=4;
      alloczero(a,j+1);
      copy(a,plbfn,j);
      plbname=a; }
    else {
      plbname=plbfn;
      prt("PLBNAME.EXT=%s does not have standard extension .plb\n\
(extension will not be removed in derived filenames)"); } }

#include "showbit.c"

  option('g')=(int)BOND;
  xwindowhints.winname=getenv("SHOWNAME");
  if (!xwindowhints.winname) {
    if (strcmp(molname,plbname))
      xwindowhints.winname=string("show %s %s",molname,plbname);
    else
      xwindowhints.winname=string("show %s",molname); }
  xwindowhints.iconame=xwindowhints.winname;
  xwindowhints.icowidth=showbit_width;
  xwindowhints.icoheight=showbit_height;
  xwindowhints.icobits=showbit_bits;

  alloc(li,256);
  readMOL(Fn("mol"));
  if (option('b')) {
    if (abs(option('b'))==1) PDBbackbone();
    else backbone(); }

  percent=strchr(plbname,'%')!=NULL;
  if (percent) {
    iplb=option('h');
    newplb(0); }
  else
    plb=fopen(plbfn,"rb");
  if (option('v')) prt("playback file %s",plbfn);
  if (!plb) ERROR(("no playback file"))

  del=option('d');

  loop (j,0,option('k')) {
    if (1!=fread(hdr+j,sizeof(float),1,plb)) ERROR(("playback file too short"))
    if (!j && hdr[j]!=ns) {
      prt("%g sites in playback file, %d sites in mol-file: will be %s",
          hdr[0],ns,ns<hdr[0]?"replicated":"truncated");
      mol2mol((int)hdr[0]); } }
  if (hdr[1]<0) {
    varL=1;
    if (hdr[1]!=-3) prt("WARNING: unsupported box flag %g",hdr[1]); }
  else {
    L[0]=L[1]=L[2]=hdr[1];
    if (option('v')) prt("cube (old format) L=%g",L[0]); }

  alloc(cfg,(long)ns*sizeof(fvector));
  readnmf(ns);

  fseek(plb,0,SEEK_END);
  maxpos=(ftell(plb)-option('k')*sizeof(float))/(sizeof(fvector)*(ns+varL));
  if (maxpos<1) ERROR(("playback file too short"))

  loop (i,0,nc) bond[i][2]=(bond[i][2]&15)*0x111;
  readcfg();

  /* find min/max for scaling */
  loop (j,0,3)
    /* in box if it exists */
    if (L[j]>0) minr[j]=0,maxr[j]=L[j];
    else maxr[j]=-(minr[j]=9e9);

  /*  09/2003: (PROBLEMS HERE - not for backbone)
      && site[i].r!=7 ensures that sites not shown are not used in centering */
  if (!SXY) loop (i,0,ns) if (cfg[i][0]<99999.0
                    /*.....&& site[i].c!=BLACK*/
                    && site[i].r!=7)
    loop (j,0,3) {
      Min(minr[j],cfg[i][j])
      Max(maxr[j],cfg[i][j]) }

  readGOL();
  /* setgmode(<=BALL) needed to set up colors - or the following: */
#if 1
  setgmode(BALL); /* the default mode */
#else
  /* works, but not faster over slow network */
  startgraph(MODE);
  setmypalette();
  setgmode(BOND2);
#endif

  if (SXY) {
    /* option -S: abs. scaling */
    loop (i,0,3) if (L[i]) center[i]=Lh[i]; }
  else if (option('s')) {
    if (option('s')<0) {
      option('s')=abs(option('s'));
      fprintf(stderr,"WARNING: negative -s ignored, use -S instead\n"); }
    typicalsize=0;
    loop (i,0,3) {
      center[i]=(maxr[i]+minr[i])/2;
      typicalsize=fmax(typicalsize,1.3+maxr[i]-minr[i]); }
    SXY=option('s')*0.009/typicalsize;
    typicalsize*=0.4; }


  zmin=(minr[0]+minr[1]+minr[2]-maxr[0]-maxr[1]-maxr[2])/3;
  zmax=-zmin+2;
#ifdef DEPTH
  zmind=zmin; zmaxd=zmax;
#endif /*# DEPTH */

  bondcolor=YELLOW;
  hbonds();
  redrawmenu=1;

  for (;;) { /* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

    if (redrawmenu) {
      if (menu) makemenu(menu);
      redrawmenu=0; }

    if (redraw) {

      ids[0]=maxx*0.5;
      ids[1]=maxy*0.5;
      Scale[0]= maxy*SXY;
      Scale[1]=-maxy*SXY;

      if (eye[2]<-2.99e33) eye[2]=4*maxy/fabs(Scale[1]);

      if (curmode<=BALL) {
        Scale[0]*=ASPECT;
        Scale[2]=ZBUFDEPTH/(0.5f+zmax-zmin);
        if (!trace) memset(zbuf,0,SCRLEN*sizeof(zbuf[0]));

        zbufmin=9e99; zbufmax=-9e99;
#ifdef DEPTH
        /* reassigning the cue according to the new z depth - slowly */
        ids[2]=(zmaxd-zmind)*0.1;
        zmind+=ids[2]; zmaxd-=ids[2];
#endif /*# DEPTH */
        ids[2]=-zmin*Scale[2];

        if (!trace)
          memset(scrbuf,palcolor(BG),SCRLEN);
        else {
          int ntrace=SCRLEN*(trace*trace-1)/512;

          loop (i,0,ntrace) {
            j=random()%SCRLEN;
            scrbuf[j]=BLACK;
            zbuf[j]=0; } }

        if (dumpmode==NFF || dumpmode==POV) opendump(); }

      ftime(&t0);

      if (curmode<=BALL) {
        /* show box and walls; see also BOND modes */
        int iw,nx=0,ny=0,i,j;
        fvector a,b;
        float grad0=grad;

        if (gradscale>0) grad=gradscale;
        else grad*=-gradscale;

        if (option('|')) {
          /* show the box */

          if (SUM(L)) {
            loop (i,0,2) loop (j,0,2) {
              a[1]=b[1]=i*L[1];
              a[2]=b[2]=j*L[2];
              a[0]=0; b[0]=L[0]; /* x */
              if (i+j==0) showbar2(a,b,1+COLLEN*2,boxpalcolor); /* R */
              else showbar(a,b,boxpalcolor);

              a[2]=b[2]=i*L[2];
              a[0]=b[0]=j*L[0];
              a[1]=0; b[1]=L[1]; /* y */
              if (i+j==0) showbar2(a,b,1+COLLEN*3,boxpalcolor); /* G */
              else showbar(a,b,boxpalcolor);

              a[0]=b[0]=i*L[0];
              a[1]=b[1]=j*L[1];
              a[2]=0; b[2]=L[2]; /* z */
              showsphere(NULL,a,grad,boxpalcolor);
              showsphere(NULL,b,grad,boxpalcolor);
              if (i+j==0) showbar2(a,b,1+COLLEN*4,boxpalcolor); /* B */
              else showbar(a,b,boxpalcolor); } }
          else {
            VO(a,=0)
            showsphere(NULL,a,grad,boxpalcolor);
            loop (i,0,3) {
              VO(b,=0)
              b[i]=10;
              //              showsphere(NULL,b,grad,boxpalcolor);
              showbar(a,b,1+COLLEN*(2+i)); } }

          /* show the walls - undefined/wrong without periodic b.c. */
          nx=L[0]/fabs(wallatomr*0.01)+0.5;
          ny=L[0]/fabs(wallatomr*0.01)+0.5;

          loop (iw,0,2) if (wallz[iw]>-9) {
            /* active of option -Z only */

            a[2]=b[2]=wallz[iw]*L[2];

            loop (i,0,nx) {
              a[0]=b[0]=(i+0.5)*(L[0]/nx);
              if (wallatomr>0)
                loop (j,0,ny) {
                  a[1]=(j+0.5)*(L[1]/ny);
                  showsphere(NULL,a,rad*wallatomr*0.7,boxpalcolor); }
              else {
                a[1]=0; b[1]=L[1];
                showbar(a,b,boxpalcolor); } }

            loop (j,0,ny) {
              a[1]=b[1]=(j+0.5)*(L[1]/ny);
              if (wallatomr>0)
                loop (i,0,nx) {
                  a[0]=(i+0.5)*(L[0]/nx);
                  showsphere(NULL,a,rad*wallatomr*0.7,boxpalcolor); }
              else {
                a[0]=0; b[0]=L[0];
                showbar(a,b,boxpalcolor); } }
          } }
        grad=grad0;
      }

      /*** SHOW molecules ***/
      switch (curmode) {

        case DUMBELL: case DUMBELL2:
          /* show atoms as small spheres */
          if (curmode==DUMBELL || curmode==DUMBELL2) {
            loop (i,0,ns) SHELLTEST(i)
              showsphere(&forclick[i],cfg[i],grad*cylsph,colormark(i)); }
          goto SHOWBAR;

        case BAR: case BAR2:
          /* show marked atoms only or record coordinates for clicking */
          loop (i,0,ns) SHELLTEST(i) {
            if (site[i].mark&4)
              showsphere(&forclick[i],cfg[i],grad*cylsph,colormark(i));
            else
              voidsphere(&forclick[i],cfg[i],grad,colormark(i)); }

        SHOWBAR:
          loop (b,0,nc) {
            i=bond[b][0]; j=bond[b][1];
            if (i==j) {
              SHELLTEST(i)
              showsphere(&forclick[i],cfg[i],site[i].r*rad,colormark(i)); }
            else if (curmode==BAR || curmode==DUMBELL) {
              SHELLTEST(i) SHELLTEST(j)
                showbar(cfg[i],cfg[j],turbo2pal[bond[b][2]&15]); }
            else { /* BAR2, DUMBELL2 & not single atom */
              SHELLTEST(i) SHELLTEST(j) { /* ? order of if's ? */
                if (b>=ncfix)
                  showbar(cfg[i],cfg[j],turbo2pal[bond[b][2]&15]);
                else {
                  if (bond[b][2]&15)
                    showbar2(cfg[i],cfg[j],site[i].c,site[j].c); } } } }

      /* show first -o sites as SPHERES anyway */
          if (option('o')) {
            int from=0,to=option('o');
            if (to<0) from=-to,to=ns;
            loop (i,from,to) showsphere(NULL,cfg[i],site[i].r*rad,colormark(i)); }
          break;

        case BALL:
          loop (i,0,ns) SHELLTEST(i)
            showsphere(&forclick[i],cfg[i],site[i].r*rad,colormark(i));
          loop (b,ncfix,nc) {
            i=bond[b][0]; j=bond[b][1];
            SHELLTEST(i) SHELLTEST(j)
            showbar(cfg[i],cfg[j],turbo2pal[bond[b][2]&15]); }
          break;

        case BOND: case BOND2:
          loop (i,0,ns) {
#ifdef SHELL
            if (!(site[i].shell&4))
              xy[i][0]=-9999;
            else
#endif /*# SHELL */
              project(xy[i],cfg[i]); }

          /* calculate projections to show box and walls */
          {
            int n=24; /* for sure: always have space to show the box */
            int iw,i,j,nx=0,ny=0;
            fvector a,b;

            nx=L[0]/fabs(wallatomr*0.01)+0.5;
            ny=L[1]/fabs(wallatomr*0.01)+0.5;
            n+=(nx+ny)*(2*nwalls);

            if (n!=nxygold) {
              if (nxygold) free(xygold);
              nxygold=n;
              allocarray(xygold,nxygold); }

            n=0;

            if (option('|')) {
              /* box */
              if (SUM(L)) {
                loop (i,0,2) loop (j,0,2) {
                  a[1]=b[1]=i*L[1];
                  a[2]=b[2]=j*L[2];
                  a[0]=0; b[0]=L[0]; /* x */
                  project(xygold[n++],a); project(xygold[n++],b);

                  a[2]=b[2]=i*L[2];
                  a[0]=b[0]=j*L[0];
                  a[1]=0; b[1]=L[1]; /* y */
                  project(xygold[n++],a); project(xygold[n++],b);

                  a[0]=b[0]=i*L[0];
                  a[1]=b[1]=j*L[1];
                  a[2]=0; b[2]=L[2]; /* z */
                  project(xygold[n++],a); project(xygold[n++],b); } }
              else {
                VO(a,=0)
                loop (i,0,3) {
                  VO(b,=0)
                  b[i]=10;
                  project(xygold[n++],a); project(xygold[n++],b); }
                loop (i,n,24) project(xygold[n++],a); /* zeros */ }

              /* walls (if any) */
              loop (iw,0,2) if (wallz[iw]>-9) {
                a[2]=b[2]=wallz[iw]*L[2];

                loop (i,0,nx) {
                  a[0]=b[0]=(i+0.5)*(L[0]/nx);
                  a[1]=0; b[1]=L[1];
                  project(xygold[n++],a); project(xygold[n++],b); }

                loop (j,0,ny) {
                  a[1]=b[1]=(j+0.5)*(L[1]/ny);
                  a[0]=0; b[0]=L[0];
                  project(xygold[n++],a); project(xygold[n++],b); }
              }
              if (n!=nxygold) ERROR(("internal n=%d nxygold=%d",n,nxygold)) } }

        default:; }

      if (percent) {
        if (nmf && iplb>=0 && iplb<3*ns) sprintf(nr,"nmf[%d]=%.2f/cm",iplb,nmf[iplb]);
        else sprintf(nr,"plb=%d",iplb); }
      else
        sprintf(nr,"%d/%d%c%d",lastpos,maxpos," >< o  ~~"[3*swing+(step>0)+2*(step<0)],abs(step));

      sliderpos=(double)(lastpos-1)/(maxpos-1);

      /* show now ! */
      if (curmode<=BALL) {
        int x,y;
        setwritemode(0);
#if XIMAGE
        loopto (y,0,maxy) loop (x,0,maxxn)
          XPutPixel(ximage,x,y,xcoltab[scrbuf[y*maxxn+x]]);
        XPutImage(display,win,gc,ximage,0,0,0,0,maxxn,maxyn);
#else /*# XIMAGE */
        loopto (y,0,maxy) loop (x,0,maxxn) putpixel(x,y,scrbuf[y*maxxn+x]);
#endif /*#!XIMAGE */
      }

      if (curmode>=BOND) {

        setviewport(0,0,maxxn-1,maxyn-1,1);
        if (dumpmode==BONDPS) {
          if (!dumpprepared) fprintf(dump,"showpage\n\
90 rotate 54 -550 translate\n\
0 setgray 1 setlinejoin 1 setlinewidth\n");
          opendump();
          fprintf(dump,"%% frame %d\n",pos); }

        setlinestyle(0,0,1);
        setwritemode(!trace);
        if (trace) {
          int i,j=random()%SCRLEN;
          int ntrace=SCRLEN*(trace*trace-1)/256;

          loop (i,0,ntrace) {
            j=random()%SCRLEN;
            putpixel(j/maxyn,j%maxyn,BLACK); } }

        /*** draw bonds - old ones are redrawn (with xor or if trace by BLUE) ***/

        /* draw fixed bonds (note: bond=bond0 for b<ns) */
        loop (b,0,ncfix) {
          i=bond[b][0]; j=bond[b][1];

          if ( (bigpoint=(i==j && option('f'))) ) bigpoint=abs(option('f'));

          if (curmode==BOND) {
            if (!first) myline(xy0,i,j,palcolor(trace?BLUE:bond0[b][2]&15));
            myline(xy,i,j,palcolor(bond[b][2]&15));
            if (dumpmode==BONDPS)
              fprintf(dump,"newpath %d %d moveto %d %d lineto stroke\n",
                      xy[i][0],xy[i][1],xy[j][0],xy[j][1]); }
          else /* BOND2 */ {
            if (!first) {
              if (trace) myline2(xy0,i,j,palcolor(BLUE));
              else { if (bond0[b][2]&15) myline2(xy0,i,j,0); } }
            if (bond[b][2]&15) myline2(xy,i,j,0); } }

        /* draw hydrogen bonds (and other variable bondings) */
        /* drawing and xor-erasing is done at once to avoid flickering */

        bigpoint=0;
        ncmax=max(nc,nc0);
        loop (b,ncfix,ncmax) {
          if (b<nc0) {
            i0=bond0[b][0]; j0=bond0[b][1];
            if (!first) myline(xy0,i0,j0,palcolor(trace?BLUE:bond0[b][2]&15)); }
          if (b<nc) {
            i=bond[b][0]; j=bond[b][1];
            myline(xy,i,j,palcolor(bond[b][2]&15));
            if (dumpmode==BONDPS)
              fprintf(dump,"newpath %d %d moveto %d %d lineto closepath stroke\n",
                      xy[i][0],xy[i][1],xy[j][0],xy[j][1]); } }

        setcolor(MODE-2); /* ORANGE */
        loop (i,0,ns) if (site[i].mark&4) {
          if (!first) circle(xy0[i][0],xy0[i][1],2);
          circle(xy[i][0],xy[i][1],2); }

        bigpoint=0;

        if (option('|')) {
        if (!first) {
          for (i=0; i<nxygold0; i+=2) {
            b=palcolor(trace?BLUE:
                       i==0?LIGHTRED:
                       i==2?LIGHTGREEN:
                       i==4?LIGHTBLUE:boxturbocolor);
            myline(xygold0,i,i+1,b); } }

        for (i=0; i<nxygold; i+=2) {
          b=palcolor(i==0?LIGHTRED:
                     i==2?LIGHTGREEN:
                     i==4?LIGHTBLUE:boxturbocolor);
          myline(xygold,i,i+1,b); }
        }

        i=nxygold0; nxygold0=nxygold; nxygold=i;
        swap=xygold0; xygold0=xygold; xygold=swap;
        swap=xy0; xy0=xy; xy=swap;

        if (hbdist || option('l')) {
          if (bondeqbond0) {
            bondeqbond0=0;
            /* hydrogen bonds: make a physical copy */
            allocarray(bond0,nbond);
            copy(bond0,bond,nbond*sizeof(bond[0]));
            nbond0=nbond; }
          if (nbond0!=nbond)
            /* bond has been enlarged */
            nbond0=enlarge(&bond0,nbond0,nbond);
          if (option('l')) copyarray(bond0,bond,nc);
          else copyarray(bond0+ncfix,bond+ncfix,nc-ncfix);
          nc0=nc; }
        first=0;
        setviewport(0,0,maxxn-1,maxyn-1,0); }

      if (curmode<BOND) if (slicez) {
        char sl[12];

        sprintf(sl,"%.3f",(slicez-ids[2])/Scale[2]+center[2]);
        setcolor(palcolor(BLACK));
        outtextxy(2,maxy-xfont.height,sl);
        setcolor(palcolor(WHITE));
        outtextxy(3,maxy-xfont.height,sl);
        setcolor(palcolor(bondcolor)); }

      setcolor(MODE-2); /* ORANGE */
      setwritemode(0);
      setfillstyle(0,0); /* = BLACK in both modes */
      if (movesel) {
        bar(1,maxy-xfont.height+2,9*xfont.width+2,maxy);
        outtextxy(1,maxy-xfont.height+3,"SELECTION"); }
      setcolor(palcolor(bondcolor));
      setwritemode(curmode==BALL?0:1);
      XFlush(display);

      if (number) {
        static int last;

        if (!menu) makemenu(0);
        setcolor(MODE-2); /* ORANGE */
        setwritemode(0);
        setfillstyle(0,BLACK); /* = BLACK in both modes */
        j=getmaxx()-strlen(nr)*xfont.width-2;
        if (last) bar(last,maxy-xfont.height+1,getmaxx()-1,maxy);
        outtextxy(j+1,maxy-xfont.height+4,nr);
        last=j;
        setcolor(palcolor(bondcolor));
        setwritemode(curmode==BALL?0:1);
        XFlush(display); }

      if (del) delay(del);

      if (dumpmode==PS || dumpmode==COLORPS || dumpmode==EPS) {
        unsigned u;

        if (!dumpprepared) {
          if (dumpmode==EPS)
            fprintf(dump,"showpage\n"); /* merging problematic - should be banned */
          else
            fprintf(dump,"showpage\n\
90 rotate 54 -550 translate %.2f 480 scale\n",
                    480.0/ASPECT/maxyn*maxxn); }

        opendump();
        if (curmode>BALL) ERROR(("internal"))

        fprintf(dump,"\
%% TEST\n\
%d %d 8 [ 1 0 0 -1 0 %d ]\n\
{ currentfile DS readhexstring pop } %simage",
                maxxn,maxyn,maxyn,
                dumpmode==PS?"":"false 3 color");

        if (dumpmode==PS)
          loop (u,0,SCRLEN)
            fprintf(dump,"%c%02x",u%20?' ':'\n',psgray[scrbuf[u]]);
        else { /* COLORPS */
          int k;
          unsigned u;

          if (option('w')) loop (k,0,3) pal[0][k]=255; /* white background */
          loop (u,0,SCRLEN) {
            fputc(u%10?' ':'\n',dump);
            loop (k,0,3)
              fprintf(dump,"%02x",
                      255-(255-pal[scrbuf[u]][k])*wscale100/100); }
          if (option('w')) loop (k,0,3) pal[0][k]=0; /* back black background */
        }
      } /* PS + COLORPS */

      if (dumpmode==PPM) {
        int k;
        unsigned u;
        char rgb[3];

        opendump();
        fprintf(dump,"P6\n%d %d\n255\n", maxxn,maxyn);

        if (option('w')) loop (k,0,3) pal[0][k]=255; /* white background */
        loop (u,0,SCRLEN) {
          loop (k,0,3)
            rgb[k]=255-(255-pal[scrbuf[u]][k])*wscale100/100;
          fwrite(rgb,1,3,dump); }
        if (option('w')) loop (k,0,3) pal[0][k]=0; /* back black background */
      } /* PPM */

      if (dumpmode==ZBUF) {
        unsigned char *zzz=(unsigned char *)zbuf;
        int i;
        double zrg=(zbufmax-zbufmin)*Scale[0]/255.9999;

        loop (i,0,SCRLEN) zzz[i]=zbuf[i]?(zbuf[i]/Scale[2]+(zmin-zbufmin))/(zbufmax-zbufmin)*255.999:0;
        opendump();
        fprintf(dump,"P5\n# %g (multiply by to get height in pixels)\n%d %d\n255\n",zrg, maxxn,maxyn);
        fwrite(zbuf,SCRLEN,1,dump); }

      if (dumpmode) {
        if (dumpstat==START) dumpstat=CONT;
        if (dumpstat==ONE || dumpstat==SERIES) closedump(0); }

    } /* redraw */

    redraw=0; reread=0;

    if (option('[')<=0) fastmove=(long)maxpos*abs(option('['))/100;
    else fastmove=option('[');
    if (fastmove<=0) fastmove=1;

    drawslider(0);

    if (kbhit()) switch (redraw++,key=readkey()) {

        /* mouse: drag left button */
        case LEFTDRAG:
          if (isinslider()) {
            drawslider(0);
            pos=(maxpos-1)*sliderpos+1.5;
            reread=1;
            redrawmenu=step;
            break; }
          else if (!menu || (dosmouse.x>=0 && dosmouse.y>=0 && dosmouse.x<maxxn && dosmouse.y<maxyn)) {
            if (movesel) {
              if (dosmouse.dx) rotateselection(1,dosmouse.dx*mouserot);
              if (dosmouse.dy) rotateselection(0,dosmouse.dy*mouserot); }
            else {
              if (dosmouse.dx) rotate(1,dosmouse.dx*mouserot);
              if (dosmouse.dy) rotate(0,dosmouse.dy*mouserot); } }
          break;

        /* mouse: drag mid button */
        case MIDDRAG:
          if (movesel)
            moveselection((double)dosmouse.dx/Scale[0],(double)dosmouse.dy/Scale[1]);
          else {
            eye[0]-=(double)dosmouse.dx/Scale[0];
            eye[1]-=(double)dosmouse.dy/Scale[1]; }
          break;

        /* mouse: drag right button */
        case RIGHTDRAG:
          if (dosmouse.dx || dosmouse.dy) {
            double rx=dosmouse.x-maxx/2.;
            double ry=dosmouse.y-maxy/2.;
            double a=Sqr(rx)+Sqr(ry);
            if (a>30) {
              a=(dosmouse.dx*ry-dosmouse.dy*rx)/a;
              if (movesel)
                rotateselection(2,a);
              else
                rotate(2,a); } }
          break;

        /* single click (no move) */
        case RIGHTCLICK:
        case MIDCLICK:
        case LEFTCLICK:

          if (key==LEFTCLICK && isinslider()) {
            drawslider(0);
            pos=(maxpos-1)*sliderpos+1.5;
            reread=1;
            redrawmenu=step; }
          else  {
            int rrmax=((unsigned char*)"\66\12\66\377")[clickmode];
            /* clickmode=0: mark ONE atom (isel) and print info:
                 wire modes: the closest atom to the click point (but max sqrt(rrmax) from it)
                 ball modes: the visible atom
               clickmode>0:
                 mark ALL atoms within radius sqrt(rrmax), print # of atoms */
            float minsel=-9e9; /* zbuf */
            int iminsel=0x7fffffff; /* wire: dist from click */
            int isel=-1;

            if (curmode<BOND) {
              /* not passed from showing */
              loop (i,0,ns) project(xy0[i],cfg[i]);
              if (!clickmode) rrmax=0x7fffffff; }

            loop (i,0,ns) {
              /* long long to avoid overflow if too much zoomed */
              long long unsigned ii=xy0[i][0]-dosmouse.x;
              long long unsigned jj=xy0[i][1]-dosmouse.y;
              long long unsigned irr=ii*ii+jj*jj;

              if (curmode<BOND && !clickmode) rrmax=forclick[i].rr;

              if (irr<=rrmax) {
                if (clickmode)
                  markatom(key,i); /* all clicked atoms marked */
                else if (curmode<BOND) {
                  /* closest to the screen will be marked */
                  if (forclick[i].z>minsel) minsel=forclick[i].z,isel=i; }
                else {
                  /* closest to the click x,y will be marked */
                  if (irr<iminsel) iminsel=irr,isel=i; } }

            } /* i */

            if (!clickmode && isel>=0) {
              markatom(key,isel);
              printmarkedinfo(isel); }
          }
          break;

        case 'I': {
          setwritemode(!trace);
          setcolor(MODE-2);
          loop (i,0,ns) {
            site[i].mark^=4;
            if (curmode>=BOND) circle(xy0[i][0],xy0[i][1],2); }
          setwritemode(0);
          break; }

        case 'U'&31:
          setwritemode(1);
          setcolor(MODE-2);
          loop (i,0,ns) if (site[i].mark&4) {
            site[i].mark-=4;
            if (curmode>=BOND) circle(xy0[i][0],xy0[i][1],2); }
          setwritemode(0);
          break;

        case 'w':
          swing=(swing+1)%3;
          if (!swing) step=0;
          if (swing==2 && step<0) step=abs(step);
          if (swing && step==0) step=1;
          redrawmenu=1;
          break;

        case 'o':
          option('o')-=option('o')>-ns; break;

        case 'O':
          option('o')+=option('o')<=ns-1; break;

        case 'F'&31:
          bodyframe=!bodyframe;
          {
            int a,b;

            loop (a,0,3) loop (b,0,3) rot[a][b]=drot[a][b]=a==b;
          }
          redrawmenu=1;
          break;

        case 'I'&31: {
          int a,b;

          if (drot[0][0]+drot[1][1]+drot[2][2]==3)
            loop (a,0,3) loop (b,0,3) rot[a][b]=drot[a][b]=a==(b+1)%3;
          else if (drot[1][0]+drot[2][1]+drot[0][2]==3)
            loop (a,0,3) loop (b,0,3) rot[a][b]=drot[a][b]=a==(b+2)%3;
          else
            loop (a,0,3) loop (b,0,3) rot[a][b]=drot[a][b]=a==b;
          }
          break;

        reinit:
        case 'i':
          fseek(plb,sizeof(float)*option('k'),0);
          pos=0; step=1;
          redrawmenu=1;
        break;

#ifdef SHELL
        case 'J':
          distq*=1.44;
        case 'j':
          distq/=1.2; if (distq==0) distq=9;
          lastpos=-1; /* forces reread even if pos not changed */
          reread=1;
          break;
#endif /*# SHELL */
        case 'K':
          clickmode=(clickmode+1)%4;
          if (clickmode)
            fprintf(stderr,"clickmode=%d: all atoms within %s radius marked\n",
                    clickmode,"small\0\0medium\0large"+clickmode*7-7);
          else
            fprintf(stderr,"clickmode=0: one best matching atom marked on click\n");
          goto ctrll;

        case 'e':
          closedump(1);
          redrawmenu=1;
          break;

        case 'W':
          dumpmode=PS;
          if (curmode==BOND) dumpmode=BONDPS;
          else goto testball;
          goto startdump;
        case 'C':
          dumpmode=COLORPS;
          goto testball;
        case 'E':
          dumpmode=EPS;
          goto testball;
        case 'L':
          dumpmode=PLB;
          goto testball;
        case 'A':
          dumpmode=ATM;
          goto testball;
        case 'P':
          dumpmode=PPM;
          goto testball;
        case 'B':
          dumpmode=ZBUF;
          goto testball;
        case 'V':
          dumpmode=POV;
          goto testball;
        case 'N':
          dumpmode=NFF;
        testball:
          if (curmode>=BOND) {
            dumpmode=NONE;
            goto badkey; }
        startdump:
          if (dump) {
            ERROR(("dump opened - closing"))
            closedump(1);
            goto badkey; }
          fprintf(stderr,
                  "dumping %s: select in graphic window:\n\
(o)ne file, (O)ne+render, (s)eries%s%s\n",
                  dumpext[dumpmode],
                  dumpmode==PPM ? ", erase old pic + (a)nimated GIF" : "",
                  dumpmode!=PPM && dumpmode!=ZBUF ? ", (m)erge, (M)erge+render" : "");
          makemenu(2);
          redrawmenu=1;
          render=animated.gif=0;
          switch (readkey()) {
            case 'O': render++;
            case 'o': dumpstat=ONE; break;
            case 's': dumpstat=SERIES; break;
            case 'a': dumpstat=SERIES; animated.gif++;
              system(string("rm %s.????.%s",percentfn,dumpext[dumpmode]));
              break;
            case 'M': render++;
            case 'm': dumpstat=START;
                      if (dumpmode!=PPM && dumpmode!=ZBUF) break;
            default: dumpmode=NONE; goto badkey; }

          preparedump(dumpmode==PLB && movesel);
          redraw=1; /* unnecessarily in most cases, but fool proof */
          break;

        case 'r':
          if (curmode>=BAR && curmode<=BALL) {
            rad/=1.0821439f;
            Max(rad,0.001f) }
          if (curmode>=BAR && curmode<=DUMBELL2) {
            grad/=1.0821439f;
            Max(grad,0.01f) }
          break;
        case 'R':
          if (curmode>=BAR && curmode<=BALL) {
            rad*=1.05f;
            if (option('u')) Min(rad,0.03f) }
          if (curmode>=BAR && curmode<=DUMBELL2) {
            grad*=1.05f;
            if (option('u')) Min(grad,2.f) }
          break;

        case 'd':
          if (curmode==DUMBELL || curmode==DUMBELL2) {
            rad=0.01f; grad=.1f; cylsph=3.f; }
          break;
        case 'D':
          if (curmode==DUMBELL || curmode==DUMBELL2) {
            rad=0.01f; grad=.4f; cylsph=1.f; }
          break;
        case 'U':
          if (curmode==DUMBELL || curmode==DUMBELL2) cylsph*=1.05;
          break;
        case 'u':
          if (curmode==DUMBELL || curmode==DUMBELL2) {
            cylsph/=1.0821439f;
            if (option('u')) Max(cylsph,1.f); }
          break;

#ifdef SLICES
        case '(':
          slicez-=ZBUFDEPTH/233; break;
        case ')':
          slicez+=ZBUFDEPTH/144; break;
        case '<':
          slicez-=ZBUFDEPTH/16; break;
        case '>':
          slicez+=ZBUFDEPTH/16; break;
        case '_':
          if (slicez) slicez=0;
          else slicez=ZBUFDEPTH/2;
          redrawmenu=1; redraw=1;
          break;
        case 'l':
          if (slicecolor) slicecolor=0;
          else slicecolor=1;
          break;
#endif /*# SLICES */

          /*        case 3: exit(0);  ctrl-c */

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

        case 'n':
          number=!number;
          goto ctrll;
        case 'L'&31:
#if 0
          if (step) {
            del=0;
            step=sign(step); // why this?
            break; }
#endif
         ctrll:
          first=1;
          if (curmode>=BOND)
            clearviewport();
          else {
            arrayzero(zbuf,SCRLEN);
            memset(scrbuf,0,SCRLEN); }
          redrawmenu=1;
          break;

        case ' ':
          if (step) step=0;
          else step=laststep;
          redraw=step;
          redrawmenu=1;
          break;

        case 'f': /* forward (NEW) */
          swing=0;
          if (!step) step=laststep;
          if (!step) step=1;
          step=abs(step);
          redrawmenu=1;
          break;

        case 'b': /* backward (NEW) */
          swing=0;
          if (!step) step=laststep;
          if (!step) step=1;
          step=-abs(step);
          redrawmenu=1;
          break;

        case 'B'&31: /* toggle backward-forward (former 'b') */
          step=-step;
          if (!step) step=1;
          redrawmenu=1;
          break;

        case 'O'&31:
        case F6:
          font=(font+2)%4; // to extend to font 4 : %6
          /* must totally redraw */
        case X_REDRAW_ME:
          newmode=curmode; curmode=UNDEF;
          goto setnewmode;
        case '!': newmode=BAR; goto setnewmode;
        case '@': newmode=BAR2; goto setnewmode;
        case '#': newmode=DUMBELL; goto setnewmode;
        case '$': newmode=DUMBELL2; goto setnewmode;
        case '%': newmode=BALL; goto setnewmode;
        case '^': newmode=BOND; goto setnewmode;
        case '&': newmode=BOND2; goto setnewmode;
        case 'g':
          newmode=(enum mode_e)(1+((int)curmode%(int)BOND2));
          goto setnewmode;
        case 'G':
          newmode=(enum mode_e)(1+(((int)(BOND2)-2+(int)curmode)%(int)BOND2));
        setnewmode:
          redrawmenu=1;
          setgmode(newmode);
          first=1;
          if (curmode>=BOND) clearviewport();
          break;

        case F2:
        case 'S'&31:
          dumponeplb(0);
          break;

        case F3:
        case 'R'&31:
          dumponeplb(1);
          break;

        case 's':
          if (abs(step==1)) {
            del+=del/2+8;
            if (del>1000) del=1000; /* max. 1 s delay */ }
          else step=(step+sign(step))/2;
          goto showdel;
        case 'S':
          if (del) del -= (del+2)/3;
          else step=(step*3+sign(step))/2;
        showdel:
          if (del || step) fprintf(stderr,"delay = %d ms, step=%d\n",del,step);
          break;

        case 't': trace=!trace; redrawmenu=1; break;
        case 'T': trace+=2; redrawmenu=1; break;

        case WHEELFORWARD:
        case '+': SXY*=1.05f; break;
        case WHEELBACKWARD:
        case '-': SXY/=1.0821439f; break;
        case '*': an*=2; if (an>2) an=(2*PI)/option('a'); redraw=0; break;
        case '/': an/=(1+(an>1e-4)); redraw=0; break;
        case 'x':
        case 'y':
        case 'z':
          if (movesel) rotateselection(key-'x',an);
          else rotate(key-'x',an);
          break;
        case 'X':
        case 'Y':
        case 'Z':
          if (movesel) rotateselection(key-'X',-an);
          else rotate(key-'X',-an);
          break;

        case '\'':
          match(cfg,ns,1e-4); break;
        case '\"':
          match(cfg,ns,1e-6); break;

        case 'm':
          movesel=!movesel;
          redrawmenu=1;
          break;

        case F10:
#if XIMAGE==1
          fprintf(stderr,"menu not available with #define XIMAGE 1 - recompile!\n");
#else /*# XIMAGE==1 */
          menu=!menu;
          if (!menu) erasebuttons();
          redraw=1;
          redrawmenu=1;
          newmode=curmode; curmode=UNDEF;
          goto setnewmode;
#endif /*#!XIMAGE==1 */

        case UP:
          if (movesel) moveselection(0,-8/Scale[1]);
          else eye[1]+=8/Scale[1];
          break;
        case DOWN:
          if (movesel) moveselection(0,7/Scale[1]);
          else eye[1]-=5/Scale[1];
          break;
        case RIGHT:
          if (movesel) moveselection(8/Scale[0],0);
          else eye[0]-=8/Scale[0];
          break;
        case LEFT:
          if (movesel) moveselection(-7/Scale[0],0);
          else eye[0]+=5/Scale[0];
          break;

        case 'A'&31: /* ctrl-A */
        case HOME:
          if (central) eye[2]*=0.9;
          //          putv(eye)
          break;
        case 'E'&31: /* ctrl-E */
        case END:
          if (central) eye[2]*=1.185870252567268;
          //          putv(eye)
          break;

        case F7:
          if (central) { SXY*=1.05f; eye[2]=1.05; }
          break;
        case F8:
          if (central) { SXY/=1.0821439f; eye[2]*=1.0821439f; }
          break;

        case ',':
          pos=0; reread=1; /* reread=1 sets step=0 ! */
          redrawmenu=step;
          break;
        case '.':
          fseek(plb,0,2);
          maxpos=(ftell(plb)-option('k')*sizeof(float))/(sizeof(fvector)*(ns+varL));
          pos=maxpos; reread=1;
          redrawmenu=step;
          break;

        case F11:  {
            int i,x=0,y=0;

            outbgtextxy(4,4,2,palcolor(BLACK),palcolor(WHITE),"ball-mode palette:");

            loop (i,0,MODE) {
              if (i) {
                x=((i-1)%COLLEN)*20+4;
                y=((i-1)/COLLEN)*20+30; }
              setwritemode(0);
              setfillstyle(0,i);
              bar(x,y,x+19,y+19); }
            outbgtextxy(4,y+4,2,palcolor(BLACK),palcolor(WHITE),"wire-mode palette:");
            y+=30;
            loop (i,0,17) {
              x=i*20+4;
              setfillstyle(0,palcolor(i));
              bar(x,y,x+19,y+19); }

            y+=26;
            outbgtextxy(4,y,2,palcolor(RED),palcolor(YELLOW),"press ESC to quit");
            while (readkey()!=ESC);
            redraw=1;
            redrawmenu=1; }
          break;

        case F4:
          if (presenter) goto reinit;
          { /* F4 */
          char line[64];
          int aty=2;

          aty=outbgtextxy(2,aty,font,palcolor(RED),palcolor(YELLOW),"WINDOW INACTIVE!    ");
          aty=outbgtextxy(2,aty,font,palcolor(RED),palcolor(YELLOW),"Respond in the      ");
          aty=outbgtextxy(2,aty,font,palcolor(RED),palcolor(YELLOW),"controlling terminal");

          printf("frame to show=");
          getsbufsize=64;
          gets(line);
          redraw=1;
          pos=atoi(line); reread=1; step=0;
          goto ctrll; }
        case '{': pos-=fastmove; reread=1; redrawmenu=step; break;
        case '}': pos+=fastmove; reread=1; redrawmenu=step; break;
        case PGUP:
        case '[':
          if (percent) { newplb(-1); goto reinit; }
          pos--; reread=1;
          redrawmenu=step;
          break;
        case PGDN:
        case ']':
          if (percent) { newplb(1);  goto reinit; }
          pos++; reread=1;
          redrawmenu=step;
          break;

        case ':': bookmark=pos; break;
        case ';': i=pos,pos=bookmark,bookmark=i; reread=1; break;

        case 'c':
          eye[0]=eye[1]=0;
          eye[2]=4*maxy/fabs(Scale[1]);
          break;

        case F5:
          if (presenter) goto reinit;
        case 'T'&31:
          /* very simple! */
          fprintf(stderr,"to find the n-th largest cluster, type (in the X-window) number 1..9 or Fn\n");
          i=readkey();
          if (isascii(i)) {
            if (isdigit(i)) findclusters(i-'0'+10*(i=='0'));
            else if (isalpha(i))  findclusters(toupper(i)-'@'); }
          else if (i>=F1) findclusters(i-F1+1);
          break;

        case '\\':
          option('|') = !option('|');
          if (curmode>=BOND) goto ctrll; /* redraw */
          break;

        case '|': if (curmode<=BALL) wallatomr=-wallatomr;
          break;

        case '=':
          central=!central;
          fprintf(stderr,"%s projection\n",central?"central":"parallel");
          break;

        case 'C'&31:
          ihbcolor=(ihbcolor+1)%NOCOL;
          goto ctrlh;

        case 'H'&31:
          xhbdist();
        ctrlh:
          hbcolor=colortab[ihbcolor].turbocolor*0x11;
          reread=1; goto ctrll;
        case 'h':
          if (!hbdist) xhbdist();
          hbdist-=0.2;
          if (hbdist<0.8) hbdist=0.8;
          fprintf(stderr,"H-bond %s-%s distance limit:=%g\n",donor,acceptor,hbdist);
          reread=1; goto ctrll;
        case 'H':
          if (!hbdist) xhbdist();
          hbdist+=0.2*1.618034;
          if (hbdist>6) hbdist=6;
          fprintf(stderr,"H-bond %s-%s distance limit:=%g\n",donor,acceptor,hbdist);
          reread=1; goto ctrll;

        case F1:
        case '?':
          lastmode=curmode;
          help();
          redrawmenu=1;
          redraw=1;
          newmode=lastmode; curmode=UNDEF;
          goto setnewmode;

        case '0':
          del=readkey()-'0';
          del=10*del+(readkey()-'0');
          step = step<0 ? -del : del;
          del=0;
          break;

        default:;
          if (key>='1' && key<='9') {
            del=0;
            step = step<0 ? '0'-key : key-'0';
            break; }

        badkey: ;
          delay(100); }

    if (step || reread) {
      if (reread) step=0;
      if (step) if (swing) {
        if (swing==2)
          if ( (pos<=1 && step<0) || (pos>=maxpos && step>0)) step=-step;
        if (swing==1) {
          if (step<0) step=abs(step);
          if (pos>=maxpos) pos=0; } }
      pos+=step;
      Max(pos,1)
      if (pos>=maxpos) {
        fseek(plb,0,2);
        maxpos=(ftell(plb)-option('k')*sizeof(float))/(sizeof(fvector)*(ns+varL)); }
      Min(pos,maxpos)

      nc=ncfix; /* ? */
      loop (i,0,nc) bond[i][2]=(bond[i][2]&15)*0x111; /* ? */
      redraw |= readcfg();
      hbonds(); }

    if (step) laststep=step;

    } /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */

 TheEnd:;
  fclose(plb);
  closedump(1);
/*.....if (curmode!=BALL) ?DOS*/
  closegraph();

  if (option('v')) {
    prt("%d fixed bonds, %d hydrogen bonds",ncfix,nc-ncfix);
    prt("maxr: %10.6f %10.6f %10.6f",maxr[0],maxr[1],maxr[2]);
    prt("minr: %10.6f %10.6f %10.6f",minr[0],minr[1],minr[2]);
    prt("size: %10.6f %10.6f %10.6f",
        maxr[0]-minr[0],
        maxr[1]-minr[1],
        maxr[2]-minr[2]); }

  return 0;
}
