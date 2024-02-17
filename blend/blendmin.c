/*
  inefficiency in blendmin.c:
  evaluation of Upot() is sometimes unnecessarily duplicated
  should store individual terms (dihedrals, ...) and print them separately
  arrays r,f,... unnecessarily allocated
  backcopy from rotated (by show) site[].r to r0 inefficient and cumbersome
*/

#define xCONTOURPLOT
/*
special patch for contour plot */

/* debug forces and cutoff */
#define xBLENDDEBUG

/* keep-keep omission, not supported under !CUTOFF */
/*.....#define OMITKEPT*/

#ifdef POLAR
#  ifdef CUTOFF
/* WARNING CUTOFF turned off because of POLAR */
#    undef CUTOFF
#  endif /*# CUTOFF */
#else /*# POLAR */
#  define CUTOFF /* cutting off the potential enabled */
#endif /*#!POLAR */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLENDMIN.C

This module evaluates energy and forces for a given configuration,
performs energy minimization, and shows the molecule graphically.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#define TRACE(X)
/*.....#define TRACE(X) prtc(X);*/

#define eyedist 0.35 /* stereo: ratio of the distance of eyes/screen height */

#include "ground.h"
#include "options.h"
#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "blendmed.h"
#include "blendmin.h"
#include "rndgen.h"
#include "units.h"


#define fullangle 0 /* 0:dihedral potential in (-180,180); 1: (0,360) */
typedef vector fvector;
#include "dihangle.c"

/* legacy patch: also in blendu.c */
int menu=1,redrawmenu=1;

#include "xdraw.h"

#include <time.h>

#ifdef POLAR
#  include "log1.c"
static vector *mu; /* dipole moments */
#endif /*# POLAR */

static double Upot(vector *r,vector *f,species_t *spec,int toprint);

extern int measure;

static int by;
static int testout=1; /* used by Upot and virial */

volatile int sig=0;
#include <signal.h>

#include "jacobi.c"

void asksig(int signr) /********************************************* asksig */
{
  if (sig++) exit(-1); /* ^C typed twice */

  (*signal)(SIGINT,asksig);
  fprintf(stderr,
    "\nINTERRUPT(%d) ACCEPTED, wait to finish cycle or type ^C again to exit immediately\n",
    signr);
}

bond_t *b0;
struct angle_s *a0;
torsion_t *d0,*i0,*ar0,*cis0;
fixdih_t *fd0;
int constraints; /* set for models with constraints (bonds/angles/fixed sites) */
int aty_n;

/* WARNING e^2 / (4 PI eps0) * NA (kcal/based electrostatics):
   5/2019 (defined NA,k): 332.063713299455 AA kcal/mol
   CODATA 2014:           332.063710541    AA kcal/mol
   CODATA 2006:           332.063709       AA kcal/mol
   cf. sim/units.h */
#define obfuscate 332.063713299455

static double ULJ,Uel,signpairs;
static int verbose2,maxpairs,checkpairs;
static maxpair_t *maxpair;
static site_t *site;

#ifdef CUTOFF
#  include "blendss1.c"
#else /*# CUTOFF */
#  ifdef GAUSSIANCHARGES
#    error NOT IMPLEMENTED - in future, blendss3.c
#  else /*# GAUSSIANCHARGES */
#    include "blendss2.c"
#  endif /*#!GAUSSIANCHARGES */
#endif /*#!CUTOFF */

#include "blenddm.c" /* dipolemoment() */

static double zeye,eye;
static int cross=3; /* can be in 0..3, but only 1 and 3 used */

static enum label_e {
  NONE, ID, TYPE, NUMBER, CHARGE,
  POLARIZABILITY, SATURATION,
  /* DIPOLE_eAA, DIPOLE_D must be at the end */
  DIPOLE_eAA, DIPOLE_D
} label;

//#  define NBONDCOL 4
//static int bondcol=0,bondcols[NBONDCOL]={0,YELLOW,WHITE,DARKGRAY};
// simplified:
#define NBONDCOL 2
static int bondcol=0,bondcols[NBONDCOL]={0,WHITE/*,YELLOW*/};
static int col,grmode=2;
struct {
  double q,dq; /* charge to find, to be in <q-dq,q+dq> */
  int i;       /* atom # to find or found */
  int x,y;     /* x,y to find (pick by mouse) */
  int l;       /* effective length of s[] */
  char s[16];  /* id or type to find */
  vector r;    /* (old) position of site[i].r */
} find = { 1e3,1e-5,-1,-1 };

#define palcolor(COL) COL
int font=2,moving;
#include "xdrawhelp.c"
#include "blendhelp.c"
#include "blendmenu.c"

static char pickinfo[256];
static unsigned ch;
static int lasti=-1,ns,redraw,center=-1,marklabeled,found;

static void draw3D(double x,double y,double z) /********************* draw3D */
{
  double q=zeye/(zeye-z);

  draw(x*q,y*q);
}

static int match3D(double x,double y,double z,int ix,int iy) /****** match3D */
{
  double q=zeye/(zeye-z);
  int rr=5;

  //if (grmode==2) rr="\36\5\13"[cross];
  //  if (grmode==3) rr="\52\13\36"[cross];
  rr="\132\14\36\70"[cross];

  ix-=SX(x*q); iy-=SY(y*q);

  return ix*ix+iy*iy<=rr;
}

/* anaglyph */
static void draw3Ds(double x,double y,double z,int stereo) /******** draw3Ds */
{
  double q=zeye/(zeye-z);

  draw(x*q+(1-q)*stereo*eye,y*q);
}

static void floodmark(int center) /******************************* floodmark */
{
  int n,i;

  site[center].count=1;

  loop (n,0,site[center].nnbr)
    if (!site[i=site[center].nbr[n]].count) floodmark(i);
}

static void makelabel(char *Us,double nr) /*********************** makelabel */
{
  char *c;

  sprintf(Us,outtext_f,nr);
  c=strend(Us);
  while (*(--c)=='0') *c=0;
  if (*c=='.') *c=0;
}

static void showlabel(int i) /************************************ showlabel */
{
  site_t *s=site+i;

  if (label) {
#ifdef POLAR
    static char Us[40],*c;
#else /*# POLAR */
    static char Us[16],*c;
#endif /*#!POLAR */

    settextstyle(0,0,1); /*?*/
    settextjustify(1,1);
    if (grmode==1) {
      setcolor((bondcols[bondcol]&7)==CYAN ? YELLOW : LIGHTCYAN); }
    else
      setcolor(col<=BLUE ? YELLOW : col==LIGHTGRAY ? CYAN : BLUE);
    up();
    draw3D(s->r[0],s->r[1],s->r[2]);
    switch (label) {
    case CHARGE:
      if (find.q<999 && fabs(find.q-s->charge)>find.dq) goto pickup;
      makelabel(Us,s->charge);
      outtext(Us);
      if (find.q<999) find.i=i;
      break;
#ifdef POLAR
    case DIPOLE_eAA:
    case DIPOLE_D:
      if (s->polar) {
        static float dipoleunit[2]={1,4.8032068};
        makelabel(Us,s->mu*dipoleunit[label==DIPOLE_D]);
/*.....      strcat(Us,label==DIPOLE_D?"D":"eAA");*/
        outtext(Us); }
      break;
    case SATURATION:
      if (s->polar && s->sat!=0) {
        makelabel(Us,s->sat);
/*.....      strcat(Us,label==DIPOLE_D?"D":"eAA");*/
        outtext(Us); }
      break;
    case POLARIZABILITY:
      if (s->polar) {
        if (!s->pol) strcpy(Us,"-");
        else if (s->polar&POL_ISO)
          makelabel(Us,((isotropicparm_t*)s->pol)->alpha);
        else {
          makelabel(Us,((axial_t*)s->pol)->parm.alpha);
          strcat(Us,",");
          makelabel(strend(Us),((axial_t*)s->pol)->parm.alphazz); }
        outtext(Us); }
      break;
#endif /*# POLAR */
    case NUMBER:
      if (find.i>=0 && find.i!=i) goto pickup;
      sprintf(Us,"%d",i); outtext(Us);
      break;
    default:
      c = label==ID ? s->id : atom[s->type].name;
      if (find.s[0] && memcmp(find.s,c,find.l)) goto pickup;
      outtext(c);
      if (find.s[0]) find.i=i;  }
    if (marklabeled) site[i].count=1,found++; }

 pickup:

  if (find.x>=0 && match3D(s->r[0],s->r[1],s->r[2],find.x,find.y)) {
    static vector last;
    static int llasti=-1,lllasti=-1;
    double rr;
    static vector dummyf;
    static angleparm_t ap; /* zeros OK */

    if (ch==RIGHTCLICK) {
      site[i].count=0;
      redraw=1; }
    else {
      double psi=-999;

      if (lasti<0) { VV(last,=s->r) lasti=i; }
      rr=SQRD(last,s->r);
      phi=-1;
      measure=2;
      if (i!=lasti && llasti!=lasti && lasti>=0 && llasti>=0) {
        anglepot(s->r,last,site[llasti].r,dummyf,dummyf,dummyf,&ap);
        if (lllasti>=0 && lllasti!=llasti && lllasti!=lasti && lllasti!=i)
          psi=dihedral(s->r,last,site[llasti].r,site[lllasti].r); }
      sprintf(pickinfo,
              infoxyz_f,
              i,s->id,atom[s->type].name,s->charge,s->r[0],s->r[1],s->r[2]);
      if (rr>0) sprintf(strend(pickinfo),infor_f,lasti,sqrt(rr));
      if (phi>=0) sprintf(strend(pickinfo),infoa_f,llasti,phi*(180/PI));
      if (psi>-999) sprintf(strend(pickinfo),infod_f,lllasti,psi);

      if (strlen(pickinfo)>79) {
        char *b=strchr(pickinfo,'(');
        char *e=strchr(pickinfo,')');

        if (b && e && b<e) memmove(b+1,e,strlen(e)+1); }

      VV(last,=s->r)
      lllasti=llasti; llasti=lasti; lasti=find.i=i;
      if (ch==LEFTCLICK || center<0) center=i;
      switch (ch) {
        case MIDCLICK: floodmark(i);
        case LEFTCLICK: site[i].count=redraw=1; }
      fprintf(stderr,"%s\n",pickinfo); }
    }
} /* showlabel */

static int darkcolor(int col) /*********************************** darkcolor */
/* make dark color */
{
  if (col==LIGHTGRAY) col=DARKGRAY;
  else if (col==WHITE) col=LIGHTGRAY;
  else if ( (col-=8)<0 ) col=BLACK;

  return col;
}

/* support for functions (w)rite and W=restore */
static struct remember_s {
  struct remember_s *next;
  int key;
  double x,y;
} *remember0;

static void remember(int key,double x,double y) /****************** remember */
{
  struct remember_s *rem;

  alloc(rem,sizeof(struct remember_s));
  rem->next=remember0;
  remember0=rem;
  rem->key=key;
  rem->x=x; rem->y=y;
}

static void forget(int back) /*************************************** forget */
{
  struct remember_s *rem,*next;
  int i,ii,jj;
  double z;

  for (rem=remember0; rem; rem=next) {
    if (back) {
      if (rem->key<0) {
        loop (i,0,ns) {
          site[i].r[0]-=rem->x;
          site[i].r[1]-=rem->y; } }
      else {
        ii=rem->key; jj=(ii+1)%3;
        loop (i,0,ns) {
          z=site[i].r[ii];
          site[i].r[ii] =  rem->x*z - rem->y*site[i].r[jj];
          site[i].r[jj] =  rem->y*z + rem->x*site[i].r[jj]; } } }

    next=rem->next;
    free(rem); }

  remember0=NULL;
}

double hotkeyE(species_t *spec,int E,double U0) /******************* hotkeyE */
{
  double tc,dipm;
  vector *r,*f;
  vector Qdiag;
  FILE *geo;
  int n,i,j,k;
  char line[128];

  ns=spec->ns;
  strcpy(spec->ext,".geo");
  geo=fopen(spec->fn,"rt");

  alloc(r,ns*sizeof(r[0]));
  alloczero(f,ns*sizeof(f[0]));
  loop (i,0,ns) VV(r[i],=site[i].r)

  tc=totalcharge(spec,-1);
  dipm=dipolemoment(spec,tc,Qdiag);
  if (E)
    U0=Upot(r,f,spec,1);
  else {
    if (out!=stdout) fprintf(stderr,"U0 = %g\n",U0);
    prt("! U0=%.12g kcal/mol = %.12g J/mol = %.12g K kB",
        U0,U0*kcal,U0*(kcal/Eunit)); }

  prt("! charge=%.4f e  dipole moment = %.9f e AA = %.9f D",
      tc,dipm,dipm*4.8032068);
  prt("! quadrupole moment [e AA^2] = diag(%.9f %.9f %.9f)",
      Qdiag[0],Qdiag[1],Qdiag[2]);
  prt("! quadrupole moment [C m^2] = diag(%g %g %g)",
      Qdiag[0]*1.60217733e-39,Qdiag[1]*1.60217733e-39,Qdiag[2]*1.60217733e-39);
  prt("! quadrupole moment [D AA = Buckingham] = diag(%g %g %g)",
      Qdiag[0]*4.8032068,Qdiag[1]*4.8032068,Qdiag[2]*4.8032068);

#ifdef POLAR
  totalpolarizability(spec);
#endif /*# POLAR */

#ifdef CONTOURPLOT
  {
    int i,j;
    FILE *cd=fopen("contour.xyz","wb");

    loopto (i,-100,100) {
      r[0][0]=i*0.05;
      fprintf(stderr,"%5.2f\r",r[0][0]);
      loopto (j,-100,100) {
        float U0,UU;
        double q;

        r[0][1]=j*0.05;
        q=spec->site[0].charge;
        spec->site[0].charge=0;
        U0=Upot(r,f,spec,0);
        spec->site[0].charge=q;
        UU=Upot(r,f,spec,0)-U0;
        if (!(UU>-1e4 && UU<1e4)) UU=1e4;
        fwrite(&UU,1,4,cd); } }
    fclose(cd); }
#endif /*# CONTOURPLOT */

  if (geo) {
    while (fgets(line,128,geo)) if (!strchr("!#",line[0])) {
      char *info=strchr(line,'!');
      if (!info) info=strchr(line,'#');
      if (!info) info="\n";

      n=sscanf(line,"%d%d%d",&i,&j,&k);
      if (i>=0 && i<ns && j>=0 && j<ns) {
        double ij=SQRD(r[i],r[j]);
        if (n==2) prt_("! dist %d-%d = %.11f %s",i,j,sqrt(ij),info);
        if (n==3 && k>=0 && k<ns) {
          double jk=SQRD(r[j],r[k]);
          prt_("! angle %d-%d-%d = %.9f %s",
              i,j,k,
              acos((ij+jk-SQRD(r[i],r[k]))/(2*sqrt(ij*jk)))*(180/PI),
              info); } } }
    fclose(geo); }

  *spec->ext=0;
  free(f);
  free(r);
  
  return U0;
}

/* #include "dumpnextppm.c" */

static int autogr; /* graphics on/off control switch (cumbersome) */
static int sdcg;   /* max(sd,cg) */

static void show(double U,species_t *spec,int gr) /******************* show */
/***
    3D picture of molecule *spec.
    Based on the HTH project (in Pascal), (c) J.Kolafa 1990

    gr=0: as in last step
    gr=1: wire mode: show bonds only
    gr=2: balls-and-sticks model, draw bonds and atom-spheres
          ordered by z coordinate of their centres
    gr=3: as above and shaded spheres
    gr=4: stereo (anaglyph)
    gr=-1: during minimization by steepest descent (gr determines color)
    gr=-2: during minimization by conjugate gradients or MC (gr determines color)
***/
{
  bond_t *b;
  int i,j,ii=0,jj,sphere;
  static int oldstrlen=20,Ucolor;
  double minz,z,dx,dy;
  static species_t *oldspec;
  static double scalex,scaley;
  static float rscale; /* of RvdW = radii */
  static int erase=0, grid=0;
  static double zzeye=6; /* 3x screen height: static now */
  static double rpert=0.3;
  double
    dd=PI/32, /* move step */
    mouserot;
  vector rot;
  char *c,Us[24];

  ns=spec->ns;
  lasti=-1;
  if (gr>0) grmode=gr;

  /* ? is count used later ? */
  loop (i,0,ns) site[i].count=0;

  if (!autogr) {
    if (spec!=oldspec) {
      /* next species */
      scaley=sqrt(ns)+0.5; /* empirical formula for automatic scale */
      oldspec=spec;
      rscale=0.25;
      bondcol=0; }

    selectfont(font);

    if (display) {
      /* X display already on */
      setmenusize();
      setviewport(0,0,maxxn-1,maxyn-1,1);
      clearviewport();
      setviewport(0,0,maxxn-1,maxyn-1,0);
      redrawmenu=1; }
    else {
      /* start graphics */
      if (startgraph(-9)<0) exit(-1);
      maxxn=getmaxx()+1;
      maxyn=getmaxy()+1;
      if (maxyn<480) font=0; else font=2;
      selectfont(font);
      menu=getenv("GUI")==NULL || !!strchr(getenv("GUI"),'b');
      setmenusize(); }
  } else {
    /* called again (from minimization) */
    switch (erase) {
      /* case 2: just draw over */
      case 1:
        /* erase previous image */
        setviewport(0,0,maxxn-1,maxyn-1,1);
        clearviewport();
        setviewport(0,0,maxxn-1,maxyn-1,0);
        break;
      case 0: {
        /* random shading */
        int i,n=(long)maxxn*maxyn/60L;

        loop (i,0,n) putpixel(irnd(maxxn),irnd(maxyn),BLACK); }
    } }

  /*  dumpnextppm(); */

  mouserot=2*PI/(maxxn+maxyn);

  selectfont(font);
  setcolor(LIGHTCYAN);
  c="showing - [F1] for help";
  if (gr<0) {
    Ucolor=13-gr;
    setcolor(Ucolor);
    c=string("minimizing(%d) - [ESC] to interrupt",sdcg);
    autogr=1; }

  settextstyle(0,0,1);
  settextjustify(2,2);
  outtextxy(maxxn-3,2,c);

  for (;;) {

    if (redrawmenu) {
      makemenu(menu);
      redrawmenu=0; }

    setviewport(0,0,maxxn-1,maxyn-1,1);

    if (pickinfo[0]==0) {
      setcolor(Ucolor);
      settextjustify(0,2);
      sprintf(Us,"U=%.12g",U);
      if (erase!=1) {
        setfillstyle(0/*EMPTY_FILL*/,BLACK);
        bar(0,0,xfont.width*(1+oldstrlen),xfont.height); }
      outtextxy(1,1,Us);
      oldstrlen=strlen(Us); }

    scalex=scaley*maxxn/maxyn;
    lwindow(0,maxxn,0,maxyn); /* within getmax */
    scale(-scalex,scalex,-scaley,scaley,0); /* in fact scale=1 would be more systematic... */
    zeye=zzeye*scaley; /* viewing point zeye*screen height from the screen */
    eye=scaley*eyedist;

    if (ch!=12) pickinfo[0]=0;
    found=0;

    if (grmode==1) {
      /* show bonds directly */
      setlinestyle(0,0,1);
      for (b=b0; b; b=b->next) {
        i=b->indx[0]; j=b->indx[1];
/*was:lline(site[i].r[0],site[i].r[1],site[j].r[0],site[j].r[1]); */
/*.....      setcolor(site[i].count && site[j].count ? LIGHTRED : bondcol);*/

        up();
        if (!bondcol) setcolor(atomcolor(site[i].type));
        else setcolor(bondcols[bondcol]);
        draw3D(site[i].r[0],site[i].r[1],site[i].r[2]);
        if (!bondcol) {
          draw3D((site[j].r[0]+site[i].r[0])/2,(site[j].r[1]+site[i].r[1])/2,(site[j].r[2]+site[i].r[2])/2);
          setcolor(atomcolor(site[j].type)); }
        draw3D(site[j].r[0],site[j].r[1],site[j].r[2]); }

      loop (i,0,ns) {
        showlabel(i);
        if (site[i].count || site[i].keep || site[i].nnbr==0) {
          float qq=zeye/(zeye-site[i].r[2]);
          int ix=SX(site[i].r[0]*qq),iy=SY(site[i].r[1]*qq);

          setlinestyle(0,0,1);
          if (site[i].nnbr==0) {
            // setcolor(WHITE);
            setcolor(atomcolor(site[i].type));
            line(ix-3,iy-3,ix+4,iy+4);
            line(ix-3,iy+3,ix+4,iy-4); }
          if (site[i].keep) {
            setcolor(i==center?LIGHTGREEN:LIGHTRED);
            circle(ix,iy,3); }
          setlinestyle(0,0,3);
          if (site[i].count) {
            setcolor(i==center?LIGHTGREEN:LIGHTRED);
            circle(ix,iy,3); } }
        setlinestyle(0,0,1); } }

    else if (grmode==4) {
      /* stereo (anaglyph): no circle over atom picked up */
      int stereo;

      setwritemode(1);
      setlinestyle(0,0,1);
      for (stereo=-1; stereo<2; stereo+=2) {
        setcolor(stereo<0 ? LIGHTRED : LIGHTCYAN);
        for (b=b0; b; b=b->next) {
          i=b->indx[0]; j=b->indx[1];
          up();
          draw3Ds(site[i].r[0],site[i].r[1],site[i].r[2],stereo);
          draw3Ds(site[j].r[0],site[j].r[1],site[j].r[2],stereo); }

        loop (i,0,ns) if (site[i].nnbr==0) {
          up();
          draw3Ds(site[i].r[0]+.1,site[i].r[1],site[i].r[2],stereo);
          draw3Ds(site[i].r[0]-.1,site[i].r[1],site[i].r[2],stereo);
          up();
          draw3Ds(site[i].r[0],site[i].r[1]+.1,site[i].r[2],stereo);
          draw3Ds(site[i].r[0],site[i].r[1]-.1,site[i].r[2],stereo);
          up();
          draw3Ds(site[i].r[0],site[i].r[1],site[i].r[2]+.1,stereo);
          draw3Ds(site[i].r[0],site[i].r[1],site[i].r[2]-.1,stereo); } }

      setwritemode(0); }

    else { /* grmode=2,3 */
      int *bonddrawn,*sitedrawn;
      int ic,nc=0;
      bond_t *bdraw=NULL;

      double maxz=-9e9;

      settextstyle(0,0,1);
      settextjustify(1,1);

      for (b=b0; b; b=b->next) nc++;
      alloczero(bonddrawn,(nc+(nc==0))*sizeof(bonddrawn[0]));
      alloczero(sitedrawn,ns*sizeof(sitedrawn[0]));

      for (;;) {

        minz=9e9;
        sphere=-1;

        /* bonds */
        for (ic=0,b=b0; b; b=b->next,ic++) if (!bonddrawn[ic]) {
          i=b->indx[0]; j=b->indx[1];
          z=site[i].r[2]+site[j].r[2];
          if (z<=minz) {
            bdraw=b; ii=ic; minz=z;
            sphere=0; } }

        /* atoms */
        loop (i,0,ns) if (!sitedrawn[i]) {
          z=site[i].r[2]*2;
          if (z<=minz) {
            ii=i; minz=z;
            sphere=1; } }

        /* Angstrom grid */
        if (grid) if (maxz*minz<0) {
          int i,s;

          s=(int)(scalex*grid);

          setlinestyle(0,0,1);
          loopto (i,-s,s) {
            setcolor(i%10 ? i%5 ? DARKGRAY : LIGHTGRAY : WHITE);
            up(); draw3D(-scalex,(double)i/grid,0); draw3D(scalex,(double)i/grid,0);
            up(); draw3D((double)i/grid,-scalex,0); draw3D((double)i/grid,scalex,0); } }

        if (sphere<0) break; /* everything drawn */

        maxz=minz;
        if (maxz>=zeye*0.95) {
          fputc('\a',stderr); zzeye=(zeye=maxz*1.1)/scaley; break; }

        if (sphere) {
          float qq=zeye/(zeye-site[ii].r[2]);
          /* this is aproximate - should be ellipses ! */
          float q=atom[site[ii].type].LJ[0].RvdW*rscale*qq;
          unsigned rx=q*scaling.x,ry=q*scaling.y;
          unsigned i,n;
          int cx=SX(site[ii].r[0]*qq),cy=SY(site[ii].r[1]*qq); /* not exact! */
          float hcx=cx+0.5,hcy=cy+0.5,co=-0.7373689,si=-0.6754903,x,y,xnew;

          sitedrawn[ii]++;
          setlinestyle(0,0,1);
          col=atomcolor(site[ii].type);
          setfillstyle(1 /*SOLID_FILL*/,col);
          if (site[ii].count) {
            setcolor(col==LIGHTRED?LIGHTGREEN:LIGHTRED);
            if (ii==center) setlinestyle(0,0,3); }
          else
            setcolor(darkcolor(col));
          fillellipse(cx,cy,rx,ry);

          col=darkcolor(col);

          if (label==NONE && cross) {
            int crx=cross+1,cry=cross+1;
            if (site[ii].count) crx=rx,cry=ry;
            if (ii==center) setlinestyle(0,0,3);
            setcolor(col);
            line(cx-crx+1,cy,cx+crx,cy);
            line(cx,cy-cry+1,cx,cy+cry); }

          if (grmode==3) {
            /* shaded spheres */
            n=rx*ry;
            x=rx-0.5; y=0;
            for (i=n; i>2; i--) {
              xnew=x*co-y*si; y=y*co+x*si; x=xnew;
              q=1-0.5/i/i;
              si*=q; co*=q;
              putpixel(hcx+x,hcy+y,col); } }

          showlabel(ii); }

        else {
          float ri[3],rj[3],dr[3],r=0,qi,qj;

          bonddrawn[ii]++;

          ii=bdraw->indx[0];
          jj=bdraw->indx[1];

          loop (j,0,3) {
            dr[j] = (rj[j]=site[jj].r[j]) - (ri[j]=site[ii].r[j]) ;
            r+=Sqr(dr[j]); }

          r=rscale/sqrt(r);
          qi=atom[site[ii].type].LJ[0].RvdW*r;
          qj=atom[site[jj].type].LJ[0].RvdW*r;
          if (qi+qj>1) {
            fputc('\a',stderr); rscale/=qi+qj+0.02; }

          setlinestyle(0,0,3);
          setcolor(DARKGRAY);
          up();
          draw3D(ri[0]+=qi*dr[0],ri[1]+=qi*dr[1],ri[2]+=qi*dr[2]);
          draw3D(rj[0]-=qj*dr[0],rj[1]-=qj*dr[1],rj[2]-=qj*dr[2]);

          setlinestyle(0,0,1);
          if (!bondcol) setcolor(atomcolor(site[ii].type));
          else setcolor(bondcols[bondcol]);
/*.....        if (site[ii].count && site[jj].count) { setcolor(LIGHTRED); setlinestyle(0,0,3); }*/
          up();
          draw3D(ri[0],ri[1],ri[2]);
          if (!bondcol) {
            draw3D((rj[0]+ri[0])/2,(rj[1]+ri[1])/2,(rj[2]+ri[2])/2);
            setcolor(atomcolor(site[jj].type)); }
          draw3D(rj[0],rj[1],rj[2]); }

        } /* for(;;) */

      free(sitedrawn);
      free(bonddrawn); } /* gr>1 */

    if (marklabeled) {
      fprintf(stderr,"%d atom%s found and marked\n",found,"s"+(found==1));
      strcpy(pickinfo,string("%d marked",found)); }

    if (pickinfo[0]) {
      //      int n=maxxn/xfont.width;
      setcolor(WHITE);
      setfillstyle(1,BLACK);
      bar(0,0,maxxn-xfont.width-4,xfont.height);
      settextstyle(0,0,1);
      settextjustify(0,2);
      outtextxy(1,0,pickinfo); }

    if (autogr)
      /* this happens during minimization */
      if (!kbhit()) goto ret;

  newkey:
    dx=dy=0;
    find.x=-1;

    if (moving) {
      char s[48];
      int sum=0;

      loop (i,0,ns) if (site[i].count) sum++;

      moving=1;
      if (center>=0) {
        sprintf(s,"moving %d atom%s, center=%d",sum,"s"+(sum==1),center);
        if (sum>1) moving=2; }
      else
        sprintf(s,"moving %d atom%s, no center",sum,"s"+(sum==1));

      setfillstyle(1,BLACK);
      setwritemode(0);
      settextstyle(0,0,1);
      settextjustify(0,2);
      bar(0,getmaxy()-xfont.height,xfont.height*strlen(s),getmaxy());
      setcolor(YELLOW);
      outtextxy(1,getmaxy()-xfont.height+1,s); }

    if (redraw) {
      /* force redraw (after click will draw marking) */
      redraw=0;
      ch=12; }
    else {
      setviewport(0,0,maxxn-1,maxyn-1,0);
      ch=readkey();
    }

    memset(rot,0,sizeof(rot));
    marklabeled=0;

    rpert=rpert*0.8+0.06;

    switch (ch) {
    case X_REDRAW_ME:
      lwindow(0,maxxn-1,0,maxyn-1);
      redrawmenu=1;
      if (kbhit()) goto newkey;
      break;
    case LEFTDRAG:
      if (dosmouse.dy) rot[0]=-dosmouse.dy*mouserot;
      if (dosmouse.dx) rot[1]=-dosmouse.dx*mouserot;
      break;

    case MIDDRAG:
      if (dosmouse.dx) dx=(double)dosmouse.dx*scalex*2/maxxn;
      if (dosmouse.dy) dy=-(double)dosmouse.dy*scaley*2/maxyn;
      break;

    case RIGHTDRAG:
      //      if (dosmouse.dy) scaley*=exp((double)dosmouse.dy/getmaxy());
      //      if (dosmouse.dx) rot[2]=-dosmouse.dx*mouserot;
      if (dosmouse.dx || dosmouse.dy) {
        double rx=dosmouse.x-maxxn/2.;
        double ry=dosmouse.y-maxyn/2.;
        double a=Sqr(rx)+Sqr(ry);
        if (a>30) /* opposite sign than in show (why?) */
          rot[2]=(dosmouse.dy*rx-dosmouse.dx*ry)/a; }
      break;

    /* single click (no move) by any button */
    case LEFTCLICK:
    case RIGHTCLICK:
    case MIDCLICK:
      find.x=dosmouse.x; find.y=dosmouse.y;
    case 12: /* force redraw (redraw=1: see above) */
      break;

    case '@':
      switch (readkey()) {
        case 'k': writeMARKED(spec,1); break;
        case 'm': writeMARKED(spec,0); break;
        case 'K': readMARKED(spec,1); break;
        case 'M': readMARKED(spec,0); break; }
      break;

    case 'P'&31:
    case 'P': {
      static char fn[16]="blend000",fn0[16];
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

      makemenu(2);
      redrawmenu=1;

      fprintf(stderr,"select one of: pPoOcCeE\n");
      ch=readkey();

      if (!strchr("pPoOcCeE",ch)) break;
      wbg=islower(ch);
      ch=tolower(ch);

      strcpy(fn+8,".ppm");

      fprintf(stderr,"writing %s\n",fn);
      f=fopen(fn,"wb");
      if (wbg) {
        col[0][0]=col[0][1]=col[0][2]=255;
        col[15][0]=col[15][1]=col[15][2]=240; }

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
        fprintf(stderr,"ERROR: Conversion to %s failed\nCheck whether \'ppm2ps\' is installed\n",fn);
      if (ch!='p') unlink(fn0);

      if (ximage) XDestroyImage(ximage);
      sprintf(fn+5,"%03d",atoi(fn+5)+1);
      break; }

    case 'e': {
      static char eraseinfo[16]="fade\0none\0over";
      erase=(erase+1)%3;
      if (menu) {
        setfillstyle(1,LIGHTGRAY);
        bar(getmaxx()-4*xfont.width-1,aty_n-1,getmaxx(),aty_n+xfont.height);
        setcolor(RED);
        settextstyle(0,0,1);
        settextjustify(2,2);
        outtextxy(getmaxx()-2,aty_n,eraseinfo+erase*5); }
      else
        fprintf(stderr,"erase mode=%s\n",eraseinfo+erase*5); }
      break;

    case 'j':
      by*=2;
      if (by>sdcg) by=3;
      if (menu) {
        setfillstyle(1,LIGHTGRAY);
        bar(getmaxx()-4*xfont.width-1,aty_n-1,getmaxx(),aty_n+xfont.height);
        setcolor(BLUE);
        settextstyle(0,0,1);
        settextjustify(2,2);
        outtextxy(getmaxx()-2,aty_n,string("%d",by)); }
      else
        fprintf(stderr,"show by %d minimization steps\n",by);
      break;

    case RIGHT: dx=1.618034/scalex; break;
    case LEFT: dx=-1/scalex; break;
    case UP: dy=1.618034/scalex; break;
    case DOWN: dy=-1/scalex; break;
    case 'A'&31:
    case HOME: zzeye/=exp(dd); break;
    case 'E'&31:
    case END: zzeye*=exp(dd); break;
    case 'x': rot[0]=-dd; break;
    case 'X': rot[0]=dd; break;
    case 'y': rot[1]=-dd; break;
    case 'Y': rot[1]=dd; break;
    case 'z': rot[2]=-dd; break;
    case 'Z': rot[2]=dd; break;
    case '>': scaley/=2; break;
    case WHEELFORWARD:
    case '+': scaley*=0.86; break;
    case '<': scaley*=2; break;
    case WHEELBACKWARD:
    case '-': scaley*=1.1; break;
    case '*': if (dd<1) dd*=2; goto newkey;
    case '/': if (dd>1e-4) dd*=0.5; goto newkey;
    case 'b': bondcol+=2;
    case 'B': bondcol=(bondcol-1)%NBONDCOL;
              redrawmenu=1;
              break;
    case '#':
    case '=': if ((grmode|1)==3)
                grid="\1\12\0\0\0\0\0\0\0\0"[grid]; /* 0->1->10->0 */
              else
                fprintf(stderr,"grid shown in graphic modes 2 and 3 only\n");
              break;
    case '\'': cross=(cross+2)&3; break;
    case 'r': rscale*=0.9; break;
    case 'R': rscale*=1.2; break;
    case 'g': grmode=(grmode&3)+1; break;
    case 'G': grmode=(grmode+3)&3; break;
    case 'c': if (colors[7]=='C') colors[7]=' ',colors[8]='C';
              else if (colors[8]=='C') colors[8]=' ',colors[11]='C';
              else colors[7]='C';
              redrawmenu=1; break;
    case 'm': moving=!moving;
              redrawmenu=1; break;
      //    case 'M': if (center>=0) floodmark(center); break;
    case 'L': marklabeled=1; redraw=1; break;
    case ' ': label=NONE; break;
    case 'i': label=ID; break;
    case 't': label=TYPE; break;
    case 'n': label=NUMBER; break;
    case 'q': label=CHARGE; break;
#ifdef POLAR
    case 'a': label=POLARIZABILITY; break;
    case 'D': label=DIPOLE_eAA; break;
    case 'd': label=DIPOLE_D; break;
    case 's': label=SATURATION; break;
#endif /*# POLAR */
    case 'F': find.s[0]=0; find.i=-1; find.q=1e3;
    case 'u'&31:
    case 'u': loop (i,0,ns) site[i].count=0;
              center=-1;
              break;
    case 'I': loop (i,0,ns) site[i].count=!site[i].count;
              center=-1;
              break;
    case 'K': loop (i,0,ns) site[i].count=site[i].keep;
              break;
    case 'k': loop (i,0,ns) site[i].keep=(enum keep_e)(!!site[i].count);
              {
                int s=0;
                loop (i,0,ns) s+=site[i].keep;
                prt("\
%d marked atoms will be kept fixed and %d left free in minimization", s,ns-s);
                }
              if (!spec->opt_k) {
                fprintf(stderr,"keep=1 set (the same as option -k1)\n");
                spec->opt_k=1; }
              break;
              //    case 'l': label=(enum label_e)((label+1)%5); break;
    /* very ugly...*/
    case 'E': U=hotkeyE(spec,1,0);
              break;
    case 'f':
      if (label>=ID && label<=CHARGE) {
        closegraph();
        switch (label) {
          case ID:
            prts_("enter atom id");
            goto sss;
          case TYPE:
            prts_("enter atom type");
          sss:
            prts(" to find, may end by * :");
            if (1!=scanf("%15s",find.s)) break;
            find.l=strlen(find.s);
            if (find.s[find.l-1]=='*') find.l--; else find.l++;
            break;
          case NUMBER:
            prts("enter atom # to find :");
            if (1!=scanf("%d",&find.i)) break;
            break;
          case CHARGE:
            prts("enter charge to find and precision:");
            if (1!=scanf("%lf%lf",&find.q,&find.dq)) break;
            find.dq+=1e-14;
            break;
          default:
            /* suppress warnings */; }
        if (startgraph(-9)<0) exit(-1);
        if (menu) makemenu(1); }
      else {
        fprintf(stderr,"don\'t know what to find: use one of ID, type, charge, number\n");
        setcolor(LIGHTRED);
        settextjustify(2,0);
        outtextxy(getmaxx(),getmaxy(),"what to find?"); }
      break;
    case 'W': forget(0); break;
    case 'w': forget(1); break;
    case 'A': averagecharges(spec,&center); break;
    case '1': case '2': case '3': case '4':
      grmode=ch-'0';
      redrawmenu=1;
      break;
    case 'Q':
    case 'Q'&31:
      closegraph();
      exit(0);
    case 'S'&31:
      if (spec->wrmol) writeMOL(spec);
      spec->wrmol=0;
      if (!spec->opt_p && !spec->opt_w) WARNING(("neither -p nor -w => nothing saved"))
      if (spec->opt_p) write3D(spec,3);
      if (spec->opt_w) write3D(spec,0);
      break;

    case ':': perturb(spec,NULL,0,rpert*=1.4);
    case ';':
    case ',':
    case '.': /* finish */
    case ESC: /* only interrupt minimization */
      moving=0;
      if (ch==ESC && !autogr) {
        if (readkey()==ESC) { closegraph(); exit(0); }
        break; }
      if (ch=='.' && !autogr) closegraph();
      sig=autogr;
      goto ret;

    case F12:
      if (readkey()==F12) { closegraph(); exit(0); }
      break;

    case F5:
    case 'O'&31:
      font=(font+2)%4;
      goto F10font;

    case F10:
      menu=!menu;
      erasebuttons();
    F10font:
      selectfont(font);
      redrawmenu=1;
      redraw=1;
      setmenusize();
      break;

    case '?': case F1: case 'h': case 'H':
      maxxn=getmaxx()+1;
      setviewport(0,0,maxxn-1,maxyn-1,1);
      help();
      maxxn=getmaxx()+1;
      if (menu) maxxn-=10*xfont.width+86;
      maxyn=getmaxy()+1;
      redraw=1;
      redrawmenu=1;
      break;

      case ']':
        if (!sdcg) sdcg=1;
        else sdcg*=2;
        if (sdcg>9999) sdcg=9999;
        goto draw_n;
        break;
      case '[':
        if (sdcg>1) sdcg/=2;
      draw_n:
        if (menu) {
          setfillstyle(1,LIGHTGRAY);
          bar(getmaxx()-4*xfont.width-1,aty_n-1,getmaxx(),aty_n+xfont.height);
          setcolor(BLACK);
          settextstyle(0,0,1);
          settextjustify(2,2);
          outtextxy(getmaxx()-2,aty_n,string("%d",sdcg)); }
        else
          fprintf(stderr,"%d minimization steps\n",sdcg);
        break;

    default:
      goto newkey;
    } /* switch (ch) */

    if (dx!=0 || dy!=0) {
      if (moving) {
        loop (i,0,ns)
         if (site[i].count) {
           site[i].r[0]+=dx;
           site[i].r[1]+=dy; } }
      else {
        remember(-1,dx,dy);
      loop (i,0,ns) {
        site[i].r[0]+=dx;
        site[i].r[1]+=dy; }
      lasti=-1; } }

    if (moving && center>=0) {
      loop (j,0,3) if (rot[j]!=0) {
        dx=cos(rot[j]); dy=sin(rot[j]);
        ii=(j+1)%3; jj=(j+2)%3;
        loop (i,0,ns) if (site[i].count) {
          double zz=site[i].r[jj]-site[center].r[jj];
          z=site[i].r[ii]-site[center].r[ii];
          site[i].r[ii] = site[center].r[ii] + dx*z + dy*zz;
          site[i].r[jj] = site[center].r[jj] - dy*z + dx*zz; }
        } }
    else loop (j,0,3) if (rot[j]!=0) {
      dx=cos(rot[j]); dy=sin(rot[j]);
      ii=(j+1)%3; jj=(j+2)%3;
      remember(ii,dx,dy);
      loop (i,0,ns) {
        z=site[i].r[ii];
        site[i].r[ii] =  dx*z + dy*site[i].r[jj];
        site[i].r[jj] = -dy*z + dx*site[i].r[jj]; }
      lasti=-1; }

      setviewport(0,0,maxxn-1,maxyn-1,1);
      clearviewport();
      if (!menu) makemenu(0);
      setviewport(0,0,maxxn-1,maxyn-1,0); }

 ret:
  forget(0);
} /* show */

static int nUpot;
static void times(int key) /****************************************** times */
{
  static time_t t0,t1;
  if (key) {
    time(&t1);
    prt("! duration=%us/%d=%gs",
        (unsigned int)(t1-t0), nUpot, (double)(t1-t0)/nUpot);
    t0=t1; }
  else
    time(&t0);
  nUpot=0;
}

#include "blendpw.c"

#ifdef POLAR
#  define NRATES 29
static double polerr=0,polerr0,rate[NRATES];
static int polnit=0;
#endif /*# POLAR */

static struct sym_s {
  struct sym_s *next;
  int axis;
  int i,j;
} *sym0;

static void makesym(vector *r) /************************************ makesym */
{
  struct sym_s *sym;

  for (sym=sym0; sym; sym=sym->next)
    if (sym->i==sym->j) r[sym->i][sym->axis]=0;
    else {
      double x=(r[sym->i][sym->axis]-r[sym->j][sym->axis])/2;
      VVVV(r[sym->i],=r[sym->j],=0.5*r[sym->i],+0.5*r[sym->j])
      r[sym->i][sym->axis]=x;
      r[sym->j][sym->axis]=-x; }
}

#include "blendu.c"

double minimize(species_t *spec,int print,enum keep_e mask) /***** minimize */
/***
    Random initial configuration (if option -r)
    Minimize energy for molecule spec.
    Greedy (steepest descent) or conjugate gradients method used
    If compiled with OMITKEPT:
    mask=FREE: all site-site interactions included
               (even if there are non-movable atoms INJAIL)
    mask=INJAIL: non-movable sites omitted
***/
{
  int ns=spec->ns,size=ns*sizeof(vector),i,sd,cg;
  unsigned it;
  vector *r,*f,*f0,*r0,*G,*H;
  double h=1e-6,U,U0;
  int shown;

  if (print) times(0);

  spec->opt_m *= ns>1; /* no minimization for 1 atom */

  sdcg=abs(spec->opt_m);

  site=spec->site;
#ifdef OMITKEPT
  keepmask=mask;
#endif /*# OMITKEPT */

  /* enable ^C */
  sig=0;
  if (spec->probe.ns==0) {
      (*signal)(SIGINT,asksig);
      fprintf(stderr,"%s energy of %s: interrupt by ^C\n",
              sdcg ? "minimizing" : "calculating",
              spec->fn); }

  /* how often print U */
  if (option('g')>0) by=sdcg/(option('g')+1);
  else if (option('g')<0) by=abs(option('g'));
  else by = sdcg>440 ? 20 : sdcg/22+1;
  if (by==0) by=1;

  ralloc(r0,size);
  ralloc(f0,size);
  {
    FILE *f;
    char line[128];

    strcpy(spec->ext,".sym");
    f=fopen(spec->fn,"rt");

    if (f) {
      while (fgets(line,128,f)) if (!strchr("!#",line[0])) {
        unsigned char axis;
        int i,j;
        int n=sscanf(line,"%c%d%d",&axis,&i,&j);
        struct sym_s *sym;
        char *d=strchr(line,'-');

        if (axis=='.' || axis=='q') break;
        if (d) j=-atoi(d);
        if (n<2) ERROR(("%s: %s : bad line",spec->fn,line))
        if (n==2) j=i;
        if (axis>='x') axis-='x';
        else if (axis>='X') axis-='X';
        if (axis>2)
          ERROR(("%s: %s : bad axis",spec->fn,line))
        if (i<0 || i>=ns || j<0 || j>=ns)
          ERROR(("%s: %s : bad site",spec->fn,line))
        do {
          ralloc(sym,sizeof(struct sym_s));
          sym->next=sym0;
          sym0=sym;
          sym->axis=axis;
          sym->i=i;
          if (d) sym->j=i;
          else sym->j=j; }
        while (d && ++i<=j); }
      prt("! symmetry on");
      fclose(f); }
    *spec->ext=0;
  }

  /*** initial configuration from the molecule ***/
  loop (i,0,ns) VV(r0[i],=site[i].r)
  makesym(r0);

#ifdef BLENDDEBUG
  prt("#### debug! site i<%d (-1=break) coordinate j={0,1,2}",spec->ns);
  for (;;) {
    static int i=0,j=0;
    static double h=1e-5;
    double e;

    for (;;) {
      getdata get(i) get(j) get(h) checkdata enddata
      if (i<0) goto debugdone;
      if (i<ns && j>=0 && j<3) break;
      prt("bad i or j"); }
    r0[i][j]-=h;
    U0=Upot(r0,f0,spec,0);
    r0[i][j]+=h+h;
    U=Upot(r0,f0,spec,0);
    r0[i][j]-=h;
    Upot(r0,f0,spec,i+1);
    e=f0[i][j]+(U-U0)/h/2;
    prt("i=%d j=%d h=%.1e f=%g -dU/dr=%g err=%g %g",
         i,   j,   h,    f0[i][j], (U0-U)/h/2, e, e/obfuscate);
  }
 debugdone:;
#endif /*# BLENDDEBUG */

  if (spec->probe.ns==3) {
    vector HH[2];
    int i,j;

    /* OBSOLETE: special patch for probing by TIP3P water */
    i=spec->probe.i;
    U0=3e33;
    loop (j,0,12) {
      makeHH(&r0[i],j);
      U=Upot(r0,f0,spec,0);
      if (U<U0) {
        U0=U;
        copy(HH,&r0[i+1],sizeof(HH)); } }
    copy(&r0[i+1],HH,sizeof(HH));
    loop (i,0,ns) VV(site[i].r,=r0[i]) }
  else

    U0=Upot(r0,f0,spec,print);

  if (U0>1e33) goto finish;
  if (print) times(1);

  if (out!=stdout) if (print) fprintf(stderr,"U = %.12g\n",U0);

  if (spec->opt_m) {
    /* minimization requested */
    constraints=option('a')==0 || option('b')==0 || fd0;
    if (constraints) {
      /* negative option -m = SD (with constraints) */
      if (spec->opt_m<0) sd=sdcg,cg=0;
      else cg=sdcg,sd=0; }
    else {
      cg=sdcg;
      /* negative option -m = SD + CG */
      if (spec->opt_m<0) sd=cg/4+2;
      else sd=0; } }
  else {
    /* no minimization */
    if (option('g')) {
      show(U0,spec,1);
      loop (i,0,ns) VV(r0[i],=site[i].r) }
    goto finish; }

  ralloc(r,size); copy(r,r0,size);
  ralloc(f,size);
  ralloc(G,size); /* conjugate grad only */
  ralloc(H,size); /* conjugate grad only */

  for (;;) {

    h*=3;
    U0=Upot(r0,f0,spec,0); /* needed because show() may move atoms */

    /* re-consider step to prevent large moves when overlap etc. */
    U=0;
    loop (i,0,ns) Max(U,SQR(f0[i]))
    U=sqrt(U); /* =max force on 1 atom */
    if (U*h>0.5) h=0.5/U; /* max move in the 1st step = 0.5 AA */
    if (U*h<0.001) h=0.001/U; /* min move in the 1st step = 0.001 AA */

    // put3(sd,cg,h)

//  #include "blendmd.c" (not finished)

    shown=0;

    if (sd) {
      /*** steepest descent method, both with and without constraints ***/
      if (print) prt("! %d steps of steepest descent%s:",sd,constraints?" (with constraints)":"");

      for (it=1; !sig; it++) {

        do {
          h*=0.55; /* magic value: works fine in a narrow valley */

          /* converg. criterion changed from 1e-12 to 1e-11 in V 2.0o */
          if (h<1e-11 || it>=sd) goto sddone;

          loop (i,0,ns) VVV(r[i],=r0[i],+h*f0[i])
          U=Upot(r,f,spec,0);
        } while (U>U0);

        copy(r0,r,size);
        copy(f0,f,size);
        U0=U;

        if (option('g')) {
          if (it && it%by==0) {
            loop (i,0,ns) VV(site[i].r,=r0[i])
            show(U0,spec,-1);
            shown++;
            /* does NOT back copy rotated r here */ } }
        else
          if (print) prtU(it,by,U,h);

        h*=2.5; /* magic: h should not shrink to zero */ } }
    put(h)
sddone:

    if (cg) {
      if (constraints) {
        /*** constraints present: MC minimization ***
             inefficient patch: constraint forces not needed but calculated */
        double qacc=1+2./ns,qrej=1-0.2/ns;
        if (print) prt("! %d steps of MC minimization (with constraints):",sd+cg);
        h*=4;

        //put2(qacc,qrej)

        if (cg) qrej=1-0.05/ns;

        for (it=1; !sig; it++) {

          if (h<1e-6 || it>=cg*2) goto cgdone;

          loop (i,0,ns) VVO(r[i],=r0[i],+h*rndgauss())

          U=Upot(r,f,spec,0);

          if (U<U0) {
            copy(r0,r,size);
            copy(f0,f,size);
            U0=U;
            h*=qacc; }
          else
            h*=qrej;

          if (option('g')) {
            if (it && it%by==0) {
              loop (i,0,ns) VV(site[i].r,=r0[i])
              show(U0,spec,-2);
              shown++;
              /* does NOT back copy rotated r here */ } }
          else
            if (print) prtU(it,by,U,h);
        }
      } else {
        /*** conjugate gradients ***
             rewritten from Numerical Recipes
             no constraints */
        double dgg,gg,h0;
        /* brrrr! U0=fp r0=p  f0=xi */
        int j,maxj=24;

        if (print) prt("! %d steps of conjugate gradients:",cg);
        /* from init, known r0 f0 U0 */
        h*=3;
        copy(G,f0,size);
        copy(H,f0,size);

        for (it=1; !sig && it<=cg; it++) {
          /***
              now we are looking for h>0 so that r=r0+h*f0 is minimum
              in most cases close to minimum we need only two
              evaluations of forces in one step
          ***/
          gg=0;
          loop (i,0,ns) gg+=SQR(f0[i]);
          if (gg==0) break;

          h0=h;
          for (j=0; ;j++) {
            loop (i,0,ns) VVV(r[i],=r0[i],+h*f0[i])
            U=Upot(r,f,spec,0);
            if (U<U0) break;
            TRACE('/')
            h*=0.5;
            if (j>=maxj) {
              /* unsuccessful - reset CG (should occur rarely!) */
              U0=Upot(r0,f0,spec,0);
              if (maxj==20) goto cgdone; /* unsuccessful again! */
              h=h0;
              maxj=20;
              copy(G,f0,size);
              copy(H,f0,size);
              TRACE('r')
              goto OK; } }

          maxj=5;

          dgg=0;
          loop (i,0,ns) dgg+=SCAL(f[i],f0[i]);
          dgg/=gg;
          if (dgg<0.8) {
            /* try one step of the secant method - efficient close to minimum */
            h*=1/(1-dgg);
            loop (i,0,ns) VV(r0[i],+=h*f0[i])
            U0=Upot(r0,f0,spec,0);
            if (U0<U) goto OK; }

          /***
              not better or too extrapolated - we must be happy with r,f,U
              this occurs once a while
          ***/

          h*=4;
          maxj=8; /* might be caused by too low h */
          TRACE('-')
          U0=U;
          copy(f0,f,size);
          copy(r0,r,size);

        OK:
          if (option('g')) {
            if (it && (it-1)%by==0) {
              loop (i,0,ns) VV(site[i].r,=r0[i])
              show(U0,spec,-2);
              shown++;
              /* does NOT back copy rotated r here */ } }
          else
            if (print) prtU(it,by,U,h);

          gg=dgg=0;
          loop (j,0,ns) {
            gg+=SQR(G[j]);
            dgg+=(f0[j][0]-G[j][0])*f0[j][0]
              +(f0[j][1]-G[j][1])*f0[j][1]
              +(f0[j][2]-G[j][2])*f0[j][2]; }
          dgg/=gg;

          copy(G,f0,size);
          loop (j,0,ns) VVV(H[j],=G[j],+dgg*H[j])
          copy(f0,H,size); } } }

  cgdone:
    //put3(h,U0,U)
    U0=Upot(r0,f0,spec,0);
    loop (i,0,ns) VV(site[i].r,=r0[i])
    if (option('g'))
      if (autogr || !shown) {
        autogr=0;
        show(U0,spec,0);
        loop (i,0,ns) VV(r0[i],=site[i].r) }
    if (print) hotkeyE(spec,0,U0);

    if (print) times(1);

    if (strchr("\e.",ch)) break;
    else if (strchr(",;:",ch)) {
      /* continue minimization */
      if (constraints) switch (ch) {
          case ',':
          case ':': cg=sdcg; sd=0; break;
          case ';': sd=sdcg; cg=0; }
      else switch (ch) {
          case ',': cg=sdcg; sd=sdcg/4+2; break;
          case ':': cg=sdcg; sd=0; break;
          case ';': sd=sdcg; cg=0; } }

    if (!option('g')) break;
  }

  /* sometimes unnecessary, but it is difficult to trace all changes in r */
  U0=Upot(r0,f,spec,print);
  if (print) times(1);

finish:
  release(r0);
  sym0=NULL;
/*.....prt("! symmetry off");*/
  return spec->Emin=U0;
}

#include "blendar.c"
#include "blendnm.c"
#include "blendess.c"
#include "blendine.c"
#include "blendvir.c"
