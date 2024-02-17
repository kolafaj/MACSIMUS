/* for plotcalc and plotfile */
#define EXTPSSTRING /* strings allow _Index ^Exponent \Symbol */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "xdraw.h"
#include "mydraw.h"

struct rg_s my_rg;
struct my_style_s my_style = {
  NULL,"Helvetica",14, 
  0, /* enc */
  1, /* minus */
  8*72,6*72,
  50,6,50,6,0,
  1,
  {6,5,4,3,2},
  {5,10,20,50,100},
  {5,10,20,50,100},
  1,0.5,2,
  'E', /* changed to EPS */
  "","", { /* colors are inverted */
    { "1 1 1", "", -1,"fullcircle",-1},
    { "1 1 0.5", "", -1,"fullcircle",-1},
    { "1 0.5 1", "", -1,"fullcircle",-1},
    { "1 0.5 0.5", "", -1,"fullcircle",-1},
    { "0.5 1 1", "", -1,"fullcircle",-1},
    { "0.5 1 0.5", "", -1,"fullcircle",-1},
    { "0.5 0.5 1", "", -1,"fullcircle",-1},
    { "0.3 0.3 0.3", "", -1,"fullcircle",-1},
    { "0.6 0.6 0.6", "", -1,"fullcircle",-1},
    { "1 1 0", "", -1,"fullcircle",-1},
    { "1 0 1", "", -1,"fullcircle",-1},
    { "1 0 0", "", -1,"fullcircle",-1},
    { "0 1 1", "", -1,"fullcircle",-1},
    { "0 1 0", "", -1,"fullcircle",-1},
    { "0 0 1", "", -1,"fullcircle",-1},
    { "0 0 0", "", -1,"fullcircle",-1}}
};

static int my_down;

double mysx(double x) /************************************************ mysx */
{
  return (x-my_rg.x)*(my_style.xsize/(my_rg.X-my_rg.x));
}

double mysy(double y) /************************************************ mysy */
{
  return (y-my_rg.y)*(my_style.ysize/(my_rg.Y-my_rg.y));
}

void myup(void) /****************************************************** myup */
{
  if (my_style.ps) {
    if (my_down) fprintf(my_style.ps," stroke\n");
    my_down=0; }
  up();
}

void mydraw(double x,double y) /************************************* mydraw */
{
  if (my_style.ps) {
    if (!my_down) fprintf(my_style.ps,"newpath ");
    fprintf(my_style.ps,"%.2f %.2f ",mysx(x),mysy(y));
    fprintf(my_style.ps,my_down?"lineto\n":"moveto ");
    my_down=1; }
  draw(x,y);
}

#ifdef EXTPSSTRING

static void changestat(int newstat,int stat)
{
  int shift;

  fprintf(my_style.ps,"/%s findfont %.2f scalefont setfont\n",
	 newstat&4?"Symbol":my_style.font,
	 my_style.fontsize*(newstat&3?3:4)/4);
  /* exponents are 3/2 times raised as indices lowered */
  shift = ((newstat&2)-(stat&2))*3 - ((newstat&1)-(stat&1))*4;
  if (shift)
    fprintf(my_style.ps,"0 %.2f rmoveto ",my_style.fontsize*0.06*shift);
}

static double ps_x,ps_y;
static int ps_col; /* from plot - passed after recoloring */

void mypssetxy(double x,double y,int col) /************************ pssetxy */
/* for $ key only... */
{
  ps_x=x;
  ps_y=y;
  ps_col=col;
  if (col<0) col=0;
}

void psstring(char *s,int align) /******************************** psstring */
/* write extended strings
 elements of string s (@ stands for any character but special):
   @ = normal size (= given by my_style.fontsize)
   $#- = info on line number # (at beginning of string)
   $#* = info on point number # (at beginning of string)
   \@ = Symbol normal size
   _@ = subscript, my_style.font
   ^@ = superscript, my_style.font
   _\@ = subscript, Symbol
   ^\@ = superscript, Symbol
   \+ \- are mathematical + -
 align: -1=left, 0=center, 1=right
*/
{
  int pass;

  if (*s=='$') {
    s++;
    while (*s && isdigit(*s)) s++;
    if (my_style.ps) {
      fprintf(my_style.ps,"gsave\n");
      mysetcolor(ps_col); /* BUG, crash if (by wrong l X Y $) not-number passed */
      my_style.pointsize=2.39; /* HARDWIRED, interacts with mysetcolor - to be fixed! */
      if (*s=='-')
        fprintf(my_style.ps,"newpath %.2f %.2f moveto %.2f %.2f lineto stroke\n",
                ps_x-my_style.fontsize*2.3, ps_y+my_style.fontsize*0.4,
                ps_x-my_style.fontsize*.3, ps_y+my_style.fontsize*0.4);
      else if (*s=='.')
        fprintf(my_style.ps,"newpath %.2f %.2f %.2f S\n",
                ps_x-my_style.fontsize*1.3,
                ps_y+my_style.fontsize*0.4,
                my_style.pointsize);
      else if (*s=='o')
        fprintf(my_style.ps,"newpath %.2f %.2f %.2f opencircle\n",
                ps_x-my_style.fontsize*1.3,
                ps_y+my_style.fontsize*0.4,
                my_style.pointsize);
      else if (*s=='=') {
        fprintf(my_style.ps,"newpath %.2f %.2f  %.2f %.2f errybar\n",
                ps_x-my_style.fontsize*1.3,
                ps_y+my_style.fontsize*0.4-my_style.pointsize*2,
                ps_y+my_style.fontsize*0.4+my_style.pointsize*2,
                my_style.pointsize);
        fprintf(my_style.ps,"newpath %.2f %.2f %.2f S\n",
                ps_x-my_style.fontsize*1.3,
                ps_y+my_style.fontsize*0.4,
                my_style.pointsize); }
      else if (*s=='O') {
        fprintf(my_style.ps,"newpath %.2f %.2f  %.2f %.2f errybar\n",
                ps_x-my_style.fontsize*1.3,
                ps_y+my_style.fontsize*0.4-my_style.pointsize*2,
                ps_y+my_style.fontsize*0.4+my_style.pointsize*2,
                my_style.pointsize);
        fprintf(my_style.ps,"newpath %.2f %.2f %.2f opencircle\n",
                ps_x-my_style.fontsize*1.3,
                ps_y+my_style.fontsize*0.4,
                my_style.pointsize); }
      else {
       fprintf(my_style.ps,"newpath %.2f %.2f moveto %.2f %.2f lineto stroke\n",
          ps_x-my_style.fontsize*2.3, ps_y+my_style.fontsize*0.4,
          ps_x-my_style.fontsize*.3, ps_y+my_style.fontsize*0.4);
        fprintf(my_style.ps,"newpath %.2f %.2f %.2f S\n",
          ps_x-my_style.fontsize*1.3, ps_y+my_style.fontsize*0.4,
          my_style.pointsize); }
      fprintf(my_style.ps,"grestore\n"); }
    s++; }

  for (pass=align<0?1:0; pass<2; pass++) {
    int stat=0; /* 1:index, 2:exponent, 4:symbol */
    char *c;
    char *what=pass?"show":"stringwidth pop";

    fprintf(my_style.ps,"(");

    for (c=s; *c; c++) {
      int newstat=0;
      char ch;

      if (*c=='_') newstat|=1,c++;
      if (*c=='^') newstat|=2,c++;
      if (*c=='\\') newstat|=4,c++;
      if ( !(ch=*c) ) c--,ch=' ';

      if (stat!=newstat) {
        fprintf(my_style.ps,") %s \n",what);
        if (!pass) what="stringwidth pop add";
        changestat(newstat,stat);
        fputc('(',my_style.ps); }

      if (strchr("()%",ch)) fputc('\\',my_style.ps);
      fputc(ch,my_style.ps);

      stat=newstat; }

    fprintf(my_style.ps,") %s \n",what);
    changestat(0,stat);
    if (!pass)
      fprintf(my_style.ps,"%d div 0 rmoveto ",align-2);
  }
}
#endif


void mysetcolor(int col) /*************************************** mysetcolor */
{
  myup();

  setcolor(col);
  if (my_style.ps) {
    float t=my_style.line[col].thick;
    float s=my_style.line[col].size;

    if (t<0) t=my_style.thick;
    if (s>0) my_style.pointsize=s; else my_style.pointsize=pointsize;
    fprintf(my_style.ps,"%s setrgbcolor [%s] 0 setdash %.2f setlinewidth /S { %s } bind def\n",
            my_style.line[col].rgb,my_style.line[col].dash,t,my_style.line[col].symbol); }
}

void mysetlinestyle(int style,int pattern,int thick) /******* mysetlinestyle */
{
  myup();
  setlinestyle(style,pattern,thick);
  if (my_style.ps)
    fprintf(my_style.ps,"[%s] 0 setdash %.2f setlinewidth\n",
            style?"1 3":"",thick*my_style.thick);
}

void mysmartdot(double x,double y) /***************************** mysmartdot */
{
  if (my_style.ps)
    fprintf(my_style.ps,"newpath %.2f %.2f %.2f S\n",
            mysx(x),mysy(y),my_style.pointsize);
  smartdot(x,y);
}

void myerrxbar(double x1,double x2,double y) /******************** myerrxbar */
{
  if (my_style.ps)
    fprintf(my_style.ps,"%.2f %.2f %.2f %.2f errxbar\n",
            mysx(x1),mysx(x2),mysy(y),pointsize);
  up(); draw(x1,y); draw(x2,y);
  point(x1,y,'|');
  point(x2,y,'|');
}

void myerrybar(double x,double y1,double y2) /******************** myerrybar */
{
  if (my_style.ps)
    fprintf(my_style.ps,"%.2f %.2f %.2f %.2f errybar\n",
            mysx(x),mysy(y1),mysy(y2),pointsize);
  up(); draw(x,y1); draw(x,y2);
  point(x,y1,'-');
  point(x,y2,'-');
}

void myaxis(int it,int n,double tick) /****************************** myaxis */
/*
  it = index in my_style.tick[]; labeling for it=0 nly
  n = axis (+1=x, -1=y)
  tick = sign of tick (1=up|right, -1=top|left)
*/
{
  double dz=1e35,d,zz,z,Z;
  double x0=0,y0=0;
  double X0=my_style.xsize,Y0=my_style.ysize;
  char xnr[32];
  char xnr0[32];
  int ex=0;

  tick*=my_style.tick[it];
  if (n>0) {
    Z=my_rg.X; z=my_rg.x;
    if (my_style.xticks[it]>=0) {
      n*=(int)my_style.xticks[it];
      d=(Z-z)/n; }
    else {
      ex++;
      d=-my_style.xticks[it]; }
    if (tick>0) Y0=0;
    else y0=my_style.ysize; }
  else {
    Z=my_rg.Y; z=my_rg.y;
    if (my_style.xticks[it]>=0) {
      n*=(int)my_style.yticks[it];
      d=(z-Z)/n; }
    else {
      ex++;
      d=-my_style.yticks[it]; }
    if (tick>0) X0=0;
    else x0=my_style.xsize; }

  fprintf(my_style.ps,"%.2f setlinewidth newpath %.2f %.2f moveto %.2f %.2f lineto stroke\n",
          my_style.framethick,x0,y0,X0,Y0);

  if (!n) return;

  if (ex) dz=d;
  else for (;;) {
    if (d>dz) break;
    dz*=0.5;
    if (d>dz) break;
    dz*=0.4;
    if (d>dz) break;
    dz*=0.5; }

  zz=((int)(z/dz)-1)*dz;
  while ((zz-z)/dz<-1e-4) zz+=dz;

  xnr0[0]=0;

  while ((zz-Z)/dz<1e-4) {
    double x,y;

    if (fabs(zz/dz)<1e-8) zz=0;
    if (n>0) x=mysx(zz),y=y0;
    else x=x0,y=mysy(zz);

    fprintf(my_style.ps,"newpath %.2f %.2f moveto %.2f %.2f rlineto stroke\n",
            x,y, tick*(n<0),tick*(n>0));
    
    if (tick>0 && it==0) {
      if (n>0) {
#ifdef EXTPSSTRING
        fprintf(my_style.ps,"newpath %.2f %.2f moveto ",x, y-my_style.fontsize*1.2);
        if (my_style.minus) sprintf(xnr,"\\%g"+(zz>=0),zz);
        else sprintf(xnr,"%g",zz);
        psstring(xnr,0);
#else
        fprintf(my_style.ps,"newpath %.2f %.2f moveto (%g) dup stringwidth pop -2 div %.2f -1.2 mul rmoveto show\n",
                x, y, zz, my_style.fontsize);
#endif
      } else {
#ifdef EXTPSSTRING
        if (my_style.minus) sprintf(xnr,"\\%g"+(zz>=0),zz);
        else sprintf(xnr,"%g",zz);
        if (strlen(xnr)>strlen(xnr0)) strcpy(xnr0,xnr);
        fprintf(my_style.ps,"newpath %.2f %.2f moveto ",
                x-my_style.fontsize*0.33,
                y-my_style.fontsize*0.33);
        psstring(xnr,1);
#else
        sprintf(xnr,"%g",zz);
        if (strlen(xnr)>strlen(xnr0)) strcpy(xnr0,xnr);
        fprintf(my_style.ps,"newpath %.2f %.2f moveto (%s) dup stringwidth pop neg %.2f -.33 mul add %.2f -0.33 mul rmoveto show\n",
                x,y, xnr,my_style.fontsize,my_style.fontsize);
#endif
      } }
    zz+=dz; }
  
  if (tick>0 && it==0) {
    if (n>0)
#ifdef EXTPSSTRING
      {
        fprintf(my_style.ps,"%% x-axis label\n\
newpath %.2f %.2f moveto\n",(x0+X0)/2, y0-3*my_style.fontsize);
        psstring(my_style.xlabel,0);
      } else {
        if (my_style.yrot) {
          fprintf(my_style.ps,"%% y-axis label\n\
gsave\n\
%.2f %.2f translate %.2f rotate\n\
0 %.2f moveto\n",
                  x0,(y0+Y0)/2,my_style.yrot,3*my_style.fontsize);
          psstring(my_style.ylabel,0);
          fprintf(my_style.ps,"grestore\n"); }
        else {
          fprintf(my_style.ps,"%% y-axis label\n\
newpath %.2f %.2f moveto (%s) stringwidth pop neg %.2f -1.3 mul add %.2f -0.33 mul rmoveto\n", x0,(y0+Y0)/2,xnr0,my_style.fontsize,my_style.fontsize);
          psstring(my_style.ylabel,1); }
    }
#else
      fprintf(my_style.ps,"%% x-axis label\n\
newpath %.2f %.2f moveto (%s) dup stringwidth pop -2 div %.2f -3 mul rmoveto show\n",
              (x0+X0)/2, y0,my_style.xlabel,my_style.fontsize);
    else
      fprintf(my_style.ps,"%% y-axis label\n\
newpath %.2f %.2f moveto (%s) dup stringwidth pop (%s) stringwidth pop add neg %.2f -1.3 mul add %.2f -0.33 mul rmoveto show\n",
              x0,(y0+Y0)/2, my_style.ylabel,xnr0,my_style.fontsize,my_style.fontsize);
#endif
  }
  
  fprintf(my_style.ps,"%.2f setlinewidth\n",my_style.thick);
}

void myopenplotps(char *fn,char *title) /********************** myopenplotps */
{
  time_t t;
  static char macros[]={"\n\
% X Y1 Y2 XTICK errybar, X1 X2 Y YTICK errxbar\n\
/errybar { newpath 3 index 1 index sub 3 index moveto 3 index 1 index add 3 index lineto 3 index 1 index sub 2 index moveto 3 index add 1 index lineto 2 index exch moveto lineto stroke } def\n\
/errxbar { newpath 3 index 2 index 2 index sub moveto 3 index 2 index 2 index add lineto 2 index 2 index 2 index sub moveto 2 index 2 index 2 index add lineto 3 index 2 index moveto 2 index 2 index lineto stroke } def\n\
% X Y RADIUS point\n\
/newxy { newpath 2 index 2 index moveto } def\n\
/square { gsave [] 0 setdash 0.886 mul newxy dup dup dup dup rmoveto -2 mul 0 rlineto 0 exch -2 mul rlineto 2 mul 0 rlineto pop pop closepath stroke grestore } def\n\
/fullsquare { 0.886 mul newxy dup dup dup dup rmoveto -2 mul 0 rlineto 0 exch -2 mul rlineto 2 mul 0 rlineto pop pop closepath fill } def\n\
/opensquare { 2 index 2 index 2 index gsave 1 setgray fullsquare grestore square } bind def\n\
/circle { gsave [] 0 setdash newpath 0 360 arc stroke grestore } def\n\
/fullcircle { newpath 0 360 arc fill } def\n\
/opencircle { 2 index 2 index 2 index gsave 1 setgray fullcircle grestore circle } bind def\n\
/diamond { gsave [] 0 setdash 1.253 mul newxy dup 0 rmoveto dup -1 mul dup dup rlineto 1 index rlineto dup rlineto pop pop closepath stroke grestore } def\n\
/fulldiamond { 1.253 mul newxy dup 0 rmoveto dup -1 mul dup dup rlineto 1 index rlineto dup rlineto pop pop closepath fill } def\n\
/opendiamond { 2 index 2 index 2 index gsave 1 setgray fulldiamond grestore diamond } bind def\n\
/up { gsave [] 0 setdash 1.34 mul newxy dup dup -0.577 mul rmoveto dup -2 mul 0 rlineto dup 1.732 mul rlineto pop pop closepath stroke grestore } def\n\
/fullup { 1.34 mul newxy dup dup -0.577 mul rmoveto dup -2 mul 0 rlineto dup 1.732 mul rlineto pop pop closepath fill } def\n\
/openup { 2 index 2 index 2 index gsave 1 setgray fullup grestore up } bind def\n\
/down { gsave [] 0 setdash -1.34 mul newxy dup dup -0.577 mul rmoveto dup -2 mul 0 rlineto dup 1.732 mul rlineto pop pop closepath stroke grestore } def\n\
/fulldown { -1.34 mul newxy dup dup -0.577 mul rmoveto dup -2 mul 0 rlineto dup 1.732 mul rlineto pop pop closepath fill } def\n\
/opendown { 2 index 2 index 2 index gsave 1 setgray fulldown grestore down } bind def\n\
/plus { gsave [] 0 setdash 1.253 mul newxy dup 0 rmoveto dup -2 mul 0 rlineto closepath stroke newxy dup 0 exch rmoveto -2 mul 0 exch rlineto pop pop closepath stroke grestore } def\n\
/cross { gsave [] 0 setdash 0.886 mul newxy dup dup rmoveto dup -2 mul dup rlineto closepath stroke newxy dup dup -1 mul rmoveto dup -2 mul exch 2 mul rlineto pop pop closepath stroke grestore } def\n\
/star { gsave [] 0 setdash 1.1 mul newxy dup 0 rmoveto dup -2 mul 0 rlineto closepath stroke newxy dup -0.5 mul dup -1.732 mul rmoveto dup dup -1.732 mul rlineto closepath stroke newxy dup -0.5 mul dup 1.732 mul rmoveto dup dup 1.732 mul rlineto pop pop closepath stroke grestore } def\n\
/dotcircle { 2 index 2 index 2 index opencircle 0.25 mul fullcircle } def\n\
"};

  time(&t);

 /* %%Page: 1 1 after 1st command is needed for mpage to work corectly */

  my_style.ps=fopen(fn,"wt");
  switch (my_style.mode) {
    case 'P': fprintf(my_style.ps,"%%!PS-Adobe-2.0\n\
%%%%Creator: MACSIMUS/plot\n\
/%s findfont %.2f scalefont setfont\n\
%%%%Page: 1 1\n\
%s\n\
gsave %.0f %.0f translate 0 setgray 1 setlinejoin\n",
                      my_style.font,my_style.fontsize,
                      macros,
                      my_style.bx,my_style.by);
              break;
    case 'L': fprintf(my_style.ps,"%%!PS-Adobe-2.0\n\
%%%%Creator: MACSIMUS/plot\n\
%%%%Orientation: Landscape\n\
/%s findfont %.2f scalefont setfont\n\
%%%%Page: 1 1\n\
%s\n\
gsave 90 rotate %.0f %.0f translate 0 setgray 1 setlinejoin\n",
                      my_style.font,my_style.fontsize,
                      macros,
                      my_style.by,-595+my_style.bx);
              break;
/*..... %!PS-Adobe 3.0 */
    case 'E': fprintf(my_style.ps,"%%!PS-Adobe-3.0 EPSF-3.0\n\
%%%%BoundingBox: 0 0 %.0f %.0f\n\
%%%%Creator: MACSIMUS/plot\n\
%%%%Title: %s\n\
%%%%CreationDate: %s\
%%%%EndComments\n\
/Helvetica findfont %.2f scalefont setfont\n\
%s\n\
gsave %.0f %.0f translate 0 setgray 1 setlinejoin\n",
                      my_style.xsize+my_style.bx+my_style.bX,
                      my_style.ysize+my_style.by+my_style.bY,
                      title,ctime(&t),
                      my_style.fontsize,
                      macros,
                      my_style.bx,my_style.by);
    default: ; }

  /* patch by T.Trnka to enable latin-1 symbol for Angstrom */
  /* ???, does not seem to work ... */
  if (my_style.enc) {
    fprintf(my_style.ps, "/latinfont\n\
  << /%s findfont {} forall >>\n\
  begin\n\
    /Encoding ISOLatin%dEncoding 256 array copy def currentdict\n\
  end\n\
definefont pop\n\
/latinfont %.2f selectfont\n", my_style.font, my_style.enc,my_style.fontsize);
    my_style.font="latinfont"; }
}

void mycloseplotps(void) /************************************ mycloseplotps */
{
  fprintf(my_style.ps,"showpage\n");
  if (my_style.mode=='E') fprintf(my_style.ps,"%%%%EOF\n");
  fclose(my_style.ps);
}

void myclip(void) /************************************************** myclip */
{
  fprintf(my_style.ps,"gsave\n\
newpath 0 0 moveto 0 %.2f lineto %.2f %.2f lineto %.2f 0 lineto closepath clip\n",
          my_style.ysize, my_style.xsize,my_style.ysize, my_style.xsize);
}
