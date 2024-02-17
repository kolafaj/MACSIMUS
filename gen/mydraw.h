/* for plotcalc and plotfile */

extern FILE *fig_dat;
extern struct my_style_s {
  FILE *ps;          /* output ps or eps file */
  char *font;        /* default font */
  float fontsize;    /* PostScript font size in pt */
  int enc;           /* encoding; 0=ASCII, 1=iso-8859-1 (contains angstrom (?)), etc. */
  int minus;         /* 0: use hyphen - for minus, 1: use Symbol - (math)*/
  float xsize,ysize; /* graph size (w/o labeling) in pt */
  float bx,bX,by,bY; /* added as margins to define the bounding box */
  float yrot;        /* rot angle of y label (normally 0 or -90) */
  int nticks;
  float tick[5];     /* [nticks] tick length in pt */
  float xticks[5],yticks[5]; /* >0: approx. # of ticks/axis, otherwise the tick */
  float thick,framethick; /* in pt */
  float pointsize;   /* PostScript point size in pt (of circle of eq. area) */
  char mode;         /* P=portrait, L=landscape, E=eps, other=no PS headers */
  char *xlabel,*ylabel;
  struct my_style_line_s { char *rgb; char *dash; float thick; char *symbol; float size; } line[16];
} my_style;

extern struct rg_s {
  double x,X,y,Y;
} my_rg;

void myup(void);
double mysx(double x);
double mysy(double y);
void mydraw(double x,double y);
void mysetcolor(int col);
void mysetlinestyle(int style,int i,int j);
void myerrxbar(double x,double y1,double y2);
void myerrybar(double x1,double x2,double y);
void mysmartdot(double x,double y);
void myaxis(int it,int n,double tick);
void myopenplotps(char *fn,char *title);
void mycloseplotps(void);
void myclip(void);

void psstring(char *s,int align);
void mypssetxy(double x,double y,int col);
