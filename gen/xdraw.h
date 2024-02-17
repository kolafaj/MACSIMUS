/*
  "new" MACSIMUS graphics interface
  see also xdraw.c

  to change the old interface (draw.h + mouselib.h) into this new interface:
  * in metamake:
    - replace draw + mouselib by xdraw
  * in c headers:
    - replace #include "draw.h" + "mouselib.h" by "xdraw.h"
  * in the code:
    - replace getch() by readkey() (cf. kbhit())
    - remove Imworking()

*/

#include "prec.h"

#include "loop.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>

/*
   WARNING: not portable because sizeof(KeySym)=8 (on my machine)
   empirically determined that KeySym<0xffff, so we spare a bit
*/
typedef unsigned short myKeySym;

/* number of fonts; see also asignments in selectfont() in xdraw.c */
/* fonts are numbered 0..NXFONTS-1 */
#define NXFONTS 5 /* NEWFONT: advance the number */

#define PALLEV 256 /* color levels in palette - true colors */

#define X_REDRAW_ME   28
#define X_UNKNOWN_KEY 29

#define LEFTDRAG      128
#define MIDDRAG       129
#define RIGHTDRAG     130
#define LEFTCLICK     131
#define MIDCLICK      132
#define RIGHTCLICK    133
#define WHEELFORWARD  134
#define WHEELBACKWARD 135

#define UP    XK_Up
#define DOWN  XK_Down
#define LEFT  XK_Left
#define RIGHT XK_Right
#define INS   XK_Insert
#define DEL   XK_Delete
#define HOME  XK_Home
#define END   XK_End
#define PGUP  XK_Prior // also XK_Page_Up
#define PGDN  XK_Next // also XK_Page_Down
#define ESC   27
#define ENTER 13
#define TAB   9
#define BACKSPACE XK_BackSpace

#define F1 XK_F1
#define F2 XK_F2
#define F3 XK_F3
#define F4 XK_F4
#define F5 XK_F5
#define F6 XK_F6
#define F7 XK_F7
#define F8 XK_F8
#define F9 XK_F9
#define F10 XK_F10
#define F11 XK_F11
#define F12 XK_F12

#include <time.h>

extern struct xwindowhints_s {
  int icowidth,icoheight; /* X icon size */
/* unsigned ? */
  char *icobits;    /* icon data */
  char *geometry;   /* standard x-window geometry (passed from -g option) */
  char *background; /* background color (passed from -bg option) */
  char *winname;    /* X window name */
  char *iconame;    /* X icon name */
  int xmin,ymin;    /* minimum window size */
  int fullscreen;   /* current state is fullscreen */
  int xfs,yfs;      /* fullscreen size, see -R in jkv */
  char lastgeometry[16]; /* geometry before fullscreen (BUG: size only, not pos) */
} xwindowhints;

/* string to be processed instead of normal input; e.g., option -I */
void addkbdstring(const unsigned char *s);

void putpixel(int x,int y,int c);

extern int x_xsize,x_ysize,x_color,x_fillcolor,x_verbose;
#define getmaxx() (x_xsize-1)
#define getmaxy() (x_ysize-1)
#define getcolor() x_color
void line(int x0,int y0,int x1,int y1);
void bar(int x0,int y0,int x1,int y1);
void rectangle(int x0,int y0,int x1,int y1);

extern GC gc; // REMOVE AFTER PORTING

extern Display *display;
extern Window win;

#define getpixel(X,Y) 0 /* NOT SUPPORTED !!! */

/* Turbo colors + more */
enum COLORS {
    BLACK,		    /* dark colors */
    BLUE,
    GREEN,
    CYAN,
    RED,
    MAGENTA,
    BROWN,
    LIGHTGRAY,
    DARKGRAY,		    /* light colors */
    LIGHTBLUE,
    LIGHTGREEN,
    LIGHTCYAN,
    LIGHTRED,
    LIGHTMAGENTA,
    YELLOW,
    WHITE
};

extern struct dosmouse_s {
  int button[5];/* 0..3: colors from black to white for buttons:
                   2=background (should be panel background, too)
                   0=text (black)
                   4=color of button underlined (key) char
                   (all colors are in the selected mapping) */
  int x,y;      /* current pointer position */
  int x0,y0;    /* start of DRAG position */
  int xl,yl;    /* previous DRAG position */
  int dx,dy;    /* difference from the last DRAG position (DRAG returned) */
  int del;      /* delay before an attempt to read another event */
  unsigned key; /* key(button) code (KeySym) */
  int drag;     /* 0=1st of drag, 1=during DRAG, 2=last of DRAG */
  int prefix;   /* repeat prefix character; 0=no repetition */
  int nbuf;     /* buffer (local macro/repeat queue) size */
  int pos;      /* position in the buffer (number of chars) */
  myKeySym *buf;/* buffer for expanded macros (myKeySym not portable!)*/
  const unsigned char *kbdstring;
                /* see addkbdstring(), replaced by NULL after read */
  int kbdstringpos; /* next char in kbdstring */
} dosmouse;
void closegraph(void);
void setcolor(int c);
void setviewport(int x0,int y0,int x1,int y1,int clip);
void clearviewport(void);
unsigned kbhit(void);
void repeatprefix(int key);
int readkey(void);
void delay(unsigned msec);
void getimage(int x0,int y0,int x1,int y1,void *buf);
void putimage(int x,int y,void *buf,int mode);
void settextstyle(int a,int b,int size);
int gettextsize(void);
void settextjustify(int x,int y);

/* fixed font control info */
extern struct xfont_s {
  int width;     /* font spacing (distance of consecutive characters) */
  int xwidth;    /* max character width */
  int height;    /* pixel range incl. accents; max. 32 */
  int top;       /* line (topmost line incl. accents =0) of top of digits */
  int base;      /* baseline (of digits) */
  int from;      /* code of the 1st character */
  int to;        /* 1+code of the last character */
  unsigned *data;/* [(code-from)*width], columns of pixel data (black=1 on white=0) */
} xfont;
void selectfont(int font);
int outtextwidth(const char *s);
int outtextppm(unsigned char (*rgb)[3],int xsize,int ysize,
               unsigned char RGB[3],
               int x,int y,const char *s);
int outtextxy(int x,int y,const char *s);
int outtextxymax(int x,int y,const char *s,int m);
int outtext(const char *s);
int outbgtextxy(int x,int y,int font,int fg,int bg,const char *s);

void setlinestyle(int mode,int pattern,int thick);
void setfillstyle(int mode,int color);
void setwritemode(int mode);
void fillellipse(int x,int y,int rx,int ry);

void circle(int x,int y,int r);
XArc *startcircles(int n);
void circles(XArc *c,int x,int y,int r);
void showcircles(XArc *c);

void setpal(int col,int r,int g,int b,int maxbright);

/* common functions */

extern struct scaling_s { REAL fromx,x,toy,y; } scaling;
//extern int pen;

int startgraph(int mode);
void up(void);
void lwindow(int x0,int x1,int y0,int y1);
void scale(REAL x0,REAL x1,REAL y0,REAL y1,int aspect);
int SX(REAL x);
int SY(REAL y);
REAL iSX(int ix);
REAL iSY(int iy);
void draw(REAL x1,REAL y1);
void lline(REAL x1,REAL y1,REAL x2,REAL y2);
void dot(REAL x1,REAL y1);
void bigdot(REAL x1,REAL y1);
extern REAL pointsize;
extern int insidecolor;
void point(REAL x1,REAL y1,char type);
#define smartdot(X,Y) point(X,Y,'o'|128)

void voidbutton(int x,int y,const char *text);
int makebutton(int x, int y, unsigned key,const char *text, const char *help);
void erasebuttons(void);
void makeslider(int x,int y,int len,double *pos);
void drawslider(int forceredraw);
int isinslider(void);
void hline(int x,int y,int dx);
#define NATX 8
extern int ATX[NATX];
void planline(int atx,char *L);
void fullscreentoggle(void);
