/* cc -g -DX11 -DDEBUG -c xdraw.c; cc -g -o outtextx -DX11 -DDEBUG outtextx.c xdraw.o -lm -lX11
   cc -g -DX11 -DDEBUG -c xdraw.c
  "NEW" MACSIMUS graphics interface
  (merges old-versions of: draw, drawx11, mouselib, mousekey; buttons)

  [ Long time ago, MACSIMUS used the Turbo C graphics library and DOS    ]
  [ full screen graphics...                                              ]
  [ Then I wrote a Turbo C graphics emulator under X11.                  ]
  [ Then I added some other layers, like DOS mouse (could be on/off) and ]
  [ extended keys; e.g., a single function key or cursor was transformed ]
  [ to a DOS 0+KEY code, and then translated to a single key again; mouse]
  [ clicks were detected, but re-analyzed later from push/release data;  ]
  [ all keys were buffered, but not mouse events.                        ]
  [ In the NSK and old simolant projects, buttons were implemented.      ]
  [ DOS is dead! This is an attempt to clean the mess a bit, without too ]
  [ many changes in the existing MACSIMUS code.                          ]

  MACSIMUS basic user interface functions (in graphics window)

  Typical usage:

    for (;;) { // hotkey and mouse event loop
      if (kbhit()) // with this line:
                   //   repeat show(), change parameters during showing
                   // without this line:
                   //   wait for a key, process it, then show()...
        switch (readkey()) {
          case 0: // needed!
            break;
          case 'A'&31: // hot key Ctrl-A
            break;
          ...
          default: // unsupported key
            warning("no such key");
        }
      show() // show something
    }

  Any event is transformed into a short unsigned int; normal characters are
  just ASCII codes.
  No support for UTF-8/UNICODE!
  ESC=Ctrl-[, TAB=Ctrl-I, etc., as in terminal or DOS

  Mouse events are transformed into (used to be in mousekey, now here):
    LEFTCLICK = left mouse pressed, left mouse released, no x,y change
    LEFTDRAG = left mouse down, x,y change
    MIDCLICK,MIDDRAG,LEFTCLICK,LEFTDRAG = similarly
  MotionNotify without a mouse pressed does not generate any user-event
  MotionNotify with a key pressed can skip some events if slow

  See also below:
    kbhit()
    readkey()
    macro[][]
*/

#define Sqr(X) ((X)*(X))

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h> /* fabs needed */
#include "xdraw.h"

#define max(A,B) ((A)<(B)?(B):(A))

/*
#define DEBUG
  #define DEBUG: detail debug info
  [ #define GUI : REMOVED (replaced by in-house buttons) ]
  [ environment variable LAZYX11 : REMOVED ]
  environment variable GUI:
    p g s - affects plot, blend, show
    v = verbose
    #HEXCOLOR - background (wire/Turbo legacy modes)
  0: verbose output
*/

int x_color,x_fillcolor; /* as set by setcolor, setfillstyle */
int x_verbose=0;
static int currentcolor=-1; /* as color of gc's */

#include "drawbit.c"

#define Isdigit(X) ((unsigned)(X)<128 && isdigit(X))
#define Islower(X) ((unsigned)(X)<128 && islower(X))
#define Isupper(X) ((unsigned)(X)<128 && isupper(X))

struct xwindowhints_s xwindowhints = {
  drawbit_width,drawbit_height,/* X icon size */
  drawbit_bits,   /* icon data */
  NULL,           /* standard x-window geometry (passed from -g option) */
  "#000000",      /* background color (passed from -bg option) */
  "JK X-draw",    /* X window name */
  "X-draw",       /* X icon name */
  128,96,         /* min window size */
  0,              /* fullscreen status */
  0,0,            /* fullscreen size */
  ""              /* lastgeometry */
};

/*                             button colors: */
struct dosmouse_s dosmouse = { {BLACK,DARKGRAY,LIGHTGRAY,WHITE, LIGHTCYAN} };

Display *display;
int screen_num,display_width,display_height,no_colors;
Window win;
XSizeHints size_hints;
XEvent event;
XColor xcolor,hwcolor;
Colormap cmap;
GC gc;

/* button remap, prepared to be general, but currently 2<->3 active via `TAB */
int ButtonMap[6]={0,Button1,Button2,Button3,Button4,Button5};

double mytime(void) /************************************************ mytime */
  /* see also cputime.c (cf. clock_gettime etc.) */
{
  struct timeval tv;
  struct timezone tz;

  gettimeofday(&tv,&tz);

  return (unsigned)tv.tv_usec/1e6+tv.tv_sec;
}

inline
static void mysleep(double dt) /*********************************** mysleep */
{
  if (dt<0) return;
  if (dt>1) dt=1; /* max sleeping time 1 s */
  usleep(dt*1000000);
}

#define XERR(X,R) { fprintf(stderr,"xdraw error: %s\n",X); return R; }
#define XERE(X) { fprintf(stderr,"xdraw ERROR: %s\n",X); exit(-3); }

#include "xbuttons.c"

long unsigned *xcoltab;

#define CLICK (LEFTCLICK-1)
#define DRAG (LEFTDRAG-1)

static struct pending_s {
  int x,y;    /* pointer position */
  int button; /* mouse button pressed; 1,2,3 or 0=none */
  int state;  /* one of CLICK, DRAG, and 0 */
} pending;

int x_xsize,x_ysize;
static int NEWdel=200000;

void addkbdstring(const unsigned char *s) /******************** addkbdstring */
{
  if (s && !*s) s=NULL;
  dosmouse.kbdstring=s;
  dosmouse.kbdstringpos=0;
}

static int getXevent(int block) /********************************* getXevent */
/*
  the low-level event (mouse and key) interface
  no macro handling here, no delays
  LEFT buttons are translated into characters
  RIGHT buttons show the associated help, a key is expected
  block=1: wait for an event (block execution) if none ready in the X-system
           when an event comes, read it (= store as pending in dosmouse.*)
  block=0: check for an event, if there is one read it as if block=1 and
           return nonzero, otherwise return 0 (do not block)
*/
{
  int aux=1;


  if (dosmouse.kbdstring) {
    dosmouse.key=dosmouse.kbdstring[dosmouse.kbdstringpos++];
    if (!dosmouse.kbdstring[dosmouse.kbdstringpos]) dosmouse.kbdstring=NULL;
    return 1; }

  if (!display) return 0;

  dosmouse.del=10; // ?

 again:

  if (block) XNextEvent(display,&event);
  else aux=XCheckMaskEvent(display,0xffffffff,&event);

  if (!aux) {
    /* no event */
    dosmouse.del=10000;
    return 0; }

#ifdef DEBUG
  fprintf(stderr,"getXevent: event came...");
#endif /*# DEBUG */

  event.type &= 0x7f; /* mask sent events */

#define NEW
#ifdef NEW

  /* ignore any ConfigureNotify without size change (just move assumed) */
  if (event.type==ConfigureNotify) {
    if (x_xsize==event.xconfigure.width && x_ysize==event.xconfigure.height) {
      if (block) goto again;
      dosmouse.del=10000;
      return 0; } }

  /* return only one X_REDRAW_ME for a series of ConfigureNotify */
  if (event.type==ConfigureNotify || event.type==Expose) {
    int del=NEWdel;

    for (;;) {
      if (event.type==ConfigureNotify) {
        x_xsize=event.xconfigure.width;
        x_ysize=event.xconfigure.height;
#  ifdef DEBUG
        fprintf(stderr,"ConfigureNotify: x-size=%d y-size=%d\n",x_xsize,x_ysize);
#  endif /*# DEBUG */
      }
#  ifdef DEBUG
      else
        fprintf(stderr,"Expose\n");
#  endif /*# DEBUG */

      usleep(del);
      del=1000;
      aux=XCheckMaskEvent(display,0xffffffff,&event);

      if (!aux) {
        return dosmouse.key=X_REDRAW_ME; }
      else {
        event.type &= 0x7f;
        if (!(event.type==ConfigureNotify || event.type==Expose)) {
          XPutBackEvent(display,&event);
          return dosmouse.key=X_REDRAW_ME;  } } } }
#endif /*# NEW */

  switch (event.type) {

    /* does not work (catch close window ..) */
    case DestroyNotify:
      exit(0);

#ifdef xDEBUG
    /* ???: when a keyboard map is changed,
       keeps generating MappingNotify forever */
    case MappingNotify:
      fprintf(stderr,"MappingNotify\n");
      return dosmouse.key=X_REDRAW_ME;
#endif /*? DEBUG */ /*# xDEBUG */

#ifndef NEW
      /* OLD VERSION, every Expose/ConfigureNotify causes redraw */
    case Expose:
#  ifdef DEBUG
      fprintf(stderr,"Expose\n");
#  endif /*# DEBUG */
      //      if (event.xexpose.count!=0) break;
      return dosmouse.key=X_REDRAW_ME;

    case ConfigureNotify:
      x_xsize=event.xconfigure.width;
      x_ysize=event.xconfigure.height;
#  ifdef DEBUG
      fprintf(stderr,"ConfigureNotify: x-size=%d y-size=%d\n",x_xsize,x_ysize);
#  endif /*# DEBUG */
      return dosmouse.key=X_REDRAW_ME;
#endif /*# NEW */

    case ButtonPress:
#ifdef DEBUG
      fprintf(stderr,"ButtonPress: last x=%d y=%d  button=%d\n",
              dosmouse.x,dosmouse.y,event.xbutton.button);
#endif /*# DEBUG */
      if (event.xbutton.button==ButtonMap[4]) return dosmouse.key=WHEELFORWARD;
      if (event.xbutton.button==ButtonMap[5]) return dosmouse.key=WHEELBACKWARD;
      pending.x=event.xbutton.x;
      pending.y=event.xbutton.y;

      if (event.xbutton.button==ButtonMap[1]) {
        struct buttons_s *b=getbutton();
        if (b) {
          do XNextEvent(display,&event); while (event.type != ButtonRelease);
          pending.button=0;
          releasebutton(b);
          return dosmouse.key=b->key; }
        pending.button=1; }
      else if (event.xbutton.button==ButtonMap[3]) {
        struct buttons_s *b=getbutton();
        if (b) {
          outbgtextxy(2,2,2,dosmouse.button[0],dosmouse.button[3],b->help);
          pending.button=0;
          do XNextEvent(display,&event); while (event.type != ButtonRelease);
          do XNextEvent(display,&event); while (event.type == MotionNotify);
          return dosmouse.key=X_REDRAW_ME; }
        pending.button=3; }
      else if (event.xbutton.button==ButtonMap[2]) pending.button=2;

      pending.state=CLICK;
      /* no key returned:
         waiting until pointer motion or button release */
      //      return dosmouse.key=0;
      goto again;

    case MotionNotify:
#ifdef DEBUG
      fprintf(stderr,"MotionNotify: x=%d y=%d  pending.button=%d\n",
              dosmouse.x,dosmouse.y,pending.button);
#endif /*# DEBUG */
      if (pending.button) {
        /* motion with button pressed -> drag event */

        /* skiping all MotionNotify's and accepting the final one only */
        while (XCheckMaskEvent(display,ButtonMotionMask,&event));

        dosmouse.xl=pending.x;
        dosmouse.yl=pending.y;

#if 0
        /* ? */
        if (pending.state==CLICK) {
          /* 1st move after click: return the initial position */
          dosmouse.dx=dosmouse.dy=0;
          pending.state=DRAG;
          dosmouse.drag=0;
          dosmouse.x=dosmouse.x0=pending.x; pending.x=event.xbutton.x;
          dosmouse.y=dosmouse.y0=pending.y; pending.y=event.xbutton.y; }
        else {
          dosmouse.drag=1;
          dosmouse.dx=event.xbutton.x-pending.x;
          dosmouse.dy=event.xbutton.y-pending.y;
          /* following drag: returning the true position */
          dosmouse.x=pending.x=event.xbutton.x;
          dosmouse.y=pending.y=event.xbutton.y; }
#else /*# 0 */
        if (pending.state==CLICK) {
          pending.state=DRAG;
          dosmouse.drag=0;
          dosmouse.x0=pending.x;
          dosmouse.y0=pending.y;
          dosmouse.dx=event.xbutton.x-pending.x;
          dosmouse.dy=event.xbutton.y-pending.y;
          dosmouse.x=pending.x=event.xbutton.x;
          dosmouse.y=pending.y=event.xbutton.y; }
        else {
          dosmouse.drag=1;
          dosmouse.dx=event.xbutton.x-pending.x;
          dosmouse.dy=event.xbutton.y-pending.y;
          /* following drag: returning the true position */
          dosmouse.x=pending.x=event.xbutton.x;
          dosmouse.y=pending.y=event.xbutton.y; }
#endif /*#!0 */

        pending.state=DRAG;
        return dosmouse.key=LEFTDRAG-1+pending.button; }

      /* motion without button pressed -> no user-event */
      dosmouse.x=event.xbutton.x; /* not needed */
      dosmouse.y=event.xbutton.y; /* not needed */
      dosmouse.key=0;
      goto again;

    case ButtonRelease:
#ifdef DEBUG
      fprintf(stderr,"ButtonRelease: x=%d y=%d\n",
              event.xbutton.x,event.xbutton.y);
#endif /*# DEBUG */
      dosmouse.drag=2; /* end of drag */

      /* release of a key that has not been pressed...? */
      if (!pending.button) goto again;

      dosmouse.x=event.xbutton.x;
      dosmouse.y=event.xbutton.y;
      if (pending.state==LEFTCLICK-1) {
        /* release after press -> click */
        aux=pending.button;
        pending.button=0;
        dosmouse.drag=3; /* not needed */
        return dosmouse.key=LEFTCLICK-1+aux; }

      /* end of drag: return the LAST position; the current is ignored (almost always is the same) */
      dosmouse.dx=event.xbutton.x-pending.x;
      dosmouse.dy=event.xbutton.y-pending.y;
      dosmouse.x=pending.x;
      dosmouse.y=pending.y;
      aux=pending.button;
      pending.button=0; /* no button pressed (multiple buttons pressed ignored) */
      return dosmouse.key=LEFTDRAG-1+aux;
      break;

    case KeyPress: {
      int n;
      char buf[8];
      KeySym key;
      XComposeStatus compose;

      n=XLookupString((XKeyEvent*)&event,buf,8,&key,&compose);

#ifdef DEBUG
      fprintf(stderr,"KeyPress: keysym=%d=\"%c\" n=%d\n",
              (int)key,(int)key,n);
      if (key>0xffff) fprintf(stderr,"WARNING: keysym>0xffff\n");
#endif /*# DEBUG */

      if (n>1) fprintf(stderr,"XRebind‐Keysym not supported, only 1st of %d chars processed\n",n);

      switch (key) {
        case UP:   case DOWN:
        case LEFT: case RIGHT:
        case INS:  case DEL:
        case HOME: case END:
        case PGUP: case PGDN:
        case BACKSPACE: // added 2023-06-25 (before, Ctrl-H = Backspace)

        case F1: case F2: case F3: case F4:  case F5:  case F6:
        case F7: case F8: case F9: case F10: case F11: case F12:
          return dosmouse.key=key;

        case XK_Shift_L: case XK_Shift_R:
        case XK_Control_L: case XK_Control_R:
        case XK_Alt_L: case XK_Alt_R:
          goto again; /* no user-event (used to be return dosmouse.key=0) */

        default:
          if (buf[0]) return dosmouse.key=buf[0]; /* normal ASCII character */
          else return dosmouse.key=X_UNKNOWN_KEY; /* any control/function key
                                                     not in the above table */
      }
    }
  }
  //  return dosmouse.key=0;
#ifdef DEBUG
  fprintf(stderr,"getXevent: ERROR, end reached and 0 returned\n");
#endif /*# DEBUG */
  goto again;
} /* getXevent() */

/* delay between kbhit() if called too often, in CLOCKS_PER_SEC (assumed 1000000) */
#define KBHITDELAY 0.01

unsigned kbhit(void) /************************************************ kbhit */
/*
  checks for an X event if not already (in dosmouse.key)
  if ESC then erases the queue
  if the queue is not empty or an event is ready, 1 is returned, otherwise 0
  waits upto 0.01 s if called too often
*/
{
  static double t0;
  double t=mytime(),dt=t-t0;
  t0=t;

  if (!display) return 0;

  mysleep(KBHITDELAY-dt);

#ifdef xDEBUG
  fprintf(stderr,"kbhit: key=\'%c\'=%x (pos=%d)\n", dosmouse.key,dosmouse.key,dosmouse.pos);
#endif /*? DEBUG */ /*# xDEBUG */

  if (!dosmouse.key) getXevent(0);

  if (dosmouse.key==ESC) {
    if (dosmouse.pos) { // buffer erased due to ESC
      dosmouse.pos=dosmouse.key=0;
      return getXevent(0); } } // ESC swallowed, check next event

  if (dosmouse.pos) return 1;
  if (dosmouse.key) return 1;

  return 0;
}

void repeatprefix(int key) /*********************************** repeatprefix */
/*
  key which works like ESC in emacs, usually '`'
*/
{
  dosmouse.prefix=key;
  dosmouse.pos=0;
  dosmouse.nbuf=10240;
  dosmouse.buf=malloc(dosmouse.nbuf*sizeof(dosmouse.buf[0]));
  if (!dosmouse.buf) {
    fprintf(stderr,"no heap - exit(8)\n");
    exit(8); }
}

static void extendbuf(int nbuf) /********************************* extendbuf */
/* resize the key buffer if needed */
{
  if (nbuf>16777216) {
    fprintf(stderr,"xdraw: too long key buffer requested - exit(4)\n");
    exit(4); }

  if (nbuf>=dosmouse.nbuf) {
    myKeySym *newbuf;

    dosmouse.nbuf=nbuf;
    newbuf=malloc(dosmouse.nbuf*sizeof(dosmouse.buf[0]));
    if (!newbuf) {
      fprintf(stderr,"no heap - exit(8)\n");
      exit(8); }
    memcpy(newbuf,dosmouse.buf,dosmouse.pos);
    free(dosmouse.buf);
    dosmouse.buf=newbuf; }
}

#define MACLEN 32
static struct macro_s {
  int n; /* length of macro (max. MACLEN) */
  myKeySym k[MACLEN];
} macro[26];
/*
  for dosmouse.prefix=`
  `A@@@@` = define macro "A" @@@@ = max 31 characters (more ignored)
            no "`" allowed inside macro (no nesting)
            26 macros "A".."Z" are available
  `A` = erase macro "A" (=replace by empty string)
  `ctrl-A = print macro "A" to stderr
  `a = invoke macro "A"
  `####@ = repeat character @ ####-times (####=decimal, max 9999)
  `####`a = repeat macro "A" ####-times (####=decimal, max 9999)
  `/ = swap mid <-> right mouse buttons
*/

int readkey(void) /***************************************************** readkey */
/*
  Read a key/event - user layer; mouse position is in mousekey.x,mousekey.y
  - 0 returned if macro interrupted (by ESC) and similar
  - an ASCII letter returns this letter, cast to (int); e.g., 'a' returns 97
  - shift-letter returns simply the uppercase letter; 'B' returns 66
  - ctrl-letter returns the ctrl code; Ctrl-C returns 3
  - other key returns its keysym
  - a mouse button returns CLICK and DRAG
  - mouse motion (without click) is suppressed (no event)
  - only shift, control, alt are suppressed (no event)
  see kbhit() for more info!
*/
{
  unsigned aux;

  if (dosmouse.prefix) { // buffer, macro, repeat support

    for (;;) {

      if (kbhit()) {

      eventcame:

#ifdef DEBUG
        fprintf(stderr,"readkey: char=\"%c\"=%x pos=%d\n",dosmouse.key,dosmouse.key,dosmouse.pos);
#endif /*# DEBUG */

        if (dosmouse.key==dosmouse.prefix) { // `
          int i,n=1,nrep,mac=0;

          getXevent(1);

          if (Isdigit(dosmouse.key)) { // `#
            for (n=0; Isdigit(dosmouse.key); getXevent(1))
              n=(n*10+dosmouse.key-'0')%10000;

#ifdef DEBUG
            fprintf(stderr,"readkey: repeat n=%d, new key=%x\n",n,dosmouse.key);
#endif /*# DEBUG */

            if (dosmouse.key==dosmouse.prefix) { // `####`
              getXevent(1);
              if (Islower(dosmouse.key)) { // `###`a
              exemac:
                mac=dosmouse.key-'a';
                nrep=macro[mac].n;
                extendbuf(dosmouse.pos+n*nrep);
                while (n--)
                  for (i=0; i<nrep; i++)
                    dosmouse.buf[dosmouse.pos++]=macro[mac].k[i];
                // on empty macro and buffer, 0 will be returned
                goto returnkey; }
              else { // `####`X (not lc => cancel)
                getXevent(1);
                goto eventcame; } }
            else { // `####@
              extendbuf(dosmouse.pos+n);
              while (n--) dosmouse.buf[dosmouse.pos++]=dosmouse.key;
              goto returnkey; } }
          else if (Islower(dosmouse.key)) { // `a
            n=1;
            goto exemac; }
          else if (Isupper(dosmouse.key)) { // `A
            mac=dosmouse.key-'A';
            getXevent(1);
            for (n=0; dosmouse.key!='`'; getXevent(1))
              if (n<MACLEN) macro[mac].k[n++]=dosmouse.key;
            macro[mac].n=n;
            if (x_verbose) fprintf(stderr,"macro %c defined\n",mac+'A');
            getXevent(1);
            goto eventcame;  }
#if 0
          /* erase macro */
          else if (dosmouse.key>=1 && dosmouse.key<=26) { // `Ctrl-A
            mac=dosmouse.key-1;
            macro[mac].n=0;
            getXevent(1);
            goto eventcame;  }
#else /*# 0 */
          /* print macro */
          else if (dosmouse.key>=1 && dosmouse.key<=26) { // `Ctrl-A
            int i;
            unsigned char c[2];

            mac=dosmouse.key-1;
            fprintf(stderr,"macro %c = |",'A'+mac);
            loop (i,0,macro[mac].n) {
              if (macro[mac].k[i]>255) c[0]='~';
              else c[0]=' ';
              c[1]=macro[mac].k[i];
              if (c[1]<32) c[0]='^',c[1]+='@';
              fprintf(stderr,"%c%c",c[0],c[1]); }
            fprintf(stderr,"|\n");
            getXevent(1);
            goto eventcame;  }
#endif /*#!0 */
          else if (dosmouse.key=='/') { // `/
            if (ButtonMap[2]==Button2) {
              ButtonMap[2]=Button3;
              ButtonMap[3]=Button2;
              if (x_verbose) fprintf(stderr,"MID AND RIGHT MOUSE BUTTONS SWAPPED\n"); }
            else {
              ButtonMap[2]=Button2;
              ButtonMap[3]=Button3;
              if (x_verbose) fprintf(stderr,"mid and right mouse buttons normal\n"); } }
          else {
            // `<anything else>
            if (x_verbose) fprintf(stderr,"macro syntax - %d ignored\n",dosmouse.key);
            getXevent(1);
            goto eventcame; } }
        else if (dosmouse.key!=0) { // any other key or event but 0 ` ESC)
          extendbuf(dosmouse.pos+64);
          dosmouse.buf[dosmouse.pos++]=dosmouse.key;
          goto returnkey; }
      } /* kbhit() */
      else if (!dosmouse.pos) { // no event
        getXevent(1);
        goto eventcame; }

    returnkey:

#ifdef DEBUG
      {
        int i;
        fprintf(stderr,"buf(%d)=",dosmouse.pos);
        loop (i,0,dosmouse.pos) fprintf(stderr,"%d%c",dosmouse.buf[i],",\n"[i==dosmouse.pos-1]);
      }
#endif /*# DEBUG */

      aux=0;
      if (dosmouse.pos) {
        aux=dosmouse.buf[0];
        memmove(dosmouse.buf,dosmouse.buf+1,dosmouse.pos*sizeof(dosmouse.buf[0]));
        dosmouse.pos--; }

      if (display) usleep(dosmouse.del); // ?

      dosmouse.key=0;

      if (aux) return aux;

      // otherwise continue with the loop (check again kbhit) ?
    } /* for(;;) */
  }
  else { // !dosmouse.prefix
    if (!dosmouse.key) getXevent(1); // else already read
    aux=dosmouse.key;
    dosmouse.key=0;
    return aux; }
}

void delay(unsigned msec) /******************************************* delay */
/* not 100% portable ... */
{
  if (!display) return; /* no delay needed in batch mode */

  usleep(msec*1000); /* cf. also nanosleep () */
}

/* one "viewport", as in the old Turbo graphics */
int isclipped;
static XRectangle XR; /* x, y, width, height */

void setviewport(int x0,int y0,int x1,int y1,int clip)
/* WARNING: cf. order of parameters in lwindow */
{
  if (!display) return;

  if (clip) {
    isclipped=1;
    XR.x=x0; XR.y=y0;
    XR.width=x1-x0+1;
    XR.height=y1-y0+1;
    XSetClipRectangles(display, gc, 0,0, &XR,1, Unsorted); }
  else {
    isclipped=0;
    XSetClipMask(display, gc, None); }
}

void clearviewport(void) /************************************ clearviewport */
{
  if (!display) return;

  if (isclipped)
    XClearArea(display,win,XR.x,XR.y,XR.width,XR.height,0);
  else
    XClearWindow(display,win);
}

long unsigned winbackground;
static int xormode; /* 0=overwrite, 1=xor */

int startgraph(int mode) /*************************************** startgraph */
{
  int atx=1,aty=1,i;
  char *gui=getenv("GUI"),*c;

  if (mode==0) return 0; /* no graphics */

  /* create X-window and set all necessary things */

  if (display) {
    /* already created - clear screen and raise instead */
    XRaiseWindow(display,win);
    XMapWindow(display,win);
    /*
       It looks like that XRaiseWindow is not sufficient, but why
       XMapWindow works, is a puzzle
    */
    clearviewport();

    {
      static XEvent dummy;
      dummy.xexpose.type=Expose;
      dummy.xexpose.display=display;
      dummy.xexpose.window=win;
      XSendEvent(display,win,False,ExposureMask,&dummy);
    }

    /*
       to eat all window exposure events (?!- at least one must come)
       (in fact, I don't know how to pop a window WITHOUT any event
    */
    i=readkey();
    if (i!=X_REDRAW_ME) fprintf(stderr,"unexpected event/key %c%c = %d\n",
                                " ^"[i<' '], i>=' '?i:i+'@', i);
    while (kbhit()==X_REDRAW_ME) (void)readkey();

    return abs(mode); }

  if (mode>16) no_colors=mode;
  else if (mode>0 && (mode&1)==0) no_colors=256; /* DOS legacy ? */
  else no_colors=16;

  xcoltab=(long unsigned*)malloc(no_colors*sizeof(xcoltab[0]));
  if (xcoltab==NULL) XERE("no heap")
  memset(xcoltab,0,no_colors*sizeof(xcoltab[0]));

  display=XOpenDisplay(getenv("DISPLAY"));
  if (!display) XERR("XOpenDisplay",-1)
  screen_num=DefaultScreen(display);
  display_width=DisplayWidth(display,screen_num);
  display_height=DisplayHeight(display,screen_num);
  //  fprintf(stderr,"xdraw: DISPLAY=%s screen=%d %dx%d\n",getenv("DISPLAY"),screen_num,display_width,display_height);

#ifndef min
#  define min(A,B) ((A)<(B)?(A):(B))
#endif /*# min */
  /* increase from 640x480? */
  //  x_xsize=min(800,display_width*2/3);
  //  x_ysize=min(600,display_height*2/3);

  /* the default if no xwindowhints.geometry */
  x_xsize=min(824,display_width*2/3);
  x_ysize=min(480,display_height*2/3);

  //  fprintf(stderr,"MIN %d %d %s\n",x_xsize,x_ysize,xwindowhints.geometry);

  if (xwindowhints.geometry) {
    /* usually -g option
       extended to shortcuts
       f or BIGxBIG = fullscreen */
    char *ch;

    switch (xwindowhints.geometry[0]) {
      /* shortcuts; f-> BIGxBIG = fullscreen */
      case 'b': xwindowhints.geometry="1024x768"; break;
      case 'g': xwindowhints.geometry="800x600"; break;
      case 'v': xwindowhints.geometry="640x480"; break;
      case 'n': xwindowhints.geometry="1680x1050"; break;
      case 'r': xwindowhints.geometry="1280x720"; break;
      case 'h': xwindowhints.geometry="1920x1080"; break;
      case 's': xwindowhints.geometry="600x600+0+0"; break;
      case 'S': xwindowhints.geometry="1024x1024"; break;
      case 'p': xwindowhints.geometry="768x1024"; break;
      case 'P': xwindowhints.geometry="960x1050"; break;
      case 'f': xwindowhints.geometry="99999x99999"; break; }

    /* parse geometry; extension: -gSIZE -> -gSIZExSIZE */
    x_xsize=x_ysize=atoi(xwindowhints.geometry);
    if ( (ch=strchr(xwindowhints.geometry,'x')) ) {
      x_ysize=atoi(ch+1);
      if ( (ch=strpbrk(ch,"+-")) ) {
	atx=atoi(ch);
        if (*ch=='-') atx--;
	if ( (ch=strpbrk(ch+1,"+-")) ) {
          aty=atoi(ch);
          if (*ch=='-') aty--; } } } }

  /* 32x32 if too small */
  if (x_xsize<32) x_xsize=32;
  if (x_xsize>display_width) x_xsize=display_width;

  if (x_ysize<32) x_ysize=32;
  if (x_ysize>display_height) x_ysize=display_height;

  /* does not take into acount window borders... */
  if (atx<0) atx+=display_width-x_xsize;
  if (aty<0) aty+=display_height-x_ysize;

  if (gui) {
    if (strchr(gui,'v')) {
      fprintf(stderr,"window geometry=%dx%d%+d%+d\n",x_xsize,x_ysize,atx,aty);
      x_verbose=1; }
    if ( (c=strchr(gui,'d')) ) NEWdel=1e6*atof(c+1)+0.5;
    if (NEWdel<1000) NEWdel=1000;
    if (NEWdel>2000000) NEWdel=2000000;
    if ( (c=strchr(gui,'#')) ) xwindowhints.background=c; }

  cmap=DefaultColormap(display,screen_num);

  /*.....winbackground=BlackPixel(display,screen_num);*/
  /* does not work generally for xor */

  /* on some systems, named colors do not work - #HEXHEX must be used */
  if (!xwindowhints.background) xwindowhints.background="#000000";

  if (XAllocNamedColor(display,cmap,xwindowhints.background,&xcolor,&hwcolor))
    winbackground=xcolor.pixel;

  /* WARNINGS (or I do not understand something here):
     Color #000000 is not drawn as black, but as background.
     On some systems, window background cannot be changed later by
     XSetBackground - is set forever by XCreateSimpleWindow.
     For sure, both ways are used here.  */

  win=XCreateSimpleWindow(display,
                          RootWindow(display,screen_num), // DefaultRootWindow(display)
			  atx,aty,
			  x_xsize,x_ysize,
			  0,
			  WhitePixel(display,screen_num),
			  winbackground);

  xwindowhints.fullscreen=0;

  if (x_xsize==display_width && x_ysize==display_height) {
    /* starting fullscreen, see also fullscreentoggle()
       see https://stackoverflow.com/questions/9065669/x11-glx-fullscreen-mode
       see also: xfullscreen.c, xfullscreenV2.c
       bug/feature: the Ubuntu panel is inaccessible */
    Atom wm_state   = XInternAtom (display, "_NET_WM_STATE", True );
    Atom wm_fullscreen = XInternAtom (display, "_NET_WM_STATE_FULLSCREEN", True );

    xwindowhints.fullscreen=1;

    XChangeProperty(display, win, wm_state, XA_ATOM, 32,
                    PropModeReplace, (unsigned char *)&wm_fullscreen, 1);

    /* ScreenCount(display) returns always 1 even if there are 2 monitors
       and the following ScreenOfDisplay returns the combined size */
    if (xwindowhints.xfs) x_xsize=xwindowhints.xfs;
    if (xwindowhints.yfs) x_ysize=xwindowhints.yfs;
    fprintf(stderr,"fullscreen: %d %dx%d\n",ScreenCount(display),x_xsize,x_ysize);
  }

#ifdef DEBUG
  fprintf(stderr,"Window ID: %d\n", (int)win);
#endif /*# DEBUG */

#if 0
  if (x_xsize>=DisplayHeight(display,screen_num) ||
      x_ysize>=DisplayHeight(display,screen_num)) {
    /* old cheep fullscreen  */
    size_hints.flags=PPosition | PSize | PMinSize | PMaxSize | PWinGravity;
    size_hints.x=0;
    size_hints.y=0;
    size_hints.min_width=x_xsize;
    size_hints.max_width=x_xsize;
    size_hints.width=x_xsize;
    size_hints.min_height=x_ysize;
    size_hints.max_height=x_ysize;
    size_hints.height=x_ysize;
    size_hints.win_gravity=CenterGravity; }
  else {
    /* window */
    size_hints.flags=PPosition | PSize | PMinSize;
    size_hints.min_width=x_xsize;
    size_hints.min_height=x_ysize; }
#else
  size_hints.flags=PPosition | PSize | PMinSize;
  size_hints.min_width=xwindowhints.xmin;
  size_hints.min_height=xwindowhints.ymin;
#endif

  {
    static XWMHints wm_hints;
    static XClassHint class_hint;
    static int xnarg; /*?*/
    static char **xarg; /*?*/

    XTextProperty winname,iconame;

    Pixmap icon_pixmap=XCreateBitmapFromData(display,win,
      xwindowhints.icobits,xwindowhints.icowidth,xwindowhints.icoheight);

    XStringListToTextProperty(&xwindowhints.winname,1,&winname);
    XStringListToTextProperty(&xwindowhints.iconame,1,&iconame);

    wm_hints.initial_state=NormalState;
    wm_hints.input=True;
    wm_hints.flags=StateHint | InputHint | IconPixmapHint;

    wm_hints.icon_pixmap=icon_pixmap;

    XSetWMProperties(display,win,&winname,&iconame,xarg,xnarg,&size_hints,&wm_hints,&class_hint);
  }

  XSelectInput(display,win,
	       ExposureMask | StructureNotifyMask
	       | KeyPressMask
	       | ButtonPressMask | ButtonReleaseMask | PointerMotionMask );

  {
    XGCValues gcvalues;

#if 1
    static char *turbo[16]={
      "#000000",
      "#000080",
      "#008000",
      "#008080",
      "#800000",
      "#800080",
      "#804000",
      "#c7c7c7", // gamma=2.2 between #ffffff and #707070
      "#707070",
      "#0000ff",
      "#00ff00",
      "#00ffff",
      "#ff0000",
      "#ff00ff",
      "#ffff00",
      "#ffffff",
};
#else /*# 1 */
    /* named colors do not work on some systems ! */
    static char *turbo[16]={
      "black",
      "darkblue",
      "darkgreen",
      "darkcyan",
      "darkred",
      "darkmagenta",
      "chocolate",
      "gray70",
      "gray40",
      "blue",
      "green",
      "cyan",
      "red",
      "magenta",
      "yellow",
      "white"};
#endif /*#!1 */
    int i;

    gcvalues.function=GXcopy;
    gc=XCreateGC(display,win,GCFunction,&gcvalues);

    if (no_colors==16) {
      loop (i,1,16)
	if (!XAllocNamedColor(display,cmap,turbo[i],&xcolor,&hwcolor))
	  fprintf(stderr,"draw/X11: %s is unknown color\n",turbo[i]);
	else xcoltab[i]=xcolor.pixel;
      xcoltab[0]=winbackground;
      setcolor(WHITE); } }

  /* this is required by show on older systems */
  XSetBackground(display,gc,winbackground);
  XMapWindow(display,win);

  /* wait until window is actually created, i.e., the first Expose */
  do {
    XNextEvent(display,&event);
#ifdef DEBUG
    fprintf(stderr,"startgraph: %s\n",
	    event.type==Expose?"Expose":"Other event");
#endif /*# DEBUG */
  } while (event.type!=Expose);

  return abs(mode);
}

void closegraph(void) /****************************************** closegraph */
/*
   graphics shut down ignored
   window remains created but inactive for future
   WARNING: colors wrong if changed palette, as in show
*/
{
  if (!display) return;
  setwritemode(0);
  setfillstyle(1,LIGHTRED);
  setcolor(YELLOW);
  selectfont(2);
  bar(0,0,53*xfont.width,xfont.height+4);
  settextstyle(0,0,1);
  settextjustify(0,2);
  outtextxy(7,5,"WINDOW INACTIVE, respond in the controlling terminal");
  setcolor(LIGHTRED);
  rectangle(0,0,getmaxx(),getmaxy());
  rectangle(1,1,getmaxx()-1,getmaxy()-1);
#if 0
  XFreeGC(display,gc);
  XCloseDisplay(display);
  display=NULL;
#endif /*# 0 */
  XFlush(display);
}

void setpal(int col,int r,int g,int b,int maxbright) /*************** setpal */
/* obsolete: maxbright<0 requests freeing colors, used only in showppm.c */
{
  int ret;

  if (!display) return;

  if (col<0 || col>no_colors)
    XERR("color number out of range",)
  else {
    if (maxbright<0) {
      if (xcoltab[col]) XFreeColors(display,cmap,xcoltab+col,1,0);
      maxbright=-maxbright; }
    if (r+g+b==0) xcoltab[col]=winbackground;
    else {
      xcolor.red=r*0x10000/maxbright;
      xcolor.green=g*0x10000/maxbright;
      xcolor.blue=b*0x10000/maxbright;
      ret=XAllocColor(display,cmap,&xcolor);
      xcoltab[col]=xcolor.pixel;
      if (ret==0) XERR("cannot allocate color",) } }
}

void setcolor(int color) /***************************************** setcolor */
{
  x_color=color;
}

void setgccolor(int color) /************************************* setgccolor */
{
  if (!display) return;

  if (color<0 || color>=no_colors)
    XERR("color",)
      if (currentcolor!=color) {
    if (xormode)
      XSetForeground(display,gc,xcoltab[currentcolor=color]^winbackground);
    else
      XSetForeground(display,gc,xcoltab[currentcolor=color]); }
}

void setfillstyle(int mode,int color) /************************ setfillstyle */
{
  /* only mode=1 (normal fill) supported */
  x_fillcolor=color;
}

void setlinestyle(int mode,int pattern,int thick)  /*********** setlinestyle */
/* extended 11/2016:
   mode=0  = normal solid
   mode>0  = pattern length, max 4; odd mode means that the pattern is repeated inverted
   pattern = number in the 256 base
   Examples:
     setlinestyle(1,0x02,1);    |  ##  ##  ##  ##  ##  ## ## ##  |
     setlinestyle(2,0x0103,1);  | ### ### ### ### ### ### ### ###|
     setlinestyle(4,0x02010203);|  #  ###  #  ###  #  ###  #  ###|
 */
{
  XGCValues gcvalues;
  unsigned long mask=GCLineStyle|GCLineWidth;

  if (!display) return;

  if (mode==0)
    gcvalues.line_style=LineSolid;
   else {
     mask|=GCDashList;
     gcvalues.line_style=LineOnOffDash; }

  gcvalues.line_width=thick;
  XChangeGC(display,gc,mask,&gcvalues);

  if (mode) XSetDashes(display,gc,0,(char*)&pattern,mode);
}

void setwritemode(int mode) /********************************** setwritemode */
{
  XGCValues gcvalues;

  if (!display) return;

  if (mode) gcvalues.function=GXxor;
  else gcvalues.function=GXcopy;
  XChangeGC(display,gc,GCFunction,&gcvalues);
  xormode=mode;
  currentcolor=-1;
}

void line(int x,int y,int xx,int yy)  /******************************** line */
{
  if (!display) return;

  setgccolor(x_color);
  XDrawLine(display,win,gc,x,y,xx,yy);
  /* was connected to the scaled graph and pen-like drawing (why?) - removed
     XDrawLine(display,win,gc,x,y,ggx=xx,ggy=yy);
  */
}

void bar(int x0,int y0,int x1,int y1) /********************************* bar */
{
  /* warning: good only if x1>x0, y1>y0 !!! */
  if (!display) return;

  setgccolor(x_fillcolor);
  XFillRectangle(display,win,gc,x0,y0,x1-x0+1,y1-y0+1);
}

void rectangle(int x0,int y0,int x1,int y1) /********************* rectangle */
{
  if (!display) return;

  setgccolor(x_color);
  XDrawRectangle(display,win,gc,min(x0,x1),min(y0,y1),abs(x1-x0),abs(y1-y0));
}

void putpixel(int x,int y,int color) /***************************** putpixel */
{
  if (!display) return;

  setgccolor(color);
  XDrawPoint(display,win,gc,x,y);
}

void fillellipse(int x,int y,int rx,int ry) /******************* fillellipse */
{
  if (!display) return;

  setgccolor(x_fillcolor);
  XFillArc(display,win,gc,x-rx,y-rx,rx*2,ry*2,0,360*64);
  setgccolor(x_color);
  XDrawArc(display,win,gc,x-rx,y-rx,rx*2,ry*2,0,360*64);
}

void circle(int x,int y,int r) /************************************* circle */
{
  if (!display) return;

  setgccolor(x_color);
  XDrawArc(display,win,gc,x-r,y-r,r*2,r*2,0,360*64);
}

/*****************************************************************************
                             multiple circles
******************************************************************************
(only the old simolant is using this)
  XArc *xcirc;
  xcirc=startcircles(n);
  loop (i,0,n) circles(c,x,y,r);
  setcolor(col); // etc.
  showcircles(xcirc);
  free(xcirc);
*/
int ncircles,icircles;
XArc *startcircles(int n)
{
  ncircles=n;
  icircles=0;
  return malloc(n*sizeof(XArc));
}

void circles(XArc *c,int x,int y,int r)
{
  if (icircles>=ncircles) return;
  c[icircles].x=x-r;
  c[icircles].y=y-r;
  c[icircles].width=r*2;
  c[icircles].height=r*2;
  c[icircles].angle1=0;
  c[icircles].angle2=360*64;
  icircles++;
}

void showcircles(XArc *c)
{
  if (!display) return;

  setgccolor(x_color);
  if (c) XDrawArcs(display,win,gc,c,ncircles);
}


/*****************************************************************************
        scaled graphics in x,y Cartesian coordinates, up+draw style
 *****************************************************************************/

/* extern: */
struct scaling_s scaling;
int pen;

static int ggx,ggy;
static int wx0,wx1=639,wy0,wy1=479; /* VGA 16col as default */

/* WARNING: cf. order of parameters in setviewport() */

void lwindow(int x0,int x1,int y0,int y1) /************************** window */
/* renamed because of window() in conio.h */
/* normally, x0<x1 and y0<y1: y0 is thus TOP and y1 BOTTOM */
{
  if (x1<=x0) x1=getmaxx();
  if (y1<=y0) y1=getmaxy();
  wx0=x0; wx1=x1;
  wy0=y0; wy1=y1;
}

void scale(REAL x0,REAL x1,REAL y0,REAL y1,int aspect) /************** scale */
{
  if (aspect) {
    scaling.x=(wx1-wx0)/(x1-x0);
    scaling.y=(wy1-wy0)/(y1-y0);
    scaling.x=scaling.y=fmin(scaling.x,scaling.y);
    scaling.fromx=x0-wx0/scaling.x;
    scaling.toy=y1+wy0/scaling.y; }
  else {
    scaling.x=(wx1-wx0)/(x1-x0);
    scaling.fromx=x0-wx0/scaling.x;
    scaling.y=(wy1-wy0)/(y1-y0);
    scaling.toy=y1+wy0/scaling.y; }
}

int SX(REAL x) /********************************************************* SX */
{
  int i=Int(scaling.x*(x-scaling.fromx)+16383.5)-16383;

  if (i<-32767) i=-32767;
  if (i>32767) i=32767;

  return i;
}

REAL iSX(int ix)
{
  return ix/scaling.x+scaling.fromx;
}

int SY(REAL y) /********************************************************* SY */
{
  int i=Int(scaling.y*(scaling.toy-y)+16383.5)-16383;

  if (i<-32767) i=-32767;
  if (i>32767) i=32767;

  return i;
}

REAL iSY(int iy)
{
  return scaling.toy-iy/scaling.y;
}

void draw(REAL x1,REAL y1) /******************************************* draw */
{
  int gx=SX(x1), gy=SY(y1);

  if (!display) return;

  setgccolor(x_color);
  if (pen) XDrawLine(display,win,gc,ggx,ggy,gx,gy);
  ggx=gx;
  ggy=gy;
  pen=1;
}

void up(void) /********************************************************** up */
{
  pen=0;
}

void lline(REAL x1,REAL y1,REAL x2,REAL y2) /************************* lline */
{
  up(); draw(x1,y1); draw(x2,y2);
}

void dot(REAL x1,REAL y1) /********************************************* dot */
{
  putpixel(SX(x1),SY(y1),getcolor());
}

/* radius of put pixel area:
   2pix/dot on average; min 0.707 to have 1 pix sure */
REAL pointsize=0.798;
int insidecolor=BLACK; /* for open points */

void point(REAL x1,REAL y1,char type) /******************************* point */
{
  REAL size=pointsize;
  int open=(unsigned)type<128,iopen;

  type=tolower(type&127);
  if (strchr("+x",type)) open=0;

  if (strchr("+dv^",type)) size*=1.414;

  if (size<=0)
    dot(x1,y1);
  else {
    REAL sx=scaling.x*(x1-scaling.fromx);
    REAL sy=scaling.y*(scaling.toy-y1);

    if (fabs(sx)<4096.0 && fabs(sy)<4096.0) loopto (iopen,0,open) {

      switch (type) {

        case 'x': case '*': {
          int
            x0=Int(sx+0.5-size),
            x9=Int(sx+0.5+size),
            y0=Int(sy+0.5-size),
            y9=Int(sy+0.5+size);

          line(x0,y0,x9,y9); line(x0,y9,x9,y0); }

          if (type=='x') break;
          else size*=1.414;

        case '-': case '+': {
          int
            x=Int(sx+0.5),
            y=Int(sy+0.5);

          if (strchr("*+-",type)) line(Int(sx+0.5-size),y,Int(sx+0.5+size),y);
          if (strchr("*+|",type)) line(x,Int(sy+0.5-size),x,Int(sy+0.5+size));
          break; }

        default: {
          int
            x =Int(sx+1-size),
            x9=Int(sx+1+size),
            y,
            y0=Int(sy+1-size),
            y9=Int(sy+1+size);
          REAL ss=Sqr(size),sv=size/2;
          unsigned col=iopen ? insidecolor : getcolor();

          for (; x<x9; x++)
            for (y=y0; y<y9; y++)
              if ((type=='o' && Sqr(x-sx)+Sqr(y-sy)<ss)
               || (type=='d' && fabs(x-sx)+fabs(y-sy)<size)
               || (type=='s' && fabs(x-sx)<size && fabs(y-sy)<size)
               || (type=='^' && y-sy<=sv && 1.7320508*fabs(x-sx)-(y-sy)<size)
               || (type=='v' && sy-y<=sv && 1.7320508*fabs(x-sx)-(sy-y)<size))
                putpixel(x,y,col);
          break; } }

      size=size*0.75-1; } }
}

void bigdot(REAL x1,REAL y1) /*************************************** bigdot */
{
  int i=SX(x1),j=SY(y1),c=getcolor();

  putpixel(i,j,c);
  putpixel(i+1,j,c); putpixel(i-1,j,c);
  putpixel(i,j+1,c); putpixel(i,j-1,c);
}

/*****************************************************************************
                          font and text interface
 *****************************************************************************/

static int textsize=1;

void settextstyle(int a,int b,int size)
{
  textsize=size;
}

int gettextsize(void)
{
  return textsize;
}

static int xjust,yjust=2;

void settextjustify(int x,int y)
{
  xjust=x; yjust=y;
}

/* To install a new font, consult pgm2xfont.c and generate a C-include file.
   Then you need to change:
   - NXFONTS in xdraw.h
   - below:
     #define XFONT xfont<NUMBER>
     #include "<NEWFONT>.c"
     #undef XFONT
   - in selectfont():
     xfonts[<NUMBER>]=xfont<NUMBER>;
   - check by: outtextx.c

NB: do not change numbering!
- even fonts used by MACSIMUS
- all fonts available in jkv
*/

#define XFONT xfont0
/* small font, former xfont7 */
#include "xfont0.c"
#undef XFONT

#define XFONT xfont1
/* 8x8, originally BGI */                                           
#include "xfbgi.c"
#undef XFONT

#define XFONT xfont2
/* medium, former xfont10 */
#include "xfont2.c"
#undef XFONT

#define XFONT xfont3
/* bold big, former xfont15.c */
#include "xfont3.c"
#undef XFONT

#define XFONT xfont4
/* new, biggest but not so bold */
#include "xfont4.c"
#undef XFONT

/* NEWFONT: add similar 3 lines as above */

                     /* width,xwidth,height,top,base,from,to,data */
struct xfont_s xfont = {8,8,9,0,7,0,256, xfont1+7}; /* default = BGI */
static int selectedfont;

void selectfont(int font) /************************************** selectfont */
{
  static unsigned *xfonts[NXFONTS];
  xfonts[0]=xfont0;
  xfonts[1]=xfont1;
  xfonts[2]=xfont2;
  xfonts[3]=xfont3;
  xfonts[4]=xfont4;
  /* NEWFONT: add a line here */

  if (font<0 || font>=NXFONTS) font=1;
  selectedfont=font;

  xfont.width=xfonts[font][0];
  xfont.xwidth=xfonts[font][1];
  xfont.height=xfonts[font][2];
  xfont.top=xfonts[font][3];
  xfont.base=xfonts[font][4];
  xfont.from=xfonts[font][5];
  xfont.to=xfonts[font][6];
  xfont.data=xfonts[font]+7;
}

static int maxlen=0x7fffff;

int outtextwidth(const char *s) /****************************** outtextwidth */
/* length in pixels returned */
{
  return strlen(s)*xfont.width*textsize;
}

int outtextppm(unsigned char (*rgb)[3],int xsize,int ysize, /*** outtext2ppm */
               unsigned char RGB[3],
               int x,int y,const char *s)
{
  unsigned X,Y,B,l=0,i=0;
  unsigned ch;

  while (*s) {
    ch=*(unsigned char*)s;

    if (ch>=xfont.from && ch<xfont.to) loop (X,0,xfont.xwidth) {
      B=xfont.data[(ch-xfont.from)*xfont.xwidth+X];

      loop (Y,0,xfont.height) {
        if (B&1) {
          int ii,jj;
          loop (ii,X*textsize,X*textsize+textsize)
            loop (jj,Y*textsize,Y*textsize+textsize)
              if (x+ii>=0 && x+ii<=xsize
               && y+jj>=0 && y+jj<=ysize) memcpy(rgb[x+ii+xsize*(y+jj)],RGB,3); }
            B>>=1; } }

    l+=xfont.width*textsize;
    x+=xfont.width*textsize;
    s++;
    if (i++>=maxlen) {
      maxlen=0x7fffff; break; }
  }

  return l;
}

int outtextxy(int x,int y,const char *s) /************************ outtextxy */
/* length in pixels returned */
{
  unsigned X,Y,B,l=0,i=0;
  unsigned ch;

  if (!display) return 0;

  setgccolor(getcolor());

  x-=xjust*(strlen(s)*xfont.width-1)*textsize/2;
  switch (yjust) {
    case 0: y-=xfont.base*textsize; break;
    case 1: y-=(xfont.base+xfont.top)/2*textsize; break;
    case 2: y-=xfont.top*textsize; }

  while (*s) {
    ch=*(unsigned char*)s;

    if (ch>=xfont.from && ch<xfont.to) loop (X,0,xfont.xwidth) {
      B=xfont.data[(ch-xfont.from)*xfont.xwidth+X];

      if (textsize==1)
        loop (Y,0,xfont.height) {
          if (B&1) XDrawPoint(display,win,gc,x+X,y+Y);
          B>>=1; }
      else
        loop (Y,0,xfont.height) {
          if (B&1) {
            int ii,jj;
            loop (ii,X*textsize,X*textsize+textsize)
              loop (jj,Y*textsize,Y*textsize+textsize)
                XDrawPoint(display,win,gc,x+ii,y+jj); }
          B>>=1; } }

    l+=xfont.width*textsize;
    x+=xfont.width*textsize;
    s++;
    if (i++>=maxlen) {
      maxlen=0x7fffff; break; }
  }

  return l;
}

int outtextxymax(int x,int y,const char *s,int ml) /*********** outtextxymax */
{
  maxlen=ml;
  return outtextxy(x,y,s);
}

int outtext(const char *s) /**************************************** outtext */
{
  return outtextxy(ggx,ggy,s);
}

int outbgtextxy(int x,int y,int font,int fg,int bg,
                const char *s)                                /* outbwtextxy */
/*
   output text as fg on bg using font, restore colors/font
   bottom y position returned
*/
{
  int oldf=selectedfont,oldc=x_color,oldbg=x_fillcolor,oldts=textsize;
  unsigned char *locs=(unsigned char*)strdup(s),*c;
  int aty=y+xfont.top,l=0,nl=0;
  unsigned char end;

  if (!display) return 0;

  setwritemode(0); /* PROBLEMATIC - should store and restore */
  selectfont(font);
  setcolor(fg);
  settextstyle(0,0,1);
  settextjustify(0,2);
  setfillstyle(1,bg);

  end=' '; c=locs;
  while (end) {
    unsigned char *e=c;

    while (!strchr("\n\r",*e)) e++;
    end=*e; if (l<e-c) l=e-c;
    c=e+1; }

  end=' '; c=locs;
  while (end) {
    unsigned char *e=c;

    while (!strchr("\n\r",*e)) e++;
    end=*e; *e=0;
    bar(x-1, aty-xfont.top-1, x+l*xfont.width+1, aty+xfont.height+2);
    outtextxy(x,aty,(char*)c); nl++;
    c=e+1;
    if (end!='\r') aty+=xfont.height+2; }

  XFlush(display);
  free(locs);
  selectfont(oldf);
  setcolor(oldc);
  setfillstyle(1,oldbg);
  settextstyle(0,0,oldts);

  return aty;
}

void fullscreentoggle(void) /****************************** fullscreentoggle */
{
  if (xwindowhints.fullscreen)
    xwindowhints.geometry=xwindowhints.lastgeometry;
  else {
    sprintf(xwindowhints.lastgeometry,"%dx%d",x_xsize,x_ysize);
    xwindowhints.geometry="f"; }
  XFreeGC(display,gc);
  XCloseDisplay(display);
  display=NULL;
  startgraph(-9);
}
