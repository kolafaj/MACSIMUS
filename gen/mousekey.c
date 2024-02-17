/*
  simple X11/BGI interface to mouse and hot keys 
  (see also drawplus.c/drawplus.h used by nsk and old version of plot)

  extended codes
  Cursors:   200 ^  208 v   205 ->  203 <-  199 Home  207 End
  mouse buttons DRAG:  left=128  middle=130  right=129
  mouse buttons CLICK: left=131  middle=133  right=132
  127=backspace
  (see drawplus.h for more codes)
*/
#define LEFTCLICK 131
#define RIGHTCLICK 132
#define MIDCLICK 133
#define LEFTDRAG 128
#define RIGHTDRAG 129
#define MIDDRAG 130
#define WHEELFORWARD 134
#define WHEELBACKWARD 135

static int mousedx,mousedy,xchar;

#ifdef DOS
#include <conio.h>
/* PATCH -- UNIFY WITH X11 !!! */
struct dosmouse_s { int x,y; } dosmouse;
void readxymouse(void)
{
  dosmouse.x=getMouseX();
  dosmouse.y=getMouseY();
}
#endif

#ifdef EMX
/* patch -- no mouse in EMX */
struct dosmouse_s { int x,y; } dosmouse;
void readxymouse(void)
{
  dosmouse.x=getMouseX();
  dosmouse.y=getMouseY();
}
#endif


int xkbhit(int mouse) /********************************************** xkbhit */
/* 
   returns 1 on key/mouse button hit, char pressed = xchar
*/
{
  static enum buttonState last[3]={buttonUp,buttonUp,buttonUp};
  static int x0,x00=-32000,y0,y00;
  enum buttonState butt;
  int i;

  if (xchar) return 1;

#ifndef EMX
  if (mouse) {

    hit(0); 
    if (dosmouse.wheel) {
      xchar=WHEELFORWARD+(dosmouse.wheel<0);
      /* this is dirty ... */
      dosmouse.wheel=0;
      return 1; }

    loop (i,0,3) {

    tryagain:
    
      butt=getButton(1<<i);
#ifdef DOS
      readxymouse();
#endif

      if (butt==buttonUp) dosmouse.ex=0; /* algorithm problem: LEFTDRAG only */

      if (last[i]==buttonUp)
	if (butt==buttonDown) {
	  x00=x0=dosmouse.x; y00=y0=dosmouse.y; }
	else ;

      else /* last=buttonDown */
	if (butt==buttonDown)
	  if (dosmouse.x==x0 && dosmouse.y==y0) goto tryagain;
	  else {
	    x00=-32000;
	    mousedx=dosmouse.x-x0; x0=dosmouse.x;
	    mousedy=dosmouse.y-y0; y0=dosmouse.y; xchar=128+i; return 1; }
	else {
	  if (dosmouse.x==x00 && dosmouse.y==y00) {
	    last[i]=butt;
	    xchar=131+i;
	    return 1; } }
      
      last[i]=butt; } }
#endif

  if (kbhit()) return 1;
  else return 0;
}

unsigned char extch(void) /******************************************* extch */
/*
  extended codes | 128
*/
{
  unsigned char c;

  while (!xkbhit(1));
  if (xchar) { int i=xchar; xchar=0; return i; }

  c=getch();
#ifdef DEBUG
  fprintf(stderr,"getch=%d=0x%x\n",c,c);
#endif
  if (c==0) return getch()|128;
  else return c;
}

