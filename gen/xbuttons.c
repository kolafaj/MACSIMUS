/* 
   buttons with text:
   - this code originates in the NSK project
   - later used in "old" (X11) simolant
   - now in MACSIMUS GUI
*/

/* best with font=2, well supported also 0,4 */
static struct buttons_s {
  struct buttons_s *next;
  int x0,x1,y0,y1;  // inside area (frame = +dy)
  unsigned key;     // returned key (unsigned char)
  const char *help; // help text
} *buttonshead = NULL;

static double buttontime;

/* slider, now only one is needed, to be extended to buttons.. */

static struct slider_s {
  int x0,x1,y0,y1;  // slider ideal range (x in [0,1])
  int X0,X1,Y0,Y1;  // slider detection range
  int yc;
  double *pos;
} slider;

#define XRG 8 // sensitive range around slider, x
#define YRG 5 // sensitive range around slider, x

#if 0
// 3D
void hline(int x,int y,int dx)
{
  setlinestyle(0,0,1);
  setcolor(dosmouse.button[3]);
  line(x,y,x+dx,y);
  line(x,y+1,x+dx,y+1);
  setcolor(dosmouse.button[1]);
  line(x,y+2,x+dx,y+2);
  setcolor(dosmouse.button[0]);
  line(x,y+3,x+dx,y+3);
}
#else
void hline(int x,int y,int dx)
{
  setlinestyle(0,0,1);
  setcolor(dosmouse.button[3]);
  line(x,y+1,x+dx,y+1);
}
#endif

int ATX[NATX];
void planline(int atx,char *L)
/* make a line of equally-spaced buttons:
   atx = left x position of the 1st button
   getmaxx() = left x position of the last button
   L = lengths of button texts, as signed char; 
       OR something PROPORTIONAL to the lengths (can create spaces, etc.)
 example 1:
   planline(maxxn+6,"\2\1\10")   
   makebutton(ATX[0],aty,'a'," a","sample button a\nhotkey=[a]");
   makebutton(ATX[1],aty,'b',"b","sample button b\nhotkey=[b]");
   makebutton(ATX[2],aty,'c',"button c","sample button c\nhotkey=[c]");
 example 2:
   planline(maxxn+6,"\3\2\1\10")   
   makebutton(ATX[1],aty,'a'," a","sample button a\nhotkey=[a]");
   makebutton(ATX[2],aty,'b',"b","sample button b\nhotkey=[b]");
   makebutton(ATX[3],aty,'c',"button c","sample button c\nhotkey=[c]");
*/
{
  int n=strlen(L),i,s=0,l;

  if (n>NATX) XERE("too many buttons in line");

  loop (i,0,n) s+=L[i];
  l=getmaxx()-atx-s*xfont.width-5; /* 5=right border */
  s=0;
  loop (i,0,n) {
    ATX[i]=atx+l*i/(n-1)+s*xfont.width;
    s+=L[i]; }
}

static void buttonframe(struct buttons_s *b,int pressed)
{
  setlinestyle(0,0,1);

  setcolor(pressed?dosmouse.button[1]:dosmouse.button[2]); // +-----
  line(b->x0-2,b->y0-2,b->x1+2,b->y0-2);                   // |
  line(b->x0-2,b->y0-1,b->x0-2,b->y1+2);                   // | ::::

  setcolor(pressed?dosmouse.button[0]:dosmouse.button[3]); //
  line(b->x0-1,b->y0-1,b->x1+1,b->y0-1);                   //  +----
  line(b->x0-1,b->y0  ,b->x0-1,b->y1+1);                   //  |::::

  setcolor(pressed?dosmouse.button[2]:dosmouse.button[0]); // :::: |
  line(b->x0-2,b->y1+2,b->x1+2,b->y1+2);                   //      |
  line(b->x1+2,b->y0-2,b->x1+2,b->y1+3);                   // -----+

  setcolor(pressed?dosmouse.button[3]:dosmouse.button[1]); // ::::|
  line(b->x0-1,b->y1+1,b->x1+1,b->y1+1);                   // ----+
  line(b->x1+1,b->y0-1,b->x1+1,b->y1+2);                   //
}

void voidbutton(int x,int y,const char *text)
{
  static struct buttons_s b;
  b.x0=x;
  b.x1=x+strlen(text)*xfont.width;
  b.y0=y-1;
  b.y1=y+xfont.height-3; /* not so high */
  setfillstyle(1,dosmouse.button[2]);
  bar(b.x0,b.y0,b.x1,b.y1);
  setcolor(dosmouse.button[0]);
  if (text[0]==' ' && text[1]) outtextxy(x+xfont.width/2,y,text+1);
  else outtextxy(x,y,text);
  buttonframe(&b,0);
}

int makebutton(int x,int y,unsigned key,const char *text,const char *help)
/*
   x,y is the TEXT start (to be compatible with outtextxy)
   note that some fonts have the first pixel column empty (except arrows)
   this was because the old SIMOLANT - to be changed?
   NEW: text starting by SPACE = shift left by half width
*/
{
  struct buttons_s *b;
  int l,hi,dy;

  setlinestyle(0,0,1);
  setwritemode(0);
  settextjustify(0,2);
  settextstyle(0,0,1);

  looplist(b,buttonshead) if (key==b->key) goto found;

  b=malloc(sizeof(struct buttons_s));
  if (!b) {
    fprintf(stderr,"no heap for button [%s]\n",text);
    sleep(1);
    return 0; }
  b->next=buttonshead; buttonshead=b;
  b->key=key;

 found:
  for (hi=0; text[hi]; hi++) if (text[hi]==key) break;

  //  fprintf(stderr,"hotkey in text=%s hi=%d key=%c\n",text,hi,key);

  l=strlen(text);
  b->x0=x-1;
  b->x1=x+l*xfont.width+1;
  //  b->y0=y-3;
  dy=xfont.height/5;
  if (dy<3) dy=3;
  b->y0=y-dy;
  b->y1=b->y0+xfont.height+2;
  b->help=help;
  setfillstyle(1,dosmouse.button[2]);
  bar(b->x0,b->y0,b->x1,b->y1);
  if (text[hi]) {
    /* underline hotkey */
    setfillstyle(1,dosmouse.button[4]);
    if (text[0]==' ' && text[1] && key!=' ') {
      /* for button text centered w. half-space around; do not underline space as hotkey */
      bar(x+(2*hi-1)*xfont.width/2+1,y+xfont.base+1,
          x+(2*hi+1)*xfont.width/2,  y+xfont.base+2+(xfont.height>14)); }
    else if (text[0]!=' ')
      /* for standard button text */
      bar(x+hi*xfont.width,    y+xfont.base+1,
          x+(hi+1)*xfont.width,y+xfont.base+2+(xfont.height>14));
  }
  setcolor(dosmouse.button[0]);
  if (text[0]==' ' && text[1])
    /* text begins with space: button text centered w. half-space around */
    outtextxy(x+xfont.width/2,y,text+1);
  else
    /* standard button text */
    outtextxy(x,y,text);

  buttonframe(b,0);

  return l*xfont.width+9;
}

void erasebuttons(void) /************************************** erasebuttons */
/*
   erase and free all buttons
*/
{
  struct buttons_s *b,*n;

  slider.pos=NULL;

  for (b=buttonshead; b; b=n) {
    n=b->next;
    free(b); }
  buttonshead=NULL;
}

static struct buttons_s *getbutton(void) /*********************** getbutton */
{
  struct buttons_s *b;

#ifdef DEBUG
  fprintf(stderr,"getbutton: %d %d\n",dosmouse.x,dosmouse.y);
#endif

  setwritemode(0);

  looplist (b,buttonshead)
    if (dosmouse.x>=b->x0-1 && dosmouse.x<=b->x1+1 && dosmouse.y>=b->y0-1 && dosmouse.y<=b->y1+1) {
      buttonframe(b,1);
      buttontime=mytime();

#ifdef DEBUG
      fprintf(stderr,"button %d=%c pressed\n",b->key,b->key);
      //      fprintf(stderr,"button %d=%c help=%s\n",b->key,b->key,b->help);
#endif
      return b;
    }

  return NULL;
}

void releasebutton(struct buttons_s *b) /********************* releasebutton */
/* redraw the button after ButtonRelease */
{

  double dt=mytime()-buttontime;

#ifdef DEBUG
  fprintf(stderr,"releasebutton: %d %d\n",dosmouse.x,dosmouse.y);
#endif
  mysleep(0.15-dt);
  setwritemode(0);
  buttonframe(b,0);
}

void drawslider(int forceredraw) /******************************* drawslider */
{
  int x,k;
  static double last=-1;

  if (!slider.pos) return;

  if (*slider.pos<0) *slider.pos=0;
  if (*slider.pos>1) *slider.pos=1;

  if (forceredraw || last!=*slider.pos) {

    x=slider.x0+*slider.pos*(slider.x1-slider.x0+0.9999);
    setlinestyle(0,0,1);
    setwritemode(0);

    setfillstyle(0,dosmouse.button[2]);
    bar(slider.x0-1,slider.y0,slider.x1+1,slider.y1);

    loopto (k,-2,2) {
      setcolor(dosmouse.button[(int)("\3\2\2\1"[k+2])]);
      line(slider.x0,slider.yc+k,slider.x1,slider.yc+k); }
    setlinestyle(0,0,3);
    setcolor(dosmouse.button[0]);
    line(x,slider.y0,x,slider.y1); }

  last=*slider.pos;
}

void makeslider(int x,int y,int len,double *pos) /*************** makeslider */
{
  slider.x0=x;
  slider.X0=x-XRG;
  slider.x1=x+len;
  slider.X1=slider.x1+XRG;
  slider.y0=y;
  slider.Y0=y-YRG;
  slider.y1=y+((xfont.height+1)/3)*2+1;
  slider.Y1=slider.y1+YRG;
  slider.yc=y+(xfont.height+1)/3;
  slider.pos=pos;

  drawslider(1);
}

int isinslider(void) /******************************************* isinslider */
{
  if (!slider.pos) return 0;
  if (slider.X0<0) return 0;
  if (dosmouse.x<slider.X0) return 0;
  if (dosmouse.x>slider.X1) return 0;
  if (dosmouse.y<slider.Y0) return 0;
  if (dosmouse.y>slider.Y1) return 0;
  *slider.pos=(double)(dosmouse.x-slider.x0)/(slider.x1-slider.x0);
  if (*slider.pos<0) *slider.pos=0;
  if (*slider.pos>1) *slider.pos=1;

  return 1;
}
