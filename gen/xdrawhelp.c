static int athelp,dathelp;

#define SHOWBUTT 64
#ifndef ORANGE
#  define ORANGE YELLOW
#endif

void H(char *k,char *msg)
/*
  H("YELLOW","");
  H("=== YELLOW ===","WHITE");
  H("GREEN;BUTT","WHITE");
  H(NULL,      "      ORANGE"); (no linefeed - ORANGE will be inside next line)
  H("",       "white        white");
*/

{
  if (!k) {
    /* yellow/orange word before next line */
    setcolor(palcolor(ORANGE)); outtextxy(xfont.width*12,athelp,msg);
    return; }
  if (*msg) {
    if (memcmp(k,"=== ",4)) {
      //   H("GREEN;BUTT","WHITE");
      char cutk[SHOWBUTT],*text;
      if ( (text=strchr(k+1,';')) && strlen(k)<SHOWBUTT && text[1]) {
        memcpy(cutk,k,SHOWBUTT);
        *(strchr(cutk+1,';'))=0;
        voidbutton(6+xfont.width*(strlen(cutk)+1),athelp,text+1);
        setcolor(palcolor(LIGHTGREEN)); outtextxy(4,athelp,cutk); }
      else {
        setcolor(palcolor(LIGHTGREEN)); outtextxy(4,athelp,k); }
      setcolor(palcolor(WHITE)); outtextxy(xfont.width*12,athelp,msg); }
    else {
      // H("=== YELLOW ===","WHITE");
      setcolor(palcolor(YELLOW)); outtextxy(4,athelp,k);
      setcolor(palcolor(WHITE)); outtextxy(xfont.width*12,athelp,msg); } }
  else {
    // H("YELLOW","");
    setcolor(palcolor(YELLOW)); outtextxy(4,athelp,k); }

  athelp+=dathelp;
}

int nextpage(void)
{
  athelp+=xfont.height/2;
  setcolor(palcolor(LIGHTCYAN));
  outtextxy(4,athelp,"ESC=exit help   scroll by wheel or standard hotkeys");
  athelp=5;

  /* NB: setmenusize unknown here */
  setviewport(0,0,getmaxx(),getmaxy(),0);

  erasebuttons();
  makebutton(getmaxx()-xfont.width*9-3,5,F1,"exit help","return back to main menu\nhotkey=[F1] or [ESC]..");

  switch (readkey()) {
    case 'p': case 'u': case 8: case PGUP: case WHEELFORWARD: case LEFTCLICK: case LEFT: case UP:
      clearviewport();
      return -2;
    case 'n': case 'd': case ' ': case PGDN: case WHEELBACKWARD: case RIGHTCLICK: case RIGHT: case DOWN:
      clearviewport();
      return 1;
    case 'Q'&31: case 'Q':
      closegraph();
      exit(0);
    case 'q': case F1:
    case ESC: case MIDCLICK:
      clearviewport();
      return 0;
    default: clearviewport();
      return -1; }
}

#define NEXTPAGE(LUP,THIS) switch (nextpage()) { case -2: goto LUP; case 0: return; case -1: goto THIS; }
