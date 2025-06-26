/* make plot
 */

int maxxn,maxyn;
int menu=1,redrawmenu=1;
int KILL=0; /* allow kill all plots, cf. getenv("GUI") */
int parent=-1; /* kill all spawn plots */

void setmenusize(void) /**************************************** setmenusize */
/* note that menu size 10*xfont.width+86 = 186 for font=2
   differs from show,blend because of 10 parameters
   thus, plot should not be used with very large font
*/
{
  maxxn=getmaxx()+1;
  if (menu) maxxn-=10*xfont.width+86-10*(font==0);
 // if (menu) maxxn-=15*xfont.width+36;

  maxyn=getmaxy()+1;
}

void makemenu(int which) /***************************************** makemenu */
{
  int atx,aty=0,dybutt,dyhdr1,dyhdr2,atx0,atxb0,i;

  if (!display) return;

  erasebuttons();
  setmenusize();

  setviewport(0,0,maxxn-1,maxyn-1,0);

  dybutt=xfont.height+10; // between rows of buttons
  dyhdr1=xfont.height+8; // before header
  dyhdr2=xfont.height+5; // after header
  atx0=maxxn+2; // text left pos
  atxb0=maxxn+6; // button text left pos

  setwritemode(0);
  settextjustify(0,2);
  setfillstyle(0,LIGHTGRAY);
  setlinestyle(0,0,1);
  if (which) bar(maxxn,0,getmaxx(),maxyn);

  aty=8;
  atx=atxb0;

  if (which==2) {
    setcolor(WHITE);
    atx=atx0;
    outtextxy(atx,aty,"Print screen:");
    atx=atxb0; aty+=dyhdr1;
    makebutton(atx,aty,'P',"PPM on black","Export screen in the Portable PixMap (P6) format\nwith black background\nhotkey=[P]");
    aty+=dybutt;
    makebutton(atx,aty,'p',"PPM on white","Export screen in the Portable PixMap (P6) format\nwith white background\nhotkey=[p]");
    aty+=dybutt;
    makebutton(atx,aty,'E',"EPS on black","Export screen in the EPS format\nwith black background\nutility ppm2ps must be installed\nhotkey=[E]");
    aty+=dybutt;
    makebutton(atx,aty,'e',"EPS on white","Export screen in the EPS format\nwith white background\nutility ppm2ps must be installed\nhotkey=[e]");
    aty+=dybutt;
    makebutton(atx,aty,'C',"PS on black","Export screen in the PS format, A4 paper\nwith black background\nutility ppm2ps must be installed\nhotkey=[C]");
    aty+=dybutt;
    makebutton(atx,aty,'c',"PS on white","Export screen in the PS format, A4 paper\nwith white background\nutility ppm2ps must be installed\nhotkey=[c]");
    aty+=dybutt;
    makebutton(atx,aty,'O',"gray PS on black","Export screen in gray PS format, A4 paper\nwith black background\nutility ppm2ps must be installed\nhotkey=[O]");
    aty+=dybutt;
    makebutton(atx,aty,'o',"gray PS on white","Export screen in gray PS format, A4 paper\nwith white background\nutility ppm2ps must be installed\nhotkey=[o]");
    aty+=dybutt;
    makebutton(getmaxx()-6*xfont.width-9,aty,ESC,"cancel","Cancel printscreen\nhotkey=[any other]");
    XFlush(display);
    return ; }

  planline(atxb0,"\4\4\4\1");
  if (menu) {
    setfillstyle(0,palcolor(CYAN));
    bar(maxxn,0,getmaxx(),9+xfont.height+xfont.top);

    makebutton(ATX[0],aty,F1,"help","get full help\nhotkeys=[?] [F1] [h]");
    makebutton(ATX[1],aty,'O'&31,"font","toggle font size\nhotkey=[F5] or [Ctrl-O]");
    makebutton(ATX[2],aty,'Q'&31,"quit","quit \'plot\'\nhotkey=[Q] or [Ctrl-Q]\nalso repeated ESC or [q]");
    makebutton(getmaxx()-xfont.width-5,aty,F10,"\22","menu on/off\nhotkey=[F10]");
  }
  else {
    /* menu toggle button - displaced from the corner, dark */
    dosmouse.button[1]=BLACK;
    dosmouse.button[2]=DARKGRAY;
    dosmouse.button[3]=LIGHTGRAY;
    makebutton(getmaxx()-xfont.width-7-font,aty+2+font,F10,"\22","menu on/off\nhotkey=[F10]");
    dosmouse.button[1]=DARKGRAY;
    dosmouse.button[2]=LIGHTGRAY;
    dosmouse.button[3]=WHITE;
    return; }

  aty+=3;
  if (KILL || parent>=0) {
    aty+=dybutt;
    makebutton(getmaxx()-8*xfont.width-5,aty,'K',parent>=0?"kill all":"KILL ALL","kill all running or spawn instances of \'plot\'\n[KILL ALL] = use command killall plot\n[kill all] = kill only plots spawned by PARENT\n             specified by option -p;\n             e.g., started by showcp or rdfg\nhotkey=[K]");
  }

  setcolor(WHITE);
  atx=atx0; aty+=dybutt;
  atx+=outtextxy(atx,aty,"print:")+3;
  atx+=makebutton(atx,aty,'#',"EPS","Export in the EPS (or PS) format\nuses ps.def to define format EPS/PS, labeling, etc.\nhotkey=[#], [F3], and quit=[F4]");
  atx+=makebutton(atx,aty,'P'&31,"PrtScr","PrintsSreen in various formats\nhotkey=[Ctrl-P] or [F2]");

  setcolor(WHITE);
  atx=atx0; aty+=dybutt;
  i=atx+=outtextxy(atx,aty,"zoom:")+3;
  atx+=makebutton(atx,aty,'X'&31,"cx", "Center around x=0\nhotkey=[Ctrl-X]\nUse keys [x] and [X] to zoom in/out,\nkeeping x=0 at center");
  makebutton(atx,aty,'Y'&31,"cy", "Center around y=0\nhotkey=[Ctrl-Y]\nUse mouse wheel or keys [y] and [Y]\nto zoom in/out, keeping y=0 at center");
  makebutton(getmaxx()-4*xfont.width-5,aty,'u',"undo","Undo the last zoom,\ncan repeat\nthere is no redo\nhotkey=[u] or [Ctrl-Z]");

  atx=i; aty+=dybutt;
  makebutton(atx,aty,'U'&31,"recalc","Recalculate scaling\nhotkey=[F9] or [Ctrl-U]");
  makebutton(getmaxx()-4*xfont.width-5,aty,'k',"init","Redraw with the initial scaling (no zoom)\nhotkey=[k]");

  atx=atx0; aty+=dybutt;
  setcolor(WHITE);

  i=atx+=outtextxy(atx,aty,"draw:")+3;
  atx+=makebutton(atx,aty,' ',"style","Cycle: draw by lines, points, lines-and-points\nActive only if the style is not specified in the command line\nhotkey=[SPACE]")+8;
  makebutton(atx+22-2*xfont.width,aty,'z',"grid","Cycle grid style: sparse/dense/none\nhotkey=[z]");

  atx=i; aty+=dybutt;
  atx+=makebutton(atx,aty,'.',".","Show by 1-pixel points (dots),\nactive only if the style is not specified in the command line\nhotkey=[.]")+2;
  atx+=makebutton(atx,aty,'p',"-","smaller points\nhotkey=[p]");
  atx+=makebutton(atx,aty,'P',"+","larger points\nhotkey=[P]");

  atx+=10;
  atx+=makebutton(atx,aty,'\\',"\\","Show by 1-pixel lines,\nactive only if the style is not specified in the command line\nhotkey=[\\]")+2;
  atx+=makebutton(atx,aty,'l',"-","thinner lines; thinnest=dotted\nhotkey=[l]");
  atx+=makebutton(atx,aty,'L',"+","thicker lines\nhotkey=[L]");

  atx=atx0; aty+=dybutt;
  setcolor(WHITE);
  if (perc) {
    atx+=outtextxy(atx,aty,"file:")+3;
    atx+=makebutton(atx,aty,'{',"-","Decrease file number:\nonly for file argds containing integer format\nhotkey=[}] or PgUp");
    atx+=makebutton(atx,aty,'}',"+","Increase file number:\nonly for file args containing integer format\nhotkey=[}] or PgUp");
    setcolor(RED); outtextxy(atx+3,aty,string("%d",percnr)); }
  else {
    atx+=outtextxy(atx,aty,"col:")+4;
    atx+=makebutton(atx,aty,'[',"-","Decrease Y-column number:\napplies to commands as \'plot FILE:1:3\' or \'plot FILE:A:C\'\nMore precisely: if the 1st character in the y-column is A..Z,\nit is decremented: plot -:A:D \32 plot -:A:C\nhotkey=[[]");
    atx+=makebutton(atx,aty,']',"+","Increase Y-column number:\napplies to commands as \'plot FILE:1:3\' or \'plot FILE:A:C\'\nMore precisely: if the 1st character in the y-column is A..Z,\nit is incremented: plot -:A:B \32 plot -:A:C\nhotkey=[]]");
    atx+=makebutton(atx,aty,'U',"Un","Reset column to the original value\nhotkey=[U]");
    if (column) { setcolor(RED); outtextxy(atx+3,aty,string("%c",column)); } }
  makebutton(getmaxx()-6*xfont.width-5,aty,'^',"errbar","Toggle showing error bars (if given)\nhotkeys=[^] [F11]");
#if PARM
  if (aty>maxyn-6*dyhdr1-3*dybutt) goto playback;
  atx=atx0; aty+=dyhdr2;
  hline(atx0,aty,getmaxx()-maxxn-5);
  aty+=12;
  setcolor(WHITE);
  atx+=outtextxy(atx,aty,"fit:")+3;
  planline(atx,"\3\4\4");
  makebutton(ATX[0],aty,'t',"std","Fit data to function, long double:\nrecommended precision\nhotkey=[t]");
  makebutton(ATX[1],aty,'T',"high","Fit data to function, emulated high-precision:\nslow, use in case of numerical problems\nhotkey=[T]");
  makebutton(ATX[2],aty,'T'&31,"fast","Fit data to function, double precision:\nonly a bit faster than long double\nhotkey=[Ctrl-T]");
  setcolor(WHITE);
  atx=atx0; aty+=dybutt;
  atx+=outtextxy(atx,aty,"nerr:")+3;
  planline(atx,"\1\2\3\3");
  makebutton(ATX[0],aty,136,"0","Do not calculate error estimates\nhotkey=[n](repeated)");
  makebutton(ATX[1],aty,137,"10","Use 10 sampling points to calculate error estimates\nalso hotkeys=[n] and [N]");
  makebutton(ATX[2],aty,138,"100","Use 100 sampling points to calculate error estimates\nalso hotkeys=[n] and [N]");
  makebutton(ATX[3],aty,139,"999","Use 999 sampling points to calculate error estimates\nalso hotkeys=[n] and [N]");

  atx=atx0; aty+=dyhdr1;
  setcolor(WHITE);
  atx+=outtextxy(atx,aty,"--- parameters ---");
  atx=atxb0; aty+=dyhdr2;
  {
    static char *help[PARM];
    int d=1+(font==0);

    loop (i,0,PARM) {
      if (!help[i]) help[i]=strdup(string("Set/change parameter %c:\nalso hot keys [%c] [%c] [Ctrl-%c]",i+'a',i+'a',i+'A',i+'A'));
      dosmouse.button[0]=parm[i].used?BLACK:DARKGRAY;
      if (i==lastparm) dosmouse.button[2]=LIGHTCYAN;
      atx+=makebutton(atx,aty,141+i,
                      string("%c",i+'a'),
                      help[i])-d+(i==4)*2; 
      dosmouse.button[2]=LIGHTGRAY;
    } }
  dosmouse.button[0]=BLACK;

  if (lastparm>=0) {
    aty+=dybutt;
    planline(atxb0,"\1\1\1\1\1\2\4");
    makebutton(ATX[0],aty,'<',"<","Decrease value of the selected parameter by step and\nrecenter the slider to [value-step,value+step]\nhotkey=[<] or [a],[b]..");
    makebutton(ATX[1]-1,aty,'>',">","Increase value of the selected parameter by step and\nrecenter the slider to [value-step,value+step]\nhotkey=[>] or [A],[B]..");
    makebutton(ATX[2]+1,aty,'*',"*","Double the step and slide range and\nrecenter the slider to [value-step,value+step]\nhotkey=[*]");
    makebutton(ATX[3],aty,'|',"|","Recenter the slider to [value-step,value+step]\nhotkey=[|]");
    makebutton(ATX[4]-1,aty,'/',"/","Halve the step and slide range and\nrecenter the slider to [value-step,value+step]\nhotkey=[/]");
    makebutton(ATX[5]+1,aty,'o',"=0","Set selected parameter to 0\nhotkey=[o] or [Ctrl-A], [Ctrl-B]..");
    makebutton(ATX[6],aty,F8,"all0","Set all parameters to 0\nand remove them from the list\nhotkey=[F8]");

    atx=atxb0+2; aty+=dyhdr1;
    sliderpos=(*parm[lastparm].parm-parm[lastparm].min)/(parm[lastparm].max-parm[lastparm].min);
    makeslider(atx+2,aty,getmaxx()-maxxn-20,&sliderpos);
  }
#endif

 playback:
  XFlush(display);
}
