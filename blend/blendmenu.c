/* make blend
 */

int maxxn,maxyn;

void setmenusize(void) /**************************************** setmenusize */
{
  maxxn=getmaxx()+1;
  if (menu) maxxn-=15*xfont.width+36;
  maxyn=getmaxy()+1;
}

void makemenu(int which) /***************************************** makemenu */
{
  int atx,aty=0,dybutt,dyhdr1,dyhdr2,atx0,atxb0;
  static char *info[5]={ "?", ": wire",":standard",": shaded",":anaglyph" };

  erasebuttons();
  //  if (!menu) return;

  setviewport(0,0,maxxn-1,maxyn-1,0);

  /* repeated: */
  selectfont(font);
  setmenusize();

  dybutt=xfont.height+10; // between rows of buttons
  dyhdr1=xfont.height+8; // before header
  dyhdr2=xfont.height+5; // after header
  atx0=maxxn+3; // text left pos
  atxb0=maxxn+6; // button text left pos

  setwritemode(0);
  settextjustify(0,2);
  setlinestyle(0,0,1);
  setfillstyle(0,LIGHTGRAY);
  if (which) bar(maxxn,0,getmaxx(),maxyn);

  aty=8;
  atx=atxb0;

  if (which==2) {
    setcolor(WHITE);
    atx=atx0;
    outtextxy(atx,aty,"Print screen:");
    atx=atxb0; aty+=dyhdr1;
    makebutton(atx,aty,'p',"PPM on white","export screen in the Portable PixMap (P6) format\nwith white background\nhotkey=[p]");
    aty+=dybutt;
    makebutton(atx,aty,'P',"PPM on black","export screen in the Portable PixMap (P6) format\nwith black background\nhotkey=[P]");
    aty+=dybutt;
    makebutton(atx,aty,'e',"EPS on white","export screen in the EPS format\nwith white background\nutility ppm2ps must be installed\nhotkey=[e]");
    aty+=dybutt;
    makebutton(atx,aty,'E',"EPS on black","export screen in the EPS format\nwith black background\nutility ppm2ps must be installed\nhotkey=[E]");
    aty+=dybutt;
    makebutton(atx,aty,'c',"PS on white","export screen in the PS format, A4 paper\nwith white background\nutility ppm2ps must be installed\nhotkey=[c]");
    aty+=dybutt;
    makebutton(atx,aty,'C',"PS on black","export screen in the PS format, A4 paper\nwith black background\nutility ppm2ps must be installed\nhotkey=[C]");
    aty+=dybutt;
    makebutton(atx,aty,'o',"gray PS on white","export screen in gray PS format, A4 paper\nwith white background\nutility ppm2ps must be installed\nhotkey=[o]");
    aty+=dybutt;
    makebutton(atx,aty,'O',"gray PS on black","export screen in gray PS format, A4 paper\nwith black background\nutility ppm2ps must be installed\nhotkey=[O]");
    aty+=dybutt;
    makebutton(getmaxx()-6*xfont.width-7,aty,ESC,"cancel","cancel printscreen\nhotkey=[any other]");
    XFlush(display);
    return ; }

  planline(atxb0,"\4\4\4\1");
  if (menu) {
    setfillstyle(0,palcolor(CYAN));
    bar(maxxn,0,getmaxx(),9+xfont.height+xfont.top);

    makebutton(ATX[0],aty,F1,"help","get full help\nhotkeys=[?] [F1] [h]");
    makebutton(ATX[1],aty,'O'&31,"font","toggle font size\nhotkey=[F5] or [Ctrl-O]");
    makebutton(ATX[2],aty,'Q'&31,"quit","quit immediately without saving\nhotkey=[Ctrl-Q]");
  }
  makebutton(getmaxx()-xfont.width-5,aty,F10,"\22","menu on/off\nhotkey=[F10]");
  if (!menu) return;

  setcolor(WHITE);
  aty+=dyhdr1;
  atx=atx0;
  outtextxy(atx,aty,string("Draw mode%s",info[grmode]));
  aty+=dyhdr2;
  planline(atxb0,"\1\1\1\1\4\4");

  dosmouse.button[2]=grmode==1?LIGHTCYAN:LIGHTGRAY;
  makebutton(ATX[0],aty,'1',"1","wire: show bonds using lines\nnonbonded atoms using crosses\nhotkey=[1]");
  dosmouse.button[2]=grmode==2?LIGHTCYAN:LIGHTGRAY;
  makebutton(ATX[1],aty,'2',"2","standard: show bonds using split-color lines\natoms using circles\nhotkey=[2]");
  dosmouse.button[2]=grmode==3?LIGHTCYAN:LIGHTGRAY;
  makebutton(ATX[2],aty,'3',"3","shaded: show bonds using split-color lines\natoms using shaded circles\nhotkey=[3]");
  dosmouse.button[2]=grmode==4?LIGHTCYAN:LIGHTGRAY;
  makebutton(ATX[3],aty,'4',"4","anaglyph\nglasses needed\nhotkey=[4]");
  dosmouse.button[2]=LIGHTGRAY;
  makebutton(ATX[4],aty,'b',"bond","cycle through several bond styles\nhotkey=[b]\nalso [B]=backwards");
  //  atx+=6;
  //  atx+=makebutton(atx,aty,'c',"C","carbon color\nhotkey=[c]");
  makebutton(ATX[5],aty,'#',"grid","show grid by 1 AA / 0.1 AA / none\nhot key=[=] or [#]");

  atx=atx0; aty+=dybutt;
  setcolor(WHITE);
  atx+=outtextxy(atx,aty,"ball size:")+3;
  atx+=makebutton(atx,aty,'r',"-","decrease ball size\nhotkey=[r]");
  atx+=makebutton(atx,aty,'R',"+","increase ball size\nhotkey=[R]");

  atx=atx0; aty+=dybutt;
  setcolor(WHITE);
  atx+=outtextxy(atx,aty,"perspective:")+3;
  atx+=makebutton(atx,aty,HOME,"-","perspective: move the view point closer to the screen\nhot key=[Home]=[ctrl-a]");
  atx+=makebutton(atx,aty,END,"+","perspective: move the view point closer to the screen\nhot key=[End]=[ctrl-e]");

  atx=atx0; aty+=dyhdr1;
  setcolor(WHITE);
  atx+=outtextxy(atx,aty,"--- label ---")+3; // 18 letters
  aty+=dyhdr2;
  planline(atxb0,"\2\4\1\2\4");
  makebutton(ATX[0],aty,'i',"id","show atom ID (as in the mol-file)\nhotkey=[i]");
  makebutton(ATX[1],aty,'t',"type","show atom TYPE (as in force-field)\nhotkey=[t]");
  makebutton(ATX[2],aty,'q',"q","show atom charge in e\nhotkey=[q]");
  makebutton(ATX[3],aty,'n',"no","show atom number (numbered from 0)\nhotkey=[n]");
  makebutton(ATX[4],aty,' ',"hide","hide labeling (redraw with no labels)\nhotkey=[q]");

  atx=atxb0; aty+=dybutt;
  planline(atxb0,"\4\4\6");
  makebutton(ATX[0],aty,'f',"find","find atom(s) matching string or number\nthe string or number is given in the terminal\napplies to one of above labelings\nfinding atom of given number is the same as clicking (marking) it\nhotkey=[f]");
  makebutton(ATX[1],aty,'L',"mark","mark (select) all found atoms\nhotkey=[L]");
  makebutton(ATX[2],aty,'F',"unfind","clear the finding string/number\nhotkey=[F]");
#ifdef POLAR
  aty+=dybutt;
  planline(atxb0,"\5\3\3");
  //  atx+=makebutton(atx,aty,'D',"ind","show induced dipoles in eAA (AA=Angstrom)");
  makebutton(ATX[0],aty,'d',"ind/D","show induced dipoles in Debye");
  makebutton(ATX[1],aty,'a',"pol","show polarizability volumes in AA^3)");
  makebutton(ATX[2],aty,'s',"sat","show saturation (some models only)");
#endif /*# POLAR */

  atx=atx0; aty+=dyhdr1;
  setcolor(WHITE);
  outtextxy(atx,aty,"--- selection ---"); // 18 letters
  aty+=dyhdr2;
  planline(atxb0,"\4\6\3");
  dosmouse.button[2]=moving?LIGHTCYAN:LIGHTGRAY;
  makebutton(ATX[0],aty,'m',"move","Toggle whether the subsequent move/rotation will\naffect the whole configuration or the SELECTION only.\n  move=mid mouse\n  x,y-rot=left mouse\n  z-rot=right mouse\nhotkey=[m]");
  dosmouse.button[2]=LIGHTGRAY;
  makebutton(ATX[1],aty,'U'&31,"unmark","Unmark all atoms\nhotkey=[ctrl-u]");
  makebutton(ATX[2],aty,'I',"inv","Invert marking\nhotkey=[I]");

  atx=atx0; aty+=dyhdr1;
  setcolor(WHITE);
  outtextxy(atx,aty,"--- optimize ---"); // 18 letters
  aty+=dyhdr2;
  planline(atxb0,"\2\2\3\1\5");
  if (constraints) {
    makebutton(ATX[0],aty,',',"MC","Minimize energy by Monte Carlo:\nuse in standard cases\nhotkey=[,]");
    makebutton(ATX[1],aty,';',"SD","Minimize energy by steepest descent\nhotkey=[;]");
    makebutton(ATX[2],aty,':',"ran","Randomize + minimize (MC)\nuse to jump from local minima\nhotkey=[:]"); }
  else {
    makebutton(ATX[0],aty,',',"CG","Minimize energy by conjugate gradients:\nuse in standard cases\nhotkey=[,]");
    makebutton(ATX[1],aty,';',"SD","Minimize energy by steepest descent:\nuse in difficult cases (overlaps etc.)\nhotkey=[;]");
    makebutton(ATX[2],aty,':',"ran","Randomize + minimize (SD+CG):\nuse to jump from local minima\nhotkey=[:]");
  }
  makebutton(ATX[3],aty,'[',"-","Halve the number of minimization steps\nhotkey[[]");
  makebutton(ATX[4],aty,']',"+","Double the number of minimization steps\nhotkey[]]");
  aty_n=aty;

  if (aty<maxyn-3*dybutt) {
    aty+=dybutt;
    planline(atxb0,"\5\3\3\3");
    makebutton(ATX[0],aty,'e',"trace","Toggle trace mode while minimization\nhotkey=[e]");
    makebutton(ATX[1],aty,'j',"by","Show minimization by given number of steps\nhotkey=[j]");
    makebutton(ATX[2],aty,'W',"put","Remember orientation and position:\napplies to the whole configuration,\nnot individual parts moved.\nWithout this key, the initial state is remembered\nhotkey=[w]");
    makebutton(ATX[3],aty,'w',"get","Restore written (or initial) orientation:\nuseful as edited input configuration for cook\nin periodic boundary conditions\nhotkey=[W]"); }
  if (aty<maxyn-4*dybutt) {
    atx=atxb0; aty+=dybutt;
    makebutton(atx+xfont.width+3,aty,'k',"keep","marked atoms will be kept fixed\nduring minimization\nhotkey=[k]");
    makebutton(getmaxx()-10*xfont.width-9,aty,'K',"mark Kept","mark kept atoms (again or kept by * in the mol-file)"); }

  aty+=dybutt;
  planline(atxb0,"\4\4\6");
  makebutton(ATX[0],aty,'S'&31,"save","Save the edited configuration now\nhotkey=[ctrl-s]");
  makebutton(ATX[1],aty,'E',"info","Print energy, dipole moment, and quadrupole moment\nto yhe controlling console\nhotkey=[E]");
  atx+=makebutton(ATX[2],aty,'.',"finish","End minimization and save, then continue\nby the next molecule (or stop if none)\nhotkey=[.]");

  atx=atx0; aty+=dyhdr1;
  if (aty>=maxyn) goto ret;
  setcolor(WHITE);
  atx+=outtextxy(atx,aty,"--- dump ---")+3; // 18 letters
  atx=atxb0; aty+=dyhdr2;
  atx+=makebutton(atx,aty,'P',"PrtScr","printscreen in various formats\nhotkey=[p]");

 ret:
  XFlush(display);
}
