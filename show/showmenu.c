/* make show
 */

void setmenusize(void) /**************************************** setmenusize */
{
  maxx=getmaxx();
  if (menu) maxx-=15*xfont.width+36;
  maxxn=maxx+1;
  maxy=getmaxy();
  maxyn=maxy+1;
}

void makemenu(int which) /***************************************** makemenu */
{
  int atx,aty=0,dybutt,dyhdr1,dyhdr2,atx0,atxb0;
  static char *info[8]={ "?", "bar","bar2","barball","barball2","ball", "wire","wire2" };
  char rg[4]="rg";

  erasebuttons();
  setmenusize();

  dybutt=xfont.height+10; // between rows of buttons
  dyhdr1=xfont.height+8; // before header
  dyhdr2=xfont.height+5; // after header
  atx0=maxxn+3; // text left pos
  atxb0=maxxn+6; // button text left pos

  setwritemode(0);
  settextjustify(0,2);
  setlinestyle(0,0,1);
  setfillstyle(0,palcolor(LIGHTGRAY));
  if (which) bar(maxxn,0,getmaxx(),maxy-xfont.height);

  aty=8;
  atx=atxb0;

  if (which==2) {
    setcolor(palcolor(WHITE));
    atx=atx0;
    outtextxy(atx,aty,string("Export %s:",dumpinfo[dumpmode]));
    atx=atxb0; aty+=dyhdr1;
    makebutton(atx,aty,'o',"one frame","export one frame in the selected format\nhotkey=[o]");
    dosmouse.button[2]=palcolor(LIGHTGRAY);
    aty+=dybutt;
    makebutton(atx,aty,'O',"One frame+render","export one frame in the selected format\nand render or show it with\nthe software registered in .startdata\nutility \'start\' must be installed\nhotkey=[O]");
    aty+=dybutt;
    makebutton(atx,aty,'s',"series","export a whole series of frames\nuntil hotkey [e] is pressed\nor button [end series] pushed\nhotkey=[s]");
    if (dumpmode==PIC) {
      aty+=dybutt;
      makebutton(atx,aty,'a',"animated GIF","make animated GIF\nall PPM files are erased before conversion\nsee option -L for parameters\nuntil hotkey [e] is pressed\nor button [end series] pushed\nhotkey=[a]"); }

    if (ISMERGED) {
      aty+=dybutt;
      makebutton(atx,aty,'m',"merge","merge a series of frames to one output file\nuntil hotkey [e] is pressed\nor button [end series] pushed\nhotkey=[m]");
      aty+=dybutt;
      makebutton(atx,aty,'M',"Merge + render","merge a series of frames to one output file\nuntil hotkey [e] is pressed\nor button [end series] pushed\nhotkey=[M]"); }

    aty+=dybutt;
    makebutton(getmaxx()-6*xfont.width-9,aty,ESC,"cancel","cancel the export\nhotkey=[any other]");
    setcolor(palcolor(WHITE));
    aty+=xfont.height+5; outtextxy(atx,aty,"The render option");
    aty+=xfont.height+1; outtextxy(atx,aty,"requires \'start\'");
    aty+=xfont.height+1; outtextxy(atx,aty,"installed.");
    aty+=xfont.height/2;
    aty+=xfont.height+1; outtextxy(atx,aty,"E.g., NFF and");
    aty+=xfont.height+1; outtextxy(atx,aty,"one frame+render");
    aty+=xfont.height+1; outtextxy(atx,aty,"will call a ray-");
    aty+=xfont.height+1; outtextxy(atx,aty,"tracer and then");
    aty+=xfont.height+1; outtextxy(atx,aty,"ImageMagick.");
    aty+=xfont.height/2;
    aty+=xfont.height+1; outtextxy(atx,aty,"Background:");
    aty+=xfont.height+1; outtextxy(atx,aty,"show -bgRRGGBB");

    if (dumpmode==ZBUF) {
      int res=96;
      char *e=getenv("RES");
      double d;

      aty+=dyhdr1;
      if (e) res=atoi(e);
      setcolor(BLACK);

      outtextxy(atx,aty,"Check the scales:");
      aty+=xfont.height+1;

      for (atx=atx0; atx<getmaxx(); atx+=res)
        line(atx,aty-3,atx,aty+xfont.height-3);
      settextjustify(1,2);
      outtextxy(atx0+res/2,aty,"1\"");

      aty+=xfont.height+1;
      for (d=atx0; d<getmaxx(); d+=res/2.54) {
        atx=d;
        setcolor(BLACK);
        if (d-(int)d>0.5) {
          atx=d+1;
          setcolor(DARKGRAY);
          line(atx,aty-3,atx,aty+xfont.height-3);
          atx--; }
        line(atx,aty-3,atx,aty+xfont.height-3); }

      setcolor(BLACK);
      settextjustify(1,2);
      outtextxy(atx0+res/5.1,aty,"1cm");

      settextjustify(0,2);

      atx=atx0;
      aty+=xfont.height+1; outtextxy(atx,aty,string("Current RES=%d,",res));
      aty+=xfont.height+1; outtextxy(atx,aty,"change environment");
      aty+=xfont.height+1; outtextxy(atx,aty,"if needed");
    }

  return ; }

  planline(atxb0,"\4\4\4\1");
  if (menu) {
    setfillstyle(0,palcolor(CYAN));
    bar(maxxn,0,getmaxx(),9+xfont.height+xfont.top);

    makebutton(ATX[0],aty,F1,"help","get full help\nhotkeys=[?] [F1] [h]");
    makebutton(ATX[1],aty,'O'&31,"font","toggle font size\nhotkey=[F6] or [Ctrl-O]");
    makebutton(ATX[2],aty,'Q'&31,"quit","quit immediately without saving\nhotkey=[Ctrl-Q]"); }

  makebutton(getmaxx()-xfont.width-5,aty,F10,"\22","menu on/off\nhotkey=[F10]");
  if (!menu) return;

  setcolor(palcolor(WHITE));
  aty+=dyhdr1;
  outtextxy(atx0,aty,string("draw mode:%s",info[curmode]));
  aty+=dyhdr2;
  planline(atxb0+xfont.width*2+15,"\1\1\1\1\1\1\1");
  dosmouse.button[2]=curmode==BAR?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[0]-15,aty,'!',"1","show bonds using white bars\nnonbonded atoms using small balls\nhotkey=[!]");
  dosmouse.button[2]=curmode==BAR2?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[1]-15,aty,'@',"2","show bonds using split-color bars\nnonbonded atoms using small balls\nhotkey=[@]");
  dosmouse.button[2]=curmode==DUMBELL?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[2]-10,aty,'#',"3","show bonds using white bars\natoms using small balls\nhotkey=[#]");
  dosmouse.button[2]=curmode==DUMBELL2?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[3]-10,aty,'$',"4","show bonds using split-color bars\natoms using small balls\nhotkey=[$]");
  dosmouse.button[2]=curmode==BALL?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[4]-5,aty,'%',"5","show atoms using balls\ndefault size=70% van der Waals\nuse r R to change\nhotkey=[%]");
  dosmouse.button[2]=curmode==BOND?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[5],aty,'^',"6","show bonds using white lines\nnonbonded atoms using small crosses\natoms can be marked by clicking in this mode\nhotkey=[^]");
  dosmouse.button[2]=curmode==BOND2?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[6],aty,'&',"7","show bonds using split-color lines\nnonbonded atoms using small crosses\natoms can be marked by clicking in this mode\nhotkey=[^]");
  dosmouse.button[2]=palcolor(LIGHTGRAY);

  atx=atx0; aty+=dybutt;
  atx+=makebutton(atx+3,aty,'\\',L[0]+L[1]+L[2]?"box":"xyz","toggle showing of:\n- box in periodic boundary conditions\n- box+walls if option -Z\n- xyz cross in vacuum (free) boundary conditions\nhotkey=[\\]")+6;
  setcolor(palcolor(WHITE));
  atx+=outtextxy(atx,aty,"atom size:")+3;
  atx+=makebutton(atx,aty,'r',"-","decrease ball size\nin modes 1-4 also bar/dumbell size\nhotkey=[r]");
  atx+=makebutton(atx,aty,'R',"+","increase ball size\nin modes 1-4 also bar/dumbell size\nhotkey=[R]");

  if (curmode==DUMBELL || curmode==DUMBELL2) {
    atx=atx0; aty+=dybutt;
    setcolor(palcolor(WHITE));
    atx+=outtextxy(atx,aty,"ball/bar")-2;
    atx+=outtextxy(atx,aty,":")+2;
    planline(atx,"\1\1\2\2");
    makebutton(ATX[0],aty,'u',"-","decrease ball:bar size ratio\nhotkey=[u]");
    makebutton(ATX[1],aty,'U',"+","increase ball:bar size ratio\nhotkey=[U]");
    makebutton(ATX[2],aty,'D',"=1","ball:bar ratio=1 (spherocylinder)\nhotkey=[D]");
    makebutton(ATX[3],aty,'d',"=d","ball:bar ratio=default (dumbbell)\nhotkey=[d]"); }

  if (curmode<DUMBELL || curmode>DUMBELL2 || aty<maxy-13*dybutt+8) {
    atx=atxb0; aty+=dybutt;
    dosmouse.button[2]=palcolor(hbdist?colortab[ihbcolor].turbocolor:LIGHTGRAY);
    atx+=makebutton(atx,aty,'H'&31,"H-bonds","toggle show H-bonds\nsee options -H (turn on and set threshold distance and color)\n-A (set acceptor)\n-D (set donor)\nhotkey=[Ctrl-H]")+1;
    dosmouse.button[2]=palcolor(colortab[ihbcolor].turbocolor);
    atx+=makebutton(atx,aty,'C'&31,"C","change color of hydrogen bonds\nhotkey=[Ctrl-C]");
    dosmouse.button[2]=palcolor(LIGHTGRAY);
    atx+=makebutton(atx,aty,'h',"-","decrease H-bond length threshold\nhotkey=[h]");
    atx+=makebutton(atx,aty,'H',"+","increase H-bond length threshold\nhotkey=[H]")+1;
    setcolor(palcolor(WHITE));
    atx+=outtextxy(atx,aty,hbdist?string("%.2f",hbdist):"off"); }

  atx=atx0; aty+=dybutt;
  setcolor(palcolor(WHITE));
  atx+=outtextxy(atx,aty,"view:")+3;
  planline(atx,"\1\1\4\3");
  makebutton(ATX[0],aty,HOME,"-","perspective: move the view point closer to the screen\ninactive for parallel projection\nhot key=[Home]=[Ctrl-A]");
  makebutton(ATX[1],aty,END,"+","perspective: move the view point farther from the screen\ninactive for parallel projection\nhot key=[End]=[Ctrl-E]");
  makebutton(ATX[2],aty,'=',"proj","toggle parallel/central projection\nhot key=[=]");
  makebutton(ATX[3],aty,TAB,"std","Set standard orientation (z-view),\nthen toggle 3 axes\nhot key=[Tab]=[Ctrl-I]");

  aty+=dyhdr1;
  setcolor(palcolor(WHITE));
  outtextxy(atx0,aty,"-- selection -----"); // 18 letters
  aty+=dyhdr2;
  planline(atxb0,"\4\6\3\2");
  dosmouse.button[2]=movesel?palcolor(ORANGE):palcolor(LIGHTGRAY);
  makebutton(ATX[0],aty,'m',"move","Toggle whether the subsequent move, save or export will affect\nthe whole configuration or the SELECTION only.\n  move=mid mouse\n  x,y-rot=left mouse\n  z-rot=right mouse\nhotkey=[m]\n(save frame by button Save [Ctrl-S], save rotated by [Ctrl-R])");
  dosmouse.button[2]=palcolor(LIGHTGRAY);
  makebutton(ATX[1],aty,'U'&31,"unmark","Unmark all atoms\nhotkey=[Ctrl-U]");
  makebutton(ATX[2],aty,'I',"inv","Invert marking\nhotkey=[I]");
  rg[1]="g123"[clickmode];
  makebutton(ATX[3],aty,'K',rg,"Cycle through 4 modes:\n0=click will mark the shown or nearest atom\n  (based on the draw mode)\n1=click will mark all atoms within a small radius\n2=click will mark all atoms within a medium radius\n3=click will mark all atoms within a large radius\nhotkey=[K]");

  aty+=dybutt;
  planline(atxb0,"\4\7\4");
  makebutton(ATX[0],aty,'S'&31,"save","Save the frame in standard orientation;\nif SELECTION (button [move]) is on, save only the selection.\nSaved files are numbered\nhotkey=[ctrl-s]");
  makebutton(ATX[1],aty,'R'&31,"rotated","Save the rotated frame;\nif SELECTION (button [move]) is on, save only it\nThe saved files are numbered\nhotkey=[ctrl-r]");
  dosmouse.button[2]=bodyframe?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[2],aty,'F'&31,"body","use body-frame coordinates\nWARNING: not intuitive to use with mouse\nhot keys [x] [X] [y] [Y] [z] [Z] [/] [*] recommended\nhotkey=[Ctrl-F]");
  dosmouse.button[2]=palcolor(LIGHTGRAY);

  if (aty>maxy-8*dybutt) goto playback;

  aty+=dyhdr1;
  setcolor(palcolor(WHITE));
  outtextxy(atx0,aty,"-- export bg=");

  atbg=aty;

  printbg();

  aty+=dyhdr2;
  if (dumpstat==SERIES || dumpstat==START) {
    makebutton(atxb0,aty,'e',"end series","end a series of dumps\nor a merged dump\nhotkey=[e]"); }
  else {
    planline(atxb0,"\4\3\3\3");
    makebutton(ATX[0],aty,'G'&31,"w/bg","white/real background toggle\nfor exported PIC EPS NFF POV\nreal = black or as given by -bgRRGGBB\nhotkey=[Ctrl-G]");
    makebutton(ATX[1],aty,'E',"EPS","printscreen in EPS format\nhotkey=[E]");

    if (curmode<=BALL) {
      makebutton(ATX[1],aty+dybutt,'P',"PIC","printscreen in PPM (Portable PixMap, P6) format\n+ show with the app registered in .startdata\n+ dump a series of PPMs (for a movie)\n+ make animated GIF (ImageMagick needed)\nhotkey=[P]");
      makebutton(ATX[2],aty,'N',"NFF","NFF file, data for Mark VandeWettering's\nReasonable Intelligent Raytracer\n(use \'ray\' to raytrace)\nhotkey=[N]");
      makebutton(ATX[2],aty+dybutt,'V',"POV","POV file, data for Persistence of Vision\nraytracer\n(use \'povray\' to raytrace)\nhotkey=[V]");
      makebutton(ATX[3],aty+dybutt,'L',"PLB","dump shown frame/trajectory in the .plb format\nhotkey=[L]");
      if (curmode>=DUMBELL)
        makebutton(ATX[3],aty,'A',"ATM","ATM file, ASCII file of format:\n#_of_atoms\n{box size | box_x box_y box_z | empty line}\nAt x y z  (repeated)\nhotkey=[A]");
      makebutton(ATX[0],aty+dybutt,'B',"ZBUF","ZBUF file, z-buffer for stereogram\n(use \'stereo\' to render)\nthe .zbuf file is a Portable GrayMap (P5)\nhotkey=[B]"); }
    else {
      outbgtextxy(ATX[2]-xfont.width+2,aty,font,palcolor(WHITE),palcolor(LIGHTGRAY),"modes");
      outbgtextxy(ATX[2]+xfont.width*5-1,aty,font,palcolor(WHITE),palcolor(LIGHTGRAY),"3-5");
      outbgtextxy(ATX[0]+xfont.width-1,aty+dybutt-1,font,palcolor(WHITE),palcolor(LIGHTGRAY),"for all functions"); }

    aty+=dybutt; }

 playback:

  atx=maxxn+3; aty=maxyn-4*dybutt-11;
  setcolor(palcolor(WHITE));
  outtextxy(atx,aty,"-- playback ------");

  aty+=dyhdr2;
  planline(atxb0,"\2\1\2\1\2\1\1");
  makebutton(ATX[0],aty,'i',"|\21","play from start, set step to 1, do not change delay\nhotkey=[i]");
  /* \23 = ||   \20 = <|   \21 = |> */
  makebutton(ATX[1]+2,aty,'b',"\20","play backwards\ncancel swing mode if any\nhotkey=[b]");

  dosmouse.button[2]=step?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[2],aty,' ',
             step?" \23":laststep<0?" \20":" \21",
             "stop/start playback (in the selected direction)\nhotkey=[SPACE]");
  dosmouse.button[2]=palcolor(LIGHTGRAY);

  makebutton(ATX[3]-2,aty,'f',"\21","play forward\ncancel swing mode if any\nhotkey=[f]");

  dosmouse.button[2]=palcolor(swing?LIGHTCYAN:LIGHTGRAY);
  makebutton(ATX[4],aty,'w',
             swing==2?"\20\21":swing?"\21\21":" ~",
             "select swing mode (shows current status):\n ~  play from start to end (default)\n \21\21 repeat playing from start to end\n \20\21 play forward, than backward, etc.\nhotkey=[w]");

  dosmouse.button[2]=palcolor(LIGHTGRAY);
  makebutton(ATX[5]+2,aty,'s',"-","play slower\n(add delays)\nhotkey=[s]");
  makebutton(ATX[6],aty,'S',"+","play faster\n(skip frames)\nhotkey=[S]");

  aty+=dybutt;
  planline(atxb0,"\5\4\5");
  dosmouse.button[2]=trace?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[0],aty,'t',"trace","trace mode\nhot key=[t]");
  dosmouse.button[2]=trace>2?palcolor(LIGHTCYAN):palcolor(LIGHTGRAY);
  makebutton(ATX[1],aty,'T',"fade","trace mode with gradual erasing the old frames\nrrepeat for faster shading\ncancel by button trace\nhotkey=[T]");
  dosmouse.button[2]=palcolor(LIGHTGRAY);
  makebutton(ATX[2],aty,'L'&31,"reset","while playback: reset delay to zero and set step to 1\nwire mode still image: redraw the configuration\nhotkey=[Ctrl-L], redraw also [n][n]\n(to cancel tracing, use button trace or key [t])");

  atx=atxb0; aty+=dybutt;
  atx+=makebutton(atx,aty,'[',"<","go by 1 frame backward\nif FILENAME contains format, go by -1 file\nhotkey=[[] or [PgUp]\nbetween < and > buttons there is a slider\nsee also hot keys [,] [.] [{] [}]");
  makeslider(atx+3,aty,getmaxx()-2*xfont.width-maxxn-34,&sliderpos);
  atx+=makebutton(getmaxx()-xfont.width-5,aty,']',">","go by 1 frame forward\nif FILENAME contains format, go by +1 file\nhotkey=[]] or [PgDn]\nbetween < and > buttons there is a slider\nsee also hot keys [,] [.] [{] [}]");
}
