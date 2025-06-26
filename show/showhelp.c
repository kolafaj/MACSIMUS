void help(void)
{
  setgmode(BOND);
  first=1;

  if (getmaxx()<639 || getmaxy()<479) selectfont(0); else selectfont(2);
  dathelp=xfont.height+1;
  athelp=5;

  setwritemode(0);
  clearviewport();
  settextjustify(0,2);

  erasebuttons();
  makebutton(getmaxx()-xfont.width*9-3,5,F1,"exit help","return back to main menu\nhotkey=[F1] or [ESC]..");

 LMENU:
  H("=== GUI MENU ===","                     [show V"VERSION" (c) J. Kolafa]");
  H("F10;\22","GUI menu on/off");
  H("right click                         ;button","show context help on the         clicked");
  H("","- each GUI button is equivalent to a hot key");
  H("","- some hot keys are not available as buttons");

  athelp+=xfont.height/2;
  H("=== QUIT ===","");
  H("ESC q F12","  needs to be typed twice to exit \'show\'");
  H("Q Ctrl-Q;quit","  exit \'show\' immediately");

  athelp+=xfont.height/2;
  H("=== MOUSE ===","");
  H("    The mouse behavior is affected by the coordinate system","");
  H("    and the selection status","");

  H("left click","mark atom and print info");
  H("mid click","mark whole molecule");
  H("right click","unmark atom");
  athelp+=xfont.height/2;
  H("left drag", "rotate in x,y around configuration center");
  H("mid drag", "move in x,y (actually, the viepoint moves wrt screen)");
  H("right drag","rotate in z (around center)");
  athelp+=xfont.height/2;
  H("mouse wheel","zoom view (rescale all incl. viewpoint distance)");

  NEXTPAGE(LMENU,LMENU)

 LDRAWMODES:
  H("=== DRAW MODES ===","");
  H("    There are 7 draw modes. Exporting the configuration","");
  H("    is available in some of them only (see EXPORT below).","");
  H("!;1","bonds=bars, nonbonded atoms=small balls");
  athelp+=3;
  H("@;2","as above, colored (bonds are split)");
  athelp+=5;
  H("#;3","as mode ! and all atoms=small balls");
  athelp+=3;
  H("$;4","as mode @ and all atoms=small balls");
  athelp+=5;
  H("%;5","show atoms using balls of real size");
  athelp+=5;
  H("^;6","wire mode: bond=line, free atom=cross");
  athelp+=3;
  H("&;7","as above, colored (bonds are split)");
  athelp+=xfont.height/2;
  H("g G","cycle through the above 7 drawing modes");
  athelp+=xfont.height/2;
  H("|;box","in modes 1-5: toggle showing the box (and walls) or xyz-cross");
  athelp+=xfont.height/2;
  H("\\","in modes 1-5: toggle style of showing the walls ");
  athelp+=xfont.height/2;
  H("Ctrl-L;reset"," ");
  H("","while playback: set delay to 0 and step to +-1");
  H("","still image: redraw (if problem; e.g., erase traces)");

  NEXTPAGE(LMENU,LDRAWMODES)

 LVIEW:
  H("=== VIEW CONFIGURATION ===","");
  H("Tab;proj",    "cycle through standard orientations || axes, also a");
  H("x X y Y z Z","rotate around axes");
  H("cursors","move the configration");
  H("/ *",    "change the angle for the above rotations");
  H("Ctrl-F;body", "toggle and reset coordinate system (screen vs. body)");
  H("",       "WARNING: using mouse is not intuitive in body-frame");
  H("c",      "recenter");
  H("=;proj",      "central/parallel projection toggle");
  H("Home End","move eye from/to screen (in central projection)");
  H("Ctrl-A Ctrl-E","  the same as Home End");
  H("+ -",    "zoom in/out (rescale all incl. viewpoint distance)");
  H("F7 F8",  "central projection: move in/out");
  H("r R",    "resize ball or bar (in the respective modes)");
  H("o O",    "+/- # of atoms always shown as balls (option -o)");
  athelp+=8;
  H("in modes 3 and 4 (bar+ball):","");
  H("u U",    "change ball/bar radius ratio");
  athelp+=3;
  H("d;=d",   "default ball/bar size ratio (dumbell)");
  athelp+=3;
  H("D;=1",   "ball/bar=1 (fat bars, useful for raytracing)");
  NEXTPAGE(LDRAWMODES,LVIEW)

 HBONDS:
  H("=== VIEW HYDROGEN BONDS (CLUSTERS) ===","");
  H("In modes 3,4 this menu is shown only if the window is tall enough.","");
  H("Ctrl-H", " ");
  H("  ;H-bonds", "toggle view of hydrogen bonds");
  athelp+=3;
  H("Ctrl-C;C","change color of hydrogen bonds");
  H("h H",    "change H-bond length threshold");
  athelp+=3;
  H("Ctrl-T F5","then type a digit, letter (a=10, etc), or F1..F12:");
  H("","mark n-th largest cluster (as defined by bonds and H-bonds)");
  H("","hint: use [move]+[save], then 'plbcluster -l'");
  H("","in presenter mode (option -P), F5=hotkey i (play from start)");
  athelp+=3;
  H("See also option -H (turn on and set parameters)","");
  H("and options -A,-D (set donor and acceptor).","");
  H(string("Current settings: donor=\"%s\" acceptor=\"%s\"",donor,acceptor),"");

  NEXTPAGE(LVIEW,HBONDS)

 LPLAYBACK:
  H("=== PLAYBACK ===","");
  H("i;|\21", "restart playback from the 1st frame, step=1 (no stride)");
  H("s S;+-", "slower (add delays) / faster (omit frames)");
  H("1..9",   "time stride (1=no omitted frame)");
  H("0",      "followed by 2 decimal digits: time stride (000=freeze)");
  H("SPACE;\23",  "freeze(pause)/go (button will change)");
  athelp+=3;
  H("t;trace","trace mode on/off");
  athelp+=3;
  H("T;fade",  "trace + fade-away mode (repeat for fading speed)");
  athelp+=3;
  H("b;\20",   "play backwards");
  H("f;\21",   "play forward");
  athelp+=3;
  H("w;~",    "cycle through one of 3 swing modes (button will change):");
  H("",       "none / autorepeat from start / forward-backward");
  H(", .",    "goto to beginning/end of the playback file");
  H("{ }",    "go by 5% of the file (or as set by option -[)");
  H("[ PgUp;<",   "go by 1 frame (file*) backwards    *applies if PLBNAME");
  H("] PgDn  ;>", "go by 1 frame (file*) forward       contains int-format");
  H("n",      "show/hide frame number info");
  H("F4",     "goto frame (enter number in the terminal); = i in presenter mode");
  H(":",      "set bookmark position");
  H(";",      "exchange bookmark and position");
  NEXTPAGE(HBONDS,LPLAYBACK)

 LEDIT:
  H("=== EDIT ===","");
  H("m;move", "move the whole configuration/move selection toggle");
  H(NULL,      "           SELECTION");
  H("",       "watch info           in the left bottom corner");
  H("I;inv",  "invert selection");
  athelp+=3;
  H("Ctrl-U;unmark"," ");
  H("","unmark all atoms");
  athelp+=3;
  H("Ctrl-S F2;save"," ");
  H("","save frame in plb-format (if marked also mol,gol)");
  athelp+=3;
  H("Ctrl-R F3;rotated"," ");
  H("", "as above, use the current rotation");
  H("K;rg",   "cycle through 4 modes of marking sensitivity:");
  H("",       "- mark the nearest/shown atom (based on mode)");
  H("",       "- mark all atoms within small range");
  H("",       "- mark all atoms within medium range");
  H("",       "- mark all atoms within large range");
  H("",       "(left = mark 1 atom, mid = mark molecule, right = unmark)");
  NEXTPAGE(LPLAYBACK,LEDIT)

 LEXPORT:
  H("=== EXPORT (in bar and ball modes only) ===","");
  H("A;ATM",      "dump ATM file(s) (atom derived from TYPE)");
  H("N    ;NFF",      "dump NFF file(s) (use \'ray\' to raytrace)");
  H("V;POV",      "dump POV file(s) (use \'povray\' to raytrace)");
  H("B    ;ZBUF",     "dump ZBUF file(s), z-buffer (use \'stereo\' to render)");
  H("","Hint: use wide window for better stereo perception");
  H("E;EPS",  "dump EPS file");
  H("P    ;PPM","dump PPM (raw Portable PixMap, P6) image(s), also  animated GIF");
  H("C;PS",  "dump PS file(s) (A4, in color)");
  H("W   ;B/W PS",  "dump PS file(s) (A4, converted to gray)");
  H("L",      "(dump PLB - out of order, use Ctrl-S)");
  H("then select action from menu or hotkey:","");
  H(" o",      "one frame will be dumped");
  H(" O",      "as above and render/view (\'start\' needed)");
  H(" s",      "a series of files will be dumped");
  H(" a",      "animated GIF (from PPM series, cf. option -L)");
  H(" m",      "one merged file will be dumped");
  H(" M",      "as above and render/view (\'start\' needed)");
  H(" e",      "close dumping a series or merged file");
  athelp+=3;
  H("output filename depends on option -v:","");
  H("-v0:","PLBNAME.#@.ext (#=frame,@=a,b,..)");
  H("-v1:","PLBNAME.####.ext (####=0000,0001,..)");

  NEXTPAGE(LEDIT,LEXPORT)

 LMISC:
  H("=== MISCELLANEOUS ===","");
  H("Ctrl-O F6","toggle small/large font size");
#  ifdef SLICES
  H("_",      "toggle cut by z-plane (in ball and bar modes)");
  H("l",      "toggle color of cutting plane");
  H("< >",    "move cutting plane in z-direction");
  H("( )",    "fine move cutting plane in z-direction");
  athelp+=xfont.height/2;
#  endif /*# SLICES */
  H("\' \"",  "match (\"=more accurate), with option -mPLBFILE:FRAME");
  athelp+=xfont.height/2;
#  ifdef SHELL
  H("J j",    "increase/decrease HOH shell width");
  H("",       "(from next frame, see options -j -J)");
  athelp+=xfont.height/2;
#  endif /*# SHELL */
#  ifdef GOLD
  H("| \\",   "toggle wall style (GOLD only, see option -\\)");
  athelp+=xfont.height/2;
#  endif /*# GOLD */
  H("F11",    "show palette (ESC to quit)");
  NEXTPAGE(LEXPORT,LMISC)

 LMACRO:
#include "macrohelp.c"
  H("--- example ---","");
  H("`Ax]Po`","define `a as: rotate by x, show next frame, dump PPM");
  H("`60`a",  "execute macro A 60 times");
  NEXTPAGE(LMISC,LMACRO)

  goto LMACRO;
}
