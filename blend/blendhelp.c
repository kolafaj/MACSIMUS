void help(void) /****************************************************** help */
{
  if (getmaxx()<639 || getmaxy()<479) selectfont(0); else selectfont(2);
  dathelp=xfont.height+1;
  athelp=5;

  setwritemode(0);
  clearviewport();
  settextjustify(0,2);

  erasebuttons();
  makebutton(getmaxx()-xfont.width*9-3,5,F1,"exit help","return back to main menu\nhotkey=[F1] or [ESC]..");


 LMENU:
  H("=== GUI MENU ===","                    [blend V"VERSION" (c) J. Kolafa]");
  H("F10;\22","GUI menu on/off");
  H("right click                         ;button","show context help on the         clicked");
  H("","- each GUI button is equivalent to a hot hey");
  H("","- some hot keys are not available as buttons");

  athelp+=xfont.height/2;
  H("=== QUIT ===","");
  H(".;finish","save and go to the next molecule");
  athelp+=3;
  H("Ctrl-Q;quit","quit blend, unsaved changes lost");
  H("ESC","interrupt minimization, type 2x to quit");
  H("F9","interrupt minimization and enter data in the terminal");

  athelp+=xfont.height/2;
  H("=== MOUSE IN MOLECULE AREA ===","");
  H("left click","mark atom");
  H("mid click","mark whole molecule (connected by bonds)");
  H("right click","unmark atom");
  athelp+=xfont.height/2;
  H("m;move","toggle whether to move selection or molecule");
  H("left drag", "rotate in x,y");
  H("mid drag", "move in x,y");
  H("right drag","rotate in z (around center)");
  athelp+=xfont.height/2;
  H("mouse wheel","zoom view");
  NEXTPAGE(LMENU,LMENU)

 LDRAWMODES:
  H("=== DRAW MODES ===","");
  H("    There are 4 draw modes:","");
  H("1;1","wire mode, nonbonded atoms=crosses");
  athelp+=3;                                                                   
  H("2;2","wire+ball mode");
  athelp+=3;                                                                   
  H("3;3","wire+ball mode, shaded");
  athelp+=3;                                                                   
  H("4;4","wire mode, anaglyph (red/cyan glasses needed)");
  athelp+=xfont.height/2;
  H("g G","cycle through the above 4 drawing modes");
  H("b B;bond","cycle through 3 bond styles");
  H("\'",     "toggle mouse click radius to select atoms");
  H("c",      "gray/dark/cyan carbons");
  athelp+=xfont.height/2;
  H("=== VIEW CONFIGURATION ===","");
  H("x X y Y z Z","rotate around axes");
  H("+ - > <","zoom in/out");
  H("r R",    "ball size: resize balls (modes 2 and 3)");
  H("* /",    "double/halve step size for the above functions");
  H("=;grid",      "grid on/off (by 1 Angstrom, modes 2,3 only)");
  H("Home End","perspective: move eye from/to screen (in central projection)");
  H("Ctrl-A Ctrl-E","  the same as Home End");
  NEXTPAGE(LMENU,LDRAWMODES)

 LMIN:
  H("=== ENERGY MINIMIZATION ===","");
  H(",;CG", "minimize (use conjugate gradients, no constraints)");
  athelp+=3;
  H(",;MC", "minimize (use Monte Carlo, with constraints)");
  athelp+=3;
  H(";;SD","minimize (steepest descent, then conjugate gradients)");
  athelp+=3;
  H(":;ran", "randomize positions and minimize as with [;]");
  H("",  "repeated [:] increases random perturbation");
  H("[;-","decrease the number of minimization steps");
  athelp+=3;
  H("];+","increase the number of minimization steps");
  athelp+=3;
  H("e;trace","how to show minimization: fade,over,none");
  athelp+=3;
  H("j;by","how often to show the course of minimization");
  athelp+=3;
  H(".;finish","save files and go to the next molecule");
  H("ESC","interrupt minimization (2x to quit)");
   NEXTPAGE(LDRAWMODES,LMIN)

 LEDIT:
  H("=== EDIT ===","");
  H("m;move", "toggle whether move the whole configuration");
  H("",       "or the selection (marked atoms)");
  H(NULL,      "           moving ..");
  H("",       "watch info           in the left bottom corner");
  H("I;inv",  "invert selection");
  athelp+=3;
  H("u;unmark","unmark marked atoms, also [Ctrl-U]");
  athelp+=3;
  H("Ctrl-S;save", "save the configuration (in plb-format)");
  athelp+=3;
  H("k;keep", "keep marked atoms fixed while minimizing");
  athelp+=3;
  H("K;mark Kept","mark kept atoms");
  H("",       "(again or kept by * in the mol-file)");
  H("A",      "replace charges of the selected atoms by their average");
  H("E;info", "calculate and print the energy, dipole, quadrupole");
  athelp+=3;
  H("W;put", "remember orientation,position of molecule (not edited fragments)");
  athelp+=3;
  H("w;get", "restore orientation and position (remembered or initial)");
   NEXTPAGE(LMIN,LEDIT)

 LLBL:
  H("=== LABELING ===","");
  H("i;ID",   "show atom IDs (as in the mol-file)");
   athelp+=3;
  H("t;type", "show atom TYPEs (as in force-field)");
   athelp+=3;
  H("q;Q",    "show atom charges in e");
   athelp+=3;
  H("n;No",   "show atom numbers (numbered from 0)");
  H("",       "NB: clicking an atom = finding its number");
  H("SPACE;hide",  "hide labeling (redraw with no labels)");
  H("f;find",  "find atom (enter in the console)");
  H("",       "one of the above labels must be on");
  H("F;unfind","label all again");
  H("L",      "mark all found-and-labeled atoms");
#    ifdef POLAR
  H("d",       "show induced dipoles in e*AA (AA=Angstrom)");
  H("D",       "show induced dipoles in Debye");
  H("p",       "show polarizability volumes (in AA^3)");
  H("s",       "show saturation (some models only)");
#    endif /*# POLAR */


  NEXTPAGE(LEDIT,LLBL)
   LDUMP:
  H("=== FILES AND DUMPS ===","");
  H("P Ctrl-P;PrtScr","");
  H("","dump screen in various formats");
  H("then select 2nd key:","");
  H(" P",  "dump screen in the PPM (Portable PixMap, P6) format");
  H(" p",  ".. as above, white background, some colors changed");
  H(" E",  "dump screen in the EPS format");
  H(" e",  ".. as above, white background, some colors changed");
  H(" C",  "dump screen in the PS format, on A4");
  H(" c",  ".. as above, white background, some colors changed");
  H(" O",  "dump screen converted to gray in the PS format, on A4");
  H(" o",  ".. as above, white background, some colors changed");

   athelp+=xfont.height/2;
  H("@", "read/write editing info (cf. mol2mol utility):");
  H("then select 2nd key:","");
  H(" m","write marked atoms to a file");
  H(" k","write kept atoms to a file");
  H(" M","read marked atoms from a file");
  H(" K","read kept atoms from a file");
   NEXTPAGE(LLBL,LDUMP)


 LMACRO:
#include "macrohelp.c"
  H("--- example ---","");
  H("`32x","rotate 32x by x axis");
   NEXTPAGE(LDUMP,LMACRO)

  goto LMACRO;
}
