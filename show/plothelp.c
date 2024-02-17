/* make plot
 */

void help(void)
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
  H("=== GUI MENU ===","                     [plot V"VERSION" (c) J. Kolafa]");
  H("F10;\22","GUI menu on/off (darker and shifted if menu off)");
  H("right click                         ;button","show context help on the         clicked");
  H("","- each GUI button is equivalent to a hot hey");
  H("","- some hot keys are not available as buttons");

  athelp+=xfont.height/2;
  H("=== QUIT ===","");
  H("ESC q F12","type twice to quit plot");
  H("Ctrl-Q Q;quit","   quit plot immediately");
  athelp+=2;
  H("K;kill all","kill all plots spawned by PARENT (given by option -p)");
  athelp+=2;
  H("K;KILL ALL","kill all plots by \'killall plot\' (if enabled by K in GUI)");
  H("","- no info nor \"plot.env\" written when killed");

  athelp+=xfont.height/2;
  H("=== MOUSE ===","");
  H("left click","digitize: print x y, on 2nd click also dx dy");
  H("mid click","print x,y in pt, useful in \"plot.eps\",\"ps.def\"");
  H("","left or mid: if TOCLIP then copy to the clipboard");
  H("right click","show file and parameter info");
  H("left drag","rectangle to zoom in");
  H("mid drag","move the graph (panning)");
  H("right drag","rescale graph");
  H("wheel","rescale graph in the y-axis");
  NEXTPAGE(LMENU,LMENU)

 LSCALE:
  H("=== SCALING ===","");
  H("+ -","zoom in/out both axes modetately");
  H("x X","zoom in/out x-axis (10x)");
  H("y Y","zoom in/out y-axis (10x)");
  H("Ctrl-X;cx","center around x=0");
  athelp+=2;
  H("Ctrl-Y;cy","center around y=0");
  H("cursors","move viewpoint (panning, cf. mid-mouse)");
  H("R","round the axes up to int-multiples of powers of 10");
  H("Ctrl-R","round the axes to int-multiples of powers of 10");
  H("u Ctrl-Z;undo","     undo last zoom");
  H("k;init","     initial zoom");
  H("F9 Ctrl-U;recalc","     recalculate scaling");
  athelp+=xfont.height/2;
  H("=== FILES ===","");
  H("r","reread files and redraw");
  H("Ins W","show file and parameter info (also right click or [W])");
  H("Del","do not show  file and parameter info (also [w])");
  H("PgUp PgDn","cycle through files if FILE contains int format; e.g., x%03d.dat");
  H("{ }","the same as above");
  athelp+=xfont.height/2;
  H("=== REDRAW ===","");
  H("Ctrl-L","redraw");
  NEXTPAGE(LMENU,LSCALE)

 LSTYLE:
  H("=== PLOT STYLE ===","");
  H("The following functions apply only if the draw style is not","");
  H("set from the command line nor by `s' in STYLE:","");
  athelp+=2;
  H("SPACE;style","toggle 3 styles: lines, points, lines+points");
  H(": ;","draw by dotted line");
  H("\\    ;\\","draw by solid line");
  H("l L","cycle through line styles (dotted/normal/thick/thicker..)");
  H(".;.","draw points as pixels (smallest dots)");
  H("p P","change point sizes");
  H("z Z;grid","cycle through grid styles (none, dense, thin)");
  athelp+=4;
  H("^ F11;errbar","    toggle showing error bars (if any, default=on)");
  athelp+=4;
  H("F5 Ctrl-O;font","    toggle font size (small/large; not for this help)");
  NEXTPAGE(LSCALE,LSTYLE)

 LCOLUMNS:
  H("=== COLUMNS ===","");
  H("[ ]","decrement/increment the column to plot (on y-axis)");
  H("",   "Example: plot -:A:B (or plot -:1:2) \32 plot -:A:C");
  H("",   "More precisely:");
  H("",   "  If the ASCII code of the 1st char in the Y-field is in A..Z,");
  H("",   "  it is changed by -1/+1.");
  H("1 2 .. 9","select the column to plot (on y-axis)");
  H("0A..0Z","select the column to plot (y-axis; 1<=column<=26)");
  H("001 .. 026","the same as above");
  H("0x 0y 0z","the same as 0A 0B 0C");
  H("0n","select \"column 0\" = line number on y-axis");
  H("U;Un","Undo column change (reset to the original)");
  athelp+=xfont.height/2;
  H("NOTES:","Columns are recoded: 0 \32 n, 1 \32 A, 2 \32 B, .., 26 \32 Z");
  H("","columns 27..99 (\32 c27..c29) are available from CLI only");
  H("BUG(?):","undo of column>26 may be incorrect (use c27 etc.)");
  NEXTPAGE(LSTYLE,LCOLUMNS)

 LDUMP:
  H("=== OUTPUT ===","");
  H("# F3;EPS","export \"plot.eps\" or \"plot.ps\" with settings from \"ps.def\"");
  H("F4","as above and Quit");
  H("","for more info, see \".../macsimus/examples/ps.def\"");

  athelp+=xfont.height/2;
  H("Ctrl-P F2;PrtScr","        dump screen in various formats, second key:"); athelp+=4;
  H(" P;PPM on black",    "        dump screen in the PPM (Portable PixMap, P6) format"); athelp+=4;
  H(" p;PPM on white",    "        = as above, white background, some colors changed"); athelp+=4;
  H(" E;EPS on black",    "        dump screen in the EPS format"); athelp+=4;
  H(" e;EPS on white",    "        = as above, white background, some colors changed"); athelp+=4;
  H(" C;PS on black",     "        dump screen in the PS format, on A4"); athelp+=4;
  H(" c;PS on white",     "        = as above, white background, some colors changed"); athelp+=4;
  H(" O;gray PS on black","        dump screen in gray in the PS format, on A4"); athelp+=4;
  H(" o;gray PS on white","        = as above, white background, some colors changed");
  NEXTPAGE(LCOLUMNS,LDUMP)

 LPARMS:
#if PARM
  H("=== PARAMETERS ===","");
  H("There are 10 parameters denoted a-j, available in","");
  H("formulas and accessible as shell variables.","");
  H("Info on parameter not present in formulas is in gray.","");
  H("See also environment variables PLOTCMD and FIT.","");
  H("Example:","a=2 plot \"[0:99]\" \"[999]:x:sin(x+a)\"");
  H("   ;a","button to get slider for parameter a");
  athelp+=2;
  H("The following functions also recenter the slider:","");
  H("a;<","decrease the value of parameter a by step");
  H("A  ;>","increase the value of parameter a by step");
  H("*;*","double the step for changing the active parameter");
  H("/  ;/","halve the step for changing the active parameter");
  H("Ctrl-A ;=0","set parameter a/the active parameter to 0");
  H("< o >","as [a] [Ctrl-A] [A] for the active parameter");
  H("F8;all0","erase all parameters");
  athelp+=xfont.height/2;
  H("v","print x,y ranges and active variables");
  H("V","export active variables to plotenv.sh:");
  H("","can be read from (ba)sh script by \". plotenv.sh\"");
  athelp+=xfont.height/2;
  H("=","switch to terminal to enter data from keyboard (with help)");
#else
  H("not available in this version","");
#endif /*# PARM */
 NEXTPAGE(LDUMP,LPARMS)

 LFIT:
  H("=== FITTING ===","");
#if PARM
  H("Example:","plot FILE \":a+b*A\"");
  H("For more help, run plot without arguments and select t as fi(t)","");
  H("See also environment variable FIT","");
  athelp+=2;
  H("t;std","fit the data to the function, long double (recommended)");
  athelp+=2;
  H("T;high","fit the data to the function, emulated high precision");
  H("Ctrl-T;fast","fit the data to the function, double precision");
  H("n N","decrease/increase the number of sampling points");
  H("","to calculate standard uncertainties of parameters");
  H(" ;100","button to set the number of sampling points (10,100,999)");
  H("Ctrl-N;0","no uncertainties calculated");
  athelp+=2;
  H("=","switch to terminal to enter data from keyboard (with help)");
#else
  H("not available in this version","");
#endif
  NEXTPAGE(LPARMS,LFIT)

 LMACRO:
#include "macrohelp.c"
  H("--- example ---","");
  H("`AaB`","define `a as: decrease parm a, increase parm b");
  H("`60`a",  "execute `a 60x");
  NEXTPAGE(LFIT,LMACRO)

  goto LMACRO;
}
