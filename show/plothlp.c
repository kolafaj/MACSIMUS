/* make plot
 */
void prtsfill(char *s)
{
  char *c;
  int n=0;

  for (c=s; *c; c++) n+=*c=='\n';

  while (n++<22) _n

  prts(s);
}

#include "guihlp.c"
#include "geohlp.c"

static char Phelp[]="\
More help (also F1 in GUI):\n\
  (i)ntro  (o)ptions  (e)nvironment  (r)ange  (f)iles  (c)olumns  (G)eometry\n\
  (s)tyle  co(l)ors   e(x)pressions  fi(t)    (m)ore   (p)rint    (g)ui   ";

static char Pintro[]="\
PLOT V"VERSION": MACSIMUS command-line plot.  Call by:\n\
\n\
  plot [OPTIONS] [RANGE] [FILE][COLUMNS] [[FILE][COLUMNS] ..]\n\
\n\
Important hotkeys:\n\
  [F1]=help   [F10]=menu on/off   RightClick menu button=context help\n\
\n\
Examples:\n\
- plot a file of 2 columns (x,y):\n\
    plot data.dat\n\
- plot a function:\n\
    plot \'[0:6.28]\' \'[111]:x:sin(x)\'\n\
- plot several files containing 3 columns (N,VAL,ERR) as 1/N vs. VAL+-ERR,\n\
  enhance the 1/N -> 0 limit by adding a/N (adjust a by button [a] or a A / *):\n\
    plot \"[0:]\" \":1/A:B+a/A:-o:C\" data1.dat data2.dat ..\n\
- fit data.dat (three columns, with errors) to a+b*A+c*A^2:\n\
    plot data.dat:A:B:o:C \"[111]:A:a+b*A+c*A^2:-:\"\n\
\n\
See also:\n\
  tabproc tab ev evu field mergetab psstring inset.sh angstrom2eps\n\
";

static char Poptions[]="\n\
OPTIONs must precede any other argument. No space in OPTION (-d1, NOT -d 1).\n\
  -a  Aspect ratio x/y=1: draw a circle by: plot -a '[99:0:6.3]:sin(x):cos(x)'\n\
  -b  Batch mode, useful for fitting: use -I or PLOTINIT with trailing Q\n\
  -c  Blank-line-separated blocks of data are drawn by advancing colors\n\
  -d# Draw slowly: # seconds and XFlush added after each point; e.g., -d.01\n\
  -e  Export file plot.env with the final values of a,b,..j; cf. \"fit.env\"\n\
  -gXGEOM  Initial X geometry (default -g824x480+1+1, less on a small screen)\n\
           Overrides environment variable PLOTGEOMETRY\n\
  -h# Initial number for FORMAT (file name with int format, see (f)iles)\n\
  -h  `%' in file name is ordinary character, not format\n\
  -ISTRING  Initial keystroke; e.g., print EPS and quit plot by -I#Q\n\
            Overrides environment variable PLOTINIT\n\
  -kKEY  Plot only lines containing KEY, no match=blank (even if !#-commmented)\n\
  -nNAME  Name for .def,.eps,.ps files [default=ps.def,plot.eps,plot.ps]\n\
  -NPLOTNAME  Overrides environment variable PLOTNAME\n\
  -p# PID of the parent: used by button [kill all] = hotkey K\n\
  -t  In fitting: do not recompile, function and precision (tT^T) must match\n\
  -v  Verbose: print ranges and parameters at exit (also hotkey v)\n\
  -z# Seed for random numbers [default=0=use time]\n\
";

static char Prange[]="\n\
RANGE argument (see (i)ntro):\n\
\n\
  [FROMX:TOX,FROMY:TOY] or [FROMX:TOX][FROMY:TOY]\n\
  [FROMX:TOX] = shortcut for [FROMX:TOX][:]\n\
  [:] = the same as missing at all\n\
\n\
Brackets [] are not metacharacters: should be protected from the shell.\n\
See (f)iles for [] meaning equidistant points.\n\
\n\
Incomplete ranges (e.g., [0:][0:]) mean that the missing bounds will be\n\
determined from FILEs.\n\
At most one RANGE argument is accepted, between OPTIONS and FILES.\n\
\n\
Example:\n\
  plot -b -I#Q \'[0:10][0:]\' data.dat\n\
";

static char Pfiles[]="\n\
FILE arguments (see (i)ntro):\n\
  FILENAME : File to plot\n\
  - : stdin (pipe) instead of file. Example:\n\
        tail -n1000 file.dat | plot -:A:C\n\
  FORMAT : C-style format accepting integer. The value of option -h (or 0 if\n\
      not given) is substituted and the file is plotted. The integer can be\n\
      in/decremented by hotkeys { } or from menu. Example:\n\
        plot -h1 :n:E sim%04d.cpa\n\
  -NAME : do not plot (right now), yet advance colors/plotstyle as if plotted\n\
  [INT] : (incl. []) Equivalent to a file of INT+1 equidistant points in the\n\
      x-interval given by parameter RANGE.\n\
  [INT:FROMX:TOX] : As above, override x-range given in RANGE. Example:\n\
        plot '[-1:1][-1:1]' '[6:0:360]:sin(PI/180*x):cos(PI/180*x)'\n\
  @ : dummy (useful to set column and plot style)\n\
  @RESPONSEFILE : read arguments (except OPTIONs) from RESPONSEFILE\n\
FILEs to plot are text files, separators=\" \\n\\t\\r,;\".\n\
Valid numbers in files: 1.2 3e1 -2d-2\n\
Invalid data in files: 1,2 (comma is separator), 2/7 (formulas not accepted)\n\
Items in \"\" (as in CSV) are not recognized\n\
#! denotes comment: line is interrupted (points are not connected)\n\
Empty line: a line is interrupted and counter (n = #0) is reset to 0.";

static char Pcolumns[]="\n\
COLUMNS arguments (see (i)ntro):\n\
  :X:Y  = column numbers or expressions to plot by the previous/default STYLE\n\
  :Y    = shortcut for :previous_X:Y or :x:Y if no previous X given\n\
  :X:Y:STYLE = plot by given STYLE (see page (s)tyle)\n\
  :X:Y:STYLE:DY = shortcut for :X:Y:STYLE:DY:DY\n\
  :X:Y:STYLE:DYL:DYH = plot vertical error bars (yy style=default);\n\
                       low and high bounds may differ\n\
  :X:Y:STYLE:DXL:DXH = plot horizontal error bars (xx style)\n\
  :X:Y:STYLE:DX:DY = plot both vertical and horizontal error bars,\n\
                     low and high bounds must be the same (symmetric)\n\
Any of parameters X,Y,DX,DY.. may be:\n\
- column number (dec. integer), can be changed from menu or by hotkeys 0..9\n\
- column number preceded by 0 (as 01), cannot be changed\n\
- expression containing columns coded as: A,B,.. or #1,#2,.. or c1,c2,..;\n\
- x,y,z are aliases for A,B,C\n\
- column 0 = #0 = c0 = @ = n = line count (from 0), reset by blank line in file\n\
  Example (plot 1-column data): plot data.dat:0:1\n\
For parameter STYLE, see (s)tyle\n\
\n\
Missing columns are copied from previous argument, default = 1:2\n\
";

static char Pexpr[]="\n\
EXPRESSIONS use the common precedence rules. Notable features:\n\
  1/2 evaluates to 0.5 (real)   NB: 0,5 IS NOT 1/2\n\
  ^ = ** = power, % = modulo, & = bitwise-and, | = bitwise-or, $ = bitwise-xor\n\
  comparisons < > <= >= == <> return 1=true, 0=false; Heaviside(x)=(x>0)\n\
Functions:\n\
  \\=sqrt (protect \\ from shell!)  cbrt\n\
  exp, ln = natural log, log = decadic log, lb = binary log\n\
  sin,cos,tan accept arg in rad, asin,acos,atan return rad\n\
  sinh,asinh,..,erf,erfc,fac,gamma,lgamma,sgn,int,floor,frac,step = Heaviside\n\
  simple arg does not need (): \\2 = sqrt 2 = sqrt(2)\n\
Parameters:\n\
  10 parameters called a,b,..j are available and adjustable from plot,\n\
  also accessible as environment variables; to export, see option -e.\n\
  Parameters can be changed from menu by slider, or by hot keys:\n\
    increment/decrement a: hotkeys aA, set a=0: Ctrl-A, change increment: /*\n\
Fitting exceptions:\n\
  Formulas are translated via C enhanced by ^ (power) and some changes:\n\
  - 1/2 evaluates to 0 (integer) => use 1./2 to get 0.5 as in plot\n\
  - not available: xor  <>  lb  int  \\  fac  frac  ** (use ^)  pi (use PI)\n\
  - often more parentheses are needed\
";

static char Pstyle[]="\n\
STYLE = string of characters:\n\
  Lines:  - = solid   = = thick   d = dotted\n\
  Dots:   . = pixel   p = small   P = medium   o = large   O = huge\n\
\n\
Decimal numbers code colors: 1=white,2=yellow,3=cyan..\n\
Colors are inverted on white background; see co(l)ors for more\n\
\n\
More keys:\n\
    c = the same as 1 and allow advancing color for next file\n\
   cc = the same as 2 and allow advancing color for next file...\n\
    C = keep last color\n\
    x = set xx style of error bars\n\
    y = set yy style of error bars (default)\n\
   xy = set xy style of error bars\n\
    s = enable style changes by [style] or SPACE\n\
    n = toggle interpretation of blank lines in files, see option -n\n\
\n\
- missing or empty STYLE = previous STYLE (color may be advanced by 1)\n\
- default STYLE = solid line, style changeable by [style] or SPACE\n\
- fixed style cannot be changed by [style] or SPACE unless key s given\n\
";

static char Pcolors[]="\n\
COLORS are coded by decimal numbers in STYLE string:\n\
\n\
  ===============================================================\n\
    #   screen    inverted           #   screen    inverted\n\
  ---------------------------------------------------------------\n\
    1 = white     black              8 = brown     blue-cyan\n\
    2 = yellow    blue               9 = dark cyan salmon\n\
    3 = cyan      red               10 = purple    light green\n\
    4 = magenta   green             11 = gray      gray\n\
    5 = green     magenta           12 = olive     purple\n\
    6 = red       cyan              13 = dark red  light cyan\n\
    7 = blue      yellow            14 = dark blue light yellow\n\
  ---------------------------------------------------------------\n\
\n\
Colors are inverted at white background, as for EPS output.\n\
";

static char Penv[]="\n\
ENVIRONMENT variables:\n\
  GUI = common GUI setup, see page (g)ui\n\
\n\
  TOCLIP = command to put data into clipboard (after left/mid clicks)\n\
    -COMMAND : via pipe, X11 example: export TOCLIP=-xclip\n\
    COMMAND : as argument, KDE example:\n\
      export TOCLIP=\"dcop --user $USER klipper klipper setClipboardContents\"\n\
\n\
  PLOTNAME = set name of window and icon, do not confuse with option -n\n\
  PLOTGEOMETRY = window geometry, see also option -g\n\
  PLOTINIT = initial keystroke, see also option -I\n\
  PLOTCMD = sh command executed before all (re)plots\n\
            parameters a,b,..,j are exported as $a,$b,..,$j\n\
\n\
  FIT4PLOT = path to fitting support [default=/home/USER/macsimus/c/fit4plot/]\n\
  a,b,..j = parameters, changeable by hotkeys a A ^A etc.\n\
  astep,bstep,..,jstep = initial steps for hotkeys a A etc.\n\
  aname,bname,..,jname = names printed in the info panel\n\
  FIT = postprocessing of fit, see (m)ore\n\
";

static char Pmore[]="\n\
ENVIRONMENT VARIABLE FIT:\n\
\n\
Postprocessing of fit; MC-sampled uncertainties with 'nerr'. One of:\n\
  expr  X      : calculate expression f(X)\n\
  integ NINT[,XFROM,XTO] : calculate integral of f(X) in range [XFROM,XTO]\n\
                 missing XFROM = 1st x in data\n\
                 missing XTO = last x in data\n\
                 NINT=# of subintervals of Gauss quadrature (4-th order)\n\
  deriv X,DX[,[-]O] : calculate centered O-th numerical derivative f^(O)(X)\n\
                 O={1,2,3,4}=2nd-order formula, O={-1,-2,-3}=4th-order formula\n\
  solve Y[,X0] : solve equation f(X)=Y by the secant method, X0=initial value\n\
The above keywords can be shortened to 1 letter {e,s,d,i}\n\
Function f(X) is:\n\
  f(X) = y-column of 3rd [FILE][COLUMNS] argument (if 3rd column is given)\n\
  f(X) = y-column of 2nd [FILE][COLUMNS] argument (if 3rd column is not given)\n\
Example:\n\
  fit F.dat w. stderr to a*exp(-b*A), print exp(2*b) w. stderr (hot keys Nt):\n\
  FIT=\"e 2\" plot F.dat:A:B:o:C \'[99]:A:a*exp(-b*A):-:\' \':A:exp(b*A)\'\n\
Batch mode example:\n\
  FIT=\"e 2\" plot -b -INtQ F.dat:A:B:o:C \'[99]:A:a*exp(-b*A):-:\' \':A:exp(b*A)\'\n\
";

static char Pfit[]="\n\
FITTING:                         WARNING: not reentrant (consider lock.sh)\n\
\n\
Input data:\n\
  x,y in columns 1,2, optional stderr (default=1) in column 3\n\
Object function:\n\
  sigma^2 = SUM_i { [y_i-f(a,b,..;x_i)]^2 / stderr_i^2 } / (n-np)\n\
where n=# of points, np=# of parameters (max. 10).\n\
Parameters a-j must be used consecutively (cannot use c without b)\n\
Examples:\n\
  plot FILE \":a+b*A\"\n\
  plot \'[0:1]\' FILE:A:B:o:C \"[111]:A:(a+b*A^(1./4))/(1+c*A)\":-:\"\n\
Expression syntax:\n\
  f(a,b,..;x) is translated by C, see e(x)pressions -> Fitting exceptions!\n\
Hints:\n\
  Pre-adjust values of parameters a,b,c.. by hotkeys aAbBcC/* if necessary\n\
    and fit by [std] = hotkey t, in case of problems [high] = T.\n\
  Standard errors of parameters by sampling: set nerr, hotkeys = nN\n\
  Command-line front-end with more help: hotkey =\n\
Parameter values are exported to fit.env, can be read by `source fit.env'\n\
See (m)ore for integral, derivative, function of x\n\
";

static char Pprint[]="\n\
Print graphs in EPS and PS formats:\n\
  - Axes, labeling etc. are set in definition file NAME.def (see -n); example:\n\
      m eps                  ! mode EPS (also Portrait, Landscape -> PS)\n\
      w 12cm 8cm             ! size of the draw area (without labeling)\n\
      x 60 5 \\r/(kg_ m^\\-3)  ! x axis left, right border in pt, x label\n\
      rotate 90              ! rotation applies to the following text\n\
      y 60 5 E_p_o_t/(J/mol) ! y axis\n\
      l 120 120 $1-E_1       ! label line in color 1 (black) at position in pt\n\
      l +0 +20 $2=E_2        ! next line color 2 (blue), errbar\n\
      t 1.5 1                ! lines and frame thickness\n\
      s 14 Helvetica         ! font size and name\n\
      f 3 3 5                ! approx. number of ticks, tick lengths (pt)\n\
      f -1 -5 2mm            ! x,y ticks by (in graph units), tick lengths\n\
    Units: cm mm in pt (in/72.27), nothing = pt (in/72)\n\
    See \".../macsimus/examples/ps.def\" for more info\n\
  - Output is (see -n) NAME.eps (default) or NAME.ps (Landscape or Portrait)\n\
  - Use hotkey # or button [EPS] to print, F4 to print and quit\n\
  - Use the middle button to get coordinates in pt=in/72 to be pasted to\n\
    commads l in \"ps.def\", or to EPS as X Y moveto () show\n\
";
