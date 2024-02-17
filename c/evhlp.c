void prtsfill(char *s)
{
  char *c;
  int n=0;

  for (c=s; *c; c++) n+=*c=='\n';

  while (n++<22) printf("\n");

  fputs(s,stdout);
}

static char Phelp[]="\n\
MORE: (i)ntro    (o)ptions      (e)nvironment  (c)ommands  (t)oggles  (r)eturn\n\
      (n)um+var  e(x)pressions  (f)unctions    (m)acros    (u)nits    p(.)units\n\
      (s)olve    (a)lgorithms   (p)lot         (d)ata      ()interactive ";

#if CALC&1
static char Pintro[]="\n\
MACSIMUS calculator with units, (c) J. Kolafa 1991--NOW.\n\
\n\
Call by:\n\
  evu [OPTIONS] [EXPR [EXPR2 ..]]          +-----------------------------+\n\
Example:                                   | Version without units:      |\n\
  evu \"exp(-1.1[kcal/mol]/298[K]/R)\"       |   ev \"exp(-1.1*kcal/298/R)\" |\n\
Quit if interactive: quit exit Ctrl-D      +-----------------------------+\n\
\n"
#else
static char Pintro[]="\n\
MACSIMUS calculator, (c) J. Kolafa 1991--NOW.\n\
\n\
Call by:\n\
  ev [OPTIONS] [EXPR [EXPR2 ..]]      +--------------------------------------+\n\
Example:                              + Example of version with units:       |\n\
  ev \"exp(-1.1*kcal/298/R)\"           +   evu \"exp(-1.1[kcal/mol]/298[K]/R)\" |\n\
Quit if interactive: quit exit Ctrl-D +--------------------------------------+\n\
\n"
#endif

"Features/caveats/bugs:\n\
* \'-2*3\' (with - as the 1st character) SUBTRACTS 6 FROM THE PREVIOUS RESULT,\n\
  this feature can be turned off by option ---\n\
  \' -2*3\' (with leading space) or (-2*3) = -6\n\
* \'047\' = 39 (leading 0 = octal number)\n\
* \'2^3^4\' = (2^3)^4 (left-associated)\n\
  \' -2^2\' = -(2^2) (unary +,- have the same precedence as binary)\n\
  \' --2\' = 4 (unary +,- are left-associated)\n\
* leading space prevents commands to be recognized: try \'deg=1;deg; deg\'\n\
* macros cannot be used with ;-separated statements (except the last one)\n"

#if CALC&1
"* cannot add different units: \'0+1[m]\' is error, \'0[m]+1[m]\' is correct\n\
* M(CuSO4.5H2O) = M(Cu)+M(S)+4.5*M(O)+M(H2O); for blue vitriol: M(CuSO4H10O5)\n"
#endif
;

static char Poptions[]="\n\
OPTIONS:\n\
--[VERB]   filter (the same as --d with no prompt)\n\
---[VERB]  -|+ as 1st char is unary -|+ (default=+,- previous result)\n\
--c[WIDTH] prompt-oriented calculator in current terminal (ncurses editing)\n\
           if WIDTH is missing, $COLUMNS is checked, otherwise 80 applies\n\
           also ()interactive from this help\n\
--dVERB    prompt-oriented calculator on a dumb terminal (no history/editing)\n\
--fFMT     format (accepting double), LF appended; must be the last option\n\
--FFMT     as above, no LF appended; must be the last option (override env)\n\
--iDATA    read file DATA with useful constants [default=$HOME/." EV "data]\n\
           see (d)ata section\n\
--i        do not read any initial constants (except embedded pi=PI)\n\
--r[VERB]  if the last result==0 then return code=1 and vice versa, e.g.:\n\
           if " EV " --r0 \"pi>3\"; then echo TRUE; else echo FALSE; fi\n\
--t        tex-formatted output\n\
--vVERB    verbosity, sum of:    1=print results of calculations\n\
                                 2=print f=(last value) after solve etc.\n\
           --v0=no output        4=print expanded macros\n\
           --v=--v1              8=print initial constants (." EV "data)\n\
           default VERB=19      16=print comments (w. option -- only)\n\
";

static char Pret[]="\n\
RETURN CODES\n\
\n\
Default return code (without option --r):\n\
  0 if no syntax error is found (TRUE for shell)\n\
  nonzero on syntax error caught (FALSE for shell)\n\
  numerical errors do not affect the return code\n\
Examples (bash):\n\
  if " EV " \"2***3\"; then echo TRUE; else echo FALSE; fi  # -> ERROR,FALSE\n\
  if " EV " \"1/0\"; then echo TRUE; else echo FALSE; fi  # -> inf,TRUE\n\
\n\
Return code with option --r:\n\
  0 if the last result is nonzero (TRUE for shell)\n\
  1 if the last result is 0 (FALSE for shell)\n\
Example (bash):\n\
  if " EV " --r \"7>3\"; then echo TRUE; else echo FALSE; fi  # -> 1,TRUE\n\
";

static char Pnumvar[]="\n\
NUMBERS:\n\
All numbers are real (converted and interpreted; e.g., 1/4 gives 0.25, not 0)\n\
Decimal numbers (examples):\n\
  1; +2.6; pi; PI; 2e7; -1.6; .8E7\n\
Fortran style accepted (a WARNING is printed once):\n\
  3d7; .3D2\n\
Integer octal and hexadecimal numbers:\n\
  077; 0xff\n\
Sexagesimal numbers (time) are always converted to hours:\n\
  DD::hh:mm:ss\n\
  2:30=2.5, 13:06:36=13.11, 1::12=36\n\
  CAVEAT: leading 0 before digit not accepted (07:32=ERROR)\n\
\n\
VARIABLES:\n\
Leading ASCII character or underscore, case-sensitive, numbers may follow\n\
# = last result\n\
pi = PI = 3.1415926535897931 predefined\n\
undefined variable is ERROR\n\
Examples:\n\
  a; a168; c_7; _3\n\
";

static char Pexpr[]="\n\
EXPRESSIONS:\n\
Unary operators (cf. functions): + - \\(sqrt)\n\
  parentheses may be omitted after a unary operator or function\n\
Binary real operators: + - * / ^(power) **(power) @(modulo)\n\
  as the 1st character of an expression they are prepended by #(last result)\n\
Binary bitwise operators: &(and) |(or)  $(xor)\n\
Boolean operators (give 0=false, 1=true):  < > <= >= == <>\n\
Assignment and calculate+assign operators (can be chained):\n\
  = :=(the same as =) += -= *= /= @=(modulo and assign)\n\
Assign the last result to variable X (the same as X=#): =X\n\
Put one expression or assignments on a line, expressions and assignments\n\
  not containing macros can also be separated by \';\'\n\
Examples:\n\
  1+sin(5*pi); pi>5; ln 2; 2^10-0x400\n\
CAVEATS:\n\
  \'a=1;-2\'  is the same as \'a=1;#-2\'\n\
  \'a=1; -2\' is the same as \'a=1;(-2)\'\n\
  chained powers are left-associated: 2^3^4 = (2^3)^4\n\
  (-8)^(1/3) is -nan, cbrt(-8) is -2\n\
  unary minus has a lower precedence than power: (-2^2') = (-(2^2))\n\
";

static char Pcommands[]="\n\
CONTROL COMMANDS:\n\
\n\
!      Anything after ! is comment\n\
*%*    set format (C-style); e.g.: %7.9f, RES=%.4g, \'%c\', 0x%x, 0%o\n\
%      set default format \" %.8g\"\n\
%?     query format\n\
rad    set radians (default)\n\
deg    set degrees\n\
=      list all variables and macros\n\
~      remove all variables and macros\n\
~NAME  remove variable or macro NAME\n\
cd PATH  change working directory, status returned\n\
sh CMD   run shell command CMD, status returned in #\n\
exe CMD  1. one # in CMD is replaced by val(#)\n\
         2. CMD is executed by shell (see ~/.evexe)\n\
         3. the 1st word of output becomes val(#)"
#if CALC&1
", unit is copied"
#endif
"\n\
write      write history to \"history." EV "\"\n\
write FILE write history to FILE\n\
verbose #  set verbose level (see --v# and --#)\n\
           hint for scripts: use option --0 and verbose 1 before the result\n\
";

static char Pdata[]="\n\
INITIAL DATA:\n\
\n\
Unless changed by option --i, file $HOME/." EV "data with useful constants\n\
is read first.  Then, ./" EV "data is tried if exists.\n\
These constants cannot be overwritten by an assignment, but you may\n\
erase them by command ~.\n\
\n\
Some control commands do not work in the initialization.\n\
Commands def and expand do not accept the \"def a=b\" form, use \"def a (b)\".\n\
\n\
The newest $HOME/." EV "data is based on 5/2019 definitions; i.e.,\n\
all constants c,h,e,NA,k are defined.\n\
Particularly, eps0, mu0 are derived from the experimental aplha.\n\
";

static char Pfunc[]="\n\
FUNCTIONS of one variable, cf. unary operators:\n\
  sqrt(also \\) cbrt abs exp ln(natural) log(common) lb(binary)\n\
  sin cos tan asin acos atan (see also commands: rad deg)\n\
  sh sinh ch cosh th tanh ash asinh ach acosh ath atanh\n\
  gamma lgamma(ln of gamma) erf erfc fact=fac(factorial)\n\
  int(integer towards zero) floor(integer towards minus infinity)\n\
  frac(towards zero) sgn step (Heaviside, step(0)=1/2)\n"
#if CALC&1
"The above functions are dimensionless and accept dimensionless arg. Examples:\n\
  \\2; ln x/ln 3; sh(asinh(2))\n\
Unit-based: unit val (using SI units; val(1[g])=0.001).\n\
M(FORMULA) returns molar mass (BUG: no parentheses accepted in FORMULA):\n"
#else
"Examples:\n\
  \\2; ln x/ln 3; sh(asinh(2))\n\
Function M() returns molar mass in g/mol:\n"
#endif
"  M(Ar) M(H2SO4) ! correct\n\
  M(CuSO4.5H2O)  ! M(Cu)+M(S)+4.5*M(O)+M(H2O)  (NOT WHAT YOU MAY EXPECT)\n\
  M(Ca(CN)2)     ! ERROR\n\
rnd(X) returns a random number according to integer argument X:\n\
  X<-1: set the seed\n\
  X=-1: return a uniformly distributed random number in [-1,1)\n\
  X=0:  return a uniformly distributed random number in [0,1)\n\
  X=1:  return a normalized Gaussian number\n\
  X>=2  return a random integer in [0..N-1] uniformly\n\
";

static char Psolve[]="\n\
SOLVE, MINIMUM, MAXIMUM:\n\
\n\
Find ONE (nearest) root or minimum or maximum of EXPR:\n\
  [RES=]{solve|min|max} X[=X0] EXPR\n\
Secant (solve) or Newton-Raphson-like (min,max) algorithm is used.\n\
Find ALL roots (minima,maxima) using interval search, then as above:\n\
  [RES=]{solve|min|max} X=FROM,TO[,BY] EXPR\n\
Where variables are:\n\
  X = independent variable\n\
  RES = variable to assign the (last) result to\n\
and expressions are:\n\
  EXPR = function of X\n\
  X0 = initial value for iterations; if missing, the current value applies\n\
  FROM,TO = range\n\
  BY = stride (Delta X) for search\n\
NOTE:\n\
  There must not be a space before keywords solve,min,max!\n\
Example:\n\
  x=solve x=1 cos x-x\n\
";

static char Pplot[]="\n\
PLOT function:\n\
  plot X=FROM,TO[,BY] EXPR\n\
Where:\n\
  X = independent variable\n\
  EXPR = function of X\n\
  FROM,TO = range of X to plot\n\
  BY = stride (Delta X)\n\
  FROM,TO,BY may be numbers or expressions\n\
\n\
MACSIMUS plot must be installed, run plot for help\n\
No space before keyword \"plot\"!\n\
\n\
Example (no space before \"plot\"):\n\
  plot x=-10*pi,10*pi sin(x)/x\n\
";


static char Palg[]="\n\
ALGORITHMS:\n\
No space allowed before a keyword!\n\
Calculate definite integral: [RES=]integ X=FROM,TO[,N] EXPR\n\
Calculate numerical derivative: [RES=]deriv X[=X0[,DX]] EXPR  ! also diff\n\
Calculate sum, product: [RES=]{sum,prod} X=FROM,TO[,BY] EXPR\n\
Repeat calculation X=EXPR N-times: [RES=]repeat X=FROM,N EXPR\n\
Iterate until precision |X-Xlast|<EPS achieved: [RES=]iter X=FROM[,EPS] EXPR\n\
Evaluate continuous fraction:\n\
  [RES=]contfrac X=FROM,TO[,BY] [BEFORE++]NUM//DENOM\n\
  Example: x=contfrac n=1,100 1++1//2; x^2\n\
Where variables are:\n\
  X = independent variable\n\
  RES = variable to assign the (last) result to\n\
and expressions are:\n\
  EXPR = function of X\n\
  X0 = value; if missing, the current value of X applies\n\
  FROM,TO = range (of plot, sum, etc.)\n\
  BY = stride (Delta X) for loop, sample points, etc.\n\
  DX = step for 4rd order numerical derivative\n\
  N = count (cast to integer) or number of subintervals (-> 8th order Gauss)\n\
";

static char Pmac[]="\n\
MACROS:\n\
Keywords (as \'def\') must not be preceded by a space!\n\
BUG: macros cannot be used with ;-separated statements (except the last one)\n\
Define macro (X will be replaced by STRING):\n\
  def X STRING\n\
  def X=EXPR  ! the same as \"def X (EXPR)\", not allowed in initial data\n\
STRING and EXPR may contain nested macros which are expanded WHEN X IS USED.\n\
Define macro as above except that the nested macros are expanded NOW:\n\
  expand X STRING\n\
  expand X=EXPR\n\
Undefine macro X:\n\
  undef X\n\
Example (without leading spaces):\n\
  verbose 5\n\
  x=1\n\
  def x1=x+1\n\
  def dxq=x1^2\n\
  expand exq=x1^2\n\
  def x1=x-1\n\
  dxq\n\
  exq\
";

static char Ptoggles[]="\n\
TOGGLES:\n\
\n\
comma   accept \',\' as the decimal point (\'.\' still valid)\n\
        use \',,\' when normally \',\' needed\n\
        prompt changed to r,>\n\
\n\
space   accept single spaces in numbers (e.g., 2*1 000 = 2000)\n\
        use double spaces when normally one needed\n\
        prompt changed to r#\n\
\n\
prompt  turn prompt on/off\n\
\n\
clip    copy the results to clipboard\n\
        command xclip.sh is used, cf. (e)nvironment\n\
        prompt changed to R>\n\
\n\
deg     set degrees, prompt changed to d>\n\
rad     set radians (default), prompt changed to r>\n\
";

#if CALC&1
static char Punits[]="\n\
UNITS:\n\
Units are written in []; the last ] may be omitted; []=previous unit:\n\
  2[kg]; 3[mol/L; 1[m/s]+2[]; 1[g/cm3]+2[\n\
Basic SI units supported: kg m s K A mol B(byte->8); more (see gen/calcu.c):\n\
  t J W Ohm Bq(1/s) Mho(S) L(liter) AA(Angstrom) ha(hectare) cal D(Debye)\n\
  at(technical) atm PSI min d(day) 1[Da]=1[u]=1[g/mol]/NA\n\
  a(tropical year 2000) a_T(tropical mean) a_J(Julian) a_G(Gregorian)\n\
All SI prefixes are understood except micro which is u; in addition:\n\
  da=10, ki=1024, Mi=1024^2, Gi=1024^3, Ti=1024^4; 5[kiB]=40960 (bit=1)\n\
Multiple prefixes are not allowed except for kg ([mkg]=[g], [mks]=ERROR)\n\
Abreviations have precedence over composed units (ha IS NOT hecto-annus)\n\
Clashes: a=annus (not are), ar=are, ha=har=hectare\n\
Any number incl. simple fraction means exponent: [m1/2]=[m0.5] means m^(1/2)\n\
Units are separated by space or dot: [W s] or [W.s] means W*s (=Joule)\n\
Slash: [kg/m s]=[kg/m/s] means kg/(m*s), [/mol]=[mol-1] means 1/mol\n\
Set preferred conversion: to UNIT; to [UNIT]\n\
Remove preference: to ~UNIT; remove all preferences: to ~\n\
\n\
Example (mass of 3 L of N2 at 300 K; note missing space before to):\n\
  p=1[bar]; V=3[dm3]; T=300[K]; m=M(N2)*p*V/R/T;to g\n\
";
static char Ppu[]="\n\
SUPPORT FOR COOK PROGRAM UNITS\n\
\n\
NUM[+UNIT] will convert the quantity given to MACSIMUS cook* program unit\n\
NUM[-UNIT] is the cook program unit of quantity with unit UNIT\n\
The units are interpreted per a particle, not per mole.\n\
WARNING: MACSIMUS uses CGS-like electric quantities,\n\
         additional factor (4*pi*eps0) may apply.\n\
Examples:\n\
  1[+eV]              ! -> 11604.518 p.u. = 1 eV in p.u. (K*k_B)\n\
  8.3144626[+J/mol K] ! -> 1 p.u. (energy p.u. = 1 K*R for mol-based units)\n\
  v=1[-m/s]           ! -> 100 m/s because lengthunit=AA, timeunit=ps\n\
  to D; mu=125[-D]    ! -> 1.4687628 D (dipole moment of 125 p.u. in D)\n\
  to D; mu=125[-C m]  ! the same as above\n\
  1[-J]               ! energy program unit\n\
  1[-mol]             ! 1/NA (makes sense in g/mol etc.)\n\
  33[-g/mol]          ! atom with mass 33 p.u. will be given in g/mol\n\
  33[-g]*NA           ! the same as above\n\
";
#else
static char Punits[]="\n\
Units are not supported by this version (ev), use evu.\n\
";
static char Ppu[]="\n\
Units are not supported by this version (ev), use evu.\n\
";
#endif

static char Penv[]="\n\
ENVIRONMENT VARIABLES:\n\
\n\
COLUMNS # of columns in the terminal\n\
\n\
RNDSEED seed for random numbers (default=0=time)\n\
\n\
HOME    home directory (used to locate ." EV "data)\n\
\n\
TOCLIP  command to copy argument to the clipboard\n\
        with prefix -: command to copy stdin to clipboard (e.g., TOCLIP=-xclip)\n\
        examples: TOCLIP=xclip.sh TOCLIP=-xclip\n\
\n\
FMT     C-style format; ' ' and \\n are stripped off for TOCLIP output\n\
EVFMT   as above, overrides FMT\n\
        NB: options -f and -F override the environment\n\
\n\
EV      string of 2 characters: will replace 1st char by the 2nd\n\
        e.g., export EV=:@ will cause that 5:2=1 (modulo, not time)\n\
";
