/* make showcp
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

static char Phelp[]="\
More help:  (i)ntro    (c)olumns      (a)ctions     (r)ange/block/lag\n\
            (o)ptions  (t)time marks  (f)ile format (e)nvironment     (g)ui ";

static char Pintro[]="\n\
SHOWCP " VERSION ": MACSIMUS convergence profile (time evolution) viewer. Call by:\n\
\n\
  showcp [OPTIONS] {SIMNAME | SIMNAME.cp | SIMNAME.cpz} [OPTIONS] [COLUMNS]\n\
\n\
INPUT FILES:\n\
  SIMNAME.cp=convergence profile file, see (f)iles for the format\n\
  SIMNAME.cpz=packed convergence profile file\n\
  SIMNAME=try SIMNAME.cp, then SIMNAME.cpz\n\
\n\
Most important keys/buttons in the plot windows:\n\
  F1=help    F10=menu   Right click on a button = help\n\
  K=[KILL ALL]=kill all plots spawn by one instance of showcp\n\
\n\
Example (SIM.cp -> SIM.cpa in blocks by 10, show Etot, temperatures, P):\n\
  SHOWCPGEOMETRY=800x400 showcp -p -b10 SIM Etot T P\n\
\n\
See also:\n\
  cpacol cp2cp cppak autocorr spectrum\n";

static char Penv[]="\n\
ENVIRONMENT\n\
\n\
  SHOWCPGEOMETRY = X-geometry used for plots\n\
  PLOTCOMMAND = command accepting ASCII data to plot [default = plot]\n\
  LINES = lines/terminal, for ASCII output\n\
  COLUMNS = columns/terminal, for ASCII output\n\
  LOCALE = to determine default ASCII output, cf. option -u\n\
  LC_ALL = as above\n";

static char Ptime[]="\n\
TIME MARKS\n\
\n\
Cook since V3.3m writes marks with simulation time and DT=h*noint to the\n\
cp-file. The cpa-file exported by options -a|-p does not contain a time column,\n\
but if the time cycle DT does not change, the graphs are plotted with the\n\
correct simulation time (in ps).\n\
\n\
The cpa-file exported by options -a2|-p2 does contain time as the first column;\n\
therefore, graphs are displayed properly. However, a range with uniform DT\n\
must be selected (by -f|-t) to get the correct statistics.\n\
\n\
If time marks are not present in the cp-file (cook prior V3.3m, other\n\
utilities), the value of DT=h*noint can be specified by option -h.\n\
Variable DT is not supported";

static char Pcolumn[]="\n\
COLUMNS\n\
\n\
  Columns can be selected by column numbers or quantity names;\n\
  if none is given, all columns are plotted/shown (based on option).\n\
\n\
  -#   select column number # (numbered from 1 in cp-file, can repeat)\n\
  NAME select column of name NAME (can repeat; always -1=Etot, -2=T)\n\
       NAMEs must be last in the argument list\n\
       NAMEs do not work with option --\n\
  \n\
  -o#[,#...] (with -p)\n\
       merge columns (given by numbers) to one plot\n\
       this option can repeat, [default = -o2,6,7 = merge Tkin,Tin,Ttr]\n\
  -o-  all plots are separate (remove default -o2,6,7)\n\
\n\
  -m#  max # of columns (must be 1st option, default=128)\n\
";

static char Pfiles[]="\n\
FORMAT OF CP-FILE:\n\
  RECORD = NCP fields by 4 bytes, 1<NCP<65536\n\
  FIRST RECORD = {CPmark,NCP,COL3,COL4,..}\n\
    CPmark = (float)-1.76398512e+37\n\
    NCP = (int)NCP (legacy: 2 most significant bytes of 4 ignored)\n\
    COL3 = (max 4 chars) name of field 3, ...\n\
           for COL1,COL2, fixed names COL1=Etot,COL2=T are assumed\n\
\n\
  DATA RECORD = NCP float numbers (columns) so that: (1st number)>CPmark\n\
\n\
  CONTROL RECORDS are detected by the 1st float, see file sim/cpmark.h\n\
    {CPmark,24 chars}:\n\
      ASCII write time (when the record was written), cf. \"date\"\n\
    {CPmarkT,int4} (for NCP=2) OR {CPmarkT,time_t=long} (for NCP>2):\n\
      write time in s since 1970-01-01 00:00:00 +0000 (UTC), cf. ctime()\n\
    {CPmarkU,float} (for NCP=2) OR {CPmarkU,double} (for NCP>2):\n\
      simulation time t [ps] of the previous record or simulation start\n\
    {CPmarkV,float} (for NCP=2) OR {CPmarkV,double} (for NCP>2):\n\
      simulation cycle DT=h*noint [ps] of the following data block\n\
    {CPmarkW,double,double} (for NCP>4)\n\
      t,DT (t=sim. time, DT=sim. cycle h*noint [ps]): preferred format\n\
";

static char Prange[]="\n\
RANGE,BlOCK,LAG:\n\
\n\
  -f#  print/export/plot from record # [default=-f1=1st record]\n\
  -t#  print/export/plot to record # (incl.), [default=0=last]\n\
  -b#  block (subaverage) for cpa-file, showing, and additional analysis\n\
       [default=0=auto-adjust to a power of 10, or screen width with -g]\n\
  -b   the same as -b1 = no blocking (showing may be slow)\n\
\n\
  -h#  legacy: simulation cycle DT=h*noint in ps (real number)\n\
       - required for the old cp format (cook < V3.3m), options -x,-c1,-c2,-a2\n\
       - optional for new cp format (WARNING printed on conflict)\n\
  -h-# legacy: as above but silently ignore any time marks in the cp-file\n\
\n\
  -l#  lag for time autocorrelation analysis [default=29]\n\
  -n#  number of consecutive blockings by 2 (not for -c1,2) [default=12]\n\
";

static char Pactions[]="\n\
ACTIONS [default=print pseudographs of selected columns]:\n\
\n\
  --   BUG, DO NOT USE (1-pass filter, only b1,-z,-f,-t,-1,-2..,no names)\n\
  -    filter: good with -a,-b etc.; consider filter | comment -# | dellines 1\n\
  -a   export ASCII file SIMNAME.cpa without time in 1st column (legacy mode)\n\
  -a2  export ASCII file SIMNAME.cpa with simulation time column prepended\n\
  -c#  sum of [default=-c0=one line of summary unless -l0; NB: -c=-c1, not -c0]:\n\
  -C#   1=write files SIMNAME.NAME.tcf with time correlation functions\n\
        2=write files SIMNAME.NAME.cov with covariances\n\
        4=detailed autocorrelation analysis (not just summary)\n\
        8=autocorrelation analysis of data blocked by option -b\n\
       16=autocorrelation analysis with removed linear T-dependence\n\
       32=more decimal digits (with 4 only)   -C = omit bad chars from names\n\
  -e   REMOVED in V2.0e\n\
  -g   print pseudographs [default if columns selected and no other action]\n\
  -|   print merged y-x pseudograph (1 line = 1 block)\n\
  -p   plot blocked selected columns with not-constant data (implies -a)\n\
  -p2  as above, with time in 1st column (implies -a2)\n\
  -x   export blocked x-y files SIMNAME.NAME.xy\n\
  -v#  with control V(t) (cook: tau.rho<0, SIMNAME.box): add (but not plot)\n\
       column of integrated and block-averaged kin. energy of flow,\n\
       where #=mol.mass of the configuration [mg/mol, int]\
";

static char Poptions[]="\n\
OPTIONS:\n\
\n\
  -@   print the # of records and stop (to be used in scripts)\n\
  -d#  delay # ms between plots launched by option -p [default = 500]\n\
  -k#  pseudograph info line: 1=range [default], 2=1st/last; 3=both lines\n\
  -r   running averages instead of time profiles (except bulk Ekin, -v)\n\
  -s#  scroll (buffer #kB, -s=-s30, obsolete)\n\
  -u#  charset for pseudographs: -u: Braille UTF-8, -u-: ASCII [df.=locale]\n\
  -z   the 1st value (or #-th of -f#) subtracted from all data before blocking\n\
  -z#,#,... as above for given columns only\n\
";
