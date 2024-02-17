  /* support of strings like init="cont" in getdata */
  static struct keys_s initkey[]={
    {"continue",0}, {"cont",0},
    {"append",1},
    {"start",2},
    {"random",3},
    {"bias",4}, {"slab",4}, // "slab" kept for compatibility */
    {"crystal",5}, {"lattice",5},
    {"asc",10},
    {"convert",999},
    {"plb",-1},
    {NULL,0}};
  static struct keys_s yesno[]={
    {"no",0}, {"false",0},
    {"yes",1},{"true",1},
    {"forever",-1}, /* for MC */
    {"debug",-1}, /* for quit: deprecated */
    {"paste",-2}, /* for quit: deprecated */
    {"shell",-4}, /* for quit: deprecated */
    {NULL,0}};
  static struct keys_s keykey[]={
    {"run",0},
    {"sort+x",1},
    {"sort+y",2},
    {"sort+z",3},
    {"sort-x",-1},
    {"sort-y",-2},
    {"sort-z",-3},
    {"cp",4},
    {"show",5},
    {"rdf",6},
    {"cn",7},
    {"shell",8},
    {"quit",9},{"exit",9},
    {NULL,0}};
  static struct keys_s pinskey[]={
    {"auto",0},
    {"sc",1},
    {"bcc",2},
    {"fcc",3},
    {NULL,0}};
  static struct keys_s thermostatkey[]={
    {"none",0},
    {"Berendsen",T_BERENDSEN}, {"friction",T_BERENDSEN},
    {"Nose",T_NOSE},{"Nose-Hoover",T_NOSE},
    {"Andersen",T_ANDERSEN},
    {"Maxwell",T_MAXWELL},
    {"AndersenCM",T_ANDERSEN_CM},
    {"MaxwellCM",T_MAXWELL_CM},
    {"Langevin",T_LANGEVIN},
    {"LangevinCM",T_LANGEVIN_CM},
    {"Bussi",T_BUSSI},
    {"BerendsenCM",T_TR},
    {"BerendsenIN",T_IN},
    {"frictions",T_FRICTIONS},
    {"MTK",T_NPT}, {"NPT",T_NPT},
    {"HooverNPT",T_HOOVERNPT}, {"Hoover",T_HOOVERNPT},
    {NULL,0}};
  static struct keys_s sortkey[]={
    {"none",0},
    {"x",1}, {"-x",-1},
    {"y",2}, {"-y",-2},
    {"z",3}, {"-z",-3},
    {NULL,0}};
  static struct keys_s loadtrkey[]={
    {"none",0},
    {"xy",1},
    {"yz",2},
    {"xz",3},{"zx",3},
    {"zyx",4},
    {"xyz",5},
    {NULL,0}};
  static struct keys_s loadzerokey[]={
    {"none",0},
    {"x",1},
    {"y",2},
    {"z",4},
    {"xy",3},
    {"xz",5},
    {"yz",6},
    {"xyz",7},
    {NULL,0}};
  static struct keys_s rescalekey[]={
    {"x",RESCALE_X},           {"xCM",RESCALE_X|RESCALE_CM},
    {"y",RESCALE_Y},           {"yCM",RESCALE_Y|RESCALE_CM},
    {"xy",RESCALE_X|RESCALE_Y},{"xyCM",RESCALE_X|RESCALE_Y|RESCALE_CM},
    {"z",RESCALE_Z},           {"zCM",RESCALE_Z|RESCALE_CM},
    {"xz",RESCALE_X|RESCALE_Z},{"xzCM",RESCALE_X|RESCALE_Z|RESCALE_CM},
    {"yz",RESCALE_Y|RESCALE_Z},{"yzCM",RESCALE_Y|RESCALE_Z|RESCALE_CM},
    {"xyz",RESCALE_XYZ},       {"xyzCM",RESCALE_XYZ|RESCALE_CM},
    {"X",RESCALE_X|RESCALE_PT},{"XCM",RESCALE_X|RESCALE_PT|RESCALE_CM},
    {"Y",RESCALE_Y|RESCALE_PT},{"YCM",RESCALE_Y|RESCALE_PT|RESCALE_CM},
    {"XY",RESCALE_X|RESCALE_Y|RESCALE_PT},{"XYCM",RESCALE_X|RESCALE_Y|RESCALE_PT|RESCALE_CM},
    {"Z",RESCALE_Z|RESCALE_PT},{"ZCM",RESCALE_Z|RESCALE_PT|RESCALE_CM},
    {"XZ",RESCALE_X|RESCALE_Z|RESCALE_PT},{"XZCM",RESCALE_X|RESCALE_Z|RESCALE_PT|RESCALE_CM},
    {"YZ",RESCALE_Y|RESCALE_Z|RESCALE_PT},{"YZCM",RESCALE_Y|RESCALE_Z|RESCALE_PT|RESCALE_CM},
    {"XYZ",RESCALE_XYZ|RESCALE_PT},       {"XYZCM",RESCALE_XYZ|RESCALE_PT|RESCALE_CM},
    {"XX",RESCALE_X|RESCALE_Y|RESCALE_PT|RESCALE_XisY},{"XXCM",RESCALE_X|RESCALE_Y|RESCALE_PT|RESCALE_XisY|RESCALE_CM},
    {"XXZ",RESCALE_XYZ|RESCALE_PT|RESCALE_XisY},{"XXZCM",RESCALE_XYZ|RESCALE_PT|RESCALE_XisY|RESCALE_CM},
    {"XXX",RESCALE_XYZ|RESCALE_PT|RESCALE_XisYisZ},{"XXXCM",RESCALE_XYZ|RESCALE_PT|RESCALE_XisYisZ|RESCALE_CM},
    {NULL,0}};
#ifdef XSECTION
  static struct keys_s xsmodekey[]={
    {"gauss",0}, {"GAUSS",8},
    {"one",1}, {"ONE",1+8},
    {"shear",2}, {"SHEAR",2+8},
    {"xyz",3}, {"XYZ",3+8},
    {"three",4}, {"THREE",4+8},
    {"x",5}, {"X",5+8},
    {"y",6}, {"Y",6+8},
    {"z",7}, {"Z",7+8},
    {NULL,0}};
#endif
#ifdef SLAB
  static struct keys_s slabgeomkey[]={
    {"drop",0}, {"droplet",0},
    {"trickle",1},
    {"slab",2},
    {"cavity",3},
    {"tunnel",4},
    {"bals",5},
    {NULL,0}};
#endif

  if (narg<2) {

static char Phelp[]="\
More help: (b)asic   (c)ontrol  (m)ethods  (v)erbosity   (p)arams  (e)nvironment\n\
           (i)nput   (o)utput   (w)rite    (a)scii dump  (d)ebug   (s)cale ";

static char Pbasic[]="\
\n\
CALL BY:\n\
\n\
  COOK [OPTIONS] SYSNAME[.ble] SIMNAME[.get] [PLBNAME[.plb]] [OPTIONS]\n\
\n\
missing SYSNAME:=SIMNAME, missing PLBNAME.plb:=SIMNAME.plb\n\
If the same OPTION appears several times, the last one applies\n\
\n\
BASIC OPTIONS\n\
\n\
-s0 batch mode [default]\n\
-s  interactive mode (terminal input/output, obsolete)\n\
\n\
-m# mode of operation (read and calculate / simulate and calculate):\n\
  -m0=read configurations SIMNAME.1, SIMNAME.2.. (recorded by -r) and process\n\
      see variables reread.from, reread.to, reread.by\n\
  -m1=read configurations from SIMNAME.plb (as recorded by -y) and process\n\
  -m2=simulate using the Verlet/SHAKE family of integrators\n\
  -m#, #>2 simulate using the GEAR/Lagrangian family of integrators\n\
\n\
-f in the playback mode (-m0,-m1): calculate energy+forces for each cfg read\n\
";

static char Pin[]="\
SELECTED INPUT AND CONTROL FILES\n\
\n\
SYSNAME.ble  force field, generated by blend\n\
\n\
SIMNAME.def  definition (default) data; if missing, SYSNAME.def is tried\n\
SIMNAME.get  input data (if -s then stdin) read after SYSNAME.def\n\
SIMNAME.stp  during run: wait for the cycle to finish then stop\n\
SIMNAME.cpi  variables recorded in convergence profiles\n\
SIMNAME.fix  see (p)arameters, option -k\n\
SIMNAME.box  density as a function of time\n\
SIMNAME.s-s  sites for radial distribution functions [default=all pairs]\n\
SIMNAME.cli  CLUSTERS: definition file\n\
\n\
PLBNAME.plb  input configuration, with init=\"plb\" (or other init<0)\n\
";

static char Pout[]="\
SELECTED OUTPUT FILES\n\
\n\
SIMNAME.cfg  configuration and restart point\n\
SIMNAME.prt  output protocol (if -s then stdout)\n\
SIMNAME.cp   convergence profile, use `showcp' to show and analyze\n\
SIMNAME.rdf  radial distribution functions, use `rdfg' to show and analyze\n\
SIMNAME.plb  trajectory, `show' to display, `plbinfo' for info and more\n\
SIMNAME.sta  statistics, use `staprt' to analyze\n\
SIMNAME.dpr  SLAB: z-density profiles (binary)\n\
SIMNAME.*.z  SLAB: z-density profiles, ascii images\n\
SIMNAME.sfr  structure factor (radial)\n\
SIMNAME.sf3d 3D structure factor\n\
SIMNAME.loc  lock-file, prevents a duplicate simulation\n\
\n\
see also (a)scii dump\n\
";

static char Pctrl[]="\
CONTROL OPTIONS\n\
\n\
-#  overrides the first no=# in SIMNAME.get\n\
-b# in batch mode [default=-15]:\n\
      -b#  check for SIMNAME.stp (stop calculation) every #-th cycle\n\
      -b-# check for SIMNAME.stp (stop calculation) every # seconds\n\
      -b0 do not check\n\
    in interactive mode (-s): set beep level (if terminal accepts \\a=alarm):\n\
      -b<=0 off  -b1 when finished  -b2 after ; in the data\n\
      -b3 as above and SIGINT, SIGTERM\n\
      -b>3: as above and print progress indicator to stderr after every step\n\
-d# check pair distances 1: when ; 2: every cycle+load 3: every step+load\n\
-d-# FREEBC: stop simulation if any atom farther than # from (0,0,0)\n\
-i# SIGINT (Ctrl-C) and SIGTERM is:\n\
    -i-1 kill (default signal handler)\n\
    -i0  ask what to do [SERIAL version interactive default]\n\
    -i1  finish cycle and read next data [PARALLEL interactive default]\n\
    -i2  finish cycle, then save all and stop [batch default]\n\
    -i5  finish step, then save all and stop (to be restarted with init=\"start\")\n\
-o0  ignore lock-file SIMNAME.loc (the same as -o-) [default=-o1]\n\
";

static char Pmethod[]="\
METHODS AND ALGORITHMS\n\
\n\
-c# constraints, sum of flags [default=9]:\n\
  1=(Gear only) correct constraints by conjugate gradients (0=SHAKE-like)\n\
  2=report constraint errors before every integration step\n\
  4=(Gear only) report constraint errors after step and before correction\n\
  8=report constraint errors after step\n\
\n\
-pKPC predictor lengths and orders (dec.digits): [default=229]\n\
  K=POLAR: ASPC (or higher order PSPC) predictor length\n\
  P=POLAR: integration of induced dipoles:\n\
    0:(Car-Parrinello without thermostatting, not recommended)\n\
    1:no prediction (i.e., last value)\n\
    2:ASPC\n\
    3:higher order\n\
  C=Verlet: length of Time Reversible Velocity Predictor [9->2]\n\
    Gear: constraint predictor order [9->Gear order]\n\
";


static char Pwrite[]="\
WRITE\n\
\n\
-e# record max # items in convergence profile SIMNAME.cp [-1=all specified]\n\
\n\
-w# write SIMNAME.cfg: 0=never\n\
                       1=after data set procesed (; in data) [default]\n\
                       2=as 1 but wait if not enough disk space\n\
\n\
-y#  record trajectory for the 1st # molecules and create data for `show'\n\
     using molcfg, output plb-file=SIMNAME.plb, cf. variable dt.plb\n\
-y-# as above except don\'t call molcfg             [default=all molecules]\n\
\n\
-n#  sum of flags:\n\
     1 = write SIMNAME.vlb (velocity playback) as for -y [0]\n\
     2 = shift SIMNAME.plb by -h/2 (in sync with SIMNAME.vlb with Verlet)\n\
     4 = center SIMNAME.plb in the box (not fully supported by utilities)\n\
     8 = POLAR: also absolute Drude positions to SIMNAME.dlb\n\
\n\
-r#  record configurations upto the (#-1)-th derivative\n\
     dt.cfg must be set, output cfg-files are SIMNAME.1,SIMNAME.2.; cf. -m0\n\
     [default=2=positions+velocities]\n\
";

static char Pascii[]="\
ASCII DUMPS\n\
\n\
-l[-]PEGFAVC  string of decimal digits:\n\
  digit 1=full precision, #>1=decimal precision\n\
  keys:\n\
    - in front of argument: configuration SIMNAME.atm [AA]\n\
    C=configuration SIMNAME.asc (MACSIMUS format) [AA]\n\
    V=velocities SIMNAME.vel [AA/ps]\n\
    A=accelerations SIMNAME.acc [AA/ps^2]\n\
    F=forces SIMNAME.for [p.u.=" STRING(forceunit) " N]\n\
    G=numerical gradients SIMNAME.gra [p.u.=" STRING(forceunit) " N]\n\
    E=F+G (should be zero) SIMNAME.err [p.u.=" STRING(forceunit) " N]\n\
    P=Drude amplitude SIMNAME.dru [AA]\n\
  The dump is performed at the and of job (\';\' in input data);\n\
  to record configurations regularly, see dt.plb, dt.cfg in the data and -n.\n\
\n\
Example:\n\
  -l-100003\n\
will dump configuration in compatible ATM format (3 dec. digits in AA)\n\
and test forces by numerical gradients (in p.u.)\n\
";

static char Pdebug[]="\
DEBUG\n\
\n\
+#  with #define CHECKHEAP==-2:\n\
    range (in B) to check array overflow and memory leak\n\
\n\
-d# check pair distances [default=do not check]:\n\
    1=once per job\n\
    2=per cycle\n\
    3=per step\n\
-d-# FREEBC: stop simulation if any atom farther than # from (0,0,0)\n\
\n\
-s# scroll buffer of # kB: type $ to activate (obsolete)\n\
\n\
See also: (a)scii dump, (s)cale, (v)erbosity\n\
";

static char Pverb[]="\
VERBOSITY\n\
\n\
-v# verbosity level, sum of [default=3]:\n\
    1=runtime info, initialization protocols\n\
    2=brief statistics, more runtime info\n\
    4=details on site-site potentials and energy terms, copy of SIMNAME.ble\n\
      constraint optimization (SHAKE), detailed statistics etc.\n\
    8 POLAR: convergence, GAUSSIAN: predicted SCF\n\
    16 POLAR: file SIMNAME.pol with summary\n\
    32 POLAR with dV: extensive info\n\
    64 momentum and angular momentum correction protocol\n\
    128 huge raw dumps, simulation stops after 10 such dumps\n\
\n\
-t  detailed time measurement [0]\n\
\n\
-@WHO@WHWRW e-mail message when finished via 'mail -s'\n\
";

static char Pparm[]="\
PARAMETERS\n\
-h#  initialization: center #-valence atoms to their neighbors [default=0=off]\n\
     e.g., -h4 will try to fix wrong pyramidal configurations of CH4\n\
-k#  keep atoms (listed in SIMNAME.fix) by harmonic force, K=# K/AA^2\n\
-k0  keep atoms exactly in place [default=-16=OFF]\n\
-k-# ANCHOR: keep atoms/molecs/axes (see SIMNAME.fix), sum of bits:\n\
     1: statistics, 2: record in CP (use SIMNAME.cpi), 4=ASCII dump, 8=err/warn\n\
-u0  all bonds are constrained (rigid) [default]\n\
-u#  force constant limit for vibrating bonds, in kcal/mol/AA2, u=K*(r-r0)^2\n\
-u-# wavenumber limit for vibrating bonds in cm-1, use -v4 for details\n\
     HINT: -u-2500 for X-H rigid,  BUGS: single bond value, no equalization\n\
-x   waters, sum of flags [default=5]:\n\
     1=use optimized code for known rigid water models (not LINKCELL)\n\
     2=keep Lennard-Jones terms H-H,O-H for TIP3P\n\
     4=doublecheck\n\
-z# random number seed [0=use time]\n\
-\\# PARALLEL: #>0: number of threads\n\
              0=use environment variable NSLOTS [default]\n\
              #<0 LINKCELL: use No.cell[0] (NOT TESTED RECENTLY)\n\
-^# POLAR: ASPC damping parameter = #/SCALE; see (s)cale\n\
    can be overriden by scf.omega in the data\n\
";

static char Pscale[]="\
SCALING FORCE FIELD\n\
\n\
Scaling factor is |#|/SCALE, where # is the option (integer).\n\
The default scalings are 1.\n\
\n\
-_SCALE set the denominator [default=100=scalings are in %]\n\
\n\
-a# POLAR: scale polarizabilities\n\
\n\
-j#  scale potential energy minima (Lennard-Jones epsilon)\n\
-j-# scale atom radii (Lennard-Jones sigma); BUG: cannot scale both\n\
\n\
-q#  scale charges (not Drude), #<0 allowed\n\
\n\
-^# scf.omega; see (p)arameters\n\
\n\
Example (run with all Lennard-Jones sigmas scaled 1.003-times)\n\
  cookew -_1000 -j-1003 SYSNAME SIMNAME \n\
";

static char Penv[]="\
ENVIRONMENT VARIABLES\n\
\n"
#ifdef PARALLEL
  "To be set by a user:\n\
  NSLOTS = number of parallel threads to be launched\n\
Example:\n\
  NSLOTS=3 cookewslcP1 ff simname\n\n"
#endif
  "Standard variables read and reported (not to be set by a user):\n\
  USER (USERNAME) HOSTNAME (HOST) PWD SHELL\n\
\n";
 
    fprintf(stderr,"%s\n\nSelected #defines:\n ",INFO);

    printdefines(stderr);

    fprintf(stderr,"\n\nThis source file = %s\n", __FILE__);
    fprintf(stderr,"This COOK = %s\n", arg[0]);

    prts("\n\
CALL BY:\n\
\n\
  COOK [OPTIONS] SYSNAME[.ble] SIMNAME[.get] [PLBNAME[.plb]] [OPTIONS]\n\
\n\
missing SYSNAME:=SIMNAME, missing PLBNAME.plb:=SIMNAME.plb\n\
\n");

    for (;;) {
      char s[8];

      prts_(Phelp);
      if (!fgets(s,8,stdin)) return 0;
      switch (s[0]) {
        case 'b': prtsfill(Pbasic); break;
        case 'i': prtsfill(Pin); break;
        case 'o': prtsfill(Pout); break;
        case 'v': prtsfill(Pverb); break;
        case 'c': prtsfill(Pctrl); break;
        case 'm': prtsfill(Pmethod); break;
        case 'w': prtsfill(Pwrite); break;
        case 'a': prtsfill(Pascii); break;
        case 'd': prtsfill(Pdebug); break;
        case 'p': prtsfill(Pparm); break;
        case 's': prtsfill(Pscale); break;
        case 'e': prtsfill(Penv); break;
        default: return 0; }
    }
  }
