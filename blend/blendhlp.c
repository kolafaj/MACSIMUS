/* make blend
 */
void prtsfill(char *s)
{
  char *c;
  int n=0;

  for (c=s; *c; c++) n+=*c=='\n';

  while (n++<21) _n

  prts(s);
}

#include "guihlp.c"
#include "geohlp.c"

static char Phelp[]="\
More help (also F1 in GUI): (i)ntro  (b)asic-opt (f)orce-field  (r)ead   (G)eom\n\
Molecules:                  fi(l)es  (o)ptions   e(x)tra-opt    (w)rite  (g)ui\n\
";

static char Pusage[]="\n\
BLEND V" VERSION ": MACSIMUS minimizer and force field interface. Call by:\n\
\n\
  blend [gen-opt] [par-opt] [ PARFILE.par | PARFILE.bin ] [PCHFILE.pch] \\\n\
    [mol1-opt] { FILE1 | FILE1.mol | FILE1.che }     \\\n\
    [[mol2-opt] { FILE2 | FILE2.mol | FILE2.che }]\n\
\n\
Numeric options are without space in DECIMAL, 0OCTAL or 0xHEX: -v259 = -v0x103\n\
\n\
Environment:\n\
  BLENDPATH = path to force field etc. data, normally MACSIMUS/blend/data\n\
  BLENDGEOMETRY = window geometry, see also option -g\n\
  BLENDINIT = initial keystroke, see also option -I\n\
  GUI = common GUI setup, see page (g)ui\n\
";

static char Pgen_opt[]="\n\
Basic options:\n\
 -o FILE[.EXT]\n\
           generate force-field file readable by cook, [default EXT=ble]\n\
 -v#       verbosity, min=0, then sum of: [default: not -o => -v0, -o => -v3]\n\
           1=list sites and their parameters\n\
           2=list bonded terms (-v3 needed for usable FILE.ble)\n\
           4=verbose (info on topological analysis, selection of terms, etc.)\n\
           8=list used pair site-site parameter terms\n"
#ifdef POLAR
"           16=self field convergence\n"
#endif
"           256=list all pair site-site parameter terms\n\
           [interactive default = -v0, -o default = -v3]\n\
 -g[#]     graphics on: minimization shown (#>0) or printed (#<0) by # steps\n\
 \"-g[#] XGEOM\" or -g[#]_XGEOM\n\
           set X-geometry (overrides BLENDGEOMETRY); e.g., \"-g5 400x300-5+5\"\n\
 -z#       random number seed [default=0=time]\n\
 -_#       output precision of force constants [default=-6]\n\
 -_-#      as above and fix small charge errors to get a multiple of e\n\
 -l#       input line buffer [default=256]\n\
 -ISTRING  initial keystroke (overrides BLENDINIT)\n\
";

static char Ppar_opt[]="\n\
Force field options and files:\n\
\n\
 PARFILE.par    force field parameter file (text format)\n\
 PARFILE.bin    its binary image (deprecated)\n\
 -b PARFILE.par create FILE.bin from FILE.par (deprecated)\n\
 -a#            scale all angle force const #% times\n\
 -a-#           set all angle force const to # kcal/mol\n\
 -a0            constrain all angles (steepest descent only)\n\
 -b#            scale all bond force const #% times\n\
 -b-#           set all bond force const to # kcal/mol/AA^2\n\
 -b0            constrain all bonds (steepest descent only)\n\
 -i#            calc. max. angle pot. energies and warn if > # kcal/mol [df.=5]\n\
 -x[RRR]EEE     Lennard-Jones scaling; R in 1/1000, E in % [default=-x1000100]\n\
 -x-DDD         density scaled x=exp(DDD/1000) times (R: x^-1/3 and E: x)\n\
 -q%            charge multiplication factor in % [default=-q100]\n\
 -t#            cutoff switch function=[#-1,#+1]; -t-12345=[12.3,12.3+4.5]\n\
 -W#            warn on (sum of): 1=missing bonds, 2=missing angles [default=3]\n\
 PCHFILE.pch    partial charge assignment file\n\
 REAFILE.rea    chemical reaction file\n\
 -f1/-f2/-f3    fabricate bonds/angles/both (if missing in PARFILE.par)\n"
#ifdef POLAR
"\nPOLAR options:\n\
\n\
 -|#   overrides polar in the par-file force field (also -\\#)\n\
 -^%   polarizability multiplication factor (df.= 100%, also -~%)\n\
 -]%   saturation energy multiplication factor (df.=100%, also -}%)\n\
 -@%   shell (aux charge for polariz.dip.) multipl. f. (df.=100%, output only)\n"
#endif
;

static char Pmol_opt[]="\n\
Molecule options:\n\
 -c#      chiralities: 0=check  1=calculate unknown  2=calculate all\n\
 -d%      add bond if dist<% of equil. dist (recommended: -d110 -- -d120)\n\
 -d-%     as above and erase all bonds first\n\
 -e       edit mol-file (not if 2D input)\n\
 -e%      2D input scaling factor [df.=100%], %<0: flat 2D, %>0: waved\n\
 -h / -h# constraint bond angles with H / light atoms<#a.u.\n\
 -h-#     as above + don\'t remove overdetermined constraints\n\
 -j       read constraint file FILE.jet\n\
 -k#      keep (freeze) while minimizing: 0=none  1=*marked  2=all but missing\n\
          3=auto 1|2 [default]    4/8=read FILE.keep/FILE.mark\n\
 -k-#     as above + nonbonded atoms always free\n\
 -m#      minimize energy by # conjugate gradient steps\n\
 -m-#     -# steepest-descent steps, then as above [default=-m-100]\n\
 -n#      export number of molecules N (negative: split mol., abs config)\n\
 -u#/-u-# print # of pairs with max/min energy\n\
 -y#      SUM: 1=center before minimization  2=after minimization [df.=3]\n\
          4=set all coord.>0   8=place mol. in box  16=do 4,8 with LJsigma=0\n\
";

static char Pread_opt[]="\n\
Read molecule options and initialize configuration:\n\
\n\
 -r0       read what is available, use -r1 if not (default)\n\
 -r1       generate random configuration\n\
 -r2       derive coordinates from FILE.che (scale by -e%)\n\
 -r3 -r13  read FILE.3db (binary) / FILE.3dt (pure xyz text) (deprecated)\n\
 -r4       read FILE.plb (1st frame)\n\
 -r4:#     read frame # from FILE.plb (1=first, -1=last)\n\
 -r4:#:#:#[:FILE] (for -AEGIR, default = 1:-1:1)\n\
           read frames #:#:#=FROM:TO:BY from FILE.plb\n\
 -r-3 -r-4 as above with reversed endian (deprecated)\n\
 -r6 -r5   generate alpha-helix peptide (-r5 = old algorithm)\n\
 -r7:PHI:PSI:OMEGA\n\
           generate peptide chain with given angles\n\
";

static char Pwrite_opt[]="\n\
Write molecule options:\n\
 -p        write FILE.plb, FILE.gol radii = 70% of vdW radii [default]\n\
 -p%       as above, atom radii scaled by given %\n\
 -p-1/-p-% as above + reverse endian (deprecated)\n\
 -p0       do not write FILE.plb, FILE.gol\n\
\n\
 -w1       write FILE.3db (obsolete; -w-1=reverse endian)\n\
 -w2..-w9  write FILE.3dt (with given # of dec. digits)\n\
 -w10:#    write FILE.pdb, #=name of H on beta C:\n\
           0:HCB1  1:HB1  2:1HCB  3:1HB  4:hydrogens omitted (default)\n\
 -w20:#    as above, derive names from atom IDs (if originate from PDB)\n\
 -w30:#    select -w10/-w20 automatically\n\
 -w40      write FILE.atm (with atom info, for Gaussian etc.)\n\
 -w80      write FILE.cfg (for cook)\n\
           (can be combined; e.g., -w123 = write FILE.3dt, FILE.atm, FILE.cfg)\n\
";

static char Pextra_mol_opt[]="\n\
Special molecule options:\n\
 -[#       (with -n-1) reorder molecule of # of sites so that masses increase\n\
           example -[3: TIP3 water = HHO\n\
 -A        generate FILE.ang of independent internal coordinates\n\
 -C#       reference atom for dipole, quadrupole, virial is: [-1]\n\
           -1=Z-centroid  -2=center-of-mass  -3=center\n\
 -D        generate FILE.msd of mean square displacements\n\
 -E@       essential dynamics, @=atoms: All Heavy Calpha Backbone\n\
 -F#       for -NE: # of frames in each plb file, #>0:cos, #<0:lin [9]\n\
 -G        eigenvalues of the gyration matrix\n\
 -H        eigenvalues of the inertia matrix (former -I)\n\
 -J#       for -NE: accuracy for Jacobi method, sqrt(#)*0.2 for cfg matching\n\
 -M#       for -NE: amplitude of motion in plb [1]; -N -M-#: ns-based guess\n\
 -N#       normal mode frequencies, #=dr for num.deriv. [1e-5 AA]\n\
 -P#       for -NE: # of plb files (most significant eigenvalues) for normal\n\
           mode and ess.dyn. visualization; #<0: least significant [0]\n\
 -R#       core size (where U->infty) for the 2nd virial coefficient, from -C\n\
 -S#       max # of SHAKE iterations (-a0,-b0)"
#ifdef POLAR
"; POLAR: of self-field iterations"
#endif
"\n -T#       temperature [K] for the 2nd virial coefficient and export to .cfg\n\
 -V#       2nd virial coeff: # of sampling points; #<0: +inversion (for dip!)\n\
";

static char Pmol_files[]="\n\
Molecule files:\n\
\n\
FILE.mol   molecule file (obtained by pdb, molcfg, blend...)\n\
FILE.che   molecule in chemical format\n\
FILE       try FILE.mol first, if does not exist then FILE.che\n\
FILE.gol   colors and sizes of atoms (used by show)\n\
FILE.plb   playback file, may contain more frames (preferred format)\n\
           (all lengths are in Angstrom)\n\
FILE.2db   2D configuration, binary\n\
FILE.2dt   2D configuration, text (x,y)\n\
FILE.3db   3D configuration, binary\n\
FILE.3dt   3D configuration, text (x,y,z)\n\
FILE.atm   3D configuration with atom info (for Gaussian etc.),\n\
           line 1 = number of atoms, line 2 = empty or box\n\
           more lines = (x,y,z) in Angstrom\n\
FILE.pdb   PDB format\n\
";
