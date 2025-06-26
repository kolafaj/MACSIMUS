/* make show
 */
void prtsfill(char *s) /******************************************* prtsfill */
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
  (i)ntro   (b)asic options  (s)how options      (c)olors      (G)eometry\n\
  (f)iles   file (o)ptions   bo(x) show options  (p)lb format  (g)ui  ";

static char Pintro[]="\n\
SHOW " VERSION ": MACSIMUS playback file viewer. Call by:\n\
\n\
  show [OPTIONS] MOLNAME[.mol|.gol|.plb] [ PLBNAME.EXT ]\n\
\n\
Important hotkeys:\n\
  [F1]=help   [F10]=menu on/off   RightClick menu button=help\n\
\n\
Environment:\n\
  GUI = common GUI setup, see page (g)ui\n\
  SHOWNAME = name of the window\n\
  SHOWNMF = name of the nmf-file to show vibrations, cf. sh/{nmf,shownmf.sh}\n\
            example: SHOWNMF=a.nmf show -Isssww a a.nm%04d.plb\n\
  SHOWPAL = number of palette entries (<=255, default=249)\n\
  RES = screen resolution in DPI (needed by `stereo', checked after [ZBUF])\n\
  SHOWGEOMETRY = window geometry, see also option -g\n\
  SHOWINIT = initial keystrokes, see also option -I\n\
";

static char Pplb[]="\n\
MACSIMUS playback file:\n\
  binary file composed of float numbers (4 bytes long)\n\
==vacuum (free) b.c.======periodic cubic box======general rectangular box======\n\
  (float)NS 0.            (float)NS L             (float)NS -3.\n\
  X[0] Y[0] Z[0]          X[0] Y[0] Z[0]          L[0] L[1] L[2]          -----\n\
  X[1] Y[1] Z[1]          X[1] Y[1] Z[1]          X[0] Y[0] Z[0]          frame\n\
     ..                      ..                      ..                     1\n\
  X[NS-1] Y[NS-1] X[NS-1] X[NS-1] Y[NS-1] X[NS-1] X[NS-1] Y[NS-1] X[NS-1] -----\n\
  X[0] Y[0] Z[0]          X[0] Y[0] Z[0]          L[0] L[1] L[2]          -----\n\
  X[1] Y[1] Z[1]          X[1] Y[1] Z[1]          X[0] Y[0] Z[0]          frame\n\
     ..                      ..                      ..                     2\n\
  X[NS-1] Y[NS-1] X[NS-1] X[NS-1] Y[NS-1] X[NS-1] X[NS-1] Y[NS-1] X[NS-1] -----\n\
  ...\n\
where\n\
  NS = number of sites (not atoms: M in TIP4P is included)\n\
  L = cubic box size, L[0] L[1] L[2] = rectangular box sizes, in Angstrom\n\
  X[i],Y[i],Z[i] = position of atom i (as in MOLNAME.mol), in Angstrom\n\
  The header is 2 floats long, frames (without the header) are repeated\n\
  Frame F starts at position 8+(F-1)*12*(NS+(L==-3))\n\
";

static char Pfiles[]="\n\
INPUT FILES:\n\
  MOLNAME.mol: mol-file, produced by `pdb', `blend', or `molcfg'\n\
  MOLNAME.gol: optional file w. colors+atom radii by `blend', ignored if -]0\n\
  PLBNAME.EXT: MACSIMUS playback file, see `(p)lb format' for definition\n\
               if missing, PLBNAME.EXT=MOLNAME.plb assumed\n\
               PLBNAME may contain an integer format; e.g., nm%04d.plb.\n\
\n\
OUTPUT FILES:\n\
  numbering: PLBNAME.0000.EXT, PLBNAME.0001.EXT, ..\n\
  numbering if option -_: PLBNAME.FRAMEa.EXT, PLBNAME.FRAMEb.EXT,..\n\
where .EXT=\n\
  .ps  : PostScript: see buttons [PS] [PS B/W] or hotkeys W C\n\
  .eps : EPS: see [EPS] or hotkey E\n\
  .nff : Neutral File Format, for ray: [NFF] or N\n\
  .pov : PovRay scene file: [POV] or V\n\
  .plb : output playback: [save] [rotated] or Ctrl-S Ctrl-R L\n\
  .zbuf: z-buffer as Portable Gray Map (P5) with `# SCALING' line: [ZBUF] or B\n\
         see utility `stereo'\n\
  .ppm : raw Portable Pixel Map (P6): [PPM] or P\n\
  .edt : editing commands for `blend -e', see option -~\n\
";

static char Pbasic[]="\n\
BASIC OPTIONS:\n\
-gOPT initial X-geometry (e.g., 320x240+100+100), overrides SHOWGEOMETRY\n\
-s#   initial scale in % (relative from default) and center [-s100]\n\
-S#   initial scale as screen width in AA\n\
-P    presenter mode: F4 or F5 = run from start (hotkey \'i\')\n\
-P#,# initial x,y-position of molecule in AA (example: -P0.5,2)\n\
-IOPT initial keystroke, overrides environment variable SHOWINIT\n\
      example: -I%$'\\t'i = ball mode, z-axis=horizontal, start playback\n\
\n\
-x# -y# -z# add (x,y,z) to coordinates, in AA (after centering by -c)\n\
-x#L -y#L -z#L as above, in the box size units\n\
-l    periodic boundary conditions wrapped (example: -l -x0.5L)\n\
\n\
-o#   # of atoms shown as balls in bar and barball modes (1-4) [default=0]\n\
-o-#  inverse: # of atoms NOT shown as balls, cf. hotkeys oO\n\
-cXYZ center molecule in x,y,z (-c111 = full center, -c0=no center)\n\
-q#   cue depth in % (far=dark), in 3D modes (not wire) [-q50]\n\
\n\
-v#   verbosity level {0,1,2} [-v1]\n\
";

static char Pcolor[]="\n\
COLOR OPTIONS:\n\
Available colors: W=white, Y=yellow, R=red, G=green, B=blue, C=cyan,\n\
                  M=magenta, O=Orange; X=omit(do not show)\n\
-@$$  atoms/bonds with $$ in ID are colored by color @:\n\
      Multiple options are scanned from right, example:\n\
        `-Rwat -X' will show `wat' in red and nothing else\n\
-p#@  recolor atom to @ if group charge > #/100 e\n\
-n#@  recolor atom to @ if group charge < #/100 e\n\
      Example: -n-80Y shows all groups < -0.8e in yellow\n\
-e0   as above, group smaller by 1 atom (neighbor)\n\
\n\
-H[dist][@] H-bond: if |A-D|<dist then colored by @ [default=off, -H=-H2.4G]\n\
-Aacceptor (with -H): atom names; wildcards ? * accepted [defaults = -AO* -DH*]\n\
-Ddonor               acceptor=donor possible\n\
-K    carbons are black (dark gray) [default=cyan]\n\
-Qopq,q NFF merge only: multiply color by opq, then opq:=(1-q)*opq+q [-Q1,0]\n\
\n\
-bgRRGGBB background for ball and bar, also exported to NFF etc [-bg000000]\n\
-w#   export with white background (also Ctrl-G [w/bg]) [default=0=off]\n\
";

static char Pfileopt[]="\n\
FILE OPTIONS:\n\
-_    change format of output file names, see (f)iles\n\
-h#   1st file # for PLBFILE containing int format\n\
-]0   ignore MOLNAME.gol and guess colors+atom sizes from MOLNAME.mol\n\
\n\
-mMATCH[:FRAME:[FIX]]\n\
      match configuration with another configuration (hot keys \'\"):\n\
      MATCH=plb-file to match, use given FRAME [default=1]\n\
      FIX=file of sites (1 per line) to match [all sites if FIX is missing]\n\
      File FIX can be obtained by blend, hotkeys @m @k\n\
\n\
-~1   append atoms clicked as `ra ATOM' to MOLNAME.edt\n\
-~2   as above, rewrite MOLNAME.edt\n\
-~-1  append atoms clicked as `rm ATOM' to MOLNAME.edt\n\
-~-2  as above, rewrite MOLNAME.edt\n\
\n\
-LOPT options for ImageMagick's command `convert' for animated GIFs, default:\n      "
      CONVERT
"\n\n\
-k#   DEPRECATED plb header size is # floats (use -k0 for .3db) [2]\n\
";

static char Pbox[]="\n\
BOX AND SLIT OPTIONS:\n\
-|@#  show box (if defined in the plb-file) in color @ and thickness #\n\
      #=real number: #>0 thickness in AA, #<0 relative to bond thickness\n\
      optional @ = color, 1 uppercase char, see (c)olors\n\
      x,y,z axes are marked by Red Green Blue\n\
      equivalent form: -\\@#, default=-|O0.1, protect | and \\ from the shell\n\
-|    do not show box\n\
-Z#,#,# show the slit pore, the 3 parameters are:\n\
      z0 = bottom wall z-position in the unit of Lz (default=none)\n\
      z1 = top wall z-position in the unit of Lz (default=none)\n\
      vdW = vdW diameter of the wall atom (0.7*vdW is shown) or grid span\n\
\n\
-[#   hot keys [ ] move by # frames\n\
-[-#  hot keys [ ] move by #% of total length file length (default=-[-5)\n\
\n\
-t1   raytracer code will contain x-z floor\n\
-t-1  as above, for a series (?)\n\
";

static char Pshow[]="\n\
SHOW OPTIONS:\n\
-d#   initial delay between frames in ms (default=0, cf. hotkeys sS)\n\
-a#   max rotation angle (xXyYzY) is 2*PI/#, initial (PI/8)/# (default=4)\n\
-f#   size of cross used to show free atoms in wire modes\n\
-f1   free atom shown by one pixel\n\
-f-#  as above, and small balls in ball modes shown as pixel [default=-2]\n\
-u    thick bars instead of dumbells (bar+ball)\n\
-u0   no limit applies to sphere/bar sizes\n\
\n\
-b1   show protein backbone only (backbone is derived from PDB)\n\
-b2   show protein backbone only (general backbone algorithm)\n\
-b-1,-b-2 as above and omit hetero-molecules\n\
-i#   set the initial number of bonds to show (default=guess, grows on demand)\n\
\n\
if compiled with #define SHELL:\n\
-j#   show water shell thick # AA around 1st molecule, cf. hotkeys jJ\n\
-J$$  name of water ($$=substring of atom ID)\n\
-Jb   use name via cook/molcfg [default=HOH]\n\
";
