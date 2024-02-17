/* \make blend

see `blendpar.h' for VERSION number !

Read the MACSIMUS manual first!

This is the main for `blend'.  There is no `blend.h'.

To compile the standard version, do
  makemake linux lj gcc
  make blend
otherwise see `metamake'

History:
2001--2024  small updates
1999        PROSIS renamed to MACSIMUS
1996        Prague V1.7
08/1996     Odense V1.6g
07/1996     Evanston V1.6f
04/1996     Guelph V1.6d
12/1995     Evanston V1.6c
11/1995     Odense V1.5
10/1993     Odense V0.663

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

#include "ground.h"

#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "blendmed.h"
#include "blendedt.h"
#include "blendgen.h"
#include "options.h"
#include "rndgen.h"

#ifdef X11
#  include "xdraw.h"
#else /*# X11 */
#  error "blend from V2.3 cannot be compiled without X11"
#endif /*#!X11 */

#ifdef POLAR
#  define THISVERSION VERSION"-POLAR"
#else /*# POLAR */
#  define THISVERSION VERSION
#endif /*#!POLAR */

int measure; /* to make intrapot.c happy */

/***
list of options:
  use x for undefined option, default value for valid option.
options are not case sensitive!
-X -X+ -X1 set option X to value 1
-X- -X0    clears option X to value 0
***/

#define x badoption
#define Q -9 /* inactive option */
int optionlist[32] =
/*   `  a   b  c d e f g h i j k   l    m n o p   q r s t u v w   x y z { |   }   ~  ? */
  {100,100,100,0,0,0,0,0,0,5,0,3,256,-100,0,0,1,100,0,x,0,0,0,0,100,3,0,0,Q,100,100,-6};
/*   @  A   B  C D E F G H I J K   L    M N O P   Q R S T U V W   X Y Z [ \   ]   ^  _ */
#undef x

static void AnotA(char *arg) /**************************************** AnotA */
{
  if (arg[2]) ERROR(("option %s : arg not allowed",arg))
}

static int iframe(char *arg) /*************************************** iframe */
{
  int i=atoi(arg);

  if (i==0) ERROR(("no frame 0: 1st frame is 1"))
  else if (i>0) i--;

  return i;
}

static void zero_opt_mwp(void) /******************************* zero_opt_mwp */
{
  if (option('w')) prt("! option -w set to zero");
  if (option('m')) prt("! option -m set to zero");
  if (option('p')) prt("! option -p set to zero");
  option('p')=option('w')=option('m')=0;
}

int main(int narg,char **arg) /*************************************** main */
{
  int i,nmol=0,che,output_file=0;
  species_t *spec=NULL;
  char *outfn=NULL;
  int plbframe=0; /* frame # for plb read  (option -r4), see also Xopt */
    /* NOTE:  1=1st, ...  -1=last, ...; 0=1st */

#include "blendhlp.c"

#if CHECKHEAP==-2
  AllocRange=64;
#endif /*# CHECKHEAP==-2 */

  {
    static char title[82]="blend";
    static char icon[32];
    int i;
    char *t=title+5,*c;

#  include "blendbit.c"

    xwindowhints.icowidth=blendbit_width;
    xwindowhints.icoheight=blendbit_height;
    xwindowhints.icobits=blendbit_bits;

    /* args pre-parsing: -o name and geometry */
    loop (i,1,narg) {
      /* default -v is -v3 with -o, otherwise -v0 */
      if (arg[i][0]=='-' && arg[i][1]=='o') option('v')=3;
      if (arg[i][0]!='-' || arg[i][1]=='o') {
        *t++=' '; c=arg[i];
        while (t<title+80 && *c) *t++=*c++; }
      else if (arg[i][1]=='g')
        if ( (c=strpbrk(arg[i]," _")) ) xwindowhints.geometry=c+1; }

    *t=0;

    if (t==title+80) copy(t-3,"..",3);
    xwindowhints.winname=title;
    copy(icon,title,32);
    if (icon[31]!=0) copy(icon+29,"..",3);
    xwindowhints.iconame=icon;
  }

  initscroll(0);
  repeatprefix('`');
  rndinit(7,0);

  /*** print help if no argument ***/
  if (narg<2) {
    prtsfill(Pusage);

    for (;;) {
      char s[8];

      prts_(Phelp);
      fgets(s,8,stdin);
      switch (s[0]) {
        case 'i': prtsfill(Pusage); break;
        case 'b': prtsfill(Pgen_opt); break;
        case 'g': prtsfill(Pgui); break;
        case 'f': prtsfill(Ppar_opt); break;
        case 'r': prtsfill(Pread_opt); break;
        case 'w': prtsfill(Pwrite_opt); break;
        case 'o': prtsfill(Pmol_opt); break;
        case 'x': prtsfill(Pextra_mol_opt); break;
        case 'm':
        case 'l': prtsfill(Pmol_files); break;
        case 'G': prtsfill(Pgeo); break;
        default: return 0; } } }

  {
    char *name;
    /* extended ... */
    if ((name=getenv("BLENDGEOMETRY"))) xwindowhints.geometry=name;
    if ((name=getenv("BLENDINIT"))) addkbdstring((unsigned char*)name);
  }

  /*** argument parsing + reading mol+cfg ***/

  loop (i,1,narg) {

    if (arg[i][0]=='-') {

      if (isupper(arg[i][1])) switch (arg[i][1]) {
        case 'A': Xopt.A=1; AnotA(arg[i]);
#if 0 /* OBSOLETE */
                  if (option('d')) {
                    WARNING(("all_dihedrals (option -d) changed to 0 because of option -A"))
                    option('d')=0; }
#endif /*# 0 */
                  zero_opt_mwp();
                  break;
        case 'C': Xopt.C=atoi(arg[i]+2); break;
        case 'D': Xopt.D=1;
                  AnotA(arg[i]);
                  zero_opt_mwp();
                  break;
        case 'E': Xopt.E=toupper(arg[i][2]);
                  if (!Xopt.E) Xopt.E=' ';
                  if (!strchr("ABCH",Xopt.E)) ERROR(("bad -E%c option",Xopt.E))
                  zero_opt_mwp();
                  break;
        case 'F': Xopt.F=atoi(arg[i]+2);
                  if (abs(Xopt.F)<2) ERROR(("-F%d too small",Xopt.F))
                  break;
        case 'G': Xopt.G=1;
                  zero_opt_mwp();
                  break;
        case 'H': Xopt.I=1;
                  zero_opt_mwp();
                  break;
        case 'I': addkbdstring((const unsigned char*)(arg[i]+2));
                  break;
        case 'J': Xopt.Jeps=fabs(atof(arg[i]+2));
                  if (Xopt.Jeps==0) Xopt.Jeps=1e-11;
                  if (option('v')&4) Xopt.Jeps=-Xopt.Jeps;
                  break;
        case 'M': Xopt.amplitude=atof(arg[i]+2);
                  if (Xopt.amplitude==0) Xopt.amplitude=1;
                  break;
        case 'N': Xopt.N=1;
                  Xopt.dr=atof(arg[i]+2);
                  if (Xopt.dr==0) Xopt.dr=1e-5;
                  break;
        case 'P': Xopt.P=atoi(arg[i]+2); break;
        case 'R': Xopt.core=atof(arg[i]+2); break;
        case 'S': Xopt.S=atoi(arg[i]+2); break;
        case 'T': Xopt.T=atof(arg[i]+2); break;
        case 'V': Xopt.V=atoi(arg[i]+2); break;
        case 'W': Xopt.W=atoi(arg[i]+2); break;
        default: ERROR(("%s is bad option",arg[i])) }
      else

      getoption(arg[i]) /* see options.h */

      switch (arg[i][1]) {
        case 'q':
          chargewarned=1;
          break;
        case 'o':
          if (output_file) ERROR(("misplaced or duplicate -o"))
          ble_file=output_file=option('o')=1;
          if (arg[i][2]) {
            alloc(outfn,strlen(arg[i])+3); /* + extension maybe */
            strcpy(outfn,arg[i]+2); }
          else if (arg[++i]) {
            alloc(outfn,strlen(arg[i])+5); /* + extension maybe */
            strcpy(outfn,arg[i]); }
          else
            ERROR(("option -o without file name"))
          break;
        case 'z':
          rndinit(7,option('z'));
          break;
        case 'r': {
          char *ch=strchr(arg[i],':');
          if (ch) {
            plbframe=iframe(++ch);
     	    Xopt.ropt=0;
            Xopt.bdihedral[0]=atof(ch);
            ch=strchr(ch,':');
            if (ch) {
	      Xopt.ropt=1;
              Xopt.toframe=iframe(++ch);
	      Xopt.bdihedral[1]=atof(ch);
              if ( (ch=strchr(ch,':')) ) {
                Xopt.byframe=Xopt.bdihedral[2]=atof(++ch);
                if (option('r')<5) if (Xopt.byframe<=0)
                  ERROR(("illegal by frame %d",Xopt.byframe))
                if ( (ch=strchr(ch,':')) ) Xopt.fn=ch+1; } } } }
          break;
        case 'w': {
          char *ch=strchr(arg[i],':');
          if (ch) PDBstyle=atoi(ch+1); } } }

    else /* no -option ==> file name */ {
      char *fn,*c;
      enum { NONE,PAR,BIN,CHE,MOL,PCH,REA } ext;

      if (option('o')) {
        /* this means that previous option was `-oNAME' or `-o NAME': open .ble */
        char *c;
        int hasext=1;
        time_t xt;

        for (c=outfn; *c; c++)
          if (*c=='/') hasext=1;     /* if the directory name contains '.' */
          else if (*c=='.') hasext=0; /* extension */
        if (hasext) strcat(outfn,".ble");
        out=fopen(outfn,"wt");
        time(&xt);
        prt_("! This file was written by MACSIMUS/BLEND V" VERSION " on %s",ctime(&xt)); }

      option('o')=0;
      if (out==NULL) ERROR(("cannot write to %s",outfn?outfn:"stdout"))
      output_file=1; /* output file MUST be opened here */

      /* making a copy: some functions alter extension,
         some functions however use arg[i] directly */
      alloc(fn,strlen(arg[i])+5);
      strcpy(fn,arg[i]);
      /* file name extension: assumed 3 letters long */
      c=fn+strlen(fn)-4;

      ext=NONE;
      if (c>fn) {
        if      (!strcmp(c,".par")) ext=PAR;
        else if (!strcmp(c,".bin")) ext=BIN;
        else if (!strcmp(c,".che")) ext=CHE;
        else if (!strcmp(c,".mol")) ext=MOL;
        else if (!strcmp(c,".pch")) ext=PCH;
        else if (!strcmp(c,".rea")) ext=REA; }

      che=0;
      switch (ext) {

        case BIN:
          if (!rwBIN(fn,NULL)) /* read bin */
            ERROR(("parameter file %s not found (check BLENDPATH)",fn))
          break;

        case PAR:
          readPAR(fn);
          if (option('b')==1)
            rwBIN(fn,c); /* make bin */
          break;

        case PCH:
          if (PCHname)
            ERROR(("%s %s multiply defined partial charge file",PCHname,arg[i]))
          PCHname=arg[i];
          break;

        case REA:
          if (REAname)
            ERROR(("%s %s multiply defined reaction file",REAname,arg[i]))
          REAname=arg[i];
          break;

        case MOL:
          spec=readMOL(fn);
          if (!spec) ERROR(("mol-file %s not found",fn))
          break;

        case NONE:
          strcat(fn,".mol");
          spec=readMOL(fn);
          if (spec) break;
          strcpy(fn,arg[i]); strcat(fn,".che");
        case CHE:
          spec=readCHE(fn);
          che++;
          if (!spec) ERROR(("che-file %s not found",fn))
          break;
        }

#ifdef POLAR
      scalepolarizabilities(option('~')*0.01,option(']')*0.01);
#endif /*# POLAR */

      scalebonds();
      scaleangles();

      if (spec) {
        int or=abs(option('r'));

        nmol++;
        spec->frame=plbframe;
        switch (or) {
          case 0: if (read3D(spec,1)) break; /* try .plb; new in V2.1a */
                  if (read3D(spec,3)) break; /* try .pla; new in V2.1a */
                  if (read3D(spec,0)) break; /* try .3db (deprecated) */
                  if (read3D(spec,2)) break; /* try .3dt (deprecated) */
                  if (che) { che=0; break; }
                  if (read2D(spec,"rb")) break; /* try .2db (deprecated) */
                  if (read2D(spec,"rt")) break; /* try .2db (deprecated) */
          case 1: randomcfg(spec); che=0; break;
          case 5: alphahelixold(spec); che=0; break;
          case 6: alphahelix(spec); che=0; break;
          case 7: makechain(spec); che=0; break;
          case 2:
          case 12: if (!che) if (!read2D(spec,or==2?"rb":"rt")) goto rerr;
                  che=0; break;
          case 3:
          case 4:
         case 13: if (read3D(spec,or==3?0:or==4?1:2)) break;
            rerr: ERROR(("no %s",spec->fn))
                  break;
          default: ERROR(("option -r%d bad value",option('r')))
          }
        if (che>0) WARNING(("%s: 3D cfg with *.che is suspicious",fn))

        autobonds(spec,option('d'));

        chiralities(spec,option('c'));
        partialcharges(spec);

        if (spec->wrmol<0) {
          if (!spec->opt_w) ERROR(("\
%s: mol-file has been renumbered and will be (re)written\n\
but no configuration will be saved because of option -w0",spec->fn))
          else if (abs(spec->opt_w)>1) WARNING(("\
%s: mol-file has been renumbered and will be (re)written\n\
but only text (.3dt) configuration will be saved: use -r13 next time",spec->fn))

          } } /* if (spec) */

      free(fn); } }


  if (option('b')!=1) {
    if (nmol<=0) WARNING(("no molecular file specified - nothing to do"))
    if (nspec<=0) WARNING(("no molecule found - nothing to do")) }
  if (nspec<=0) return 0;

  _n

  if (abs(option('_'))<4 || abs(option('_'))>12)
    WARNING(("option -_%d not in recommended range [4..12]",option('_')))
  fillprec();

  /* list all sites */
  listofsites();

  if (Xopt.V) {
    switch (nspec) {
      case 1: newemptyspec(); appendspec(spec0->next,spec0);
              *(spec0->ext)=0; prt("%s duplicated",spec0->fn);
      case 2: newemptyspec(); appendspec(spec0->next->next,spec0); appendspec(spec0->next->next,spec0->next);
              *(spec0->ext)=0; prt("%s+%s merged",spec0->fn,spec0->next->fn);
              break;
      case 3: prt("virial warning: the 3rd species specified explicitly, 3=1+2 assumed");
              break;
      default: ERROR(("virial: bad number of species")); }
    if (nspec!=3) ERROR(("nspec !=3"))
    spec0->Xopt.V=0;
    spec0->next->Xopt.V=0;
    spec0->next->next->opt_m=0; }

  /* build force field for all species (molecules) */
  for (spec=spec0; spec; spec=spec->next) {
    react(spec,REAname);
    if (spec->edit) moledit(spec);

    build(spec);

    if (option('i')) {
      prt("! max angle energy = %g kcal/mol",Uanglemax);
      if (Uanglemax>option('i')) WARNING(("too large angle energy (%g kcal/mol) around site %d\n\
*** (check, e.g., pyramidal error of tetrahedral angles)",Uanglemax,ianglemax)) }

    /* force fix charge if not rounded and option -_ */
    if (option('_')<0 && fabs(roundedcharge(spec,-1))>1e-13) spec->wrmol++;

    if (spec->wrmol) {
      if (spec->opt_p || spec->opt_w) writeMOL(spec);
      else
          if (!spec->probe.ns)
          WARNING(("%s: mol-file not updated because neither -w nor -p",spec->fn)) }
    if (spec->opt_p) write3D(spec,3);
    if (spec->opt_w) write3D(spec,0); }

  _n
  if (nspec!=newnspec) prt("! molecule(s) split: new nspec=%d",newnspec);
  fclose(out);

  /* this is really very nasty patch: */
  if (ble_file) if (nspec!=newnspec) {
    char s[256];
    FILE *oout;

    /* # of species changed because of molecule split: rewrite it in .ble file */
    if (rename(outfn,tempname()) | !(oout=fopen(tempname(),"rt")))
      ERROR(("cannot rename %s to %s",outfn,tempname()))
    out=fopen(outfn,"wt");
    if (!out) ERROR(("cannot rewrite %s",outfn))
    while (fgets(s,256,oout)) {
      /* rewrite the correct # of species */
      if (!memcmp("nspec=",s,6)) sprintf(s+6,"%d !fixed\n",newnspec);
      prts_(s); }
    fclose(out);
    unlink(tempname()); }

  return 0;
}
