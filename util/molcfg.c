/* make molcfg
 */
#include "ground.h"

#define MAXVAL 12
#define LINE 256

int main(int narg, char **arg)
{
  char fn[LINE];
  int off,ns,totns,iarg,i,rep,irep,from=0,golns;
  char line[LINE];
  char id[16],type[16],charge[16],chir[16];
  char *prefix,prefixnr[16],*suffix,suffixnr[16],*c;
  static char parset[32],oldparset[32];
  int ii,nnb,nb[MAXVAL],len;
  int pass;
  FILE *mol,*gol,*cfgmol,*cfggol=NULL;
  int isgol=1;
  int *N=NULL;

  initscroll(0);

  if (getenv("MOLCFG")) in=NULL; /* do not ask what to do in case of error */

  if (narg<3) {
    fputs("\
Writes SIMNAME.mol [optionally also SIMNAME.gol] of a whole configuration\n\
from SPECIES1.mol SPECIES2.mol ... [and SPECIES1.gol..., if available]\n\
To be used by `show' to show a trajectory recorded by cook*\n\
NB: cook* (version V2.3e and later) calls molcfg automatically\n\
ENVIRONMENT: export MOLCFG=anything causes silent rewriting old SIMNAME.mol\n\
\n\
Usage 1:\n\
  molcfg SPECIES1 SPECIES2 ... SIMNAME\n\
repeats SPECIES1.mol [and SPECIES1.gol] N[0] times, SPECIES2.mol N[1] times, ..\n\
where N[0], N[1], ... are read from SIMNAME.def\n\
N[] may be given by formulas, but only limited set of variables is understood\n\
(double: x,y,z,a,b,c,L[3],rho; int: i,j,k,n,init,no,noint); e.g. n=3 N[0]=n^3*4\n\
Example 1:\n\
  molcfg Li Al Cl I salt\n\
\n\
Usage 2:\n\
  molcfg -COUNT1[PREFIX1:SUFFIX1:FROM1] SPECIES1 \\\n\
    [-COUNT2][PREFIX2:SUFFIX2:FROM2] SPECIES2 ... SIMNAME\n\
repeats SPECIES1.mol COUNT1 times, etc.; SIMNAME.def is ignored\n\
PREFIX and SUFFIX are added to atom id and may contain format (as %d) to hold\n\
file number, starting from FROM (default FROM=0)\n\
Example 2:\n\
  molcfg -1 cyto -10w hoh -3:p%d proton -1 etoh config\n\
",stderr);
  exit(0); }

  loop (pass,0,2) {
    rep=1; off=0; prefix=suffix=""; prefixnr[0]=suffixnr[0]=0;
    if (pass) {
      strcpy(fn,arg[narg-1]); strcat(fn,".mol");
      cfgmol=fopen(fn,"rt");
      if (cfgmol) {
        WARNING(("%s exists and will be overwritten",fn))
        fclose(cfgmol); }

      cfgmol=fopen(fn,"wt");
      if (!cfgmol) ERROR(("cannot (re)write %s",fn))

      if (isgol) {
        fn[strlen(fn)-3]='g';
        cfggol=fopen(fn,"wt");
        if (!cfggol) ERROR(("cannot (re)write %s",fn))
        fprintf(cfggol,"!sphere#1\n! %s\n%d\n",fn,totns); }

      fprintf(cfgmol,"configuration:");
      loop (iarg,1,narg) fprintf(cfgmol," %s",arg[iarg]);
      fprintf(cfgmol,"\n\n");
      if (parset[0]) fprintf(cfgmol,"parameter_set = %s\n",parset);
      fprintf(cfgmol,"number_of_atoms = %d\n\n\
atoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",totns); }

    if (arg[1][0]!='-') if (pass==0) {
      int init=0,no=0,noint=0;
      double L[3]={0,0,0};
      double rho=1000;
      double x=0,y=0,z=0,a=0,b=0,c=0;
      int i=0,j=0,k=0,n=0;

      strcpy(fn,arg[narg-1]); strcat(fn,".def");
      in=fopen(fn,"rt");
      if (!in) ERROR(("open %s",fn))
      fprintf(stderr,"reading %s:",fn);
      alloczero(N,100*sizeof(N[0]));

      _Id.echo=-1;
      getdata
        getvec(N,,100)
        get(x) get(y) get(z)
        get(a) get(b) get(c)
        get(i) get(j) get(k) get(n)
        get(rho) getvec(L,,3)
        get(init) get(no) get(noint)
        enddata
      fclose(in);
      if (getenv("MOLCFG")) in=NULL; else in=stdin;
      loop (golns,0,100)
        if (N[golns]) {
          fprintf(stderr," N[%d]=%d",golns,N[golns]);
          if (golns>narg-2) WARNING(("no mol-file for nonzero N[%d]",golns)) }
      fprintf(stderr,"\n"); }

    loop (iarg,1,narg-1) {

      if (N) rep=N[iarg-1];

      if (arg[iarg][0]=='-') {
        if (N) ERROR(("cannot combine # of molecules from def-file with options"))
        if (!isdigit(arg[iarg][1])) ERROR(("%s no count",arg[iarg]))
        rep=-atoi(arg[iarg]);
        for (prefix=arg[iarg]; *prefix && strchr("-0123456789",*prefix); prefix++);
        prefix=strdup(prefix);

        if ( (suffix=strchr(prefix,':')) ) {
          *suffix++=0;
        if ( (c=strchr(suffix,':')) ) *c++=0,from=atoi(c); }
        else suffix=""; }

      else loop (irep,0,rep) {

        strcpy(fn,arg[iarg]); strcat(fn,".mol");
        if (!pass) {
          if (irep==0)
            fprintf(stderr,"molcfg: file %s, count=%d, 1st copy offset=%d, ",fn,rep,off);
          if (irep==rep-1)
          fprintf(stderr,"last copy offset=%d\n",off); }
        mol=fopen(fn,"rt");
        if (!mol) ERROR(("open %s",fn))

        fn[strlen(fn)-3]='g';
        gol=fopen(fn,"rt");
        if (!gol) {
          if (!irep && !pass) fprintf(stderr,"no %s\n",fn);
          isgol=0; }
        else {
          loop (golns,0,3)
            if (!fgets(line,LINE,gol)) ERROR(("%s too short",fn))
          golns=atoi(line); }

        for (;;) {
          if (!fgets(line,LINE,mol)) ERROR(("mol-file too short"))
          ns=0;
          if (sscanf(line,"parameter_set = %s",parset)) {
            if (oldparset[0]) if (strcmp(parset,oldparset))
              fprintf(stderr,"WARNING: parameter_set=%s and %s\n",parset,oldparset);
            strcpy(oldparset,parset); }
          if (sscanf(line,"number_of_atoms = %d",&ns) && ns>0) break; }
        for (;;) {
          if (!fgets(line,LINE,mol)) ERROR(("mol-file too short"))
          if (!memcmp(line,"atoms",5)) break; }
        for (;;) {
          if (!fgets(line,LINE,mol)) ERROR(("mol-file too short"))
          if (line[0]!='!') break; }

        if (gol) if (ns!=golns)
          ERROR(("ns=%d in mol-file, ns=%d in %s",ns,golns,fn))

        loop (i,0,ns) {
	  char *tok=strtok(line,"\n\r\t ");
	  int inb;

	  if (!tok) ERROR(("%s: line too short",fn))
          ii=atoi(tok);
          if (ii!=i) ERROR(("%s: mol-file numbering"))
	  strcpy(id,strtok(NULL,"\n\r\t "));
	  strcpy(type,strtok(NULL,"\n\r\t "));
	  strcpy(charge,strtok(NULL,"\n\r\t "));
      	  strcpy(chir,strtok(NULL,"\n\r\t "));
	  nnb=atoi(strtok(NULL,"\n\r\t "));
	  if (nnb>MAXVAL) ERROR(("%s: # of bonds exceeds MAXVAL=%d",fn,MAXVAL))

          loop (inb,0,nnb) {
	    tok=strtok(NULL,"\n\r\t ");
	    if (!tok) ERROR(("%s: missing nbr in line",fn))
	    nb[inb]=atoi(tok); }

	  if (pass) {
            sprintf(suffixnr,suffix,irep+from);
            sprintf(prefixnr,prefix,irep+from);
            len=11-strlen(suffixnr)-strlen(prefixnr)-strlen(id);
            if (len<0) len=0;
            fprintf(cfgmol,"%3d %s%s%s%*s%-4s %8s %1s %d",
                    off+i,prefixnr,id,suffixnr,-len," ",type,charge,chir,nnb);
            loop (ii,0,nnb) fprintf(cfgmol," %d",nb[ii]+off);
            fprintf(cfgmol,"\n"); }
          if (!fgets(line,LINE,mol) && i<ns-1)
            ERROR(("mol-file too short")) }

        if (gol) {
          loop (i,0,ns) {
            if (!fgets(line,LINE,gol)) ERROR(("gol-file too short"))
            if (pass && cfggol) fputs(line,cfggol); }
          fclose(gol); }

        fclose(mol);
        off+=ns; }
      } /* iarg */

    if (!pass)
      fprintf(stderr,"total %d sites\n", totns=off);
    else
      if (totns!=off) ERROR(("internal (file changed?)"))
    } /* pass */

  fclose(cfgmol);
  if (cfggol) fclose(cfggol);

  return 0;
}
