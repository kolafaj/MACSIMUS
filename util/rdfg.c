/* make rdfg
   02/2024 option -v added
   10/2022 bug (write to closed file) fixed
   05/2020 options -s, -l added
   03/2019 RDFGGEOMETRY added
   01/2017 better plot [KILL ALL] support
   04/2016 Reverse removed
   07/2012 extended to coord.n. and Kirkwood-Buff integrals
           unused options removed/simplified
   07/1999 bug fixed: coord.number was 1/2 of correct value for the same atoms
*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"
#include <unistd.h>

#include "rdf.h"

#define NFILE 4

struct outtype_s {
  int yes;
  char *ext;
  char sys[10240];
  char *info;
  char *hdr;
  char *fmt;
} outtype[NFILE] = {
  {0,".g","plot","radial distribution function","g(r)","%8.4f"},
  {0,".hist","plot","histogram (per 1 measurement)","G(r)","%9.2f"},
  {0,".cn","plot","running (cumulative) coordination number","N(r)","%9.3f"},
  {0,".kb","plot","running Kirkwood-Buff integral","K(r)","%8.4f"},
};

struct sites_s {
  struct sites_s *next;
  char *name,*nname; /* nname active for a paor only */
} *shead,*lhead,*xshead,*xlhead;

struct sites_s *getlist(char *str,char *pairsep) /****************** getlist */
/*
   str is 1 before the ,-separated list starts
   pairsep applies for pairs (option -l,-l-)
   newly allocated list returned
*/
{
  char *c=strdup(str),*x,*sep;
  struct sites_s *l,*head=NULL;

  for (;;) {
    *c++=0;
    if (!*c) break;
    x=strchr(c,',');
    if (x) *x=0;
    else x=c+strlen(c);

    allocone(l);
    l->name=c;
    if (pairsep) {
      if ( (sep=strpbrk(l->name,"-.")) ) {
        l->nname=sep+1;
        *sep=0; }
      else
        l->nname=l->name; }
    else
      l->nname=NULL;
    c=x;
    l->next=head;
    head=l; }

  return head;
}

int main(int narg, char **arg) /*************************************** main */
{
  rdf_t *rdf;
  double q,g,r,coord=-1,rmin=9e99,V,Vref=0;
  int ir,iarg,plot=0,bin=0;
  enum outtype_e { G,HIST,CN,KB } iout;
  char *rfmt="%6.3f";
  char site0[8],site1[8];
  double rho0=1,rho1=1;
  long sumhist;
  int first=1,nplots=0,nsel=0;
  int pid=getpid();
  char *geometry=getenv("RDFGGEOMETRY");

  static char fn[256],fnmin[256]="(unknown)";

  char *dot;

  if (!geometry) geometry="720x480";

  initscroll(0);

  if (narg<2) {
    fprintf(stderr,"\
Calculate and plot radial distribution functions, coordination numbers, and\n\
Kirkwood-Buff integrals from cook binary files SIMNAME.rdf. Call by:\n\
  rdfg [OPTIONs] SIMNAME[.rdf] [OPTIONs]\n\
Options:\n\
  -c[FMT] calculate running coordination numbers SIMNAME.SITE1.SITE2.cn [%s]\n\
  -g[FMT] calculate RDFs SIMNAME.SITE1.SITE2.g [set output format, df.=%s]\n\
  -h[FMT] print histogram SIMNAME.SITE1.SITE2.hist [%s]\n\
  -k[FMT] calc. running Kirkwood-Buff integrals SIMNAME.SITE1.SITE2.kb [%s]\n\
          (binary mixture: print hints to get the partial molar volumes)\n\
  -lSITE1.SITE2[,SITE1.SITE2]  select given pairs only (any order)\n\
                               SITE = SITE.SITE, SITE1-SITE2 = SITE1.SITE2\n\
  -l-SITE1.SITE2[,SITE1.SITE2] exclude pairs\n\
          if none of -s,-l (without -) given, the default is select all\n\
          -s-,-l- applies after -s,-l or if any (or excludes from all)\n\
          -s,-l,-s-,-l- cannot repeat (the last one applies)\n\
  -p      plot the generated files (command `plot' accepting -pPID needed:\n\
          if more than 1 plot then [KILL ALL] kills only them)\n\
  -rFMT   format for r [default=%s]\n\
  -sSITE[,SITE..]              select only functions with SITEs present\n\
  -s-SITE[,SITE..]             exclude SITEs\n\
  -vV     set <V> (only if <V>=0, as after free b.c. calculations)\n\
\n\
\n\
Environment Variable:\n\
  RDFGGEOMETRY = size of window exported to plot in the form XxY [720x480]\n\
Example:\n\
  RDFGGEOMETRY=600x300 rdfg simul.rdf -l-H4-H5,H4-H6 -r%%7.4f -g -p\n\
See also:\n\
  plb2rdf rdfdrop rdfconv gblock\n\
  rdfgv1 (old version) smoothg (out of order)\n",
            rfmt,
            outtype[G].fmt,
            outtype[HIST].fmt,
            outtype[CN].fmt,
            outtype[KB].fmt);
    exit(1); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') switch (arg[iarg][1]) {
      case 'c':
        nplots++;
        outtype[2].yes++;
        if (arg[iarg][2]) outtype[2].fmt=arg[iarg]+2;
        break;
      case 'g':
        nplots++;
        outtype[0].yes++;
        if (arg[iarg][2]) outtype[0].fmt=arg[iarg]+2;
        break;
      case 'h':
        nplots++;
        outtype[1].yes++;
        if (arg[iarg][2]) outtype[1].fmt=arg[iarg]+2;
        break;
      case 'k':
        nplots++;
        outtype[3].yes++;
        if (arg[iarg][2]) outtype[3].fmt=arg[iarg]+2;
        break;
      case 'l':
        if (arg[iarg][2]=='-') xlhead=getlist(arg[iarg]+2,".-");
        else lhead=getlist(arg[iarg]+1,".-");
        break;
      case 'r':
        rfmt=arg[iarg]+2;
        break;
      case 's':
        if (arg[iarg][2]=='-') xshead=getlist(arg[iarg]+2,NULL);
        else shead=getlist(arg[iarg]+1,NULL);
        break;
      case 'p':
        plot++;
        break;
      case 'v':
        Vref=atof(arg[iarg]+2);
        break;
      default:
        ERROR(("%s is unknown option",arg[iarg])) }
    else {
      if (fn[0]) ERROR(("two file names, %s and %s\n\
*** note that the syntax of rdfg changed: to calculate rdfs and plot them, use\n\
***   rdfg %s -g -p",fn,arg[iarg],fn))

      strcpy(fn,arg[iarg]); }

  if (!fn[0]) ERROR(("no RDF file name"))

  /* this is likely unsafe */
  remove(string("/tmp/plot.%d",pid));

  dot=strend(fn);
  if (strlen(fn)>4 && !strcmp(dot-4,".rdf")) dot-=4;

  strcpy(dot,".rdf");
  VarOpen(fn,"r");
  fprintf(stderr,"rdfg:");
  while (VarFile.size) {
    int4 vs=VarFile.size,i4,go;
    double d8;
    struct sites_s *l;

    if (vs==4) {
      VarRead(&i4,4);
      fprintf(stderr," nsites=%d",i4); }
    if (vs==8) {
      VarRead(&d8,8);
      fprintf(stderr," grid=%g",d8); }

    if (vs>8) {
      rdf=VarGetSds();
      if (rdf->size<sizeof(rdf_t)-2*sizeof(unsigned4) || rdf->size>8000000)
	ERROR(("rdf->size=%d is suspicious - endian?",rdf->size))

      if (rdf->nhist*sizeof(rdf->hist[0])+((char*)rdf->hist-(char*)rdf) != rdf->size)
        WARNING(("This file file was likely obtained by an old version of cook (< V2.7f)\n\
*** use  rdfconv  to convert it to the new version"))

      if (!rdf->npair)
        go=0;
      else if (shead==NULL && lhead==NULL)
        go=1;
      else {

        go=0;

        looplist (l,lhead) {
          if (!strcmp(rdf->name[0],l->name) && !strcmp(rdf->name[1],l->nname) ) go=1;
          if (!strcmp(rdf->name[0],l->nname) && !strcmp(rdf->name[1],l->name) ) go=1; }

        looplist (l,shead) {
          if (!strcmp(rdf->name[0],l->name)) go=1;
          if (!strcmp(rdf->name[1],l->name)) go=1; } }

      looplist (l,xlhead) {
        if (!strcmp(rdf->name[0],l->name) && !strcmp(rdf->name[1],l->nname) ) go=0;
        if (!strcmp(rdf->name[0],l->nname) && !strcmp(rdf->name[1],l->name) ) go=0; }

      looplist (l,xshead) {
        if (!strcmp(rdf->name[0],l->name)) go=0;
        if (!strcmp(rdf->name[1],l->name)) go=0; }

      if (go) {
        nsel++;

	loop (ir,0,rdf->nhist) if (rdf->hist[ir]) break;
        r=(ir+0.5)/rdf->grid;
        /* 1st nonzero g(r) */
        if (r<rmin) {
          rmin=r;
          sprintf(fnmin,"%s(%d)-%s(%d) hist[%d]=%d",
                  rdf->name[0],rdf->ns[0],
                  rdf->name[1],rdf->ns[1],
                  ir,rdf->hist[ir]); }

        if (outtype[KB].yes) bin++;

        loop (iout,0,NFILE) if (outtype[iout].yes) {
          strcpy(dot,"."); strcat(fn,rdf->name[0]);
          strcat(fn,".");  strcat(fn,rdf->name[1]);
          strcat(fn,outtype[iout].ext);
          if (plot) {
            if (first) {
              if (nplots>1)
                strcat(outtype[iout].sys,
                       string(" -p%d -g%s+%d+%d %s", pid,geometry,iout*12,iout*12,fn));
              else
                strcat(outtype[iout].sys,
                       string(" -g%s+%d+%d %s",geometry,iout*12,iout*12,fn));  }
            else
              strcat(outtype[iout].sys,string(" %s",fn)); }

          out=fopen(fn,"wt");

          if (rdf->V==0) {
            if (Vref==0) {
              Vref=1e3;
              WARNING(("<V>=0 and no option -v, 1e3[AA3] used")) }
            V=Vref; }
          else {
            if (Vref==0)
              V=rdf->V/rdf->nmeas;
            else {
              static int pass=0;
              V=Vref;
              if (pass==0) WARNING(("option -v replaces <V>"))
              pass=1; } }

          prt("### %s ###\n# %s:%s\n\
#%ld sites of %s[%ld]  %ld sites of %s[%ld]  %.0f pairs\n\
#%ld measurements  <V>=%f  grid=%ld\n\
#number densities [AA^-3]: rho%s=%g  rho%s=%g",
              outtype[iout].info,
              rdf->name[0],rdf->name[1],
              (long)rdf->ns[0],rdf->name[0],(long)rdf->indx[0],
              (long)rdf->ns[1],rdf->name[1],(long)rdf->indx[1],
              rdf->npair,
              (long)rdf->nmeas,V,(long)rdf->grid,
              rdf->name[0],rdf->ns[0]/V,
              rdf->name[1],rdf->ns[1]/V);

          if (iout==KB && rdf->indx[0]!=rdf->indx[1]) {
            strcpy(site0,rdf->name[0]);
            strcpy(site1,rdf->name[1]);
            rho0=rdf->ns[0]/V;
            rho1=rdf->ns[1]/V; }

          prt("# r       %s",outtype[iout].hdr);

          q = V/(4*PI*rdf->npair*(double)rdf->nmeas/Cub(rdf->grid));

          sumhist=0;

          loop (ir,0,rdf->nhist) {

            if (iout>=CN) {
              r=(ir+1)/rdf->grid;
              sumhist += rdf->hist[ir];
              g=(double)sumhist/rdf->nmeas;
              if (iout==CN) {
                /* coordination number */
                if (rdf->indx[0]==rdf->indx[1])
                  coord=g*=2./rdf->ns[0];
                else {
                  coord=g/rdf->ns[1];
                  g/=rdf->ns[0]; } }
              else {
                /* Kirkwood-Buff */
                g*=V/rdf->ns[0]/rdf->ns[1];
                if (rdf->indx[0]==rdf->indx[1]) g*=2;
                g-=4*PI/3*Cub(r); } }
            else {
              /* RDF */
              r=(ir+0.5)/rdf->grid;
              if (iout==G) g=rdf->hist[ir]*q/(ir*(ir+1.0)+1.0/3);
              else g=rdf->hist[ir]/rdf->nmeas; }

            prt_(rfmt,r);
            prt_(" ",r);
            prt_(outtype[iout].fmt,g);
            if (iout==CN) prt_(outtype[iout].fmt,coord);
            _n }
          fclose(out);
          out=stdout; }

        first=0;
        free(rdf); } } }

  if (nsel) fprintf(stderr,"%d functions selected\n",nsel);
  else Error("no function selected - nothing to do");

  if (plot) loop (iout,0,NFILE) if (outtype[iout].yes) {
    putenv(string("PLOTNAME=%s",outtype[iout].info));
    strcat(outtype[iout].sys," &");
    fprintf(stderr,"\n%s\n",outtype[iout].sys);
    if (system(outtype[iout].sys)) fprintf(stderr,"... failed!\n");
    sleep(1); }

  fprintf(stderr," minr=%g for %s\n",rmin,fnmin);

  ir=0;
  loop (iout,0,NFILE) ir+=outtype[iout].yes;
  if (!ir)
    fprintf(stderr,"no output file because none of -g -h -c -k specified (rdfg w/o args for help)\n");

  *dot=0;
  if (bin==3)
    printf("\
Hints to determine partial molar volumes of a binary mixture\n\
WARNING: tests passed with problems - bug in the r^3 term?\n\
export a=%.9g\n\
export b=%.9g\n\
mergetab \\\n\
  %s.%s.%s.kb:1:2 \\\n\
  %s.%s.%s.kb:1=1:2 \\\n\
  %s.%s.%s.kb:1=1:2 | \\\n\
  tabproc A '(1+b*(D-C))/(a+b+a*b*(B+D-2*C))' > %s.V%s\n\
mergetab \\\n\
  %s.%s.%s.kb:1:2 \\\n\
  %s.%s.%s.kb:1=1:2 \\\n\
  %s.%s.%s.kb:1=1:2 | \\\n\
  tabproc A '(1+a*(B-C))/(a+b+a*b*(B+D-2*C))' > %s.V%s\n\
plot %s.V%s:A:'B+a*A^3' %s.V%s\n",
           rho0,rho1,
           fn,site0,site0,
           fn,site0,site1,
           fn,site1,site1,
             fn,site0,
           fn,site0,site0,
           fn,site0,site1,
           fn,site1,site1,
             fn,site1,
           fn,site0,fn,site1);

  return 0;
}
