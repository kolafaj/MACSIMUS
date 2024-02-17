/* make rdfgv1
OLD VERSION for cook < V2.7f
see rdfg.c for the new version (good for large systems)
7/99 bug fixed: coord.number was 1/2 of correct value for the same atoms
*/
#define xOLDFORMAT /* old format of *.rdf files -- changes since !!! */

#include "ground.h"
#define SDS
#include "alloc.h"
#include "varfile.h"

#ifdef Reverse
#define VarReverse VarFile.reverse
#else
int VarReverse; /* dummy */
#endif

#ifdef OLDFORMAT

/* old version (June 1993) */
/* for site-site radial distribution function for one pair of sites */
typedef struct /* SDS */ {
  int size;      /* sds */
  int npair;     /* # of pairs of given type in the system */
  int nmeas;     /* # of measurements on the system */
  double grid;   /* 1/grid is step in the adopted units */
  double V;      /* sum of volumes (to work eventually for NPT) */
  int nhist;     /* length of array hist */
  int hist[1];   /* [nhist] */
  } rdf_t;

#else

#include "rdf.h"

#endif

int main(int narg, char **arg)
{
  rdf_t *rdf;
  double q,g,r,coord,rmin=9e99;
  char fnmin[64];
  int ir;
  char key=0,plot=0;
  char sys[10240];
  char *fmt="%8.4f";
  long sumhist;
  char *plotcommand;
#ifdef DOS  
  char *arg2="d";
#else  
  char *arg2="u";
#endif  

  char fn[128],*dot;

  plotcommand=getenv("PLOTCOMMAND");
  if (!plotcommand) plotcommand="plot";

  initscroll(0);

  prt("*** g(r) from *.rdf *** (c) J.Kolafa ***");
#ifdef OLDFORMAT
  prt("*** old format ***");
#endif

  if (narg<2) {
    fprintf(stderr,"\
Write (and plot) ascii *.g files from binary SIMNAME.rdf file.  Call by:\n\
  %s SIMNAME[.rdf] [KEY [FMT]]\n\
KEY = string of:\n"
#ifndef DOS
"  u  (as UNIX), separate files SIMNAME.SITE1.SITE2.g\n"
#endif
"  d  (as DOS), separate files SITE1SITE2.g\n\
  o  one merged file SIMNAME.g\n\
  m  print min. distances (hint: pipe to sort -n)\n\
  p  start `plot' (or command in env. PLOTCOMMAND) of generated *.g file(s)\n\
  r  reverse endian on input\n\
The default is KEY=u FMT=%%8.4f\n\
FMT=format for g(r) (default=%%8.4f)\n\
Example (generate files with DOS-names, prec.=5 dec.digits, plot them):\n\
  %s simul.rdf dp %%9.5f\n",arg[0],arg[0]);
  exit(1); }

  if (narg>2) arg2=arg[2];
  
  if (strchr(arg2,'o')) key=0;
  if (strchr(arg2,'d')) key='D';
  if (strchr(arg2,'u')) {
#ifdef DOS
    Error("UNIX file format not supported for 16 bit DOS application");
#endif
    key='U'; }
  if (strchr(arg2,'r')) {
#ifdef OLDFORMAT
    ERROR(("option (r)everse endian not supported with OLDFORMAT"))
#endif
#ifdef Reverse
    VarReverse=1;
#else
    ERROR(("option (r)everse endian but no compile-time Reverse"))
#endif
    }
  if (strchr(arg2,'p')) plot=1;
  
  if (narg>3) fmt=arg[3];
  
  fprintf(stderr,"g(r) format = %s",fmt);

  strcpy(fn,arg[1]); dot=strend(fn);
  if (strlen(fn)>4 && !strcmp(dot-4,".rdf")) dot-=4;

  if (plot) {
    strcpy(sys,plotcommand);
    if (!key) {
      strcpy(dot,".g");
      out=fopen(fn,"wt");
      strcat(sys," ");
      strcat(sys,fn); } }
  else {
    if (!key) {
      fprintf(stderr,"\ngnuplot\n");
      strcpy(dot,".g");
      out=fopen(fn,"wt");
      fprintf(stderr,"p \'%s\' w l\n",fn); } }

  strcpy(dot,".rdf");
  VarOpen(fn,"r");
  while (VarFile.size) {
    int4 vs=VarFile.size,i4;
    double d8;

    if (vs==4) {
      VarRead(&i4,4); ReverseOrder(&i4,4);
      fprintf(stderr," nsites=%d",i4); 
      if (!plot) if (key) fprintf(stderr,"\ngnuplot\n"); }
    if (vs==8) {
      VarRead(&d8,8); ReverseOrder(&d8,8);
      fprintf(stderr," grid=%g",d8); }

    if (vs>8) {
      rdf=VarGetSds();
      if (rdf->size<sizeof(rdf_t)-2*sizeof(unsigned4) || rdf->size>8000000)
	ERROR(("rdf->size=%d is suspicious - endian?",rdf->size))

      if (VarReverse) {
	int i;

	ReverseOrder(rdf->indx,sizeof(rdf->indx[0]));
	ReverseOrder(rdf->indx+1,sizeof(rdf->indx[0]));
	ReverseOrder(rdf->ns,sizeof(rdf->ns[0]));
	ReverseOrder(rdf->ns+1,sizeof(rdf->ns[0]));
	ReverseOrder(&rdf->npair,sizeof(rdf->npair));
	ReverseOrder(&rdf->grid,sizeof(rdf->grid));
	ReverseOrder(&rdf->V,sizeof(rdf->V));
	ReverseOrder(&rdf->nmeas,sizeof(rdf->nmeas));
	ReverseOrder(&rdf->nhist,sizeof(rdf->nhist));
	loop (i,0,rdf->nhist)
	  ReverseOrder(rdf->hist+i,sizeof(rdf->hist[0])); }


      if (rdf->npair) {
	if (key) {
	  if (key=='D') {
	    strcpy(fn,rdf->name[0]);
	    strcat(fn,rdf->name[1]); strcat(fn,".g"); }
	  else {
	    strcpy(dot,".");	strcat(fn,rdf->name[0]);
	    strcat(fn,".");	strcat(fn,rdf->name[1]);
	    strcat(fn,".g"); }
	  if (plot) {
	    strcat(sys," "); strcat(sys,fn); }
	  else
	    fprintf(stderr,"p \'%s\' w l\n",fn);
	  out=fopen(fn,"wt"); }

	prt("\n# %s:%s\n\
# %ld sites of %s[%ld]  %ld sites of %s[%ld]  %ld pairs\n\
# %ld measurements  <V>=%f  grid=%ld",
#ifdef OLDFORMAT
#error no longer valid
#else
	    rdf->name[0],rdf->name[1],
	    (long)rdf->ns[0],rdf->name[0],(long)rdf->indx[0],
	    (long)rdf->ns[1],rdf->name[1],(long)rdf->indx[1],
#endif
	    (long)rdf->npair,(long)rdf->nmeas,rdf->V/rdf->nmeas,(long)rdf->grid);

	header("# r       g(r)     hist  coord1 [coord2] #");

	q = rdf->V/(4*PI*rdf->npair*Sqr((double)rdf->nmeas)/Cub(rdf->grid));
	sumhist=0;

	loop (ir,0,rdf->nhist) if (rdf->hist[ir]) break;

	if (strchr(arg2,'m'))
	  printf("%6.3f %s-%s\n",
		 (ir+0.5)/rdf->grid,
		 rdf->name[0],rdf->name[1]);
	if ((ir+0.5)/rdf->grid<rmin) {
	  rmin=(ir+0.5)/rdf->grid;
	  sprintf(fnmin,"%s(%d sites) %s(%d sites) hist[%d]=%d",
		  rdf->name[0],rdf->ns[0],
		  rdf->name[1],rdf->ns[1],
		  ir,rdf->hist[ir]); }

	loop (ir,0,rdf->nhist) {
	  r=(ir+0.5)/rdf->grid;
	  g = rdf->hist[ir]*q/(ir*(ir+1.0)+1.0/3);
	  sumhist += rdf->hist[ir];
	  if (sumhist || rdf->hist[ir+1]) {
	    prt_("%6.3f ",r);
	    prt_(fmt,g);
	    prt_(" %7i",rdf->hist[ir]);
	    coord=(double)sumhist/rdf->nmeas;
	    if (rdf->indx[0]==rdf->indx[1]) 
	      prt_("%8.3f ",2*coord/rdf->ns[0]);
	    else {
	      prt_("%8.3f ",coord/rdf->ns[0]);
	      prt_("%8.3f ",coord/rdf->ns[1]); }
	    _n }
          }
	header("");
	_n
        if (key) fclose(out); }
      else
	fprintf(stderr,"\n%s-%s: 0 pairs",
		rdf->name[0],rdf->name[1]);

      free(rdf); } }

  if (plot) {
#ifndef DOS
    strcat(sys," &");
#endif
    fprintf(stderr,"\n%s\n",sys);
    system(sys); }
  else
    fprintf(stderr,"\n");

  fprintf(stderr,"min distance=%.3f for %s\n",rmin,fnmin);

  return 0;
}
