/* cc -O2 -o plb2nbr plb2nbr.c -lm
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  typedef double vector[3];
  vector *r,L={0,0,0};
  int i,j,N,frame,oldformat;
  double dist,R=1;
  FILE *mol,*plb;
  FILE *in;
  float f[3];
  static int hist[1000];
  char *COLORS="wROYGCBM";
  char MORECOLOR='W';

  if (narg<3) {
    fprintf(stderr,"\
Read plb file and sort sites to files according to # of neighbors.\n\
Call by:\n\
  plb2nbr PLBFILE DIST [DIAM [COLORS [MORECOLOR]]]\n\
DIST = distance to define neighbors\n\
DIAM = diameter [default=1]\n\
COLORS = string of {ROYGCBMWroygcbmw}, 1st char = no neighbor, etc.;\n\
         Red Orange Yellow Green Cyan Blue Magenta White;\n\
         lowercase = tiny sphere [default=wROYGCBM]\n\
MORECOLOR = color if there are more neighbors than in COLORS [default=W]\n\
see also: c/shownbr hbonds plb2cryst\n");
    exit(0); }

  in=fopen(arg[1],"rb");
  if (!in) Error(arg[1]);
  dist=atof(arg[2]);
  if (narg>3) R=atof(arg[3]);
  if (narg>4) COLORS=arg[4];

  if (fread(f,4,2,in)!=2) Error(".plb file too short");

  N=f[0];
  if (N<1 || N>16777216) Error("bad N in .plb file");
  oldformat=f[1]>=0;
  if (oldformat) L[0]=L[1]=L[2]=f[1];

  allocarrayzero(r,N);

  for (frame=0;;frame++) {

    loop (i,0,1000) hist[i]=0;

    if (!oldformat) {
      if (3!=fread(f,4,3,in)) Error(".plb file too short");
      L[0]=f[0]; L[1]=f[1]; L[2]=f[2]; }

    loop (i,0,N) {
      if (3!=fread(f,4,3,in)) Error(".plb file too short");
      r[i][0]=f[0]; r[i][1]=f[1]; r[i][2]=f[2]; }

    mol=fopen(string("nbr%04d.mol",frame),"wt");
    fprintf(mol,"neighbors\n\
parameter_set=unknown\n\
number_of_atoms=%d\n\
atoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",N);

    plb=fopen(string("nbr%04d.plb",frame),"wt");
    f[0]=N; f[1]=-3;
    fwrite(f,4,2,plb);
    f[0]=L[0]; f[1]=L[1]; f[2]=L[2];
    fwrite(f,4,3,plb);

    loop (i,0,N) {
      int nnbr=-1; /* i-i will count for 1 */
      char col;
      double rad=R/2;

      loop (j,0,N) {
	double dd;
	if (L[0])
	  dd=Sqr(L[0]/2-fabs(L[0]/2-fabs(r[i][0]-r[j][0])))
	    +Sqr(L[1]/2-fabs(L[1]/2-fabs(r[i][1]-r[j][1])))
	    +Sqr(L[2]/2-fabs(L[2]/2-fabs(r[i][2]-r[j][2])));
	else
	  dd=Sqr(r[i][0]-r[j][0])+Sqr(r[i][1]-r[j][1])+Sqr(r[i][2]-r[j][2]);

	nnbr+=dd<Sqr(dist); }

      if (nnbr>=1000) Error("more than 1000 neighbors");
      hist[nnbr]++;

      if (nnbr<strlen(COLORS)) col=COLORS[nnbr];
      else col=MORECOLOR;
      if (islower(col)) col=toupper(col),rad/=2;
      //      fprintf(gol,"%c %f\n",col,rad);
      fprintf(mol,"%d %c%03d/HS HS-%c 0 0 0\n",i,col,(int)(rad*100+0.5),col);
      f[0]=r[i][0];
      f[1]=r[i][1];
      f[2]=r[i][2];
      fwrite(f,4,3,plb); }

    fclose(plb);
    fclose(mol);

    printf("%d:",frame);
    loop (i,0,1000) if (hist[i]) printf(" %d=%d",i,hist[i]);
    printf("\n"); }

  return 0;
}
