/* cc -O2 -o plb2cryst plb2cryst.c -lm
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  typedef double vector[3];
  vector *r;
  int i,j,N,frame,varL;
  double dist,R=1,alpha=170,cosa;
  FILE *mol,*gol,*plb;
  FILE *in;
  float f[3],L[3];
  int NRED=0,SCALE=1;
  static int hist[1000];

  if (narg<3) {
    fprintf(stderr,"Read plb file and sort sites to files according to crystal-like structure.\n\
Call by:\n\
  plb2cryst PLBFILE DIST {ALPHA|-SMARTLIM} [DIAM [[-]NRED [SCALE]]]\n\
DIST = max. distance to define neighbors\n\
ALPHA = min. (nbr)-(i)-(nbr) angle, in deg [170]\n\
SMARTLIM = as above, algorithm with weight\n\
DIAM = diameter [default=1]\n\
NRED = # of neighbors shown RED (then, MAGENTA BLUE CYAN GREEN YELLOW WHITE)\n\
       [default=0]\n\
-NRED : nbrs < NRED are RED, nbrs > NRED+DN*6 are WHITE (default: modulo 7)\n\
SCALE : nbrs in [NRED,NRED+SCALE) are RED, etc. [default=1]\n\
see also: findcryst c/shownbr plb2nbr\n");
    exit(0); }

  in=fopen(arg[1],"rb");
  if (!in) Error(arg[1]);
  dist=atof(arg[2]);
  if (narg>3) alpha=atof(arg[3]);
  if (narg>4) R=atof(arg[4]);
  if (narg>5) NRED=atoi(arg[5]);
  if (narg>6) SCALE=atoi(arg[6]);

  if (R<=0) Error("negative R");

  cosa=cos(PI/180*alpha);
  fread(f,4,2,in);
  
  N=f[0];
  if (N<1 || N>1000000) Error("bad N");
  L[0]=f[1];  varL=L[0]<0;
  if (!varL) L[1]=L[2]=L[0];

  allocarrayzero(r,N);

  for (frame=0;;frame++) {

    loop (i,0,1000) hist[i]=0;

    if (varL) if (3!=fread(L,4,3,in)) exit(0);
    loop (i,0,N) {
      if (3!=fread(f,4,3,in)) exit(0);
      r[i][0]=f[0]; r[i][1]=f[1]; r[i][2]=f[2]; }

    mol=fopen(string("nbr%04d.mol",frame),"wt");
    fprintf(mol,"neighbors\n\
parameter_set=unknown\n\
number_of_atoms=%d\n\
atoms\n",N);

    gol=fopen(string("nbr%04d.gol",frame),"wt");
    fprintf(gol,"!\n!\n%d\n",N);

    plb=fopen(string("nbr%04d.plb",frame),"wb");
    f[0]=N; f[1]=-3;
    fwrite(f,4,2,plb);
    fwrite(L,4,3,plb);

    loop (i,0,N) {
      int npairs;
      char *col="WHITE";
      int nnbrs=0;
      int n1,n2,k;
      vector rnbr[N];
      vector Rnbr[N];
      double s;

      loop (j,0,N) {
        double dd;

        if (L[0])
          dd=Sqr(L[0]/2-fabs(L[2]/2-fabs(r[i][0]-r[j][0])))
            +Sqr(L[1]/2-fabs(L[2]/2-fabs(r[i][1]-r[j][1])))
            +Sqr(L[2]/2-fabs(L[2]/2-fabs(r[i][2]-r[j][2])));
        else
          dd=Sqr(r[i][0]-r[j][0])+Sqr(r[i][1]-r[j][1])+Sqr(r[i][2]-r[j][2]);

        if (dd<Sqr(dist)) {
          if (i!=j) {
            loop (k,0,3) {
              rnbr[nnbrs][k]=r[j][k]-r[i][k];
              if (rnbr[nnbrs][k]>L[k]/2) rnbr[nnbrs][k]-=L[k];
              if (rnbr[nnbrs][k]<-L[k]/2) rnbr[nnbrs][k]+=L[k];
	      Rnbr[nnbrs][k]=rnbr[nnbrs][k];
              rnbr[nnbrs][k]/=sqrt(dd); }
            nnbrs++; } } }
    
      npairs=0;

      if (alpha>0) 
        loop (n1,0,nnbrs)
          loop (n2,0,n1) {
            s=0;
	    loop (k,0,3) s+=rnbr[n1][k]*rnbr[n2][k];
	    if (s<cosa) npairs++; }
      else {
	double dd[N],err1,err2,err3;

        loop (n1,0,nnbrs) {
	  dd[n1]=Sqr(Rnbr[n1][0])+Sqr(Rnbr[n1][1])+Sqr(Rnbr[n1][2]);

	  loop (n2,0,n1) {
	    s=0;
	    
	    loop (k,0,3) s+=rnbr[n1][k]*rnbr[n2][k];

	    err1=(dd[n1]-1)/(Sqr(dist)-1);
	    err2=(dd[n2]-1)/(Sqr(dist)-1);
	    err3=(s+1)/(cosa+1);
	    if (Sqr(err1)+Sqr(err2)+2*err3<4) npairs++; } } }


      fprintf(mol,"%d HS-%02d HS 0 0 0 0\n",i,npairs);
      switch (npairs) {
        case -1: case -2: case -3: case -4:
        case 0: col="RED"; break;
        case 1: col="MAGENTA"; break;
        case 2: col="BLUE"; break;
        case 3: col="CYAN"; break;
        case 4: col="GREEN"; break;
        case 5: col="YELLOW"; break;
        case 6: col="WHITE"; break;
        default: col="RED"; }

      if (npairs>=1000) Error("more than 1000 neighbors/straight pairs");
      hist[npairs]++;

      npairs+=1001*SCALE-abs(NRED);
      npairs/=SCALE;
      if (NRED<0) {
        npairs-=1001;
        if (npairs<0) npairs=0;
        if (npairs>6) npairs=6; }
      else npairs%=7;
      if (npairs<0) Error("");

      col="R\0M\0B\0C\0G\0Y\0W"+npairs*2;
      fprintf(gol,"%s %f\n",col,R/2);
      f[0]=r[i][0];
      f[1]=r[i][1];
      f[2]=r[i][2];
      fwrite(f,4,3,plb); }

    fclose(plb);
    fclose(mol);
    fclose(gol); 

    printf("%d:",frame);
    loop (i,0,1000) if (hist[i]) printf(" %d=%d",i,hist[i]);
    printf("\n"); }

  return 0;
}
