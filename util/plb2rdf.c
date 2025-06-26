/* cc -O2 -Wall -o plb2rdf plb2rdf.c -lm 
 */
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"

typedef float fvector[3];

int OFF1=-1,NS1=-1,NM1=-1;
int OFF2=-1,NS2=-1,NM2=-1;
char *key1="i0",*key2="i0";
char KEY1,KEY2;         // =key[0]
#define NI 1024   // max. # of atom indices in the key or atoms in a molecule
int indx1[NI],indx2[NI];
double weight1[NI],weight2[NI];

double RANGE=15,cq;     // RDF range
double GRID=0;    // z-histogram grid (per unity)
double DR=0;      // z-histogram grid (delta z); DZ=1/GRID
int NH;           // ((int)(GRID*RANGE)) // # of z-histogram bins

double *P0,*P1;

int getindx(char *arg,int *indx,double *weight) /******************* getindx */
{
  int i=0;
  char *c;

  for (c=arg; *c; ) {
    c++;
    if (!*c) break;
    if (i>=NI) Error("to many indices in arg");
    indx[i]=(int)(weight[i]=atof(c));
    i++;
    c=strchr(c,','); 
    if (!c) break; }

  return i;
}

int main(int narg,char **arg) /**************************************** main */
{
  FILE *plb;
  char name[256]="";
  int iarg,nf=0,ns,from=1,to=0x7fffffff;
  float hdr[2];
  fvector L,*R,*r1,*r2;
  double Vav=0,cq,rr;
  int zerohist=0;
  int i,n1,n2,j,homo;
  int nindx1,nindx2;
  vector Lh,dr;

  if (narg<3) {
    fprintf(stderr,"\
RDFs from trajectory, periodic b.c.  Call by:\n\
  plb2rdf [OPTIONS] PLB-FILE [OPTIONS] > OUTPUT.g\n\
Files:\n\
  PLB-FILE input plb-file (new format with L required)\n\
  OUTPUT.g output RDF (r,g(r)) with #-info lines\n\
Options:\n\
  -dDR     r-grid, in histogram bin width (only one of -g, -d accepted)\n\
  -fFIRST  first frame in NAME.plb analyzed [default=1]\n\
  -gGRID   r-grid, in histogram bins per 1A [default=1]\n\
  -kKEY    defines the vector for analysis of angular correlations, KEY:\n\
    iA1,A2,... list of sites to make average [default=i0]\n\
               A1,A2,... are atom numbers taken from the mol-file\n\
    mM0,M1,..,MNS  weighed average; CM if M are masses\n\
  -mM      number of molecules of given kind in a configuration (frame)\n\
  -nNS     number of sites in the molecule measured\n\
  -oOFFSET where given block of M molecules start (in # of sites) [0]\n\
           (options -k,-m,-n,-o may repeat to denote the second molecule)\n\
  -rRANGE  max distance [%g]\n\
  -tLAST   last  frame in NAME.plb analyzed [default=EOF]\n\
  -vVERBOSE print also zero values of RDF (good in scripts)\n\
Example (g_OO for 3-site water, Oxygen=1):\n\
  plb2rdf water.plb -ki1 -m500 -n3 -g25 | plot -\n\
See also:\n\
  densprof rdfdrop plb2diff plbmsd plb2hist plb2cryst hbonds plb2dens plb2sqd\n\
  plbinfo plbconv plbcut plbfilt plb2asc asc2plb\n",
            RANGE);
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'd': DR=atof(arg[iarg]+2); break;
        case 'f': from=atoi(arg[iarg]+2); break;
        case 'g': GRID=atof(arg[iarg]+2); break;
        case 'k': if (key1) key2=key1; key1=arg[iarg]+2; break;
        case 'm': NM2=NM1; NM1=atoi(arg[iarg]+2); break;
        case 'n': NS2=NS1; NS1=atoi(arg[iarg]+2); break;
        case 'o': OFF2=OFF1; OFF1=atoi(arg[iarg]+2); break;
        case 'r': RANGE=atof(arg[iarg]+2); break;
        case 't': to=atoi(arg[iarg]+2); break;
        case 'v': zerohist++; break;
        default: Error("unknown option"); }
    else 
      strcpy(name,arg[iarg]);
      
  if (!name[0]) Error("no NAME");

  if (NS1<=0) Error("number of sites undefined, use option -n");
  if (NM1<=0) Error("number of molecules undefined, use option -m");
  if (OFF1<0) OFF1=0;
  if (NS2<=0) NS2=NS1;
  if (NM2<=0) NM2=NM1;
  if (OFF2<0) OFF2=OFF1;
  if (!key2) key2=key1;
  cq=Sqr(RANGE);
  homo=OFF1==OFF2 && !strcmp(key1,key2);
  if (homo && NM1!=NM2) Error("homoatomic RDF and NM1!=NM2"); 
  if (homo && NS1!=NS2) Error("homoatomic RDF and NS1!=NS2"); 

  if (from<1) Error("FROM<1");
  if (to<from) Error("empty frame range");

  if (DR!=0) {
    if (GRID!=0) Error("both -g and -z specified");
    GRID=1/DR; }
  if (!GRID) GRID=100;
  
  if (!(plb=fopen(name,"rb"))) Error(name);
  nindx1=getindx(key1,indx1,weight1);
  nindx2=getindx(key2,indx2,weight2);
  KEY1=key1[0];
  KEY2=key2[0];
 
  if (fread(hdr,4,2,plb)!=2) Error("file too short (no header)");
  ns=hdr[0];
  if (hdr[1]!=-3) Error("only new (L3) format of plb-files is accepted");
  fseek(plb,(ns+1)*(from-1)*sizeof(fvector),SEEK_CUR);
  if (OFF1+NM1*NS1>ns)
    Error("selected sites do not exist in the configuration\n\
*** check options -o,-n,-m");
  if (OFF2+NM2*NS2>ns)
    Error("selected sites (2nd molecule) do not exist in the configuration\n\
*** check options -o,-n,-m");

  allocarray(R,ns*3);

  switch (tolower(KEY1)) {
    case 'm': if (nindx1!=NS1) Error("KEY=m: bad number of weights (masses)"); break;
    case 'i': if (nindx1<1) Error("KEY=i: no site"); break;
    default: Error("bad KEY"); }

  switch (tolower(KEY2)) {
    case 'm': if (nindx2!=NS2) Error("KEY=m: bad number of weights (masses)"); break;
    case 'i': if (nindx2<1) Error("KEY=i: no site"); break;
    default: Error("bad KEY"); }

  printf("# plb2rdf");
  loop (iarg,1,narg) printf(" %s",arg[iarg]);
  printf("\n#  r     g_i%c(r)","ji"[homo]);

  NH=(int)(GRID*RANGE);
  allocarrayzero(P0,NH+1);
//allocarrayzero(P1,NH);

  for (; from<=to; from++) {
    double V;
   
    if (fread(L,sizeof(fvector),1,plb)!=1) goto end;
    if (fread(R,sizeof(fvector),ns,plb)!=ns) Error("unexpected EOF");

    V=L[0]*L[1]*L[2];
    VV(Lh,=0.5*L)
    Vav+=V;

    loop (n1,0,NM1) {
      vector c1;
      double m1=0;
      int TO=NM2;

      if (homo) TO=n1;

      r1=R+OFF1+n1*NS1;
      VO(c1,=0) \

      switch (KEY1) {
        case 'i': 
          loop (j,0,nindx1) VV(c1,+=r1[indx1[j]])
          VO(c1,/=nindx1);
          break;
        case 'm':
          loop (j,0,nindx1) {
            VV(c1,+=weight1[j]*r1[j])
            m1+=weight1[j]; }
          VO(c1,/=m1)
          break; }

      loop (n2,0,TO) {
        vector c2;
        double m2=0;

        r2=R+OFF2+n2*NS2;
        VO(c2,=0) \

        switch (KEY2) {
          case 'i': 
            loop (j,0,nindx2) VV(c2,+=r2[indx2[j]]);
            VO(c2,/=nindx2);
            break;
          case 'm':
            loop (j,0,nindx2) {
              VV(c2,+=weight2[j]*r2[j])
              m2+=weight2[j]; }
            VO(c2,/=m2)  
          break; }
        
        VVV(dr,=c2,-c1)
        loop (j,0,3) {
          while (dr[j]<-Lh[j]) dr[j]+=L[j];
          while (dr[j]>Lh[j]) dr[j]-=L[j]; }
        rr=SQR(dr);
        if (rr<cq) {
          int ind=(int)(sqrt(rr)*GRID);
          P0[ind]++; } } }

    nf++; }
 end:

  loop (i,0,NH) {
    double r=(i+0.5)/GRID;
    double DV=(Cub(i+1)-Cub(i))*4*PI/3;
    double rdf=P0[i]*Vav/(Sqr(nf)*DV*NM1*NM2)*Cub(GRID)*(homo+1);
    if (rdf||zerohist) printf("%8.5f %9.6f\n",r,rdf); }    

  printf("# %d frames analyzed, <V>=%g\n",nf,Vav/nf);

  fclose(plb);

  return 0;
}
