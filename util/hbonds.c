/* cc -O2 -o hbonds hbonds.c -lm
 */
#include "../gen/include.h"

int NS; /* sites in one water molecule */
#define NB 12 /* max # of H-bonds per molecule */
int H1,H2,O,M,LL,OFF;
int *rentab; /* [NS] */
int periodic; /* 1 for periodic b.c., 0 for vacuum b.c. */

int ren(int i) /******************************************************** ren */
/* renumber NS-site water to 3-site water */
{
  return ((i-OFF)/NS)*3+rentab[(i-OFF)%NS]+OFF;
}

typedef float vector[3];
vector *r;
float L[3];
double Lh[3];
FILE *plbin;
int varL,ns,nbonds;
struct site_s {
  int nnbr;
  int nbr[NB];
} *s;

void addbond(int i,int j) /***************************************** addbond */
{
  if (i>=ns || j>=ns) Error("index out of range (bad OFF?, bad SITES?)");
  if (s[i].nnbr>=NB) Error("too many bonds");
  s[i].nbr[s[i].nnbr++]=j;
  if (s[j].nnbr>=NB) Error("too many bonds");
  s[j].nbr[s[j].nnbr++]=i;
  nbonds++;
}

double distq(int i,int j) /******************************************* distq */
{
  int k;
  double rr=0;

  if (periodic)
    loop (k,0,3)
      rr+=Sqr(fabs(fabs(r[i][k]-r[j][k])-Lh[k])-Lh[k]);
  else
    loop (k,0,3)
      rr+=Sqr(r[i][k]-r[j][k]);

  return rr;
}

double cosangle(int O1,int H,int O2) /***************************** cosangle */
  /* cosangle=cos(PI-angle(O1-H-O2)) */
{
  int k;
  double rr1=0,rr2=0,scal=0,dr1,dr2;

  if (periodic)
    loop (k,0,3) {
      dr1=r[O1][k]-r[H][k]; if (dr1>Lh[k]) dr1-=L[k]; if (dr1<-Lh[k]) dr1+=L[k];
      dr2=r[O2][k]-r[H][k]; if (dr2>Lh[k]) dr2-=L[k]; if (dr2<-Lh[k]) dr2+=L[k];
      rr1+=Sqr(dr1);
      rr2+=Sqr(dr2);
      scal+=dr1*dr2; }
  else
    loop (k,0,3) {
      dr1=r[O1][k]-r[H][k];
      dr2=r[O2][k]-r[H][k];
      rr1+=Sqr(dr1);
      rr2+=Sqr(dr2);
      scal+=dr1*dr2; }

  // prt("%g cosalpha",-scal/sqrt(rr1*rr2));
  
  return -scal/sqrt(rr1*rr2);
}

int readplb(void) /************************************************* readplb */
{
  int k;

  if (varL) {
    if (fread(L,12,1,plbin)!=1) return 0;
    loop (k,0,3) Lh[k]=L[k]/2; }
  if (fread(r,12,ns,plbin)!=ns) return 0;

  return 1;
}

int main(int narg,char **arg) /**************************************** main */
{
  int i,j,frame,writeplb,nmol;
  FILE *plb,*mol,*gol,*hb;
  float hdr[2];
  double ho,hoq,ca;
  char *sem;
  char *name,*dot,*SITES=arg[2];
  double S=0,SQ=0;
  
  if (narg<3) {
    fprintf(stderr,"\
Generate H-bonds for a model of water. Call by:\n\
  hbonds NAME.plb [-]SITES[:OFF] [[-]O-H_DIST [ANGLE]\n\
both old and new plb formats accepted, output type=input type\n\
output = series of files NAME.####.{plb,mol} for show in 3-site water format\n\
         files NAME.####.hb water numbers in ACCEPTOR DONOR order\n\
ARGUMENTS:\n\
  SITES = string of {H,O,M,L} defining the water model\n\
  -SITES = as above, without files NAME.####.{plb,mol}\n\
  OFF = offset (# of atoms other than water at the beginning)\n\
  O-H_DIST = H-bond distance threshold (default=2.4), periodic b.c.\n\
  -O-H_DIST = as above, ignore periodic b.c.\n\
  angle = PI-angle(O-H-O) threshold, in degrees (default=no threshold)\n\
Example:\n\
  hbonds solution.plb HOMH:12 -2.4\n\
See also:\n\
  hbonds4 (old version)\n\
");
    exit(0); }

  plbin=fopen(arg[1],"rb");
  if (!plbin) Error(arg[1]);
  name=strdup(arg[1]);
  dot=strend(name)-4;
  if (strcmp(dot,".plb")) Error("extension .plb expected");
  *dot=0;

  if (*SITES=='-') writeplb=0,SITES++;
  else writeplb=1;

  if (narg>3) ho=atof(arg[3]); else ho=2.4;
  periodic=ho>0;
  ho=fabs(ho);
  hoq=ho*ho;

  if (narg>4) ca=atof(arg[4]); else ca=0;
  ca=cos(ca*PI/180);

  if (fread(hdr,4,2,plbin)!=2) Error("plb too short");
  ns=(int)hdr[0];
  if (ns<1 || ns>100000) Error("too many sites");
  L[0]=L[1]=L[2]=hdr[1];
  loop (i,0,3) Lh[i]=L[i]/2;
  varL=hdr[1]<0;

  allocarray(r,ns);
  allocarray(s,ns);

  if ( (sem=strchr(SITES,':')) ) {
    OFF=atoi(sem+1);
    NS=sem-SITES; }
  else {
    OFF=0;
    NS=strlen(SITES); }

  allocarray(rentab,NS);
  j=0;
  loop (i,0,NS) {
    rentab[i]=j;
    switch (SITES[i]) {
      case 'H': H1=H2; H2=i; j++; break;
      case 'L': LL=i; break;
      case 'M': M=i; break;
      case 'O': O=i; j++; break; } }

  //  printf("H1=%d H2=%d M=%d O=%d\n",H1,H2,M,O);
  
  for (frame=0; readplb(); frame++) {

    loop (i,0,ns) s[i].nnbr=0;
    nmol=0;
    for (i=OFF; i<ns; i+=NS) {
      nmol++;
      addbond(i+H1,i+O);
      addbond(i+H2,i+O); }

    if (writeplb) hb=fopen(string("%s.%04d.hb",name,frame),"wt");

    nbonds=0;

    for (i=OFF; i<ns; i+=NS)
      if (ca==1) {
        for (j=OFF; j<ns; j+=NS) if (i!=j) {
          if (distq(i+H1,j+O)<hoq) {
            addbond(i+H1,j+O);
            if (writeplb) fprintf(hb,"%d %d\n",j/NS,i/NS); }
          if (distq(i+H2,j+O)<hoq) {
            addbond(i+H2,j+O);
            if (writeplb) fprintf(hb,"%d %d\n",j/NS,i/NS); } } }
      else
        for (j=OFF; j<ns; j+=NS) if (i!=j) {
          if (distq(i+H1,j+O)<hoq && cosangle(j+O,i+H1,i+O)>ca) {
            addbond(i+H1,j+O);
            if (writeplb) fprintf(hb,"%d %d\n",j/NS,i/NS); }
          if (distq(i+H2,j+O)<hoq && cosangle(j+O,i+H2,i+O)>ca) {
            addbond(i+H2,j+O);
            if (writeplb) fprintf(hb,"%d %d\n",j/NS,i/NS); } }

    if (writeplb) fclose(hb);

    printf("%d %d frame nbonds\n",frame,nbonds);
    S+=nbonds;
    SQ+=Sqr(nbonds);

    if (writeplb) {
      mol=fopen(string("%s.%04d.mol",name,frame),"wt");
      //    gol=fopen(string("%s.%04d.gol",name,frame),"wt");
      plb=fopen(string("%s.%04d.plb",name,frame),"wb");
      //    fprintf(gol,"!\n!\n%d\n",(ns-OFF)/NS*3+OFF);
      fprintf(mol,"\
water HB structure\n\
parameter_set = water\n\
number_of_atoms = %d\n\
\n\
atoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",(ns-OFF)/NS*3+OFF);
      hdr[0]=(ns-OFF)/NS*3+OFF;
      fwrite(hdr,4,2,plb);
      if (varL) fwrite(L,4,3,plb);

      loop (i,0,ns) {
        if (i<OFF) {
          fprintf(mol,"%3d G124/%d-a X 0 0 0\n", i,i);
          // fprintf(gol,"GREEN 1.24\n");
          fwrite(r[i],4,3,plb); }
        else if (strchr("HO",SITES[i%NS])) {
          fprintf(mol,"%3d %s/%d-%c%c%03d %c 0 0 %d",
                  ren(i),
                  SITES[i%NS]=='H'?"W056":"R124",
                  i,SITES[i%NS],OFF?'b':'a',(ns-OFF)/NS,SITES[i%NS],s[i].nnbr);
          loop (j,0,s[i].nnbr) fprintf(mol," %d",ren(s[i].nbr[j]));
          fprintf(mol,"\n");
          //  if (SITES[i%NS]=='H') fprintf(gol,"WHITE 0.56\n");
          //  else fprintf(gol,"RED 1.24\n");
          fwrite(r[i],4,3,plb); } }
      fclose(mol);
      // fclose(gol);
      fclose(plb); }

  }

  S/=nmol;
  SQ/=Sqr(nmol);
  
  printf("# %7.4f %7.4f H-bonds_per_water stderr\n",
         S/frame,sqrt((SQ/frame-Sqr(S/frame))/(frame-1)/NS));

  return 0;
}
