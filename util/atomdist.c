/* cc -Wall -O2 -o atomdist atomdist.c -lm
 */
#include "../gen/include.h"

typedef float vector[3];
#define VVV(A,B,C) { A[0] B[0] C[0]; A[1] B[1] C[1]; A[2] B[2] C[2]; }
#define SQR(A) (A[0]*A[0]+A[1]*A[1]+A[2]*A[2])

int main(int narg,char **arg) /**************************************** main */
{
  int i,j,i1,i2,ns=-1,n=0,k;
  int varL=0;
  FILE *plb;
  float hdr[2];
  vector L,r,r1={0,0,0},r2={0,0,0};
  vector *cfg=NULL;

  double dr[3];
  double rr,sumrr=0;

  if (narg<4) {
    fprintf(stderr,"\
Atom-atom distances (periodic b.c.) from playback files. Call by:\n\
  atomdist FILE I1 I2\n\
where\n\
  FILE      MACSIMUS playback file (typically with extension .plb)\n\
  I1,I2     atom indices (see the mol-file)\n\
  I1<0      calculates the minimum distance from I2\n\
  I1<0,I2<0 calculates the minimum distance (all pairs)\n\
");
    exit(0); }

  i1=atoi(arg[2]);
  i2=atoi(arg[3]);
  if (i1<0 && i2>=0) k=i1,i1=i2,i2=k; /* from i1>0 */

  plb=fopen(arg[1],"rb");
  if (!plb) Error("atomdist: plb-file not found");
  if (fread(hdr,sizeof(hdr),1,plb)!=1) Error("atomdist: plb-file no no header");
  ns=hdr[0];
  printf("%d sites, ",ns);

  if (hdr[1]>=0)
    printf("fixed box size L=%.8g (old format)\n",L[0]=L[1]=L[2]=hdr[1]);
  else if (hdr[1]<0) {
    varL=1;
    printf("variable box (new format)\n"); }
  else {
    fprintf(stderr,"WRONG PARAMETER L=%g\n",hdr[1]);
    exit(0); }

  if (ns<=0 || ns>1e7) {
    fprintf(stderr,"%d wrong number of sites (bad endian or not a playback file)",ns);
    exit(1); }

  if (i1<0 || i2<0) {
    allocarray(cfg,ns);
    for (;;) {
      double minrr=9e9;

      if (varL) if (1!=fread(L,sizeof(L),1,plb)) goto end;
      if (ns!=fread(cfg,sizeof(r),ns,plb)) goto end;

      if (i1<0 && i2<0) {
        loop (i,0,ns)
          loop (j,0,i) {
            VVV(dr,=cfg[i],-cfg[j])
            if (L[0]+L[1]+L[2]>0)
              loop (k,0,3) {
                while (dr[k]<-L[k]/2) dr[k]+=L[k];
                while (dr[k]>L[k]/2) dr[k]-=L[k]; }
            rr=SQR(dr);
            Min(minrr,rr) }
        printf("%9.6f\n",sqrt(minrr));
        sumrr+=minrr;
        n++; }
      else {
        loop (i,0,ns) if (i-i1) {
          VVV(dr,=cfg[i1],-cfg[i])
          if (L[0]+L[1]+L[2]>0)
            loop (k,0,3) {
              while (dr[k]<-L[k]/2) dr[k]+=L[k];
              while (dr[k]>L[k]/2) dr[k]-=L[k]; }
	  rr=SQR(dr);
	  Min(minrr,rr) }
        printf("%9.6f\n",sqrt(minrr));
        sumrr+=minrr;
        n++; }
    } }
  else
    for (;;) {
      if (varL) if (1!=fread(L,sizeof(L),1,plb)) goto end;
      loop (i,0,ns) {
        if (1!=fread(r,sizeof(r),1,plb)) goto end;
        if (i==i1) memcpy(r1,r,sizeof(r));
        if (i==i2) memcpy(r2,r,sizeof(r)); }
      VVV(dr,=r1,-r2)
      if (L[0]+L[1]+L[2]>0)
        loop (k,0,3) {
          while (dr[k]<-L[k]/2) dr[k]+=L[k];
          while (dr[k]>L[k]/2) dr[k]-=L[k]; }
      rr=SQR(dr);
      printf("%9.6f\n",sqrt(rr));
      sumrr+=rr; n++; }

 end:
  printf("%d cfgs  <rr>=%g  sqrt(<rr>)=%g\n", n,sumrr/n,sqrt(sumrr/n));
}
