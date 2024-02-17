/* cc -O2 -o plbsmooth plbsmooth.c -lm

   history:
   03/2002 V1 (smoothpl)
   10/2022 V2 (variable L, polishing, endiand and no-ns removed)
*/
#include "../gen/include.h"

typedef float real;
#include "../gen/vector3d.h"

int main(int narg, char **arg) /*************************************** main */
{
  vector **r,*auxr,fr,L0;
  double dr[3],Lsm[3];
  double *window;
  float hdr[2];
  int ns,s,smooth,gauss=0,i,j,no=0;
  FILE *f,*ff;
  int tri=0,varL;
  char *plus;

  if (narg<4) {
    fprintf(stderr,"\
Smooth atomic motion. Call by:\n\
  smoothplb IN OUT WINDOW[+GAUSS]\n\
ARGUMENTS:\n\
  IN     input playback file\n\
  OUT    output playback file\n\
  WINDOW rectangular window of length WINDOW: [1/WINDOW,...,1/WINDOW]\n\
  GAUSS  number of [0.5,0.5] windows to convolute\n\
  SITES  # of sites\n\
See also:\n\
  smooth smoothpl(old version)\n\
Warning:\n\
  should support fluctuating box, but not debugged enough\n\
");
  exit(0); }

  f=fopen(arg[1],"rb");
  if (!f) {
    fprintf(stderr,"no plb-file\n");
    exit(1); }

  ff=fopen(arg[2],"wb");

  if (!strcmp(arg[1],arg[2])) {
    fprintf(stderr,"IN and OUT files identical\n");
    exit(-1); }

  /* preparing the window */
  smooth=atoi(arg[3]);
  if (arg[3][0]=='+') tri=1,smooth=atoi(arg[3]+1);

  if (smooth>0) {
    char *ch=strchr(arg[3]+1,'+');
    double sum;
    
    if (ch) gauss=abs(atoi(ch+1));

    allocarray(window,smooth+gauss);
    if (tri)
      loop (i,0,(smooth+1)/2) window[i]=window[smooth-1-i]=i+1;
    else
      loop (i,0,smooth) window[i]=1;
    sum=0;
    loop (i,0,smooth) sum+=window[i];
    loop (i,0,smooth) window[i]/=sum;

    loop (j,0,gauss) {
      window[smooth]=0;
      for (i=smooth; i>0; i--)
        window[i]=(window[i]+window[i-1])/2;
      window[0]=window[0]/2;
      smooth++; } }
  else
    Error("smoothplb: bad value of WINDOW");

  allocarray(r,smooth);

  if (2!=fread(hdr,sizeof(float),2,f)) Error("smoothplb: input file too short");
  fwrite(hdr,sizeof(float),2,ff);
  ns=hdr[0];
  varL=hdr[1]<0;
  if (!varL) L0[0]=L0[1]=L0[2]=hdr[1];

  printf("number of sites=%d\n",ns);
  if ((hdr[0]-ns)!=0 || ns<1 || hdr[0]>16777216.) Error("smoothplb: bad hdr");

  /* allocate and pre-read */
  loop (s,0,smooth) {
    allocarray(r[s],ns+1); /* r[][ns] = L */
    if (varL) {
      if (1!=fread(r[s][ns],sizeof(vector),1,f))
        Error("smoothplb: not enough data (L expected)"); }
    else
      VV(r[s][ns],=L0)

      if (ns!=fread(r[s],sizeof(vector),ns,f))
      Error("smoothplb: not enough data (r expected)"); }

  for (;;) {
    no++;

    if (varL) fwrite(r[smooth-1][ns],sizeof(vector),1,ff);

    VO(Lsm,=0)
    loop (s,0,smooth) VV(Lsm,+=window[s]*r[smooth-1][ns])

    loop (i,0,ns) {
      VO(dr,=0)
      loop (s,0,smooth) VVV(dr,+=window[s]*r[s][i],/r[smooth-1][ns])
      VVV(fr,=dr,*Lsm)
      fwrite(fr,sizeof(vector),1,ff); }

    auxr=r[0];
    loop (s,1,smooth) r[s-1]=r[s];
    r[smooth-1]=auxr;

    if (varL)
      if (fread(auxr+ns,sizeof(vector),1,f)!=1) break;
    if (fread(auxr,sizeof(vector),ns,f)!=ns) break; }

  printf("%d smoothed frames\n",no);
  fclose(ff);
  fclose(f);
  
  return 0;
}
