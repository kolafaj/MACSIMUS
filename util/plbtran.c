/* cc -O2 -Wall -o plbtran plbtran.c -lm
*/
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"

int axis=3;

void rotate(float *r,double angle) /******************************** rotate */
{
  double sa=sin(angle), ca=cos(angle);
  double o[3][3];
  float rr[3];
  int i=(axis+1)%3,j=(i+1)%3,k=(j+1)%3,a;

  o[i][i]=ca; o[i][j]=-sa; o[i][k]=0;
  o[j][i]=sa; o[j][j]= ca; o[j][k]=0;
  o[k][i]= 0; o[k][j]=  0; o[k][k]=1;

  loop (a,0,3) rr[a]=SCAL(o[a],r);

  VV(r,=rr)
}

int shear=0;

void doshear(float *r,double angle) /***************************** doshear */
{
  double sa=sin(angle);
  int i=(axis+1)%3,j=(i+1)%3;

  if (shear>0) r[j]+=sa*r[i];
  else r[i]+=sa*r[j];
}

int *trankey;

void Tran(float r[])
{
  if (trankey) {
    float old[4]={r[0],r[1],r[2],0};
    int i;

    loop (i,0,3) r[i]=old[trankey[i]]; }
}

int main(int narg,char **arg)
{
  char *fn1=NULL,*fn2=NULL;
  int iarg,i,k,f=0,ns,varL=0,periodic=0,recenter=0,mirror=0,glide=0,from=1,screw=0;
  float r[3],center[3]={0,0,0},L[3];
  FILE *plb1,*plb2;
  double angle=0,dangle=0;

#if 0
  int mplL[3]={0,0,0};

  -xX[L]   reference point or vector R=(X,Y,Z), default R=(0,0,0):\n\
  -yY[L]      center of rotation or reflection, or displacement vector\n\
  -zZ[L]      appended L means in the units of the respective box size\n\

#endif

  if (narg<3) {
    fprintf(stderr,"Rotate a plb-file. Call by:\n\
  plbtran [INPUT.plb [OUTPUT.plb]] OPTIONS\n\
OPTIONS:\n\
  -/AXIS   rotation axis or reflection plane, one of {x,y,z} or {0,1,2}\n\
           if not given, only -t may apply\n\
  -aANGLE[d] rotation angle (for 1st frame), in radians (d=degrees)\n\
  -dDANGLE the rotation angle is ANGLE+(frame-1)*DANGLE (1st frame=1)\n\
  -m       mirror reflection (of plane perpendicular to AXIS at R_AXIS)\n\
           can be combined with rotation\n\
  -g       gliding reflection (as -m + displace by the remaining two R_i)\n\
           (rotation is ignored)\n\
  -q       screw (rotation + shift in the AXIS direction)\n\
  -s       shear (similar as rotation; ANGLE should be small)\n\
  -S       as -s, swapped coordinates\n\
  -c       put the rotation center to the origin and set the box size to zero\n\
           (if -g: only the box size to zero)\n\
  -l       normalize periodically to the basic cell (void if -c)\n\
           NB: carefully normalized cfg (crystal) to the basic cell needed!\n\
  -fFRAME  first frame to transform (only -l applies to frames 1..FRAME-1)\n\
  -tXYZ    where XYZ = any permutation of {xyz0}, transformation includes box\n\
           -tzxy means that x:=z, x:=y, y:=z, 0 means 0 assigned\n\
           performed before any other transformation\n\
Example (reflection about origin):\n\
  plbtran cfg.plb rot.plb -/x -a180d\n\
Example (rotate crystal reflection about origin):\n\
  plbtran cfg.plb rot.plb -/x -a180d\n\
See also:\n\
  cook data: load.tr\n\
  plbstack (for rot. by 90 deg and axes swap)\n\
  plbrot (old deprecated version)\n");
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') switch (arg[iarg][1]) {
      case 'a': angle=atof(arg[iarg]+2);
                if (strchr(arg[iarg],'d')) angle*=PI/180; break;
      case 'd': dangle=PI/180*atof(arg[iarg]+2); break;
#if 0
      case 'x': center[0]=atof(arg[iarg]+2); if (strchr(arg[iarg],'L')) mplL[0]=1; break;
      case 'y': center[1]=atof(arg[iarg]+2); if (strchr(arg[iarg],'L')) mplL[1]=1; break;
      case 'z': center[2]=atof(arg[iarg]+2); if (strchr(arg[iarg],'L')) mplL[2]=1; break;
#endif
      case 'l': periodic++; break;
      case 'm': mirror++; break;
      case 'g': glide++; break;
      case 'q': screw++; break;
      case 's': shear++; break;
      case 'S': shear--; break;
      case 'c': recenter++; break;
      case '/': axis=(arg[iarg][2]-'0')&7; break;
      case 'f': from=atoi(arg[iarg]+2); break;
      case 't': {
        static int t[4];
        int i;

        trankey=t;
        loop (i,0,3) {
          t[i]=arg[iarg][2+i]-'x';
          if (t[i]<0 || t[i]>2) t[i]=3; }
        break; }
      default: Error("bad option"); }
    else {
      if (fn1) fn2=arg[iarg];
      else fn1=arg[iarg]; }

  if (fn1 && !strcmp(fn1,fn2)) Error("infile=outfile");
  if (fn1) plb1=fopen(fn1,"rb"); else plb1=stdin;
  if (fn2) plb2=fopen(fn2,"wb"); else plb2=stdout;

  fread(r,4,2,plb1);
  varL=r[1]<0;
  VO(L,=r[1])
  ns=r[0];
  if (ns<1 || ns>16777216) Error("wrong ns - endian?");

  r[1]=-3;
  fwrite(r,4,2,plb2);

  for (f=1;;f++) {
    double a=angle+dangle*f;

    if (varL) if (3!=fread(L,4,3,plb1)) goto end;
    Tran(L);
    if (axis<3) center[axis]=L[axis]/2;

    if (recenter) VO(L,=0) /* not periodic */
    fwrite(L,4,3,plb2);

    loop (i,0,ns) {
      if (3!=fread(r,4,3,plb1)) goto end;
      Tran(r);

      if (axis<3) {
        if (f>=from) {
          if (shear) {
            if (screw) { loop (k,0,3) if (k!=axis) r[k]-=center[k]; }
            else VV(r,-=center)
            doshear(r,a);
            if (mirror) r[axis]=-r[axis];
            if (!recenter) VV(r,+=center) }
          else if (glide) {
            loop (k,0,3)
              if (k==axis) r[k]=center[k]-r[k];
              else r[k]+=center[k]; }
          else {
            if (screw) { loop (k,0,3) if (k!=axis) r[k]-=center[k]; }
            else VV(r,-=center)
            rotate(r,a);
            if (mirror) r[axis]=-r[axis];
            if (!recenter) VV(r,+=center) } }

          if (periodic)
            loop (k,0,3) if (L[k]) {
              while (r[k]<0) r[k]+=L[k];
              while (r[k]>=L[k]) r[k]-=L[k]; } }

      if (3!=fwrite(r,4,3,plb2)) Error("cannot write"); }

  }

 end:
  fclose(plb1);
  fclose(plb2);

  return 0;
}
