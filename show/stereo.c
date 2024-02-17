/* make stereo
                      =============================
                       STEREO (c) J. Kolafa 1994-5
                        updated (ppm format) 1998
                        simplified and fixed 2016
                      (old version = util/stereo.c)
                      =============================
DESCRIPTION
^^^^^^^^^^^
This program calculates stereograms normally with the eye cross point behind
the paper, i.e., your eyes should watch an imaginary point lying behind the
paper in approximately the same distance as the eye-paper distance (this is
for the recommended default of the x-pattern size one half of the eye's
distance).  Different people may prefer different distances or even
stereograms with the cross point above the paper.

METHOD
^^^^^^
The programs works with rows of pixels in the x-direction.  It starts from a
periodically repeating `pattern' from left (the periodicity is given by the
condition that both eyes must see the same pattern) and maps the images of
pixels on the paper with their positions on the 3D scene.  The following
figure shows how one element of the pattern denoted by A repeats both on the
paper and in the scene.

            (1)    (2)
**A*****A*****      ******A****
             *      *              scene
             A**A**A*
               ..
  -A--A--A--A-A-A-A--A--A--        paper
             .  .
            .   .
           .    .
          O     O                  eyes

It is seen that the pattern shortens if the scene is closer to the eyes; if we
go from left then a piece of the pattern disappears at point (1) and another
piece of pattern must appear at point (2).  The latter case is a source of
technical problems how to create such a pattern (see  rightedge  and
rightedgep  in input data).

VERSIONS
^^^^^^^^
For some obsolete stuff (direct output to some HP printers), see the old
macsimus/util/stereo.c

INPUT SCENE
^^^^^^^^^^^
The input scene must be given in the form of z-buffer as a PGM (raw
Portable GrayMap, P5, depth=256) file, black=far, white=front.
The first 4 lines of this file must be:
P5
# SCALING (any text)
XSIZE YSIZE
255
Where SCALING is the number transforming gray level (in [0,255]) into the
height in the pixel (voxel) unit

The default extension is ".zbuf", can be obtained from show V2.0 and newer.

OUTPUT
^^^^^^
The output picture is PPM (Portable PixMap, P6, depth=256).

*/

#include "ground.h"
#include "rndgen.h"

#define inc(X) X++ /* translated from Pascal */

int xpat=384, ypat=600; /* pattern size in printer pixels */

typedef unsigned char byte;
typedef byte rgb[3];

rgb **pat; /* [xpat][ypat] */

int xpaper,ypaper;
int shrink=2;

/* hidden algorithm parameters (I no longer remeber what they are good for,
   except phase which fixes a bit the artifact of left-based scanning */
int up=1;
int phase=1,rightedge=6;
float rightedgep=0;

int mc=20;
double gm=1.4;
double contrast=0.5; /* [-4..2] */
int xzbuf,yzbuf;

rgb **sheet;

rgb rndcol;
void makerndcol(void)
{
  int k;
  double x;

  loop (k,0,3) {
    x=rnd();
    x=pow(x,1/gm);
    if (contrast>rnd()) {
      if (x<0.5) x=0; else x=0.999; }
    rndcol[k]=x*256; }
}

void rndpat(void) /* =============================================== rndpat */
/*
   random pattern, with some MC processing making patches
*/
{
  int i,j,k;
  int mcl=mc*xpat*ypat*Sqr(shrink);

  loop (i,0,xpat) loop (j,0,ypat) {
    makerndcol(); memcpy(pat[i][j],rndcol,3); }

  while (mcl--) {
    i=irnd(xpat);
    j=irnd(ypat);
    if (rnd()<0.2) loop (k,0,3) pat[(i+1)%xpat][j][k]=pat[i][j][k];
    else loop (k,0,3) pat[i][(j+1)%ypat][k]=pat[i][j][k]; }
}

void *readpnm(char *fn,int nb,int scaleup,double *ret)
/*
   reading PPM(nb=3)/PGM(nb=1) file
   allocated rgb** / unsigned char** returned
   scaleup works with ppm only
*/
{
  FILE *f;
  int i,j;
  int nx,ny;
  rgb **r;
  unsigned char **b;
  int depth,k;
  char line[256];

  f=fopen(fn,"rb");
  if (f==NULL) ERROR(("cannot open %s",fn))

  fgets(line,256,f);
  if (line[1]==5 && nb!=1) ERROR(("%s: bad file type",fn))
  if (line[1]==6 && nb!=3) ERROR(("%s: bad file type",fn))

  for (;;) {
    fgets(line,256,f);
    if (line[0]=='#') *ret=atof(line+2);
    else break; }

  sscanf(line,"%d%d",&nx,&ny);
  fgets(line,256,f);
  sscanf(line,"%d",&depth);
  if (depth!=255) ERROR(("%d bad depth",depth))

  if (nb==3) {
    xpat=nx*scaleup,ypat=ny*scaleup;
    alloc2Darray(r,xpat,ypat);
    loop (j,0,ypat) if (j%scaleup)  {
      loop (i,0,xpat) loop (k,0,3) r[i][j][k]=r[i][j-(j%scaleup)][k]; }
      else {
        loop (i,0,xpat) if (i%scaleup)  {
          loop (k,0,3) r[i][j][k]=r[i-(i%scaleup)][j][k]; }
          else { loop (k,0,3) r[i][j][k]=getc(f); } }

    fclose(f);
    return r;
  }
  else {
    alloc2Darray(b,ny,nx);
    xzbuf=nx,yzbuf=ny;
    loop (j,0,ny)
      loop (i,0,nx) b[j][i]=getc(f);
    fclose(f);
    return b;
  }
}

int res=96;
int *x,*y; /* array[-xpat..xpaper-1] */

void row(int aty, byte *z) /* ========================================= row */
/*
   basic stereo algorithm for one row of pixels
*/
{
  int i,j,k,xshift,yshift,zlast,jx;

  xshift=0;
  if (phase) loop (i,0,xpaper) if (z[i]>xshift) xshift=z[i];
  xshift=xpat*3-xshift; /* shift,phase at the max hill onset */
  yshift=aty % ypat;
  loop (i,-xpat,0) {
    x[i]=(i+xshift) % xpat;
    y[i]=yshift; }

  zlast=0;
  i=0;
  do {
    if (z[i]<zlast) {
      loop (k,i,xpaper)
        if (z[k]+k-i>=zlast) goto kfound;
      k=xpaper;
    kfound: k--;
      if (up) {
        jx=max(0,2*i-k);
        for (j=i-1; j>=jx; j--) if (z[j]<zlast-(i-j)) { k=i; break; } }
      if (k-i>rightedge) {
        /* now the pattern in [i,k] will be seen by the right eye only
           => we generate it by shifting the normal pattern in y-dir */
        xshift=xpat-z[i];
        yshift=(aty + ypat/6 + 2*ypat*i/(3*xpaper)) % ypat;
        loopto (j,i,k) {
          if (rnd()<rightedgep) x[j]=irnd(xpat);
        else x[j]=x[j-xshift];
          y[j]=yshift; }
        i=k;
        goto edged; }
    }

    k=i-xpat+z[i];
    if (k<-xpat || k>=xpaper) DISASTER(("x subrange"))
    x[i]=x[k];
    y[i]=y[k];

  edged:

    zlast=z[i];
    inc(i);
  } while (i<xpaper);

  loop (i,0,xpaper) loop (k,0,3) sheet[i][aty][k]=pat[x[i]][y[i]][k];
}

void printPPM(char *name) /* ========================================== printPPM */
{
  FILE *f=fopen(name,"wb");
  int i,j,k,ii;
  int (*row)[3],xxx,yyy;

  fprintf(f,"P6\n%d %d\n255\n",xxx=xpaper/shrink,yyy=ypaper/shrink);
  alloczero(row,xxx*sizeof(row[0]));

  loop (j,0,ypaper) {
    if (j%shrink==0) memset(row,0,xxx*sizeof(row[0]));
    loop (i,0,xpaper) loop (k,0,3)
      row[i/shrink][k]+=sheet[i][j][k];
  if (j%shrink==shrink-1) loop (ii,0,xxx) {
    unsigned char RGB[3];
    loop (k,0,3) RGB[k]=row[ii][k]/Sqr(shrink);
    fwrite(&RGB,1,3,f); }
  }
  fclose(f);
  prt("%s written",name);
}


int main(int narg,char **arg) /* ===================================== main */
{
  char *name=NULL,*patname=NULL,*env;
  byte **zbuf; /* [yzbuf][xzbuf]: note x:y order! */
  float scale=1,eye;
  int i,j,jj,iin;
  byte *z;
  long rndseed=0;
  int iarg=0;
  char *viewer=NULL;
  double zbufscaling=0,dummy=0,ratio=6;

  env=getenv("RES");
  if (env) res=atoi(env);

  if (narg<2) {
    fprintf(stderr,"\
Stereogram V 2.0, part of MACSIMUS. (c) J. Kolafa 1994-2017. Call by:\n\
  stereo OPTIONS NAME [PATNAME.ppm] OPTIONS\n\
FILES:\n\
  NAME.zbuf = input z-buffer, Portable GrayMap (P5) format, depth=255\n\
    2nd line should be `# SCALING', otherwise SCALING=1\n\
    (dark levels are multiplied by SCALING to get z-buffer in pixel units)\n\
  PATNAME.ppm = optional input pattern, Portable PixMap (P6) format, depth=255\n\
    The x-size should be XPAT=D*RES*SHRINK (for -aSHRINK) or XPAT=D*RES\n\
    (for -a-SHRINK), where D is the eye distance at background in inches (\").\n\
    Recommended D: 1.4\" < D < 2.5\".\n\
    The pattern must be periodic (not in y if it is tall enough).\n\
  NAME.ppm = output stereogram, Portable PixMap (P6) format, depth=255\n\
ENVIRONMENT:\n\
  RES = screen (or print) resolution, in DPI (integer) [96 if not given]\n\
GENERAL OPTIONS:\n\
  -rRES = (int) screen (or print) resolution, in DPI [current value=%d]\n\
     Overrides the environment variable RES\n\
     It is important to have the correct RES set!\n\
  -aSHRINK = (int) antialiasing, usually not necessary for print [default=2]\n\
     The final image is scaled 1/SHRINK-times.  If PATNAME.ppm is used,\n\
     it must be SHRINK-times larger than the final pattern.\n\
  -a-SHRINK = as above except that if PATNAME.ppm is used, it is SHRINK-times\n\
     enlarged so that the final pattern size equals the original one.\n\
     Not recommended because of worse quality than high-resolution pattern.\n\
  -eRATIO = (float) ratio (distance from screen):(eyes distance) [default=6]\n\
  -sSCALE = (int) NAME.ppm will be SCALE-times the size of FNAME.zbuf [df.=1]\n\
  -zRNDSEED = (int) random number seed [default=0=use time]\n\
  -vVIEWER = start viewer on the final image\n\
RANDOM PATTERN OPTIONS (only if PATNAME.ppm is not given):\n\
  -xXPAT = (int) x pattern size [default=384] (-> D=2\" with other defaults)\n\
  -yYPAT = (int) y pattern size [default=600]\n\
  -cCONTRAST = (float) contrast for random colors, in [0,1] [default=0.5]\n\
  -gGAMMA = (float) random pattern brightness, useful range = [0.5,3] [df.=1.4]\n\
  -mMC = (int) number of MC steps for growing patches [default=20]\n\
SEE ALSO: show (can export ZBUF of molecules)\n\
",res);
    exit(0);
  }

  loop (iarg,1,narg) if (arg[iarg][0]=='-') {
    double d=atof(arg[iarg]+2);
    switch (arg[iarg][1]) {
      case 's': scale=d; break;
      case 'a': shrink=d; break;
      case 'e': ratio=d; break;
      case 'm': mc=d; break;
      case 'c': contrast=d; break;
      case 'g': gm=d; break;
      case 'x': xpat=d; break;
      case 'y': ypat=d; break;
      case 'r': res=d; break;
      case 'v': viewer=arg[iarg]+2; break;
      case 'z': rndseed=d; break; } }
  else {
    if (name) patname=arg[iarg];
    else name=arg[iarg]; }

  initscroll(0);
  rndinit(0,rndseed);
  fprintf(stderr,"%s\n",name);

  if (patname) {
    if (shrink<0)
      pat=readpnm(patname,3,-shrink,&dummy);
    else
      pat=readpnm(patname,3,1,&dummy); }
  else {
    alloc(pat,sizeof(pat[0])*xpat);
    loop (i,0,xpat)
    alloczero(pat[i],sizeof(pat[0][0])*ypat);

    rndpat(); }

  shrink=abs(shrink);
  scale*=shrink;

  /* === z-buffer === */
  zbuf=readpnm(string("%s.zbuf",name),1,1,&zbufscaling);
  if (zbufscaling)
    prt("z-scaling = %g read from %s.zbuf",zbufscaling,name);
  else {
    prt("WARNING %s.zbuf does not contain z-scaling, 1 assumed\n",name);
    zbufscaling=1; }

  zbufscaling/=ratio;

  xpaper=xzbuf*scale+xpat+0.999;
  ypaper=yzbuf*scale;

  prt("screen/printer resolution = %d DPI",res);
  eye=(double)xpat/(res*shrink);
  prt("eye distance at background = %.2f\" = %.2f cm%s",eye,eye*2.54,eye<1.4 || eye>2.5?", WARNING: not in [1.4\",2.5\"]":"");

  prt("z-buffer will be pre-scaled %g times",scale);
  prt("final zbufscaling = %g",zbufscaling);
  prt("working image = %d x %d, the final will be scaled 1/%d-times",
      xpaper,ypaper,shrink);
  prt("picture size = %.2f cm x %.2f \" = %.2f cm x %.2f cm",
      (double)xpaper/(res*shrink),(double)ypaper/(res*shrink),
      (double)xpaper/(res*shrink)*2.54,(double)ypaper*2.54/(res*shrink)*2.54);

  alloc2Darray(sheet,xpaper,ypaper);

  allocarray(z,xpaper);
  allocarray(x,xpaper+xpat); x+=xpat;
  allocarray(y,xpaper+xpat); y+=xpat;

  /* generating the stereogram by rows of pixels */
  loop (j,0,ypaper) {
    jj=j / scale;
    if (jj>=yzbuf) ERROR(("zbuf y range"))
    loop (i,0,xpat) z[i]=0;
    loop (i,xpat,xpaper) {
      iin=(i-xpat) / scale;
      if (iin<xzbuf) z[i]=zbuf[jj][iin]*zbufscaling;
      else z[i]=0; }
    row(j,z);
  }

  printPPM(string("%s.ppm",name));

  if (viewer) system(string("%s %s.ppm &",viewer,name));

  return 0;
}
