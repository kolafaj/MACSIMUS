/* cc -O2 -o ppm2ps ppm2ps.c -lm
Converts pnm (P1-P6) to PostScript / (c) J.Kolafa 1997- */

#include "../gen/include.h"
#include <time.h>
#include "pnm.c"
#include "ddither.c"

#define PT (72/2.54) /* points per cm */

/* A4 paper as the default */
double papw=21*PT;
double paph=29.7*PT;
#define TARGETDPI 600.0 /* for -r0; should be double */

FILE *out;

int xsize,ysize,depth;
int model=0,inv=0;
char *modelname[5]={"auto","B/W","gray","RGB","CMYK"};
char *infn=NULL,*outfn=NULL;
double aspect=1;
double res=-75;
enum { AUTO,PORTRAIT,LANDSCAPE,EPS } pstype=AUTO;
int maxgray=0;
int maxrgb[3]={0,0,0};
double margin=1;

static int nrle=0,nchar=80;
static char s[12],s0[12];

void putrle(void) /************************************************** putrle */
{
  if (strcmp(s,s0) || nrle==255) {
    if (nrle) {
      if (nchar>=72) {
	fprintf(out,"\n");
	nchar=0; }
      nchar+=fprintf(out,"%02x%s",nrle,s0); }
    nrle=0;
    strcpy(s0,s); }
  nrle++;
}

int main(int narg,char **arg)
{
  int i,j,byte=0,bwthreshold=127;
  int iarg;
  char *paper="A4";
  char hdr[1024];
  double xcm=0,ycm=0,xpt,ypt;
  static double atx=0,aty=0,xx,yy;
  static double Xcm=0,Ycm=0,compress=0,iscomment=1;
  int *parm;
  int setres;
  time_t t;
  struct page_s {
    int i; /* this page, see -P */
    int n; /* number of pages, see -P */
  } page={1,1};

  out=stdout;
  nsq=0;

  if (narg<2) {
    fprintf(stderr,"\
PBM, PGM, PPM (P1-P6) to PS converter, (c) J.Kolafa 1997-2018\n\
  ppm2ps [OPTIONS] [INFILE [OUTFILE]] [OPTIONS]\n\
OPTIONS:\n\
  -#   do not copy PPM/PGM/PBM comment to PS/EPS\n\
  -0   autoselect output format [default]\n\
  -1[-#] output B/W image [default for P1,P4], #=threshold [default=127]\n\
  -1#  output dithered B/W image, #=dither square (recommended 2, useful 1..4)\n\
  -2   output gray-scale image [default for P2,P5]\n\
  -3   output RGB image [default for P3,P6]\n\
  -4   output CMYK image (some printers like this)\n\
  -A#  A# paper (#=0..10) [default=A4]\n\
  -B#  B# paper (#=0..10)\n\
  -G#  gamma correction [default=1.0=linear]\n\
  -L   Letter paper\n\
  -P#,##  this is page number # out of ## (to concatenate all such pages)\n\
  -R#[,#,#]  fix dithering by adding pseudo-random 5*5%%8 (-R32 for DeskJet)\n\
  -X#  x size in cm, overrides -r\n\
  -Y#  y size in cm, overrides -r; both -X and -Y override both -r and -a\n\
  -a#  y pixel aspect ratio (-a1.2 for VGA) [default=1]\n\
  -c   compress image (run-length encoding), use if large one-color areas\n\
  -d#[,#,#]  input depth for gray [RGB], in 0..255 [default=depth in file]\n\
  -e   output=EPS (EPSF-3.0; -x,-y,-L,-A,-B ignored)\n\
  -g   the same as -d248,248,248 (16 bit gray adjustment of white)\n\
  -i   invert colors/grayscale (processed after -d,-g,-w) [default=-i0]\n\
       -i1=-i;  -i-1=auto (invert if 1st pixel is black = darker than 248)\n\
       -i4=invert P4,P1 only;  -i-4=invert but P4,P1;  -i5,-i6 similar\n\
  -l   output=PS landscape [default=autoselect]\n\
  -m#  margins for resolution autoselect (-r-#) in cm [default=1]\n\
  -p   output=PS portrait [default=autoselect]\n\
  -r#  resolution in DPI                                    [default=-75]\n\
       #<0: minimum resolution (=small pictures small, large fit to page)\n\
       #=0: smart autoselect good for %gdpi printers\n\
  -w   the same as -d248,252,248 (16 bit TrueColor adjustment of white)\n\
  -x#  x-position in cm: [default=0=center], #>0=from left, #<0=from right\n\
  -y#  y-position in cm: [default=0=center], #>0=from top, #<0=from bottom\n\
FILES:\n\
  INFILE = Portable Pixel|Gray|B/W format, ASCII (P1-P3) or raw (P4-P6)\n\
  OUTFILE = PostScript or EPS (with -e)\n\
  missing OUTFILE = stdout, missing INFILE and OUTFILE = filter\n\
See also:\n\
  ps2pss (instead of psmerge for A4) psstring epsbg epscrop\n\
  doc2pic eps2ppm\n",TARGETDPI);
   exit(0); }

  loop (iarg,1,narg) {
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'P': {
          char *c=strchr(arg[iarg],',');
          if (!c) Error("ppm2ps: option -P: comma expected");
          page.i=atoi(arg[iarg]+2);
          page.n=atoi(c+1); }
          break;
        case 'L': 
          paper="Letter"; paph=11*72; papw=8.5*72; break;
        case 'B': xpt=atof(arg[iarg]+2)/2-0.25; goto A;
        case 'A':
          xpt=atof(arg[iarg]+2)/2;
        A:
          paper=arg[iarg]+1;
          paph=(int)(1000/pow(2,xpt-0.25))*(PT/10);
          papw=(int)(1000/pow(2,xpt+0.25))*(PT/10);
          break;
        case 'a': aspect=atof(arg[iarg]+2); break;
        case 'r': res=atof(arg[iarg]+2); break;
        case 'm': margin=atof(arg[iarg]+2); break;
        case 'x': atx=atof(arg[iarg]+2); break;
        case 'y': aty=atof(arg[iarg]+2); break;
        case 'X': Xcm=atof(arg[iarg]+2); break;
        case 'Y': Ycm=atof(arg[iarg]+2); break;
        case 'p': pstype=PORTRAIT; break;
        case 'l': pstype=LANDSCAPE; break;
        case 'i': inv=arg[iarg][2] ? atoi(arg[iarg]+2) : 1; break;
        case 'c': compress++; break;
        case 'e': pstype=EPS; break;
        case '0': model=0; break;
        case '1': model=1;
          if (arg[iarg][2]) {
            if (!atoi(arg[iarg]+2)) Error("ppm2ps: bad -1 option");
            else if (isdigit(arg[iarg][2])) makesq(atoi(arg[iarg]+2));
            else bwthreshold=atoi(arg[iarg]+3); }
          break;
        case '2': model=2; break;
        case '3': model=3; break;
        case '4': model=4; break;
        case '#': iscomment=0; break;
        case 'R': parm=dither; goto get3;
        case 'd': parm=Depth;
get3:
          {
            char *c;
            parm[0]=parm[1]=parm[2]=atoi(arg[iarg]+2);
            c=strchr(arg[iarg],',');
            if (c) {
              parm[1]=atoi(++c);
              c=strchr(c,',');
              if (c) parm[2]=atoi(++c); } }
          break;
        case 'G': Gamma=atof(arg[iarg]+2); break;
        case 'g': Depth[0]=Depth[1]=Depth[2]=248; break;
        case 'w': Depth[0]=Depth[2]=248; Depth[1]=252; break;
        case 0: break; /* - is dummy arg */
        default: Error("ppm2ps: unknown option"); }
    else if (infn)
      if (outfn) Error("ppm2ps: too many file parameters");
      else outfn=arg[iarg];
    else infn=arg[iarg]; }

  if (page.n<1 || page.n>10000)
    Error("ppm2ps: number of pages out of range, check option -P");
  if (page.i<1 || page.i>page.n)
    Error("ppm2ps: page number out of range 1..(number of pages), check option -P");

  i=openpnm(infn);
  if (model==0) model=i;

  if (outfn) out=fopen(outfn,"wt");

  fprintf(stderr,"input = %s, P%d, %d x %d\n",
          infn?infn:"(stdin)", type, xsize,ysize);
  if (comment[0]) fputs(comment,stderr);
  if (aspect<0.01 || aspect>100) Error("ppm2ps: invalid aspect");

  if (strchr("\1\2\4\5",type) && model>=3)
    fprintf(stderr,"WARNING: B/W or gray input, color output - check options\n");
  if (strchr("\1\4",type) && model==2)
    fprintf(stderr,"WARNING: B/W, gray output - check options\n");

  if (pstype==AUTO) {
    if (xsize>ysize*aspect) pstype=LANDSCAPE;
    else pstype=PORTRAIT; }
  if (pstype==PORTRAIT) xx=papw,yy=paph;
  else xx=paph,yy=papw;

  /* picture size given in cm ==> res ignored */
  setres=0;
  if (Xcm) {
    xcm=Xcm;
    setres++;
    res=xsize*2.54/xcm; }
  if (Ycm) {
    if (Xcm) {
      setres++;
      aspect=Ycm/ysize/Xcm*xsize;
      fprintf(stderr,"\npixel aspect=%g\n",aspect);
      xcm=Xcm,res=xsize*2.54/xcm;
      ycm=Ycm; }
    else {
      ycm=Ycm;
      setres++;
      res=ysize*2.54/ycm*aspect; } }

  if (!setres) {
    if (res>0) {
      /* given resolution in DPI */
      xcm=xsize*2.54/res; ycm=ysize*2.54/res*aspect; }
    else if (res<=0) {
      int autores;

      /* given minimum resolution in DPI */
      res=-res;
      if ( (autores=res==0) ) {
        res=TARGETDPI/8.; /* good dither at TARGETDPI */
        while (res<59) res*=2;
        while (res>99) res/=2; }
      xcm=xsize*2.54/res;
      ycm=ysize*2.54/res*aspect;
      if (xcm>xx/PT-2*margin || ycm>yy/PT-2*margin) {
        res=fmax(xsize*2.54/(xx/PT-2*margin),ysize*2.54/(yy/PT-2*margin));
        if (autores) {
 	 /* certain "intelligent" res good for LaserJet */
 	 if (res<TARGETDPI && res>75) res=TARGETDPI/((int)(TARGETDPI/res)); } }
      xcm=xsize*2.54/res;
      ycm=ysize*2.54/res*aspect; } }

  xpt=xsize*72/res; ypt=ysize*72/res*aspect;

  atx*=PT; aty*=-PT;
  if (atx==0) atx=(xx-xpt)/2; else if (atx<0) atx=xx-xpt+atx;
  if (aty==0) aty=(yy-ypt)/2; else if (aty<0) aty=yy-ypt+aty;

  fprintf(stderr,"output = %s, %s\n\
%.2fcm x %.2fcm = [width=%gin] [height=%gin], %g dpi\n",
          outfn?outfn:"(stdout)",
          modelname[model],
          xcm,ycm,
          xcm/2.54,ycm/2.54,
          res);

      time(&t);

  sprintf(hdr,"\
%%%%Creator: MACSIMUS/ppm2ps\n\
%%%%Title: %s\n\
%%%%CreationDate: %s\
", infn, ctime(&t));

  switch (pstype) {

    case PORTRAIT:
      if (page.i==1)
        fprintf(out,"\
%%!PS-Adobe-2.0\n\
%%%%Orientation: Portrait\n\
%%%%PaperSize: %s\n\
%%%%Pages: %d\n\
%s\
%%%%EndComments\n",
                paper,
                page.n,
                hdr);

      fprintf(out,"\
%%%%Page: %d %d\n\
%%%%BeginPageSetup\n\
%%%%PageOrientation: Portrait\n\
%%%%PageBoundingBox: 0 0 %d %d\n\
%%%%EndPageSetup\n\
gsave\n\
%.2f %.2f translate\n",
              page.i,page.n,
              (int)(papw+.999),(int)(paph+.999),
              atx,aty);

      break;

    case LANDSCAPE:
      if (page.i==1)
        fprintf(out,"\
%%!PS-Adobe-2.0\n\
%%%%Orientation: Landscape\n\
%%%%PaperSize: %s\n\
%%%%Pages: %d\n\
%s\
%%%%EndComments\n",
                paper,
                page.n,
                hdr);

      fprintf(out,"\
%%%%Page: %d %d\n\
%%%%BeginPageSetup\n\
%%%%PageOrientation: Landscape\n\
%%%%PageBoundingBox: 0 0 %d %d\n\
%%%%EndPageSetup\n\
gsave\n\
90 rotate %.2f %.2f translate\n",
              page.i,page.n,
              (int)(papw+.999),(int)(paph+.999),
              atx,aty-papw);
      break;

    case EPS: {
      fprintf(out,"\
%%!PS-Adobe-3.0 EPSF-3.0\n\
%s\
%%%%BoundingBox: 0 0 %.0f %.0f\n\
%%%%ExactBoundingBox: 0 0 %f %f\n\
%%%%Pages: 1\n\
%%%%LanguageLevel: 1\n\
gsave\n",
              hdr,
              xsize/res*72,ysize/res*72*aspect,
              xsize/res*72,ysize/res*72*aspect);

      if (model>2) fprintf(out,"%%%%Extensions: %s\n",modelname[model]); }
  }

  fprintf(out,"%% image type=%s\n",modelname[model]);
  if (iscomment && comment[0]) {
    char *c;
    for (c=comment; *c; c++) if (*c<' ') *c=' ';
    fprintf(out,"%% %s\n",comment); }

  if (compress) {
    fprintf(out,"\
%% run-length encoding\n\
72 %.3f div dup %.4f mul scale\n\
/N 0 def\n\
/C 1 string def\n\
/S %d string def\n\
%d %d %d [ 1 0 0 -1 0 %d ]\n\
{ N 0 eq {\n\
  /N currentfile C readhexstring pop 0 get def\n\
  currentfile S readhexstring pop pop }\n\
if /N N 1 sub def S } ",
            res,aspect,
            model<3?1:model,
            xsize,ysize,model==1?1:8,ysize);

    if (model<3) fprintf(out,"image");
    else fprintf(out,"false %d colorimage",model);

    loop (j,0,ysize) {
      readpnm(inv,j);
      loop (i,0,xsize) {

        switch (model) {
          case 1:
            /* inefficient repetition of the dither code in both compress/noncompress versions */
            byte<<=1;
            if (nsq) {
              if (pix2pix(i,j,gray[i])) byte|=1; }
            else {
              if (gray[i]>bwthreshold) byte|=1; }
            if ((i&7)==7) {
              sprintf(s,"%02x",byte);
              putrle();
              byte=0; }
            break;
          case 2:
            sprintf(s,"%02x",gray[i]);
            putrle();
            Max(maxgray,gray[i]);
            break;
          case 3:
            sprintf(s,"%02x%02x%02x",rgb[i][0],rgb[i][1],rgb[i][2]);
            goto getmaxrgbc;
          case 4: {
            int c,m,y,k;

            k=0;
            c=255-rgb[i][0];
            m=255-rgb[i][1];
            y=255-rgb[i][2];
            k=c; if (m<k) k=m; if (y<k) k=y;

            sprintf(s,"%02x%02x%02x%02x",c-k,m-k,y-k,k);
            getmaxrgbc:
            putrle();
            loop (k,0,3) Max(maxrgb[k],rgb[i][k])
              break; }
        } }
      if (model==1 && xsize&7) {
        byte<<=8-(xsize&7);
        sprintf(s,"%02x",byte);
        putrle();
        byte=0; } }
    s[0]=0;
    putrle(); }

  /* not compressed */
  else {
    fprintf(out,"\
72 %.3f div dup %.4f mul scale\n\
/S %d string def\n\
%d %d %d [ 1 0 0 -1 0 %d ]\n\
{ currentfile S readhexstring pop } ",
            res,aspect,
            model>2?xsize*model:model==2?xsize:(xsize+7)/8,
            xsize,ysize,model==1?1:8,ysize);

    if (model<3) fprintf(out,"image");
    else fprintf(out,"false %d colorimage",model);

    loop (j,0,ysize) {
      readpnm(inv,j);
      loop (i,0,xsize) {

        switch (model) {
          case 1:
            if (i%256==0) putc('\n',out);
            byte<<=1;
            if (nsq) {
              if (pix2pix(i,j,gray[i])) byte|=1; }
            else {
              if (gray[i]>bwthreshold) byte|=1; }
            if ((i&7)==7) {
              fprintf(out,"%02x",byte);
              byte=0; }
            break;
          case 2:
            if (i%32==0) putc('\n',out);
            fprintf(out,"%02x",gray[i]);
            Max(maxgray,gray[i]);
            break;
          case 3:
            putc(" \n"[i%10==0],out);
            fprintf(out,"%02x%02x%02x",rgb[i][0],rgb[i][1],rgb[i][2]);
            goto getmaxrgb;
          case 4: {
            int c,m,y,k;

            k=0;
            c=255-rgb[i][0];
            m=255-rgb[i][1];
            y=255-rgb[i][2];
            k=c; if (m<k) k=m; if (y<k) k=y;

            putc(" \n"[i%8==0],out);
            fprintf(out,"%02x%02x%02x%02x",c-k,m-k,y-k,k);
            getmaxrgb:
            loop (k,0,3) Max(maxrgb[k],rgb[i][k])
              break; }
        } }
      if (model==1 && xsize&7) {
        byte<<=8-(xsize&7);
        fprintf(out,"%02x",byte);
        byte=0; }
    } }

  fputs("\n\
grestore\n\
showpage\n",out);
  if (pstype==EPS) fprintf(out,"%%%%EOF\n");

  fclose(out);
  fclose(pnm);

  if (model==2) {
    if (maxgray!=255)
      fprintf(stderr,"\a\nWARNING max gray level (white)=%d, consider option -g\n\n",
              maxgray); }
  else if (model>2) {
    if (maxrgb[0]!=255 && maxrgb[1]!=255 && maxrgb[2]!=255)
      fprintf(stderr,"\a\nWARNING ");
    if (maxrgb[0]!=255 || maxrgb[1]!=255 || maxrgb[2]!=255)
      fprintf(stderr,"max RGB level (white)=%d,%d,%d, consider option -w\n\n",
              maxrgb[0],maxrgb[1],maxrgb[2]); }

  return 0;
}
