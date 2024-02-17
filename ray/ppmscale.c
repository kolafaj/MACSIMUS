/* cc -O2 -o ppmscale ppmscale.c -lm
*/
#include "../gen/include.h"

double getscale(char *c)
{
  char *s=strchr(c,'/');
  double d=atof(c);

  if (!s) return d;
  else return d/atof(s+1)+1e-14;
}

double GAMMA=2.2; /* gamma value of display */
double iGAMMA;
int isgamma;

double Gamma(int x)
/* gamma function: should contain also the initial linear part */
{
  return pow(x/255.,GAMMA);
}

double InvGamma(double x)
/* inverted gamma function */
{
  return 255.5*pow(x,iGAMMA);
}

int ppmsize(char *fn,int w)
{
  FILE *f=fopen(fn,"rb");
  char line[256];
  int x=1000,y=1000;

  if (!f) Error("ppmscale -W,-H: no such file");

  fgets(line,256,f);
  if (line[0]!='P') Error("ppmscale -W,-H: bad header");

  do {
    if (!fgets(line,256,f)) Error("ppmscale -W,-H: format");
  } while (line[0]=='#');
  sscanf(line,"%d%d",&x,&y);

  fclose(f);

  return w?y:x;
}


#define DQ (1e-5) /* quantization error if exactly between levels */

int main(int narg,char **arg)
{
  FILE *in,*out;
  char line[256];

  int xxsize,yysize; /* old picture sizes */
  int depth;         /* old picture depth; probably only 255 */

  int xsize,ysize; /* new picture sizes */
  int iarg,j;

  int type; /* P4, P5, P6 */
  int outtype; /* P5, P6 */
  int inv=0; /* input invert flag */
  int Inv=0; /* output invert flag */
  int raw=0; /* raw (don't average) flag */
  int quant=0; /* quantization */
  int width=0,height=0;
  double tabGamma[256];

  typedef unsigned char rgb_t[3];
  rgb_t **rgb;

  int rowsize;
  unsigned char *row;

  int ywin,yfrom,yto; /* y-window (sliding part of file) size and contents */

  double xscale,yscale=0,q;

  if (narg<4) {
    fprintf(stderr,"Rescale image with proper smooth pixel averaging.  Call by:\n\
  ppmscale {INFILE|-} {OUTFILE|-} [XSCALE [YSCALE]] [OPTIONS]\n\
INFILE: one of .pbm (bitmap, P4), .pgm (graymap, P5), .ppm (pixmap, P6)\n\
OUTFILE: .ppm (pixmap, P6) or .pgm (graymap, P5), depth=255\n\
  -: denotes stdin and stdout as i/o files \n\
XSCALE YSCALE are real scaling factors, simple ratios allowed\n\
  default is YSCALE=0 which means YSCALE=XSCALE\n\
OPTIONS:\n\
  -xWIDTH  output width in pixels (YSCALE=XSCALE set if one of -x, -y)\n\
  -XFILE   output width is the same as of PNM file FILE\n\
  -yHEIGHT output height in pixels\n\
  -YFILE   output height is the same as of PNM file FILE\n\
  -qQUANT  forces quantization (use -q2 for B/W output)\n\
  -gGAMMA  gamma value (of output device used for viewing), default=2.2\n\
  -g       the same as -g1 (linear)\n\
  -r       raw mode (pickup the nearest pixel), default=average pixels\n\
  -i       invert b/w on input (useful for pbm)\n\
  -I       invert b/w on output\n\
Legacy: -w=-x  -W=-X -h=-y -H=-Y\n\
Example (scale a bitmap to 40%% x 40%%, inverting black/white on input)\n\
  ppmscale test.pbm small.pgm 2/5 0.4 -i\n\
See also: ppm2ppm ppmmap gammappm jkv\n");
  exit(0); }

  if (strcmp(arg[1],"-")) in=fopen(arg[1],"rb");
  else in=stdin;
  if (!in) Error(arg[1]);

  if (strcmp(arg[2],"-")) out=fopen(arg[2],"wb");
  else out=stdout;
  if (!out) Error(arg[2]);

  xscale=getscale(arg[3]);
  if (narg>4) yscale=getscale(arg[4]);

  loop (iarg,3,narg) if (arg[iarg][0]=='-') switch (arg[iarg][1]) {
    case 'x':
    case 'w': width=atoi(arg[iarg]+2); break;
    case 'X':
    case 'W': width=ppmsize(arg[iarg]+2,0); break;
    case 'y':
    case 'h': height=atoi(arg[iarg]+2); break;
    case 'Y':
    case 'H': height=ppmsize(arg[iarg]+2,1); break;
    case 'i': inv=1; break;
    case 'I': Inv=1; break;
    case 'r': raw=1; break;
    case 'q': quant=atoi(arg[iarg]+2); break;
    case 'G': case 'g': GAMMA=atof(arg[iarg]+2); break;
    default: Error("bad option or negative scaling factor"); }

  if (yscale==0) yscale=xscale;

  if (GAMMA==0) GAMMA=1;
  isgamma=GAMMA!=1;
  iGAMMA=1./GAMMA;
  if (raw && isgamma) {
    fprintf(stderr,"WARNING: GAMMA ignored in the raw mode\n");
    isgamma=0; GAMMA=1; }

  if (isgamma) {
    int x;
    loop (x,0,256) tabGamma[x]=Gamma(x); }

  if (quant) if (quant<2 || quant>255) Error("QUANT must be in [2,255]");

  fgets(line,256,in);
  if (line[0]!='P') Error("ppmscale: bad header");
  type=atoi(line+1);
  if (type<4 || type>6) Error("ppmscale: unsupported file type\n\
supported: P4 (bitmap), P5 (graymap), P6 (color pixmap)");

  outtype=type;
  if (outtype==4) outtype=5;

  do {
    if (!fgets(line,256,in)) Error("ppmscale: format in");
  } while (line[0]=='#');
  sscanf(line,"%d%d",&xxsize,&yysize);
  if (type!=4) {
    fgets(line,256,in);
    sscanf(line,"%d",&depth); }
  else depth=255;
  if (depth!=255) Error("ppmscale: sorry, only depth=255 supported");

  /* scaling to given pixels requested */
  if (width)  yscale=xscale=width/(xxsize-1e-9);
  if (height) xscale=yscale=height/(yysize-1e-9);
  if (width*height) {
    xscale=width/(xxsize-1e-9);
    yscale=height/(yysize-1e-9); }

  switch (type) {
    case 4:
      rowsize=(xxsize+7)/8;
      alloc(row,rowsize);
      break;
    case 5:
      rowsize=xxsize;
      alloc(row,rowsize);
      break;
    case 6:
      rowsize=xxsize*3; }

  q=255.0/depth*(xscale*yscale);
  if (quant) q+=DQ;

  fprintf(out,"P%d\n%d %d\n255\n",
	  outtype,xsize=xxsize*xscale,ysize=yysize*yscale);

  outtype=2*outtype-9; /* 1 or 3 */

  fprintf(stderr,"ppmscale: %d x %d depth=%d --> %d x %d depth=255\n",
	  xxsize,yysize,depth,xsize,ysize);

  ywin=2+1.000000000001/yscale;
  alloc(rgb,ywin*sizeof(rgb[0]));
  loop (j,0,ywin) alloc(rgb[j],xxsize*sizeof(rgb[0][0]));

  yfrom=yto=0;

  /* new picture y-loop */
  loop (j,0,ysize) {
    double Yfrom=j/yscale,Yto=(j+1)/yscale;
    int jfrom=Yfrom,jto=Yto;
    int i,k;

    /* sliding input picture window */
    while (jfrom>yfrom) {
      int h;
      rgb_t *aux=rgb[0];

      loop (h,1,ywin) rgb[h-1]=rgb[h];
      rgb[ywin-1]=aux;
      yfrom++; }

    /* filling input picture window from file */
    while (yto<yfrom+ywin) {
      int ii,yy;

      if (yto>=yysize) break;
      yy=yto-jfrom;
      switch (type) {
        case 4:
	  fread(row,rowsize,1,in);
	  loop (ii,0,xxsize)
	    rgb[yy][ii][0]=!(row[ii/8]&(128>>(ii&7)))*255;
	  break;
        case 5:
	  fread(row,rowsize,1,in);
	  loop (ii,0,xxsize) rgb[yy][ii][0]=row[ii];
	  break;
        case 6:
	  fread(rgb[yy],rowsize,1,in); }
      if (inv)
	loop (ii,0,xxsize) loop (k,0,outtype) rgb[yy][ii][k]=~rgb[yy][ii][k];
      yto++; }


    /* new picture x-loop */
    loop (i,0,xsize) {
      double pix[3]={0,0,0};
      double Xfrom=i/xscale,Xto=(i+1)/xscale;
      int ifrom=Xfrom,ito=Xto,jj;

      if (raw)
	/* raw mode */
	loop (k,0,outtype)
	  pix[k]=rgb[(int)((Yfrom+Yto)/2)-jfrom][(int)((Xfrom+Xto)/2)][k]
	  /(xscale*yscale);

      else loopto (jj,jfrom,jto) {
	/* calculating output pixel (x,y) as weighted sum over jj,ii */
	int ii,djj=jj-jfrom;
	double wy;

	if (Yfrom>jj)
	  if (Yto<jj+1) wy=1./yscale; else wy=jj+1-Yfrom;
	else
	  if (Yto<jj+1) wy=Yto-jj; else wy=1;

	loopto (ii,ifrom,ito) {
	  double w;

	  if (Xfrom>ii)
	    if (Xto<ii+1) w=wy/xscale; else w=wy*(ii+1-Xfrom);
	  else
	    if (Xto<ii+1) w=wy*(Xto-ii); else w=wy;

	  if (quant && jj==(jfrom+jto)/2 && ii==(ifrom+ito)/2) w+=DQ;

	  if (isgamma) loop (k,0,outtype) pix[k]+=tabGamma[rgb[djj][ii][k]]*w;
	  else loop (k,0,outtype) pix[k]+=rgb[djj][ii][k]*w; } }

      loop (k,0,outtype) {
	unsigned char outpix;

	//	if (isgamma==2) outpix=255.5*sqrt(pix[k]*q);	else
	if (isgamma) outpix=InvGamma(pix[k]*q);
	else outpix=pix[k]*q+0.5;

	if (quant) {
	  outpix=pix[k]*q/255*quant;
	  if (outpix>=quant) outpix=quant-1;
	  outpix=outpix*255/(quant-1); }

	if (Inv) outpix=~outpix;
	putc(outpix,out); }

    } }

  fclose(out);
  fclose(in);

  return 0;
}
