/* cc -O2 -o ppm2ppm ppm2ppm.c
   cuts out pic (ppm) file
*/
#define _ANSI_C_SOURCE
#include "../gen/include.h"

FILE *in,*out;
char line[256],*a;

int xsize,ysize,depth;
#define y0 y_0
#define y1 y_1
int x0,x1,y0,y1;
int xstride=1,ystride=1;
int dith=172;

int intype,outtype,inv;
unsigned char rgb[3];
unsigned char inmask=0;
static unsigned char b,outmask=128;
int rgb2g[4]={30,59,11,100}; /* conversion RGB->GRAY, [3]=sum */

int mygetc(FILE *in) /*********************************************** mygetc */
{
  int i;

  if (fscanf(in,"%d",&i)==0) i=-1;

  return i;
}

int getrgb(FILE *in) /*********************************************** getrgb */
{
  int k,g;

  switch (intype) {

    case '1': {
      static unsigned char b;

      if (!inmask) {
        inmask=128;
        g=mygetc(in);
        if (g<0) return -1;
        b=g; }
      if (b&inmask) rgb[0]=255;
      else rgb[0]=0;
      if (!inv) rgb[0]=~rgb[0];
      rgb[1]=rgb[2]=rgb[0];
      inmask >>= 1;
      break; }

    case '4': {
      static unsigned char b;

      if (!inmask) {
        inmask=128;
        g=getc(in);
        if (g<0) return -1;
        b=g; }
      if (b&inmask) rgb[0]=255;
      else rgb[0]=0;
      if (!inv) rgb[0]=~rgb[0];
      rgb[1]=rgb[2]=rgb[0];
      inmask >>= 1;
      break; }

    case '2':
      g=mygetc(in); if (g<0) return -1;
      g=g*255/depth;
      if (inv) g=~g;
      rgb[0]=rgb[1]=rgb[2]=g;
      break;

    case '5':
      g=getc(in); if (g<0) return -1;
      if (depth<=255) g=g*255/depth;
      else g=(g+256*getc(in))*255/depth;
      if (inv) g=~g;
      rgb[0]=rgb[1]=rgb[2]=g;
      break;

    case '3':
      loop (k,0,3) {
	g=mygetc(in); if (g<0) return -1;
	g=g*255/depth;
	if (inv) g=~g;
	rgb[k]=g; }
      break;

    case '6':
      loop (k,0,3) {
	g=getc(in); if (g<0) return -1;
	if (depth<=255) g=g*255/depth;
	else g=(g+256*getc(in))*255/depth;
	if (inv) g=~g;
	rgb[k]=g; }
      break;

    default:
      return -2; }

  return 0;
}

int i,j;

int rev(int c) /******************************************************** rev */
/* reverse order of bits */
{
  return (c&128)>>7 | (c&64)>>5 | (c&32)>>3 | (c&16)>>1
         | (c&8)<<1 | (c&4)<<3 | (c&2)<<5 | (c&1)<<7;
}

char *sepP0;

void flushbyte(void) /******************************************** flushbyte */
{
  outmask=128;
  if (outtype=='4') fputc(b,out);
  else if (outtype=='1') {
    int x;
    loop (x,0,8) {
      fprintf(out,b&128?" 1":" 0");
      b<<=1; } }
  else fprintf(out,"%s%d",sepP0,rev(b));
  b=0;
  sepP0=",";
}

void putrgb(FILE *out) /********************************************* putrgb */
{
  int k;
  int gray=(rgb2g[0]*rgb[0]+rgb2g[1]*rgb[1]+rgb2g[2]*rgb[2])/rgb2g[3];

  switch (outtype) {
    case '0':
    case '1':
    case '4': {
      static int dither[15][17],init;
      int *ditherrow=&dither[0][0]; /* standard order assumed */

      if (!init) {
        /* Simple... */
        int k;

        if (dith==0) dith=172;
#if 0
        /* the paranoic compiler shouts loud */
        if (dith>0) loop (k,0,255) dither[0][k]=(k*dith)%255;
        else loop (k,0,255) dither[0][k]=-dith;
#else
        /* equivalent cleaned code */
        if (dith>0) loop (k,0,255) ditherrow[k]=(k*dith)%255;
        else loop (k,0,255) ditherrow[k]=-dith;
#endif

        init++; }

      gray += dither[i%15][j%17];

      if (gray<255) b |= outmask;
      outmask >>= 1;
      if (!outmask) flushbyte();
      break; }

    case '5':
      fputc(gray,out);
      break;

    case '6':
      loop (k,0,3) fputc(rgb[k],out);
      break;

    case '3':
      fprintf(out,"%3d %3d %3d\n",rgb[0],rgb[1],rgb[1]);
      break;

    case '2':
      fprintf(out," %3d",gray);
      break;

    default:
      fprintf(stderr,"P%c:",outtype);
      Error("unsupported output type"); }
}

int main(int narg,char **arg) /**************************************** main */
{
  int comm=0;
  int justwritten=0;
  int R=-1,G=-1,B=-1;

  if (narg<3) {
   fprintf(stderr,"\
Extract rectangle and/or convert pnm files. Call by:\n\
  ppm2ppm {INFILE|-} {OUTFILE|-}\n\
    [-xFROMX] [-yFROMY] [{-XTOX|-wWIDTH}] [{-YTOY|-hHEIGHT}] [-fCUTFRAME]\n\
    [-sXSTRIDE] [-tYSTRIDE]\n\
    [-POUTPUTTYPE] [-dDITHER|-d-THRESHOLD]\n\
    [-RRED] [-GGREEN] [-BBLUE] [-WWHITE]\n\
    [-rSETRED -gSETGREEN -bSETBLUE] [-i] [-#] -POUTTYPE\n\
ARGUMENTS:\n\
  OUTFILE is the same format as INFILE unless -P\n\
  supported input formats: P1-P6\n\
  supported output formats: P1-P6, P0=X bitmap in C\n\
    -Pa/-Pb convert to asc/bin(raw)\n\
  TOX and TOY are not included (C-style)\n\
  negative FROMX or TOX is subtracted from the X-range (-X or width), so is Y\n\
  more -x -X -y -Y -w -h are combined\n\
  -f removes given frame from all sides\n\
  -P4 only: -dDITHER = mult.factor for 15x17 simple dither rectangle [df.=172]\n\
            -d-THRESHOLD = no dithering; lower THRESHOLD = darker\n\
  RED GREEN BLUE WHITE: RGB -> GRAY conversion, default=30,59,11,100\n\
  SETRED SETGREEN SETBLUE: set RGB\n\
  -i invert\n\
  -# preserve comments (default=comments stripped off)\n\
BUG: Output P2,3,5,6 converted to depth=255, information loss may occurr\n\
See also:\n\
  ppmframe ppmcenter\n\
  ppmscale # use for for smooth (halftone) scaling\n\
");
   exit(0); }

  loop (i,3,narg) if (arg[i][0]=='-') {
    if (arg[i][1]=='#') comm++;
    else if (arg[i][1]=='P') outtype=arg[i][2]; }

  if (strcmp("-",arg[1])) in=fopen(arg[1],"rb");
  else in=stdin;
  if (!in) Error(arg[1]);

  if (strcmp("-",arg[2])) out=fopen(arg[2],"wb");
  else out=stdout;
  if (!out) Error(arg[2]);

  if (!fgets(line,256,in)) Error("ppm2ppm: INFILE too short");
  intype=line[1];
  if (!outtype) outtype=intype;
  if (strchr("rRbB",outtype)) {
    outtype=intype;
    if (intype<'4') outtype=intype+3; }
  if (strchr("aA",outtype)) {
    outtype=intype;
    if (intype>='4') outtype=intype-3; }
  if (outtype>'0') fprintf(out,"P%c\n",outtype);
  if (line[0]!='P' || intype<'1' || intype>'6')
    Error("bad header: P1--P6 expected");
  do {
    if (line[0]=='#' && comm) {
      if (outtype=='0') fprintf(out,"//");
      fputs(line,out); }
    if (!fgets(line,256,in)) Error("ppm2ppm: INFILE bad format");
    } while (line[0]=='#');
  sscanf(line,"%d%d",&xsize,&ysize);
  if (strchr("2356",intype)) {
    if (!fgets(line,256,in)) Error("ppm2ppm: INFILE too short");
    sscanf(line,"%d",&depth);
    if (depth!=255) fprintf(stderr,"WARNING depth=%d: will be converted to 255\n",depth); }

  x0=0; x1=xsize;
  y0=0; y1=ysize;

  loop (i,3,narg) {
    a=arg[i];
    if (a[0]=='-') a++;
    j=atoi(a+1);
    switch (*a) {
      case 'i': inv++; break;
      case 'R':rgb2g[0]=j; break;
      case 'G':rgb2g[1]=j; break;
      case 'B':rgb2g[2]=j; break;
      case 'W':rgb2g[3]=j; break;
      case 'r': R=j; break;
      case 'g': G=j; break;
      case 'b': B=j; break;
        /*     case 'W': R=G=B=255; break; */
      case 'x': if (j<0) x0=x1+j; else x0=j; break;
      case 'y': if (j<0) y0=x1+j; else y0=j; break;
      case 'X': if (j<0) x1+=j; else x1=j; break;
      case 'Y': if (j<0) y1+=j; else y1=j; break;
      case 'f': x0=y0=j; x1-=j; y1-=j; break;
      case 'w': x1=x0+j; break;
      case 'h': y1=y0+j; break;
      case 'd': dith=j; break;
      case 's': xstride=j; break;
      case 't': ystride=j;
      case 'P':
      case '#': break;
      default: Error("ppm2ppm: bad option"); } }

#if 0
  if (x0>=x1 || x0<0 || x1>xsize) Error("ppm2ppm: bad x range");
  if (y0>=y1 || y0<0 || y1>ysize) Error("ppm2ppm: bad y range");
#else
  if (x0>=x1) Error("ppm2ppm: bad x range");
  Max(x0,0) Min(x1,xsize)
  if (y0>=y1) Error("ppm2ppm: bad y range");
  Max(y0,0) Min(y1,ysize)
#endif
  if (xstride<1 || ystride<1) Error("ppm2ppm: bad stride(s)");

  fprintf(out,
 	 outtype=='0'
 	 ?"#define Xbitmap_width %d\n#define Xbitmap_height %d\n"
 	 :"%d %d\n",
 	 i=(x1-x0+xstride-1)/xstride,
 	 j=(y1-y0+ystride-1)/ystride);

  if (strchr("2356",outtype)) fprintf(out,"%d\n",255);

  fprintf(stderr,"P%c %d x %d --> P%c %d x %d\n",intype,xsize,ysize,outtype,i,j);

  sepP0="static unsigned char Xbitmap_bits[] = {\n";

  loop (j,0,ysize) {
    justwritten=0;

    loop (i,0,xsize) {
      if (getrgb(in)) {
        Error("ppm2ppm: unexpected eof");
        goto end; }

      if (i>=x0 && i<x1 && j>=y0 && j<y1
 	 && (i-x0)%xstride==0
 	 && (j-y0)%ystride==0) {
        if (R>=0) rgb[0]=R;
        if (G>=0) rgb[1]=G;
        if (B>=0) rgb[2]=B;
        putrgb(out);
        justwritten++; } }

    if (strchr("12",outtype)) fprintf(out,"\n");

    if (justwritten && (outtype=='4' || outtype=='0') && outmask!=128)
      flushbyte();
    sepP0="\n,";

    inmask=0; }

 end:
  if (outtype=='0') fprintf(out,"};\n");

  fclose(out);
  fclose(in);

  return 0;
}
