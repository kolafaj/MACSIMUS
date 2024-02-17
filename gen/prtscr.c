/* extended eprtscr.c : see there */

static FILE *prn;

#define MAXX 1024 /* must be a multiple of 16 */

static unsigned int swap(unsigned int w)
{
char x,*c=(char*)&w;
x=c[0]; c[0]=c[1]; c[1]=x;
return w;
}

void PrtScr(char mode)
/*
  mode:
    e 1 2 9 n  9 pin epson, 120dpi in x-direction, 72dpi in y-direction
    p 7 8 t      24 pin epson
    l 3 h        HP LaserJet
    d 6          HP DeskJet Econofast B/W
(see ldraw.pas)
*/
{
unsigned int i,j,k,n,b,off,maxx=getmaxx(),maxy=getmaxy(),iw,w,
  mapcol,mapdata,pixbits,ipix;
unsigned int *dataptr;
static unsigned char dot24[MAXX][3];
unsigned char *dot=(unsigned char*)dot24,mask;
unsigned int *hp=(unsigned int*)dot24;
/* original fix-length, MAXX=640
static struct image_s {
  unsigned int xsize,ysize;
  struct l_s {
    unsigned int b8[MAXX/16],b4[MAXX/16],b2[MAXX/16],b1[MAXX/16];
    } l[8];
  } image;
*/
static struct image_s {
  unsigned int xsize,ysize;
  unsigned int data[2*MAXX];
  } image;
struct viewporttype ViewPort;
enum {PIN9,PIN24,HP,DJ,UNKNOWN} printer=UNKNOWN; /* don't change order ! */
char *initstr[4]={
  "\33@\33A\10\r", /* 9: init + spacing 8"/72 */
  "\33@\33""3\30\r\33U1", /* 24:180 dpi, unidirectional */
  "\33E\33&l0O\33&l6D\33&l66P\33&l0E\33&k0G\33&l0L\33*t150R", /* HP LJ */
  "\33E\33&l0O\33&l6D\33&l66P\33&l0E\33&k0G\33*o-1M\33*o2D\33&l0M\33&l0L\33*t150R"};
     /* DJ B/W econofast, depletion=2, normal paper */
char *prtname=getenv("PRTSCR");

if (strchr("eE129",mode)) printer=PIN9,w=1;
if (strchr("pPtT78",mode)) printer=PIN24,w=3;
if (strchr("hHlL3",mode)) printer=HP,w=1;
if (strchr("dD6",mode)) printer=DJ,w=1;

if (printer==UNKNOWN) { fputc('\a',stderr); return; }

if (maxx>MAXX-1) maxx=MAXX-1;
mapcol=(maxx+1)/16;

getviewsettings(&ViewPort);
setviewport(0,0,maxx,maxy,0);

pixbits=0;
j=imagesize(0,0,maxx,7);
if (j==10+maxx*4) pixbits=4; /* 16color mode */
if (j==7+maxx) pixbits=1;    /* 1color mode */
                             /* can easily extend: I do not need them */

mapdata=mapcol*pixbits;

if (pixbits) {

  prn=fopen(prtname && *prtname ? prtname : "LPT1","wb");
  fputs(initstr[printer],prn);

  i=0;
  do {

    if (printer==PIN24) memset(dot24,0,sizeof(dot24));

    loop (iw,0,w) {

      getimage(0,i,maxx,i+7,&image);
      if (printer==PIN9) memset(dot,0,maxx+1);
      mask=128;
      loop (j,0,8) {
        loopto (k,0,maxx/16) {
          off=16*k;
  /* orig fix-length image:
          b=swap(image.l[j].b1[k] | image.l[j].b2[k] | image.l[j].b4[k] | image.l[j].b8[k]);
  */
          dataptr=&image.data[j*mapdata+k];
          b=0;
          for (ipix=0; ipix<mapdata; ipix+=mapcol) b |= dataptr[ipix];
          if (ipix>1) b=swap(b);
/*.....          b=swap(dataptr[0] | dataptr[mapcol] | dataptr[mapcol*2] | dataptr[mapcol*3]);*/
          if (printer<HP) loop (n,0,16) {
            if ((b & 0x8000)!=0)
              if (printer==PIN9) dot[off+n]|=mask;
              else dot24[off+n][iw]|=mask;
            b<<=1; }
          else
            hp[k]=b;
          }

        if (printer>=HP) {
          for (n=maxx/16; n>0; n--) if (hp[n]!=0) goto lhp;
          n=0;
          lhp: n++;
          fprintf(prn,"\33*b%dW",n*2);
          loop (k,0,n) hp[k]=swap(hp[k]);
          fwrite(hp,n,2,prn); }

        mask>>=1; }

      i+=8; }

    if (printer<HP) {
      if (printer==PIN9) {
        for (n=maxx; n>0; n--) if (dot[n]!=0) goto l2;
        n=0;
        l2: n++; }
      else
        n=maxx; /* pessimistic */

      fputs(printer==PIN9 ? "\33L" : "\33*\'",prn);
      fwrite(&n,2,1,prn);

      if (printer==PIN9) fwrite(dot,1,n,prn); else fwrite(dot24,3,n,prn);
      fputs("\r\n",prn); }

    } while (i<=maxy);

  if (printer>=HP) fputs("\33*rB\14",prn); /* end graphics, FF */
  else fputs("\33""2\r\n",prn);  /* 1""/6 CR LF */
  fclose(prn); }
else
  fputc('\a',stderr);
setviewport(ViewPort.left,ViewPort.top,ViewPort.right,ViewPort.bottom,ViewPort.clip);
}
