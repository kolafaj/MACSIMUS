/***********************************************************************
 * $Author: markv $
 * $Revision: 1.2 $
 * $Date: 88/10/04 14:30:44 $
 * $Log:	pic.c,v $
 * Revision 1.2  88/10/04  14:30:44  markv
 * Changed pixel writing primitives to write individual pixels rather
 * than scanlines.  Simplifies certain loops inside Screen.
 * 
 * Revision 1.1  88/09/11  11:00:41  markv
 * Initial revision
 * 
 ***********************************************************************/

#include <stdio.h>
#include "defs.h"
#include "pic.h"
#include "extern.h"

/*======================================================================*/
/* PIC.C								*/
/*                                                                      */
/* Simple routines for outputting a pixel map....                       */
/*                                                                      */
/* Mark VandeWettering, markv@cs.uoregon.edu                            */
/*======================================================================*/

Pic *
PicOpen(filename,x,y)
 char *filename ;
{
	Pic	*tmp ;
	int i;

	tmp = (Pic *) malloc(sizeof (Pic)) ;
	tmp -> filename = (char *) malloc (strlen(filename)+1) ;
	strcpy(tmp->filename, filename);

        /* JK: "w" changed to "wb" */
        if (((tmp -> filep)=fopen(filename, "wb"))==NULL) {
		perror( filename );
		exit( 1 );
	}

	tmp -> x = x ; tmp -> y = y ;

	fprintf(tmp -> filep, "P6\n#");
	for (i=0; i<mainargc; i++)  fprintf(tmp -> filep," %s",mainargv[i]);
	fprintf(tmp -> filep, "\n%d %d\n255\n", x, y);

	return(tmp);
}

#if 0
int
PicWritePixel(pic, color)
 Pic *pic ;
 Color color ;
{
	fputc((unsigned char) (255.0 * color[0]), pic -> filep) ;
	fputc((unsigned char) (255.0 * color[1]), pic -> filep) ;
 return fputc((unsigned char) (255.0 * color[2]), pic -> filep) ;
}
#else
double pow(double,double);
int
PicWritePixel(pic, color)
 Pic *pic ;
 Color color ;
{
  int i,e;
  for (i=0; i<3; i++) e=fputc((unsigned char) (255.0 * pow((double)color[i],(double)Gamma)+0.5), pic -> filep) ;
  return e;
}
#endif

int
PicClose(pic)
 Pic *pic ;
{
        return fclose(pic -> filep) ;
}


/* added by JK: */
pic_t *bg; /* background picture */

pic_t *ReadPic(char *fn)
{
pic_t *pic;
FILE *in;
char line[256];
int depth;
int i;
Vec avcol;

avcol[0]=avcol[1]=avcol[2]=0;

if (!fn) return NULL;

in=fopen(fn,"rb");
if (!in) {
  fprintf(stderr,
    "no background file %s, using uniform color %f %f %f\n",
    fn,
    BackgroundColor[0],BackgroundColor[1],BackgroundColor[2]);
  return NULL; }

pic=(pic_t*)malloc(sizeof(pic_t));
fgets(line,256,in);
do fgets(line,256,in); while (line[0]=='#');
sscanf(line,"%d%d",&pic->x,&pic->y);
fgets(line,256,in);
sscanf(line,"%d",&depth);

fprintf(stderr,"background file %s: %dx%d, depth=%d\n",fn,pic->x,pic->y,depth);
if (depth!=255) perror("bad depth");

pic->c=(void*)malloc(pic->x*pic->y*sizeof(pic->c[0]));
if (!pic->c) perror("ReadPPM: no heap");

for (i=0; i<pic->x*pic->y; i++) {
  avcol[0]+=(pic->c[i][0]=getc(in));
  avcol[1]+=(pic->c[i][1]=getc(in));
  avcol[2]+=(pic->c[i][2]=depth=getc(in));
  if (depth<0) perror("ReadPPM: unexpected eof"); }

avcol[0]*=1./255/(pic->x*pic->y);
avcol[1]*=1./255/(pic->x*pic->y);
avcol[2]*=1./255/(pic->x*pic->y);

if (abs(bgmode)==2) {
  VecCopy(avcol,BackgroundColor);
  printf("used "); }
VecPrint("averaged bg color =",avcol);

return pic;
}

void
GetPic(pic_t *pic,Flt x,Flt y,int mode,Vec c)
/*
 * returns color c
 * x,y in <-0.5,0.5> correspond to the picture
 * if x,y are out picture:
 *   mode=0: BackGround returned
 *   mode=1: the picture is tiled (periodic b.c.)
 *   mode=-1: the picture is anti-tiled (all odd are flipped)
 */
{
#define sign(X) (X<0?-1:1)
struct i_s { int x,y; } i;
int indx;
extern double rndcos();

x+=bg->xshift;
y+=bg->yshift;

#define norm(X) \
  i.X=(int)(X*pic->X+0.5*sign(X))+pic->X/2;               \
  if (i.X<0 || i.X>=pic->X) {                             \
    if (mode==0) { VecCopy(BackgroundColor, c); return; } \
    while (i.X<0) i.X+=64*pic->X;                         \
    if (mode>0) i.X=i.X%pic->X;                           \
    else {                                                \
      i.X=i.X%(2*pic->X);                                 \
      if (i.X>=pic->X) i.X=2*pic->X-1-i.X; } }

/*.....printf("%f %f %d %d\n",x,y,i.x,i.y);*/

x= -x;
norm(x)
norm(y)

indx=pic->x*i.y+i.x;
c[0]=pic->c[indx][0]/255.;
c[1]=pic->c[indx][1]/255.;
c[2]=pic->c[indx][2]/255.;
/*.....if (i.x%64<2 || i.y%64<2) c[0]=c[1]=c[2]=0;*/
}
