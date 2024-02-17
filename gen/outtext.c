static int textsize=1;

void settextstyle(int a,int b,int size)
{
textsize=size;
}

static int xjust,yjust=2;

void settextjustify(int x,int y)
{
 xjust=x; yjust=y;
}

#include "fbgi.c"
#include "fc7x13.c"
#include "fc9x15.c"

unsigned short *font_data=font_fbgi;

/* ??? DOS ??? */
int font_width=8;  /* includes spacing (after character) */
int font_height=9; /* pixel range incl. accents */
int font_top=0;    /* line (topmost line incl. accents =0) of top of digits */
int font_base=7;   /* baseline */
/* see also default: in selectfont() */

void selectfont(int font)
{
 switch (font) {
   case 3: font_data=font_fc9x15; font_width=9; font_height=15; 
           font_top=1; font_base=11;
           break;
   case 2: font_data=font_fc7x13; font_width=7; font_height=13; 
           font_top=1; font_base=10;
           break;
   default: font_data=font_fbgi; font_width=8; font_height=9; 
            font_top=0; font_base=7;
            break;
   }
}

int outtextxy(int x,int y,char *s)
{
 int X,Y,B,l=0;

 x-=xjust*(strlen(s)*font_width/2-1)*textsize;
 switch (yjust) {
   case 0: y-=font_base*textsize; break;
   case 1: y-=(font_base+font_top)/2*textsize; break;
   case 2: y-=font_top*textsize; }

 g_lock();
 while (*s) {
   loop (X,0,font_width) {
     B=font_data[(*(unsigned char*)s)*font_width+X];

     if (textsize==1) loop (Y,0,font_height) {
       if (B&1) putpixel(x+X,y+Y,getcolor());
       B>>=1; }
     else loop (Y,0,font_height) {
       if (B&1) {
	 int ii,jj;
	 loop (ii,X*textsize,X*textsize+textsize)
	   loop (jj,Y*textsize,Y*textsize+textsize)
	   putpixel(x+ii,y+jj,getcolor()); }
       B>>=1; } }

   l+=font_width*textsize;
   x+=font_width*textsize; s++; }
 g_unlock();

 return l;
}

int outtext(char *s)
{
 return outtextxy(ggx,ggy,s);
}
