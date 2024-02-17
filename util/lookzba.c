#define PURIFY

#include "ground.h"
#include "alloc.h"

int main(int narg,char **arg)
{
typedef unsigned char byte;

byte **zbuf; /* [yzbuf][xzbuf]: note x:y order! */
int xzbuf=320,yzbuf=200;
FILE *zb=fopen(arg[1],"rb");
int i,j;

initscroll(0);

getdata
  get(xzbuf) get(yzbuf)
  checkdata enddata

/* === z-buffer === */
alloc(zbuf,sizeof(zbuf[0])*yzbuf);
loop (i,0,yzbuf)
  alloc(zbuf[i],sizeof(zbuf[0][0])*xzbuf);

/* note that in the z-buf file y=0 is the top line: reflex (?) */
{
int minz=255,maxz=0,j;
loop (i,0,yzbuf) {
  if (fread(zbuf[i],xzbuf,1,zb)!=1) Error("z-buf file too short");
  loop (j,0,xzbuf) {
    Min(minz,zbuf[i][j]); 
    Max(maxz,zbuf[i][j]) } }
fclose(zb);
loop (i,0,yzbuf) loop (j,0,xzbuf) zbuf[i][j]-=minz;
put2(minz,maxz)
}

#ifdef PURIFY
/* purify */
#define rg 8
loop (i,0,yzbuf-rg) loop (j,0,xzbuf-rg) if (zbuf[i][j]) {
  int x,y,s=0;
  loop (x,0,rg) loop (y,0,rg) s+=!!zbuf[i+x][j+y];
  if (s<rg*5/2) zbuf[i][j]=0; }

zb=fopen("pur.zba","wb");
loop (i,0,yzbuf)
  if (fwrite(zbuf[i],xzbuf,1,zb)!=1) Error("pur.zba write");
fclose(zb);
#endif

for (;;) {
int x,y;

getdata get(x) get(y) checkdata enddata

for (j=19; j>=0; j--)
  loop (i,0,20) prt_("%3d%c",(int)zbuf[x+i][y+j],i==19 ? '\n' : ' ');
  
}
}
