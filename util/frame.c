/* cc -o frame frame.c
 gcc -o frame.exe frame.c
*/

#include "../gen/include.h"

int main(int narg,char **arg)
{
char fn[64],old[64],*c,*dot=NULL;
FILE *plb,*new;
long i,ns,NS=0;
long pos,frame=1;
float r[3];
char key='b';

if (narg<2) {
  fprintf(stderr,"Extracts frame from a plb-file. Call by:\n\
  frame {FILE|FILE.plb|FILE.p00|FILE.ppp} [FRAME [NS [{a|t|b}]]]\n\
FILE  means FILE.plb, if not found then FILE.p00, if not found then FILE.ppp\n\
FRAME is frame #, default=1=1st frame, -1=last frame\n\
NS    optional number of sites to truncate; 0=all\n\
a=t=ASCII file FILE.3dt created\n\
b=bin file FILE.3db (default)\n");
  exit(0); }

strcpy(fn,arg[1]);
for (c=fn; *c; c++) if (*c=='.') dot=c;
if (dot && (!strcmp(dot+1,"plb") || !strcmp(dot+1,"p00") || !strcmp(dot+1,"ppp")))
  plb=fopen(fn,"rb");
else {
  strcat(fn,".plb");
  plb=fopen(fn,"rb");
  if (!plb) {
    strcpy(fn,arg[1]); strcat(fn,".p00");
    plb=fopen(fn,"rb"); 
    if (!plb) {
      strcpy(fn,arg[1]); strcat(fn,".ppp");
      plb=fopen(fn,"rb"); } } }

if (!plb) Error("no plb-file");

printf("%s opened\n",fn);

if (narg>2) frame=atoi(arg[2]);
if (frame==0||frame<-1) Error("bad frame");
if (narg>3) NS=atoi(arg[3]);
if (narg>4) key=arg[4][0];
if (!strchr("atbATB",key)) Error("bad key, use one of a,b");

if (fread(r,sizeof(float),2,plb)!=2) Error("no header");
ns=r[0];
if (ns<=0) Error("bad header:ns");
printf("%ld sites\n",ns);
if (NS>ns || NS<0) Error("bad NS");

if (frame<0) {
  fseek(plb,0,SEEK_END);
  pos=ftell(plb);
  frame=(pos-8)/(ns*sizeof(float)*3L);
  if (pos!=8+frame*(ns*sizeof(float)*3L))
    fprintf(stderr,"WARNING: bad file length\n");  }

rewind(plb);
fseek(plb,8+(frame-1L)*(ns*sizeof(float)*3L),SEEK_CUR);

dot=NULL;
for (c=fn; *c; c++) if (*c=='.') dot=c;
if (!dot) Error("?");

strcpy(dot,".3dx");
strcpy(old,fn);
strcpy(dot,key=='b'?".3db":".3dt");
remove(old);
rename(fn,old);

new=fopen(fn,key=='b'?"wb":"wt");
if (!new) Error("open output");

if (NS==0) NS=ns;
loop (i,0,NS) {
  if (fread(r,sizeof(float),3,plb)!=3) Error("read plb or too short");
  if (key=='b') fwrite(r,sizeof(float),3,new);
  else fprintf(new,"%9.4f %9.4f %9.4f\n",r[0],r[1],r[2]); }
fclose(new);
printf("%ld sites of frame %ld written to %s\n",NS,frame,fn);
fclose(plb);
}
