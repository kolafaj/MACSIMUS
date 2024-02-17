/* cc -O2 -o ppminfo ppminfo.c
*/
#include "../gen/include.h"

FILE *in;
char line[256];

int xsize,ysize,depth;

int main(int narg,char **arg)
{
 int iarg,type,com;

 if (narg<2) {
   fprintf(stderr,"\
Call by:\n\
  ppminfo {FILE|-} (verbose)\n\
  ppminfo FILE FILE ... (brief)\n\
Return codes:\n\
  1 if PNM and depth<=255\n\
  0 if PNM and depth>255\n\
  other on error\n\
Note: in shell, 0 is treated as true:\n\
  if ppminfo x.ppm then; echo \"depth<=255\"; fi\n\
");
  exit(0); }

 loop (iarg,1,narg) {

   if (strcmp(arg[iarg],"-")) in=fopen(arg[iarg],"rb");
   else in=stdin;

   if (!in) Error(arg[iarg]);

   printf("%s: ",arg[iarg]);
   fgets(line,256,in);
   if (narg>2) *strchr(line,'\n')=0;
   fputs(line,stdout);
   if (line[0]!='P') Error("bad header: Pn expected");
   type=line[1];
   if (type<'1' || type>'6') Error("bad type: P{1,2,3,4,5,6} expected");
   do {
     if (!fgets(line,256,in)) Error("format in");
     com=line[0]=='#';
     if (!com) if (narg<=2) printf("size: ");
     if (narg<=2) fputs(line,stdout);
     } while (com); 
   sscanf(line,"%d%d",&xsize,&ysize);
   if (strchr("2356",type)) {
     fgets(line,256,in);
     sscanf(line,"%d",&depth);
     if (narg>2) printf(" %d %d %d\n",xsize,ysize,depth);
     else printf("depth: %d\n",depth); }
   else if (narg>2) printf(" %d %d\n",xsize,ysize);
   fclose(in); }

 return depth<=255;
}
