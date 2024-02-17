/* cc -O2 -o density density.c -lm 
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
 FILE *plb;
 char line[128],*ch;
 int i,sum=0,nf=0,ns,frame;
 float r[3],R[3],rad,radq;
 float hdr[2];

 if (narg<6) {
   fprintf(stderr,"Call by:\n\
  density FILE.plb[:FROMFRAME] X Y Z R\n\
calculates # of atoms in given sphere of center (X,Y,Z) and radius R\n");
   exit(0); }

 strcpy(line,arg[1]);
 if ( (ch=strchr(line,':')) ) *ch=0,frame=atoi(ch+1);
 else frame=1;
 if (frame<1) Error("frame<1");
 
 if (!(plb=fopen(line,"rb"))) Error(line);
 R[0]=atof(arg[2]);
 R[1]=atof(arg[3]);
 R[2]=atof(arg[4]);
 rad=atof(arg[5]);
 radq=Sqr(rad);
 
 fread(hdr,4,2,plb);
 ns=hdr[0];
 fseek(plb,ns*(frame-1)*sizeof(r),SEEK_CUR);

 for (;;) {
   int nin=0;

   loop (i,0,ns) {
     if (fread(r,sizeof(r),1,plb)!=1)
       if (i) Error("unexpected EOF");
       else goto end;
     nin += Sqr(r[0]-R[0])+Sqr(r[1]-R[1])+Sqr(r[2]-R[2]) < radq; }
   printf("%4d %g\n",nin,nin/(4*PI/3*rad*radq));
   nf++; sum+=nin;
   }
 end:
 fprintf(stderr,"%g %g\n",sum/(double)nf,sum/(double)nf/(4*PI/3*rad*radq));

 fclose(plb);
 return 0;
}
