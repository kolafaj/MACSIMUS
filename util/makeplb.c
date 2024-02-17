/* make makeplb
 */
#include "ground.h"
#include "alloc.h"
#include "rndgen.h"

int main(int narg,char **arg)
{
 FILE *plb;
 float hdr[2];
 float (*r)[3],(*rL)[3],cm[3];
 float L,dr; 
 int i,k,ns,iframe,nframes;

 if (narg<5) {
   fprintf(stderr,"Call by:\n\
  makeplb FILE L NS NFRAMES DR\n");
   exit(0); }

 initscroll(0);
 rndinit(0,0);

 L=atof(arg[2]); 
 ns=atoi(arg[3]); 
 nframes=atoi(arg[4]);
 dr=atof(arg[5]);
 plb=fopen(arg[1],"wb");

 hdr[0]=ns; hdr[1]=L;
 alloc(r,ns*12);
 alloc(rL,ns*12);

 fwrite(hdr,4,2,plb);

 loop (i,0,ns) loop (k,0,3) r[i][k]=rnd()*L;

 loop (iframe,0,nframes) {
   cm[0]=cm[1]=cm[2]=0;
   loop (i,0,ns) 
     loop (k,0,3) {
       if (i<ns/2)
	 r[i][k] += (rndcos()-rndcos()+rndcos()-rndcos()+rndcos()-rndcos())*dr;
       else
	 r[i][k] = r[i-ns/2][k] + (rndcos()-rndcos()+rndcos()-rndcos()+rndcos()-rndcos())*dr;
/*.....if (iframe==13 && i==7 && k==1) r[i][k]+=0.4*L;*/
       cm[k] += r[i][k]; }

   loop (k,0,3) cm[k] /= ns;

   loop (i,0,ns) 
     loop (k,0,3) {

       rL[i][k] = r[i][k] - cm[k];

       while (rL[i][k]>=L) rL[i][k]-=L;
       while (rL[i][k]<0) rL[i][k]+=L; }

   fwrite(rL,12,ns,plb); }

return 0;
}
