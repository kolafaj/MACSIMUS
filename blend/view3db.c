/* cc -O2 -o view3db view3db.c
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
 float r[3];
 FILE *f;

 if (narg<2) {
   fprintf(stderr,"View a 3db file. Call by:\n\
  view3db 3db-FILE\n\
  view3db 3db-FILE > 3dt-FILE\n");
   exit(0); }

 f=fopen(arg[1],"rb");
 if (!f) Error(arg[1]);

 _n

 while (fread(r,4,3,f)==3)
   printf("%8.5f %8.5f %8.5f\n",r[0],r[1],r[2]);
 fclose(f);

 return 0;
}
