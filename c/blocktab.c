/* cc -O2 -o blocktab blocktab.c -lm
*/
#include "../gen/include.h"

#define NL 32000
char line[NL];

int ncol;
double *data,*sum;

int nextline(void) /*********************************************** nextline */
{
  char *c;
  int i;

  do
    if (fgets(line,NL,stdin)==NULL) return 0;
  while ( strchr("!#\n",line[0]) || !(c=strtok(line," \t\n")) );

  loop (i,0,ncol)
    if (c) {
      data[i]=atof(c);
      c=strtok(NULL," \t\n"); }
    else {
      fprintf(stderr,"blocktab: not enough data in line\n");
      exit(1); }
  
  return 1;
}

int main(int narg, char **arg) /*************************************** main */
{
  int block,i,j,makesubav=1;
  char *fmt=" %.12g";

  if (narg<3) {
    fprintf(stderr,"\
Blocking (making sub-averages) of equal-grid data.  Call by:\n\
  blocktab [-]BLOCK NCOL [FORMAT] < INPUT > OUTPUT\n\
Arguments:\n\
  BLOCK  block length, make subaverage\n\
  -BLOCK block length, make sums over blocks\n\
  NCOL   # of columns\n\
  FORMAT output format, default=\"%s\"\n\
See also:\n\
  gblock smooth mergetab tabproc sum redtab datablock sparsetab\n",fmt);
    exit(0); }

  block=atoi(arg[1]);
  if (block<0) makesubav=0,block=-block;
  if (block<1) {
    fprintf(stderr,"blocktab: invalid block size\n");
    exit(1); }
  
  ncol=atoi(arg[2]);
  if (narg>3) fmt=arg[3];
  
  allocarray(data,ncol);
  allocarrayzero(sum,ncol);
  
  j=0;
  while (nextline()) {
    loop (i,0,ncol) sum[i]+=data[i];
    j++;
    if (j==block) {
      
      loop (i,0,ncol) {
        if (makesubav) sum[i]/=block;
        printf(fmt,sum[i]);
        sum[i]=0; }
      printf("\n"); 
      j=0; } }
  
  if (j)
    fprintf(stderr,"blocktab: incomplete last block (%d line%s) ignored\n",
            j,"s"+(j==1));
  
  return 0;
}

