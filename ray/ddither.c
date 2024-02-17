/* square dithering, #included from dither.c and ppm2ps.c
This version implements a square diagonally placed pattern:
. . . . . . . .
 . . . . . . .
. . . . . . . .
 . . . . . . .
. . . . . . . .
 . . . . . . .
 (see sqdither.c for a parallel pattern)
*/

int nsq=6; /* square half-size */
int **sq; /* the pattern [j=0:2*nsq-1][i=-sq+1:sq-1] */

struct shade_s {
  double x;
  int i,j; };

int compar(const void *a,const void *b)
{
  return ((struct shade_s*)a)->x < ((struct shade_s*)b)->x;
}

void makesq(int NSQ) 
/* make dither patterns */
{
  int i,j,n,k;
  struct shade_s *shade; /* [n] */

  nsq=NSQ;
  if (nsq<=0) nsq=6;
  n=nsq*nsq*2;

  allocarray(shade,n);

  allocarray(sq,nsq*2-1);
  sq+=nsq-1;
  k=0;
  loop (j,-nsq+1,nsq) {
    allocarrayzero(sq[j],nsq*2);
    loop (i,0,2*nsq) sq[j][i]=-1; /* check only */ }

  loop (i,0,2*nsq) {
    int to=i;
    if (to>=nsq) to=nsq*2-1-to;
    loopto (j,-to,to) {
      shade[k].x=-cos(PI*(i+0.5)/nsq)*cos(PI*j/nsq);
      shade[k].i=i;
      shade[k].j=j;
      k++; } }

  if (k!=n) Error("ddither: pattern size (internal error)"); 

  qsort(shade,n,sizeof(shade[0]),compar);

  loop (k,0,n) sq[shade[k].j][shade[k].i]=(255*k+127)/n;

#if 0
  {
    FILE *f=fopen("pat.ppm","wt");
    fprintf(f,"P3\n%d %d\n255\n",nsq*2,nsq*2-1);
    loop (j,-nsq+1,nsq)
      loop (i,0,nsq*2)
        if (sq[j][i]<0) 
	  fprintf(f,"255 0 0\n"); 
        else
	  fprintf(f,"%d %d %d\n",(unsigned char)sq[j][i],(unsigned char)sq[j][i],(unsigned char)sq[j][i]);
    fclose(f); }
#endif

  if (n!=nsq*nsq*2) Error("ddither: pattern");
  free(shade);
}

int pix2pix(int i,int j,int inpix)
{
  int a=(i+j)/(2*nsq), b=(i-j+200000*nsq)/(2*nsq)-100000;
  int ic=a+b, jc=a-b;
  int ii=i-ic*nsq, jj=j-jc*nsq;

  if (sq[jj][ii]<0) Error("ddither: range check (internal error)");

  return sq[jj][ii]<inpix;
}
