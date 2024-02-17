REAL Asum[N],Aqsum[N];
volatile int sig;

#define MAXn 32768
int n;
double x[MAXn],y[MAXn],dy[MAXn],dy0[MAXn],y_0[MAXn];

REAL ff(REAL *A)
{
  REAL s=0,f;
  int i;

  loop (i,0,n) {
    f=func(A,x[i]);
    s+=Sqr((f-y[i])/dy[i]); }

  return s;
}

REAL xsqrt(REAL x)
{
  if (x<0) return 0;
  else return sqrt(x);
}

void dody(double err,double sigma)
{
  int i;

  if (err==-3) loop (i,0,n) dy[i]=dy0[i]*sigma;
  else if (err<0) loop (i,0,n) dy[i]=dy0[i];
  else loop (i,0,n) dy[i]=pow(fabs(y[i]),err)*sigma;
}

int main(int narg,char **arg)
{
  char line[1024];
  FILE *f;
  REAL D=1e-5,eps=1e-8,par=1;
  int method=0,maxit=1000;
  int i,nerr=0,ierr,isstdev=0;
  double err=-2,sigma,sqsum;

  initscroll(0);
  rndinit(0,0);

  if (narg<2) {
    fprintf(stderr,"Call by (normally from plot):\n\
  fit INPUTFILE\n");
    exit(0); }

  f=fopen(arg[1],"rt");
  if (!f) ERROR(("open %s",arg[1]))

  while (fgets(line,1024,f)) {
    if (line[0]=='#') continue;
    if (n>=MAXn) ERROR(("tab ofl, increase MAXn and recompile"))
    i=sscanf(line,"%lf%lf%lf",x+n,y+n,dy0+n);
    if (i==1) ERROR(("%s bad format",line))
    if (i<=0) continue;
    if (i==2) {
      if (n) dy0[n]=dy0[n-1];
      else dy0[n]=1; }
    if (i==3) isstdev++;
    n++; }

  for (;;) {
