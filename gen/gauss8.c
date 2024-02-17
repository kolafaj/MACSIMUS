/*
  Gauss 4 point formula for numerical integration, error ~ n^-8
  8n function calls
*/

REAL Gauss8(REAL (*f)(REAL),REAL a,REAL b,int n)
{
REAL h,sum,w1,w2,h1,h2,x;
int i;
#ifdef HIGHPREC
static REAL q1,q2,w;
static int pass;
if (!pass) {
  pass++;
  q1=sqrt((REAL)120);
  q2=sqrt((15-q1)/140);
  q1=sqrt((15+q1)/140);
  w=0.25-sqrt((REAL)5/864); 
}
#else
const REAL
  q1=0.430568155797026287612,
  q2=0.169990521792428132401,
  w=0.173927422568726928687;
#endif

if (n<=0) return 0;

sum=0;
h=(b-a)/n;
w1=h*w; w2=h/2-w1;
h1=h*q1; h2=h*q2;
loop (i,0,n) {
  x=h*(i+0.5)+a;
  sum+=(f(x-h1)+f(x+h1))*w1+(f(x-h2)+f(x+h2))*w2; }
return sum;
}
