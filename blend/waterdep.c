/* cc -O2 -Wall -o waterdep waterdep.c -lm
 */
#include "../gen/include.h"

double dsin(double x)
{
  return sin(PI/180*x);
}

double dcos(double x)
{
  return cos(PI/180*x);
}

char *getsite(char *site,int n,char *from)
{
  *site++=n+'0';
  *site++=' ';
  while (*from>='0') *site++=*from++;

  //  fprintf(stderr,"%d %s\n",n,from); sleep(1);

  return from;
}

int main(int narg,char **arg)
{
  double OH,HOH,MO,LO,LOL=180;
  double xH,yH,wM,wH,zL,xL;
  double al[3]={1.4679,1.5284,1.4146};
  double AHH,ALL,AHL;
  char key;
  int sg,n,iarg;
#define XXX "? ?   "
  char L1[8]=XXX,L2[8]=XXX,H1[8]=XXX,H2[8]=XXX,O[8]=XXX,M[8]=XXX,X[8]=XXX;
  char *a;

  if (narg<6) {
    fprintf(stderr,"Calculate dependants of water models and fluctuating charge data.\n\
          H    y\n\
         /     ^              NOTES:\n\
    L.  /      |                HOH<180 deg assumed\n\
      'O-M     +-----> x        LOL<180 deg => L_x<0 (NE6 model)\n\
    L\"  \\    z\"                 LOL>180 deg => L_x>0 (L=split M)\n\
         \\\n\
          H\n\
Call by:\n\
  waterdep SITES OH HOH M MO\n\
  waterdep SITES OH HOH L LO LOL\n\
ARGUMENTS:\n\
  SITES = comma-separated site list as in ble-file (1st letter=[OHLM], 4 chars)\n\
  M = site M in the plane (M is given as is)\n\
  L = site L (lone electron) (L is given as is)\n\
  OH = O-H distance in A\n\
  HOH = H-O-H angle in degrees\n\
  MO = M-O distance in A\n\
  LO = L-O distance in A\n\
  LOL = L-O-L angle in degrees\n\
Examples):\n\
  waterdep L6NE,H6NE,O6NE,M6NE,L6NE,H6NE 0.98 108 M 0.23       # NE6:M\n\
  waterdep L6NE,H6NE,O6NE,M6NE,L6NE,H6NE 0.98 108 L 0.8892 111 # NE6:L1,2\n\
  waterdep LFQ4,HFQ4,OFQ4,LFQ4,HFQ4 0.98 108 L 0.2 330         # FQ4:L1,2\n\
  waterdep L54D,H54D,O54D,L54D,H54D 0.98 108 L 0.525 180       # POL4D:L1,2\n\
");
    exit(0); }

  for (n=0,a=arg[1]; *a; n++) {
    if (*a==',') a++;
    if (*a=='L') { memcpy(L1,L2,8); a=getsite(L2,n,a); }
    else if (*a=='H') { memcpy(H1,H2,8); a=getsite(H2,n,a); }
    else if (*a=='O') a=getsite(O,n,a);
    else if (*a=='M') a=getsite(M,n,a);
    else a=getsite(X,n,a); }

  OH=atof(arg[2]);
  HOH=atof(arg[3]);
  key=arg[4][0];
  iarg=6;

  LO=MO=atof(arg[5]);
  if (key=='L') {
    iarg=7;
    LOL=atof(arg[6]);
    MO=-LO*dcos(LOL/2); }

  if (iarg<narg) al[0]=atof(arg[iarg++]);
  if (iarg<narg) al[1]=atof(arg[iarg++]);
  if (iarg<narg) al[2]=atof(arg[iarg++]);

  xH=OH*dcos(HOH/2);
  yH=OH*dsin(HOH/2);
  wM=(xH-MO)/xH;
  wH=MO/xH/2;

  if (key=='L') printf("\n\
! %s\n\
! OH=%g HOH=%g LO=%g LOL=%g\n",arg[1],OH,HOH,LO,LOL);
  if (key=='M') printf("\n\
! %s\n\
! OH=%g HOH=%g MO=%g\n",arg[1],OH,HOH,MO);

  printf("\nbonds\n\
!i atom   i atom  K[kcal/mol]  r[AA]   calc.     Upot\n\
 %s  %s  450 %.8g  %g 0.0\n\
 %s  %s  450 %.8g  %g 0.0\n\
! bonds equivalent to constrained angles:\n\
 %s  %s  200 %.8g  %g 0.0\n",
         H1,O,OH,OH,
         H2,O,OH,OH,
         H1,H2,yH*2,yH*2);

  printf("\ndependants\n\
!i atom   #   i atom  weight...\n");

  if (key=='M')
    printf("M  %s  3  %s %11.8f  %s %11.8f  %s %11.8f e\n",
               M,     O, wM,     H1,wH,     H2, wH);

  else {
    for (sg=-1; sg<=1; sg+=2) {
      zL=sg*LO*dsin(LOL/2);
      xL=-LO*dcos(LOL/2);

      printf("L  %s  3  %s %11.8f  %s %11.8f  %s %11.8f",
             sg<0?L1:L2, O,wM,     H1,wH,     H2,wH);

      /* local x-coord system */
      printf("  %11.8f %11.8f %11.8f",-1/xH,0.5/xH,0.5/xH);
      printf("  %11.8f %11.8f %11.8f",0.0,0.5/yH,-0.5/yH);
      printf("  %11.8f",zL);
      printf("  %11.8f %11.8f %11.8f",zL/xH,-zL/2/xH,-zL/2/xH);
      printf("  %11.8f %11.8f %11.8f e\n",0.0,-zL/2/yH,zL/2/yH); }

    AHH=2*Sqr(yH)/al[1];
    ALL=2*Sqr(zL)/al[2];
    AHL=(AHH+ALL-2*Sqr(xH-xL)/al[0])/4;

    printf("\nfq4\n\
!i atom   i atom  i atom  i atom  AHH         ALL         ALH\n\
 %s  %s  %s  %s  %11.8f %11.8f %11.8f\n",
           H1,H2,L1,L2,AHH,ALL,AHL); }

  return 0;
}
