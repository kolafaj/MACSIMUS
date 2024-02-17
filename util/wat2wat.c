/* cc -O2 -Wall -o wat2wat wat2wat.c -lm
*/

#include "../gen/include.h"
#define real float
#include "../gen/vector3d.h"

char key;
int atm;
double Q=1;

void mywrite(float *r) /******************************************** mywrite */
{
  if (atm) printf("%c %9.5f %9.5f %9.5f\n",key,VARG(r));
  else fwrite(r,4,3,stdout);
}

int main(int narg,char **arg) /**************************************** main */
{
  static float header[2];
  vector *r;
  vector L;
  vector shift={0,0,0};
  int i,ns,nm,l1,i1,l2,i2,iarg,varL;
  double BETA=111;
  double OH=0.9572;
  double ALPHA=104.52;
  double MO=0.1546;
  double QO,QH;
  char *comma;

  if (narg<3) {
    fprintf(stderr,"\
Convert configurations by various water models. Call by:\n\
  wat2wat INPUT-SITES OUTPUT-SITES [OPTIONS] < INPUT-PLB > OUTPUT-PLB\n\
SITES = string of OHML; M=tip4p-like charged site, L=NE or ST2 lone pair\n\
OPTIONS:\n\
  -qQ    scale coordinates+box by Q [default=1]\n\
  -x     export in the .atm format without box size\n\
  -l     export in the .atm format with box size included\n\
  -rHO   HO distance, in A [0.9572] (note: NE=0.98)\n\
  -mMO   MO distance, in A [TIP4P/2005=0.1546] (TIP4P/Ice=0.1577, NE=0.23)\n\
  -aHOH  the HOH angle, in degrees [104.52] (NE=108)\n\
  -bBETA the LOL angle, in degrees [NE=111]\n\
  -sX[,Y[,Z]]\n\
         shift (add to all distances), before scaling (-q) [default=0]\n\
O and H are copied, M,L calculated if necessary\n\
See also:\n\
  m2m plb2plb plbmerge mol2mol ice\n");
    exit(0); }
    
  loop (iarg,3,narg) switch (arg[iarg][1]) {
    case 'q': Q=atof(arg[iarg]+2); break;
    case 'x': atm=1; break;
    case 'l': atm=2; break;
    case 'r': OH=atof(arg[iarg]+2); break;
    case 'a': ALPHA=atof(arg[iarg]+2); break;
    case 'b': BETA=atof(arg[iarg]+2); break;
    case 'm': MO=atof(arg[iarg]+2); break;
    case 's':
      shift[0]=shift[1]=shift[2]=atof(arg[iarg]+2);
      comma=strchr(arg[iarg],',');
      if (comma) {
        shift[1]=shift[2]=atof(comma+1);
        comma=strchr(comma+1,',');
        if (comma) shift[2]=atof(comma+1); }
      break;
    default: Error("wat2wat: bad option"); }

  QH=OH*cos(ALPHA/360*PI);
  QO=(QH-MO)/QH;
  QH=(1-QO)/2;

  if (fread(header,4,2,stdin)!=2) Error("no header on input");
  ns=header[0];
  varL=header[1]<0;
  if (!varL) L[0]=L[1]=L[2]=header[1]*Q;
  l1=strlen(arg[1]);
  l2=strlen(arg[2]);
  nm=ns/l1;
  if (ns%l1) Error("wat2wat: total number of sites is not an integer multiple of water sites");
  header[0]=nm*l2;
  if (atm) printf("%d\n",nm*l2);
  else {
    header[1]=-3;
    fwrite(header,4,2,stdout); }
  
  allocarray(r,ns);
  
  for (;;) {
    if (varL) {
      if (fread(L,4,3,stdin)!=3) break;
      VO(L,*=Q) }
    i=fread(r,12,ns,stdin);
    if (i==0) break;
    if (i!=ns) Error("wat2wat: plb is truncated");
    if (!atm) fwrite(L,4,3,stdout);
    if (atm==2) printf("  %9.5f %9.5f %9.5f\n",VARG(L));
    else if (atm==1) printf("\n");
    loop (i,0,ns) {
      VV(r[i],+=shift)
      VO(r[i],*=Q) }

    for (i=0; i<ns; i+=l1) {
      vector rO,rH1,rH2,rM,rL1,rL2,rX;
      int nH=0,nL=0,nO=0,nM=0,nX=0;
    
      loop (i1,0,l1) switch (arg[1][i1]) {
        case 'O': nO++; VV(rO,=r[i+i1]); break;
        case 'M': nM++; VV(rM,=r[i+i1]); break;
        case 'H': nH++; 
                  if (nH==1) VV(rH1,=r[i+i1]) 
                  else VV(rH2,=r[i+i1]) 
                  break;
        case 'L': nL++; 
                  if (nL==1) VV(rL1,=r[i+i1]) 
                  else VV(rL2,=r[i+i1]); 
        default: nX++; VV(rX,=r[i+i1]); }
        
      if (nO!=1 || nH!=2) Error("wat2wat: wrong number of O or H");  
      if (!nM)
        VVVV(rM,=QO*rO,+QH*rH1,+QH*rH2)
      if (nL!=2) {
        vector a,p,h1,h2;
        double r;
        
        VVV(h1,=rH1,-rO)
        VVV(h2,=rH2,-rO)
        VECT(p,h1,h2)
        VVV(a,=h1,+h2)
        r=sin(PI/360*BETA)/sqrt(SQR(p));
        VO(p,*=r)
        r=cos(PI/360*BETA)/sqrt(SQR(a));
        VO(a,*=r)
        
        VVVV(rL1,=rO,-a,+p)
        VVVV(rL2,=rO,-a,-p) }

      nH=0;
      loop (i2,0,l2) switch (key=arg[2][i2]) {
        case 'O': mywrite(rO); break;
        case 'M': mywrite(rM); break;
        case 'H': nH++; 
                  if (nH==1) mywrite(rH1); 
                  else mywrite(rH2); 
                  break;
        case 'L': nL++; 
                  if (nL==1) mywrite(rL1); 
                  else mywrite(rL2); 
                  break;
        default: 
          if (nX) mywrite(rX); 
          else Error("wat2wat: no X on input"); } } }

  return 0;
}
