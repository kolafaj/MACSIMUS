/* cc -O2 -o project project.c -lm
 */
#include "../gen/include.h"

double r[3],x,y;
double zeye;
double rot[3][3]={{1,0,0},{0,1,0},{0,0,1}};

void project(void) /****************************** project */
{
  double q=zeye+(rot[2][0]*r[0]+rot[2][1]*r[1]+rot[2][2]*r[2]);

  q=zeye/q;

  x=q*(rot[0][0]*r[0]+rot[0][1]*r[1]+rot[0][2]*r[2]);
  y=q*(rot[1][0]*r[0]+rot[1][1]*r[1]+rot[1][2]*r[2]);
}

void rotate(int axis,double angle) /******************************** rotate */
{
 double sa=sin(angle), ca=cos(angle);
 double o[3][3], rot2[3][3];
 int i=(axis+1)%3,j=(i+1)%3,k=(j+1)%3,a,b;

 o[i][i]=ca; o[i][j]=-sa; o[i][k]=0;
 o[j][i]=sa; o[j][j]= ca; o[j][k]=0;
 o[k][i]= 0; o[k][j]=  0; o[k][k]=1;

 loop (a,0,3) loop (b,0,3)
   rot2[a][b]=o[a][0]*rot[0][b]+o[a][1]*rot[1][b]+o[a][2]*rot[2][b];
 
 loop (a,0,3) loop (b,0,3)
   rot[a][b]=rot2[a][b];
}

int main(int narg,char **arg)
{
  char line[1024];
  double u;

  if (narg<5) {
    fprintf(stderr,"Call by:\n\
  project unit phi theta alpha fromz < xyz-file > xy-file\n");
    exit(0); }

  zeye=atof(arg[5]);
  u=atof(arg[1]);
  if (u==0) u=PI/180;
  rotate(2,u*atof(arg[2]));
  rotate(0,u*atof(arg[3]));
  rotate(2,u*atof(arg[4]));

  while (gets(line)) {
    int i=sscanf(line,"%lf%lf%lf",r,r+1,r+2);
    if (i<3) puts(line);
    else {
      project();
      printf("%.9g %.9g\n",x,y); } }

  return 0;
}
