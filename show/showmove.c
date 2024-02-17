int lasti=0; /* llasti is local */

void moveselection(double dx,double dy)
{
  int i,k;
  double d[3]={dx,dy,0},D[3];

  D[0]=drot[0][0]*d[0]+drot[1][0]*d[1]+drot[2][0]*d[2];
  D[1]=drot[0][1]*d[0]+drot[1][1]*d[1]+drot[2][1]*d[2];
  D[2]=drot[0][2]*d[0]+drot[1][2]*d[1]+drot[2][2]*d[2];

  loop (i,0,ns) if (site[i].mark&4)
    loop (k,0,3) cfg[i][k]+=D[k];
}

void rotateselection(int axis,double angle) /*************** rotateselection */
{
  double sa=sin(angle), ca=cos(angle);
  double o[3][3];
  int i=(axis+1)%3,j=(i+1)%3,k=(j+1)%3;
  double center[3]={cfg[lasti][0],cfg[lasti][1],cfg[lasti][2]};
  double d[3],D[3];
  double RR[3][3],R[3][3];

  o[i][i]=ca; o[i][j]=-sa; o[i][k]=0;
  o[j][i]=sa; o[j][j]= ca; o[j][k]=0;
  o[k][i]= 0; o[k][j]=  0; o[k][k]=1;

  loop (i,0,3) loop (j,0,3) {
    RR[i][j]=0;
    loop (k,0,3) RR[i][j]+=o[i][k]*drot[k][j]; }
  loop (i,0,3) loop (j,0,3) {
    R[i][j]=0;
    loop (k,0,3) R[i][j]+=drot[k][i]*RR[k][j]; }

  loop (i,0,ns) if (site[i].mark&4) {
    loop (k,0,3) d[k]=cfg[i][k]-center[k];

    D[0]=R[0][0]*d[0]+R[0][1]*d[1]+R[0][2]*d[2];
    D[1]=R[1][0]*d[0]+R[1][1]*d[1]+R[1][2]*d[2];
    D[2]=R[2][0]*d[0]+R[2][1]*d[1]+R[2][2]*d[2];

    loop (k,0,3) cfg[i][k]=D[k]+center[k]; }
}
