/*.....#include "cpmark.h"*/

void inertiamatrix(char *IG,species_t *spec) /**************** inertiamatrix */
/*
  IG="inertia" -> the inertia matrix I = SUM mi (ri^2 delta - ri ri)
  IG="gyration" -> the gyration matrix G = SUM mi ri ri
  delta = unit tensor
*/
{
  int ns=spec->ns,i,j,k;
  FILE *cp,*plb;
  int NCP;
  float rcp[7];
  unsigned it;
  double M,m,rr;
  vector CM,IM[3],R[3],R0[3];
  static int ord0[3]={0,1,2};
  double *passIM[3],*passR[3];

  if (sizeof(double)*3 != sizeof(vector)) ERROR(("vector bad length"))

  site=spec->site;
  if (abs(option('r'))!=4) ERROR(("-r4 or -r-4 required"))

  loop (j,0,3) {
    passIM[j]=IM[j];
    passR[j]=R[j];
    loop (k,0,3) R0[j][k]=j==k; }

  cp=fopen(string("%s.cp",IG),"wt");
  rcp[0]=CPmark;
  if (IG[0]=='g') NCP=7;
  else NCP=3;
  *(int4*)(&rcp[1])=NCP;
  /*                          3---4---5---6---7--- */
  if (IG[0]=='g') copy(rcp+2,"G_zzRg_xRg_yRg_zRgyr",(NCP-2)*4);
  else            copy(rcp+2,"I_zz",(NCP-2)*4);

  fwrite(rcp,4,NCP,cp);

  plb=fopen(string("%s.mol",IG),"wt");
  fprintf(plb,"for %s axes: see blend -%c\n\
\n\
number_of_atoms = 6\n\
\n\
atoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n\
  0 X1         N     0.0000  0    1  1\n\
  1 X2         N     0.0000  0    1  0\n\
  2 Y1         CT    0.0000  0    1  3\n\
  3 Y2         CT    0.0000  0    1  2\n\
  4 Z1         ST    0.0000  0    1  5\n\
  5 Z2         ST    0.0000  0    1  4\n\
\n",IG,toupper(IG[0]));
  fclose(plb);
  prt("! %s.mol created",IG);

  plb=fopen(string("%s.plb",IG),"wb");
  rcp[0]=6;
  rcp[1]=0;
  fwrite(rcp,4,2,plb);

  prt("! %s matrix analysis",IG);

  for (it=0;;it++) {
    if (spec->frame>spec->Xopt.toframe) break;
    if (!read3D(spec,-1)) break;
    spec->frame+=spec->Xopt.byframe;
    M=0;
    VO(CM,=0)
    memset(IM,0,sizeof(IM));
    loop (i,0,ns) {
      m=atom[site[i].type].mass;
      M+=m;
      VV(CM,+=m*site[i].r) }
    VO(CM,/=M)
    loop (i,0,ns) {
      vector dr;
      m=atom[site[i].type].mass;
      VVV(dr,=site[i].r,-CM)
      if (IG[0]=='g')
        loop (j,0,3) loop (k,0,3) IM[j][k]+=m*dr[j]*dr[k];
      else {
        rr=SQR(dr);
        loop (j,0,3) loop (k,0,3) IM[j][k]+=m*((j==k)*rr-dr[j]*dr[k]); }
    }
    m=Jacobi(3,passIM,passR,spec->Xopt.Jeps);

    /* principal axes multiplied by "partial radii of gyration"
       then, the following Cartesian cross or `molecule':

            O
            |  O
            | /
            |/
       O----+----O
           /|
          / |
         O  |
            O

       where the masses of `atoms' O are 1/2 and the lengths of
       radiusvectors are R_x, R_y, R_z, has the same tensor (matrix) of
       inertia as the original molecule */

    prt("! diagonalized %s matrix in g/mol*AA^2: %g %g %g",
        IG,
        IM[0][0],IM[1][1],IM[2][2]);
    prt_("! principal moments of inertia in kg*m^2:");
    if (IG[0]=='g')
      prt(" %g %g %g",
        (IM[2][2]+IM[1][1])*1.66053904e-47,
        (IM[0][0]+IM[2][2])*1.66053904e-47,
        (IM[0][0]+IM[1][1])*1.66053904e-47);
    else 
      prt(" %g %g %g",
        IM[0][0]*1.66053904e-47,
        IM[1][1]*1.66053904e-47,
        IM[2][2]*1.66053904e-47);

    loop (i,0,3)
      loop (j,0,3) R[i][j]*=sqrt(fmax(0,IM[i][i])/M);

    { /* this code reorders the indices so that the above mentioned
         cross is as similar as possible as the cross from the previous
         step.  Hence, the principal axes do not swap from step to step
         (if possible, of course, if two R_i become almost the same,
         this is impossible in principle)
         BUG: the sign of both vectors in each O---O `bond' is undefined.
            This does not cause any problems while `show'-ing, but might
            cause problems in further calculations */

      double s,ss=-9;
      static int ord[6][3]={{0,1,2},{0,2,1},{1,2,0},{1,0,2},{2,0,1},{2,1,0}};
      int iord=-1;

      loop (i,0,6) {
        s=0;
        loop (j,0,3) s+=fabs(SCAL(R0[ord0[j]],R[ord[i][j]]));
        if (s>ss) iord=i,ss=s; }

      if (iord<0) ERROR(("intertia matrix: numeric problem"))
      copy(ord0,ord[iord],sizeof(ord0)); }

    m=IM[0][0]+IM[1][1]+IM[2][2];
    loop (i,0,3) {
      double Mdiag=IM[ord0[i]][ord0[i]];

      if (Mdiag<0) Mdiag=0;
      rcp[i]=m-Mdiag;
      rcp[i+3]=sqrt(Mdiag/M); }
    rcp[6]=sqrt(m/M);
    fwrite(rcp,4,NCP,cp); /* for I, only part written */

    loop (i,0,3) {
      loop (j,0,3)
        rcp[j+3]=-(rcp[j]=R[ord0[i]][j]);
      fwrite(rcp,4,6,plb); }

    copy(R0,R,sizeof(R)); }

  fclose(cp);
  fclose(plb);

  prt("! %d frames analyzed, inertia.cp and inertia.plb closed",it);
}
