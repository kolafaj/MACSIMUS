typedef struct sort_s {
  double l;
  int i;
} sort_t; /* sizeof(sort_t)<=2*sizeof(double) assumed, dirty memory re-use */

static int dblcmp(const void *a,const void *b)
{
  sort_t *A=(sort_t *)a;
  sort_t *B=(sort_t *)b;
  if (A->l<B->l) return -1;
  else return A->l>B->l;
}

void writevibration(species_t *spec,double **R,int n,sort_t *sort,
  vector *r0,char *tit,int ns)
/* used to write motion of normal modes and essential dynamics */
{
if (R) {
  int from=0,to=n,k,indx,ix,j,i;
  FILE *plb;
  float r[3];
  double m,a;
  double ampl=Xopt.amplitude;

  /* negative -M: compromise between 2D and 3D molecules */
  /* prior 2.2h: ampl=-ampl*sqrt(ns) -- 2D molecules */
  if (ampl<0) ampl=-0.07*ampl*pow(ns,0.4);

  if (spec->Xopt.P>0) to=spec->Xopt.P;
  else from=n+spec->Xopt.P;
  Max(from,0)
  Min(to,n)

  if (from<to) {
    prt("! writing frames from %d to %d (incl)",from,to-1);
    loop (j,from,to) {
      indx=sort[j].i;
      sprintf(spec->ext,".%s%04d.plb",tit,j);
      plb=fopen(spec->fn,"wb");
      if (!plb) ERROR(("cannot write to %s",spec->fn))
      else {
        if (tit[0]=='n') {
          /* normal modes:
             ampl = max amplitude of moves over all atoms */
          m=a=0;
          loop (i,0,n) {
            R[indx][i]/=sqrt(atom[site[i/3].type].mass);
            a+=Sqr(R[indx][i]);
            if (i%3==2) { Max(m,a) a=0; } }
          m=ampl/sqrt(m); }
        else {
          /* essential dynamics:
             Xopt.amplitude = in units of the ess. motion = sqrt(eigenval) */
          m=0;
          loop (i,0,n) m+=Sqr(R[indx][i]);
          m=spec->Xopt.amplitude*sqrt(-sort[j].l)
            * sqrt((2+(spec->Xopt.F<0))/m);
          /*
            note: the last term 2+(opt_F<0) guarantees that the stdev
            (mean quadratic amplitude) is the same for the shown
            harmonic (opt_F>0) or linear/\/\/\/\ (opt_F<0) motion
          */
          }
        r[0]=ns; r[1]=0;
        fwrite(r,4,2,plb);
        loop (ix,0,abs(spec->Xopt.F)) {
          double x=ix/(double)(abs(spec->Xopt.F)-1);

          if (spec->Xopt.F>0) x=m*cos(x*PI);
          else x=m*(1-2*x);
          loop (i,0,ns) {
            loop (k,0,3) r[k]=r0[i][k]+R[indx][i*3+k]*x;
            fwrite(r,4,3,plb); } }
        }
      fclose(plb); }
    prt("! %s%04d.plb -- %s%04d.plb written",tit,from,tit,to-1); }
  }
}

void normalmodes(species_t *spec) /***************************** normalmodes */
{
  int ns=spec->ns,n=3*ns,size=n*sizeof(double),i,j;
  double **A,**R=NULL,*d0;
  double rr,x;
  vector *r0,*fp,*fm;
  sort_t *sort;
  FILE *nmv;

#define D(X) ((double*)(X[0]))

  if (sizeof(double)*3 != sizeof(vector)) ERROR(("vector bad length"))

  prt("! normal mode analysis, %dx%d matrix of d2U/drdr, dr=%g",n,n,spec->Xopt.dr);

  ralloc(A,n*sizeof(A[0])); /* this is the mark for release */
  loop (i,0,n) ralloczero(A[i],size);
  if (spec->Xopt.P) {
    ralloc(R,n*sizeof(R[0]));
    loop (i,0,n) ralloczero(R[i],size); }
  ralloc(r0,size); d0=D(r0);
  ralloc(fp,size*2); fm=fp+ns; /* .. because of sorting eigenvectors */

  site=spec->site;
  loop (i,0,ns) VV(r0[i],=site[i].r)

  loop (i,0,n) {
    rr=d0[i];
    d0[i]=rr+spec->Xopt.dr;
    Upot(r0,fp,spec,0);
    d0[i]=rr-spec->Xopt.dr;
    Upot(r0,fm,spec,0);
    d0[i]=rr;
    loop (j,0,n) {
      x=(D(fp)[j]-D(fm)[j])/spec->Xopt.dr
        /sqrt(atom[site[i/3].type].mass*atom[site[j/3].type].mass);
      A[i][j]-=x;
      A[j][i]-=x; }
    }
  /* ... A is four times the correct (1/2)d2U/dxdy/sqrt(mi mj) */

  rr=Jacobi(n,A,R,spec->Xopt.Jeps);
  /* ... and now:
     A=diag(omega_i^2),
     R[i][j]/sqrt(mj) = i-th normal mode eigenvector */
  prt("! Jacobi diagonalization err=%g",rr);

  sort=(struct sort_s*)fp;
  if (sizeof(struct sort_s)>2*sizeof(double)) ERROR(("internal"))
  loop (i,0,n) {
    sort[i].l=A[i][i];
    sort[i].i=i; }
  qsort(sort,n,sizeof(struct sort_s),dblcmp);

  strcpy(spec->ext,".nmf");
  nmv=fopen(spec->fn,"wt");

fprintf(nmv,"\
# normal mode frequencies\n\
# matrix n=%d   d2U/drdr=-df/dr, dr=2*%g   Jacobi diag.err.=%g\n\
#  i        [THz]          [1/cm]\n",n,spec->Xopt.dr,rr);
 rr=kcal*1e23/Sqr(4*PI); /* [A] = kcal/mol/AA^2/(g/mol) = kcal*1e23 */
  loop (i,0,n) {
    x=sort[i].l*rr;
    if (x<0) x=-sqrt(-x); else x=sqrt(x);
    fprintf(nmv,"%4d %15.8f %15.8f\n",i,x*1e-12,x/29979245800.); }

  fclose(nmv);
  prt("! %s written",spec->fn);

  writevibration(spec,R,n,sort,r0,"nm",ns);
  release(A);
}

