#define WRITEMATCH /* forces writing the matched trajectory (with
                      rotations and translations removed) */

#include "cpmark.h"
vector *match_r;

void match_center(vector *r,int ns)
{
  vector c;
  int i;

  VO(c,=0)
  loop (i,0,ns) VV(c,+=r[i])
  loop (i,0,ns) VV(r[i],-=(1./ns)*c)
}

void match_rotate(vector *r,int ns,int axis,double angle)
{
  double sa=sin(angle), ca=cos(angle);
  double o[3][3];
  int i=(axis+1)%3,j=(i+1)%3,k=(j+1)%3;
  vector dr;

  o[i][i]=ca; o[i][j]=-sa; o[i][k]=0;
  o[j][i]=sa; o[j][j]= ca; o[j][k]=0;
  o[k][i]= 0; o[k][j]=  0; o[k][k]=1;

  /* optimize later ! */
  loop (i,0,ns) {
    dr[0]=o[0][0]*r[i][0]+o[0][1]*r[i][1]+o[0][2]*r[i][2];
    dr[1]=o[1][0]*r[i][0]+o[1][1]*r[i][1]+o[1][2]*r[i][2];
    dr[2]=o[2][0]*r[i][0]+o[2][1]*r[i][1]+o[2][2]*r[i][2];
    VV(r[i],=dr) }
}

void match_displace(vector *r,int ns,int axis,double d)
{
  int i;

  loop (i,0,ns) r[i][axis]+=d;
}

double cfgdist(vector *r,int ns)
{
  double sum=0;
  int i;

  loop (i,0,ns) sum+=SQRD(match_r[i],r[i]);

  return sum;
}

void match(vector *r,int ns,double eps)
/* rotate and displace molecule to get the best match with match_r[] */
/* consts for minimization: */
#define QGROW 1.1
#define QSHRINK 0.99
{
  double dist,dist0;
  vector *aux;
  static double dr=0.01,da=0.001; /* optimize? */

  alloc(aux,sizeof(aux[0])*ns);
  copy(aux,r,sizeof(aux[0])*ns);

  dist=dist0=cfgdist(r,ns);
  printf("%g --> ",sqrt(dist/ns));

  da=sqrt(da/eps); dr=sqrt(dr/eps);

  do {
    copy(aux,r,sizeof(aux[0])*ns);
    match_rotate(aux,ns,irnd(3),da*rndgauss());
    dist=cfgdist(aux,ns);
    if (dist<dist0) {
      dist0=dist;
      copy(r,aux,sizeof(aux[0])*ns);
      da*=QGROW; }
    else
      da*=QSHRINK;

    copy(aux,r,sizeof(aux[0])*ns);
    match_displace(aux,ns,irnd(3),dr*rndgauss());
    dist=cfgdist(aux,ns);
    if (dist<dist0) {
      dist0=dist;
      copy(r,aux,sizeof(aux[0])*ns);
      dr*=QGROW; }
    else
      dr*=QSHRINK;

    } while (fabs(dr)+fabs(da*10)>eps);

  printf("%g da=%g dr=%g\n",sqrt(dist0/ns),da,dr);
  free(aux);
}

static char *atname(char *name)
/* extracts atom symbolic name as CA in TYR93aCA */
{
  char *e=strend(name)-2;

  if (e<name) return name;
  while (e>=name) {
    if (isdigit(*e) || islower(*e)) return e+1;
    e--; }

  return name;
}

static int essinclude(int mode,int i)
{
  char *c=atname(site[i].id);

  switch (mode) {
    case 'H': /* heavy atoms */
      return *c!='H';
    case 'B': /* backbone */
      return !strcmp(c,"CA") || !strcmp(c,"N") || !strcmp(c,"C");
    case 'C': /* Calpha */
      return !strcmp(c,"CA");
    default:
      ERROR(("bad -E%c option",mode))
    case 'A': /* all atoms */
      return 1; }
}

void essential(species_t *spec) /************************************ analyze */
{
  int ns=spec->ns,i,j,ness,specframe=spec->frame;
  FILE *ess,*mol=NULL;
#ifdef WRITEMATCH
  FILE *essmatch=NULL;
#endif
  unsigned it;
  int n,size;
  double **A,**R=NULL,*d0,*d;
  double rr,x;
  vector *r0,*fp;
  sort_t *sort;
  int *ren;

#define D(X) ((double*)(X[0]))

  if (sizeof(double)*3 != sizeof(vector)) ERROR(("vector bad length"))

  site=spec->site;
  ness=0;
  ralloc(ren,sizeof(int)*ns); /* this is the mark for release */
  loop (i,0,ns)
    if ( (site[i].backbone=essinclude(spec->Xopt.E,i)) ) ren[i]=ness++;
  if (abs(option('r'))!=4) ERROR(("-r4 or -r-4 required"))

  prt("! %d atoms for essential dynamics",ness);

  if (ness!=ns) {
    /* new mol-file called "ess.mol" to view ess####.plb */
    //    mol=fopen("ess.mol","wt");
    strcpy(spec->ext,".ess.mol");
    mol=fopen(spec->fn,"wt");
    if (!mol) ERROR(("cannot write to %s",spec->fn))
    fprintf(mol,"mode=%c for essential dynamics\n\
\n\
number_of_atoms = %d\n\
\n\
atoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",spec->Xopt.E,ness);
    j=0;
    loop (i,0,ns) if (site[i].backbone) {
      int *k;
      int innn,nnn=0,nnb[MAXVAL];

      fprintf(mol,"%3d %8s %4s 0 0 ", j,site[i].id,atom[site[i].type].name);
      loopnbr(k,i) {
        /* bonded atoms will be bonded also in ess.mol */
        if (site[*k].backbone) nnb[nnn++]=ren[*k];
        else if (spec->Xopt.E=='C') {
          int *kk;
          /* Calpha: upto 3rd neighbors will be bonded in ess.mol */
          loopnbr(kk,*k) if (*kk!=i) {
            if (site[*kk].backbone) nnb[nnn++]=ren[*kk];
            else {
              int *kkk;
              loopnbr(kkk,*kk) if (*kkk!=*k) {
                if (site[*kkk].backbone) nnb[nnn++]=ren[*kkk]; } } } } }
      fprintf(mol,"%d",nnn);
      loop (innn,0,nnn) fprintf(mol," %d",nnb[innn]);
      fprintf(mol,"\n");
      j++; }

    fputc('\n',mol);
    fclose(mol);
    prt("! %s written",spec->fn); }

  n=3*ness; size=n*sizeof(double);

#ifdef WRITEMATCH
  if (spec->Xopt.P) {
    float hdr[2];

    strcpy(spec->ext,"essmatch.plb");
    essmatch=fopen(spec->fn,"wb");
    hdr[0]=ness; hdr[1]=0;
    fwrite(hdr,4,2,essmatch); }
#endif

  ralloc(A,n*sizeof(A[0]));
  loop (i,0,n) ralloczero(A[i],size);
  ralloc(R,n*sizeof(R[0]));
  loop (i,0,n) ralloczero(R[i],size);
  ralloczero(r0,size); d0=D(r0);
  ralloc(fp,size*2);
  d=D(fp);

  spec->opt_p=spec->opt_w=0; /* unify! */
  prt("# WARNING: -p0 -w0 applied");

  for (it=0;;it++) {
    if (spec->frame>spec->Xopt.toframe) break;
    if (!read3D(spec,-1)) break;
    spec->frame+=spec->Xopt.byframe;
    j=0;
    loop (i,0,ns) if (site[i].backbone) { VV(fp[j],=site[i].r) j++; }
    match_center(fp,ness);
    if (it)
      match(fp,ness,sqrt(spec->Xopt.Jeps)*0.1);
    else {
      ralloc(match_r,size);
      copy(match_r,fp,size); }

    j=0;
    loop (i,0,ns) if (site[i].backbone) {
      VV(r0[j],+=fp[j])
      j++; }
    loop (i,0,n) loop (j,0,n) A[i][j]+=d[i]*d[j];

#ifdef WRITEMATCH
  /* writing matching cfg (i.e., rotated and translated to match the 1st) */
  if (essmatch) {
    float f[3];

    j=0;
    loop (i,0,ns) if (site[i].backbone) {
      VV(f,=fp[j])
      j++;
      fwrite(f,12,1,essmatch); } }
#endif
    }

#ifdef WRITEMATCH
  if (essmatch) {
    fclose(essmatch);
    fprintf(stderr,"! essmatch.plb written\n"); }
#endif

  prt("! %d frames analyzed",it);

  loop (i,0,n) loop (j,0,n) A[i][j]=(A[i][j]-d0[i]*d0[j]/it)/it;

  /* r0=average (in old version, r0 was replaced by the last frame) */
  j=0;
  loop (i,0,ns) if (site[i].backbone) {
    /*.....  VV(r0[j],=site[i].r)*/
    VVO(site[i].r,=r0[j],/=it)
    j++; }

  {
    /* writing center */
    FILE *plbr0=fopen("essr0.plb","wb");
    /*.....float hdr[2]={ness,0}; - not portable*/
    float hdr[2];

    hdr[0]=ness; hdr[1]=0;
    fwrite(hdr,4,2,plbr0);
    fwrite(r0,12,ness,plbr0);
    fclose(plbr0);
  }

  rr=Jacobi(n,A,R,spec->Xopt.Jeps);
  prt("! Jacobi diagonalization err=%g",rr);

  strcpy(spec->ext,".ess");
  ess=fopen(spec->fn,"wt");

  alloc(sort,n*sizeof(sort[0]));
  if (sizeof(struct sort_s)>2*sizeof(double)) ERROR(("internal"))
  loop (i,0,n) {
    sort[i].l=-A[i][i];
    sort[i].i=i; }

  qsort(sort,n,sizeof(struct sort_s),dblcmp);

  fprintf(ess,"\
# essential eigenvalues mode=%c\n\
# matrix n=%d (ness=%d) %d frames  Jacobi diag.err.=%g\n\
#  i        lambda         lambda/ness    ln(lambda/ness)\n",
    spec->Xopt.E,n,ness,it,rr);
  loop (i,0,n) {
    x=-sort[i].l;
    if (x<=0) break;
    fprintf(ess,"%4d %15.7g %15.7g %15.8f\n",i,x,x/ness,log(x/ness)); }

  fclose(ess);
  prt("! %s written",spec->fn);

  writevibration(spec,R,n,sort,r0,"ess",ness);

  if (R) {
    int from=0,to=n;

    spec->frame=specframe;

    if (spec->Xopt.P>0) to=spec->Xopt.P;
    else from=n+spec->Xopt.P;
    Max(from,0)
    Min(to,n)

    if (from<to) {
      FILE *cp=fopen("ess.cp","wb");
      float *rcp;
      int ncp=max(2,to-from),indx;

      alloczero(rcp,4*ncp+1);
      rcp[0]=CPmark;
      *(int4*)(&rcp[1])=ncp;
      loop (j,2,ncp) sprintf((char*)(rcp+j),"%d",j);
      fwrite(rcp,4,ncp,cp);

      for (it=0;;it++) {
        if (spec->frame>spec->Xopt.toframe) break;
        if (!read3D(spec,-1)) break;
        spec->frame+=spec->Xopt.byframe;
        j=0;
        loop (i,0,ns) if (site[i].backbone) { VV(fp[j],=site[i].r) j++; }
        match_center(fp,ness);
        if (it)
          match(fp,ness,sqrt(spec->Xopt.Jeps)*0.1);
        else
          copy(match_r,fp,size);

        j=0;
        loop (i,0,ns) if (site[i].backbone) {
          VV(fp[j],-=r0[j])
          j++; }
        loop (j,from,to) {
          indx=sort[j].i;
          rr=x=0;
          loop (i,0,n) {
            rr+=Sqr(R[indx][i]);
            x+=d[i]*R[indx][i]; }
          rr*=ness;
          /* this normalization leads to rcp=sqrt(SUM dr^2_i/ness), where dr_i
             are displacements in the direction of the eigenvector R */
          rcp[j-from]=x/sqrt(rr); }
        fwrite(rcp,4,ncp,cp); }

      fclose(cp); } }

  release(ren);
}
