/* cc -O2 -o naclcryst naclcryst.c -lm
 */
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"
typedef int intvector[3];

void swap(int *a,int *b) /********************************************* swap */
{
  int c=*a;
  *a=*b;
  *b=c;
}

int gcd(int a,int b) /************************************************** gcd */
{
  a=abs(a);
  b=abs(b);
  for (;;) {
    if (b==0) return a;
    a=a%b; swap(&a,&b); }
}

int divisor(intvector a) /****************************************** divisor */
/* gcd provided that SUM(a)/gcd is even */
{
  int g=gcd(a[0],a[1]);

  g=gcd(g,a[2]);

  if ((SUM(a)/g)&1) {
    if (g&1) Error("odd vector");
    return g/2; }

  return g;
}

void canonsign(intvector a) /************************************* canonsign */
/* 1) positive a[2]
   2) positive a[1]
   3) positive a[0]
*/
{
  if (a[2]>0) return;
  if (a[2]<0) goto inv;
  if (a[1]>0) return;
  if (a[1]<0) goto inv;
  if (a[0]>0) return;
  if (a[0]<0) goto inv;
  return; /* a=0: cannot happen */
 inv: VO(a,*=-1)
}

struct vec_s {
  struct vec_s *next;
  intvector v;
} *head;

struct found_s {
  struct found_s *next;
  int N;
  intvector a,b;
} *fhead;

int order;
double Na_Cl=2.82; /* see -a */

int eq(intvector a,intvector b) /**************************************** eq */
{
  switch (order) {
    case 0: return abs(a[0])==abs(b[0]) && abs(a[1])==abs(b[1]) && abs(a[2])==abs(b[2]);
    case 1: return abs(a[0])==abs(b[0]) && abs(a[2])==abs(b[1]) && abs(a[1])==abs(b[2]);
    case 2: return abs(a[1])==abs(b[0]) && abs(a[2])==abs(b[1]) && abs(a[0])==abs(b[2]);
    case 3: return abs(a[1])==abs(b[0]) && abs(a[0])==abs(b[1]) && abs(a[2])==abs(b[2]);
    case 4: return abs(a[2])==abs(b[0]) && abs(a[0])==abs(b[1]) && abs(a[1])==abs(b[2]);
    case 5: return abs(a[2])==abs(b[0]) && abs(a[1])==abs(b[1]) && abs(a[0])==abs(b[2]); }
  return -1;
}

int inslab(intvector a,intvector r) /******************************** inslab */
{
  double scal=SCAL(a,r);
  return scal>=0 && scal<SQR(a);
}

int main(int narg,char **arg) /**************************************** main */
{
  intvector c,b,i;
  int iarg,rg,n,Nmax=0,Nmin=0,sg,raw=0,ncut=0,nn,mn,verbose=1,ordered,canonical=0;
  struct vec_s *v,*vv;
  struct found_s *f;
  double maxaspect=2,lasp;
  vector O[3],CM;
  typedef float floatvector[3];
  floatvector *R;
  int *mark;
  double d;
  char *FORMAT=NULL;

  if (narg<4) {
    fprintf(stderr,"\
NaCl crystal matching periodic b.c. with given (cx,cy,cz) plane. Call by:\n\
  naclcryst cx cy cz [OPTIONS]\n\
WHERE c=[cx,cy,cz] is a Miller index or its integer multiple (the pattern\n\
  repeats in the z-direction), cx+cy+cz must be even (periodic Na matches Na)\n\
OPTIONS:\n\
  -a#   lattice constant (Na-Cl distance) [default=2.8203]\n\
  -x#   the aspect ratio of the cell will be in [1/#,#] [default=2]\n\
  none of -n,-m: a,b will satisfy |a|<=|b|<=|c| [default]\n\
  -n#   number of atoms in the output cell, with -m the minimum number of atoms\n\
  -m#   maximum number of atoms\n\
  -c#   cut to max # atoms closest to the center of mass, # must be even\n\
  -o[FORMAT] output as plb,mol files, FORMAT must accept integer [nacl-%%04d]\n\
  -r    raw (do not rotate to the final box)\n\
  -s    sorted to the canonical order, exclusive with -n,-m\n\
  -q    quiet\n\
ALGORITHM:\n\
  1. find all (upto max.size) integer vectors b perpendicular to c=(cx,cy,cz)\n\
  2. remove b=0 and odd bx+by+bz\n\
  3. of these, find vector a perpendicular to b\n\
  4. throw away those not satisfying the number of atoms or aspect\n\
  5. write plb,mol files\n\
Example:\n\
  naclcryst 3 4 5 -x1.5 -m600\n\
See also:\n\
  lattice ice naclcryst.sh\n");
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') {
      double d=atof(arg[iarg]+2);
      switch (arg[iarg][1]) {
        case 'n': Nmin=d; break;
        case 'm': Nmax=d; break;
        case 'c': ncut=d; break;
        case 'r': raw=1; break;
        case 'a': Na_Cl=d; break;
        case 'x': maxaspect=d; break;
        case 's': canonical=1; break;
        case 'q': verbose=0; break;
        case 'o': FORMAT=arg[iarg]+2; break; } }
    else
      c[0]=c[1],c[1]=c[2],c[2]=atoi(arg[iarg]);

  if (canonical) {
    c[0]=abs(c[0]);
    c[1]=abs(c[1]);
    c[2]=abs(c[2]);
    if (c[0]>c[1]) swap(c,c+1);
    if (c[1]>c[2]) swap(c+1,c+2);
    if (c[0]>c[1]) swap(c,c+1); }

  if (SUM(c)&1) Error("vector c=(cx,cy,cz) is odd, cx+cy+cz should be even");

  if (!FORMAT || !FORMAT[0]) FORMAT="nacl-%04d";
  if (!strchr(FORMAT,'%')) Error("-oFORMAT does not contain format");

  if (maxaspect<1) maxaspect=1/maxaspect;
  lasp=log(maxaspect)*2;

  if (Nmax==0) Nmax=Nmin;

  ordered=Nmax==0 && Nmin==0;
  if (ordered) {
    rg=sqrt(SQR(c))+0.99999;
    Nmax=SQR(c)*(sqrt(SQR(c))+0.99999); }
  else {
    if (canonical) Error("-s cannot be used with -n,-m");
    rg=sqrt(Nmax/sqrt(SQR(c)))*maxaspect; }


  if (verbose) printf("N in [%d,%d]  search range=%d  maxaspect=%g\n",Nmin,Nmax,rg,maxaspect);

  n=0;
  loopto (b[0],-rg,rg)
    loopto (b[1],-rg,rg)
      loopto (b[2],-rg,rg) if (SCAL(c,b)==0) if ((SUM(b)&1)==0) if (SQR(b)) {
        allocone(v);
        n++;
        v->next=head; head=v;
        VV(v->v,=b) }

  if (verbose) printf("%d vectors perpendicular to c generated in given range\n",n);

  n=0;
  looplist (v,head)
    looplist (vv,v->next)
      /* cumbersome - better by vector product ? */
      if (SCAL(v->v,vv->v)==0) {
        int N;
        long unsigned NN=(long unsigned)SQR(c)*(long unsigned)SQR(v->v)*(long unsigned)SQR(vv->v);

        N=(int)sqrt(NN+0.5);
        if ((long)N*N!=NN) {
          fprintf(stderr,"N=%d NN=%ld: N^2!=NN - algorithm error or overflow",N,NN);
          exit(1); }

        if (N<Nmin || N>Nmax) goto excl;

        if (lasp<40) {
          /* check aspect ratio */
          if (fabs(log((double)SQR(v->v)/SQR(vv->v)))>lasp) goto excl;
          if (fabs(log((double)SQR(c)   /SQR(vv->v)))>lasp) goto excl;
          if (fabs(log((double)SQR(c)   /SQR(v->v)))>lasp) goto excl; }

        if (ordered) {
          if (SQR(c)<SQR(v->v)) goto excl;
          if (SQR(v->v)<SQR(vv->v)) goto excl; }

        /* skip equivalent boxes */
        looplist (f,fhead) loop (order,0,6) {
          if (eq(v->v,f->b) && eq(vv->v,f->a)) goto excl;
          if (eq(v->v,f->a) && eq(vv->v,f->b)) goto excl;
          if (eq(c,f->b) && eq(v->v,f->a) && eq(vv->v,c)) goto excl;
          if (eq(c,f->a) && eq(v->v,f->b) && eq(vv->v,c)) goto excl;
          if (eq(c,f->a) && eq(v->v,c) && eq(vv->v,f->b)) goto excl;
          if (eq(c,f->b) && eq(v->v,c) && eq(vv->v,f->a)) goto excl; }

        /* append to final list */
        allocone(f);
        f->next=fhead; fhead=f;
        VV(f->b,=v->v)
        VV(f->a,=vv->v)
        f->N=N;
        n++;

        excl:; }

  if (verbose) printf("%d crystal%s found\n",n,"s"+(n==1));
  if (n>9999) Error("More than 9999 crystals to write, decrease -m and/or -x");
  else if (n==0 && verbose) Error("No crystal matches criteria.\n\
      Increase -m, multiply [cx,cy,cz] by an interger, or increase -x");

  rg=rg*3+2;

  looplist (f,fhead) {
    static char fn[16];
    static int nr=0;
    FILE *plb,*mol;
    float r[3];
    vector L,asp;
    int k;
    double masp;

    if (canonical) {
      /* c[0,1,2] is sorted, |a|,|b|,|c| is sorted */

      canonsign(f->a);
      canonsign(f->b);

      if (c[0]==c[1]) {
        if (f->b[0]>f->b[1]) { swap(f->a,f->a+1); swap(f->b,f->b+1); } }
      if (c[1]==c[2]) {
        if (f->b[1]>f->b[2]) { swap(f->a+1,f->a+2); swap(f->b+1,f->b+2); } }
      if (c[0]==c[1]) {
        if (f->b[0]>f->b[1]) { swap(f->a,f->a+1); swap(f->b,f->b+1); } }

      /* probably enough for most practical cases */ }

    if (FORMAT) {
      char *end;
      sprintf(fn,FORMAT,nr);
      end=strend(fn);
      if (verbose) printf("===== writing %s.{plb,mol}. Show by: show -I%% %s\n",fn,fn);
      strcpy(end,".plb"); plb=fopen(fn,"wb");
      strcpy(end,".mol"); mol=fopen(fn,"wt"); }
    nr++;

    r[0]=ncut?ncut:f->N;
    r[1]=-3;
    if (FORMAT) fwrite(r,4,2,plb);
    if (raw) VO(r,=0)
    else {
      r[0]=L[0]=Na_Cl*sqrt(SQR(f->a));
      r[1]=L[1]=Na_Cl*sqrt(SQR(f->b));
      r[2]=L[2]=Na_Cl*sqrt(SQR(c)); }

    if (FORMAT) {
      fwrite(r,4,3,plb);
      fprintf(mol,"\
nacl crystal\n\n\
parameter_set = sea\n\
number_of_atoms = %d\n\n\
atoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",ncut?ncut:f->N); }

    masp=1;
    asp[0]=sqrt((double)SQR(c)/SQR(f->b));   Max(masp,asp[0]) Max(masp,1/asp[0])
    asp[1]=sqrt((double)SQR(f->b)/SQR(f->a));Max(masp,asp[1]) Max(masp,1/asp[1])
    asp[2]=sqrt((double)SQR(f->a)/SQR(c));   Max(masp,asp[2]) Max(masp,1/asp[2])

    printf("%-5d maxasp=%.4f  a=[%d,%d,%d] b=[%d,%d,%d] c=[%d,%d,%d] r=%d%s",
           f->N, masp,
           f->a[0],f->a[1],f->a[2],
           f->b[0],f->b[1],f->b[2],
           c[0],c[1],c[2],
           divisor(f->a)*divisor(f->b)*divisor(c),
           canonical?" C":"");

    if (ncut) printf(" -> N=%d",ncut);
    printf("\n");

    /* rotation matrix */
    VV(O[0],=f->a) d=sqrt(SQR(O[0])); VO(O[0],/=d)
    VV(O[1],=f->b) d=sqrt(SQR(O[1])); VO(O[1],/=d)
    VV(O[2],=c)    d=sqrt(SQR(O[2])); VO(O[2],/=d)

    if (verbose) {
      printf("L = %9.6f %9.6f %9.6f\n", VARG(L));
      printf("|a| = %9.6f  |b| = %9.6f  |c| = %9.6f\n", sqrt(SQR(f->a)), sqrt(SQR(f->b)), sqrt(SQR(c)));
      printf("c/a=%.4f a/b=%.4f b/c=%.4f\n\n", VARG(asp)); }

    n=0;
    loop (sg,0,2)
      loopto (i[0],-rg,rg)
        loopto (i[1],-rg,rg)
          loopto (i[2],-rg,rg) {
            if (inslab(c,i) && inslab(f->b,i) && inslab(f->a,i) && sg+SUM(i)&1) {
              n++; } }

    nn=n;
    allocarray(R,nn);
    allocarray(mark,nn);
    loop (n,0,nn) mark[n]=1;

    n=0;
    VO(CM,=0)
    loop (sg,0,2)
      loopto (i[0],-rg,rg)
        loopto (i[1],-rg,rg)
          loopto (i[2],-rg,rg) {
            if (inslab(c,i) && inslab(f->b,i) && inslab(f->a,i) && sg+SUM(i)&1) {
              if (raw) VV(r,=Na_Cl*i)
              else loop (k,0,3) r[k]=Na_Cl*SCAL(O[k],i);
              VV(CM,+=r)
              VV(R[n],=r)
              n++; } }

    if (n!=nn) Error("");
    VO(CM,/=nn)

    if (ncut) {
      do {
        n=0;
        loop (sg,0,2) {
          double mrr=0,rr;

          mn=0;
          loopto (i[0],-rg,rg)
            loopto (i[1],-rg,rg)
              loopto (i[2],-rg,rg) {
                if (inslab(c,i) && inslab(f->b,i) && inslab(f->a,i) && sg+SUM(i)&1) {
                  if (mark[n]) {
                    rr=SQRD(R[n],CM)+n*1e-9;
                    if (rr>mrr) mrr=rr,mn=n; }
                  n++; } }
          mark[mn]=0; }
        mn=0; loop (n,0,nn) mn+=mark[n]; }
      while (mn>ncut);
      if (mn!=ncut) Error("parity (-c must be even"); }

    n=mn=0;
    loop (sg,0,2)
      loopto (i[0],-rg,rg)
        loopto (i[1],-rg,rg)
          loopto (i[2],-rg,rg) {
            if (inslab(c,i) && inslab(f->b,i) && inslab(f->a,i) && sg+SUM(i)&1) {
              if (FORMAT && mark[n]) {
                fwrite(R[n],4,3,plb);
                fprintf(mol,"%d %s/%s%d %s %d 0 0 0\n",
                        mn,
                        sg?"G170":"M090",
                        sg?"Cl":"Na",n,sg?"CL":"NA",1-2*sg);
                mn++; }
              n++; } }

    if (FORMAT) {
      fclose(plb);
      fclose(mol); }
    if (n!=f->N) Error("N does not match");
    free(mark);
    free(R); }

  return 0;
}
