/* cc -Wall -O2 -o densprof densprof.c -lm
 */
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"

inline double sqr(double x)
{
  return x*x;
}

typedef float fvector[3];

int OFFSET,NS=1,NM;
int OFFSETA,NSA,NMA;
char *key;
char KEY;         // =key[0]
#define NI 1024   // max. # of atom indices in the key or atoms in a molecule
int indx[NI];
double weight[NI];

double Lz=100;
double RANGE;     // 4*Lz dangling H
double GRID=0;    // z-histogram grid (per unity)
double DZ=0;      // z-histogram grid (delta z); DZ=1/GRID
int NH;           // ((int)(GRID*RANGE)) // # of z-histogram bins
int COSGRID=0;    // cos(theta)-histogram grid (per unity)

struct hist_s {
  int n;
  double P1,P2,P3,P4;
  int *coshist; // [2*COSGRID]
} *hist; // [NH]

double *smoothmask;

double smooth(int iz,int indx) /************************************* smooth */
{
  double x,sum=0;
  int i;

  loop (i,0,NH) {
    if (indx) x=(&(hist[(iz+i)%NH].P1))[indx-1];
    else x=hist[(iz+i)%NH].n;
    sum+=x*smoothmask[i]; }
  return sum;
}

int getindx(char *arg) /******************************************** getindx */
{
  int i=0;
  char *c;

  for (c=arg; *c; ) {
    c++;
    if (!*c) break;
    if (i>=NI) Error("to many indices in arg");
    indx[i]=(int)(weight[i]=atof(c));
    i++;
    c=strchr(c,',');
    if (!c) break; }

  return i;
}

double DIST(fvector a, fvector b,fvector L) /************************** DIST */
{
  int i;
  double rr=0;

  loop (i,0,3)
    rr+=sqr(L[i]/2-fabs(L[i]/2-fabs(a[i]-b[i])));

  return sqrt(rr);
}

int main(int narg,char **arg) /**************************************** main */
{
  FILE *plb,*out;
  char name[256],*ch;
  int iarg,i,nf=0,ns,from=1,to=0x7fffffff,n,iz;
  float hdr[2];
  fvector L,*R,*r;
  double Aav=0,Fav=0,PHASE=0;
  int nsites,zerohist=0;
  int PER=0,VERBOSE=1;
  int GAUSS=0,RECT=0,SMOOTH,prtmask=0;
  double ROHW=2.45;
  double ROH2=2.38;
  int F=0,k,simpleper=0,smoothF=0;
  double *cf,*sf,*ampl;

  if (narg<3) {
    fprintf(stderr,"\
Selected density profiles and angular correlations from trajectory. Call by:\n\
  densprof [OPTIONS] NAME [OPTIONS]\n\
FILES:\n\
  NAME.plb             input plb-file (new format with L required)\n\
  NAME+OFFSET-KEY.lp   output Pn(cos(theta)) distributions (Pn=Legendre polyn.)\n\
  NAME+OFFSET-KEY.z=#.##.cos  \n\
                       output cos(theta) distributions (if option -c)\n\
OPTIONS:\n\
  -cGRID   slab-based histograms of cos(theta), GRID per unity [default=%d]\n\
  -dDZ     z-grid, in histogram bin width (no -p; only one of -g, -d)\n\
  -fFROM   first frame in NAME.plb analyzed [default=1]\n\
  -FNF     make Fourier analysis of every z-profile, use with -kq\n\
  -F-NF    as above, smooth the coefficients (does not work very well :()\n\
  -gGRID   z-grid, in histogram bins per 1A [default=1] or total (if -p)\n\
  -hPHASE  with -p, add phase PHASE in [0,1] to z/Lz\n\
  -kKEY    defines the vector for analysis of its cos(theta), KEY:\n\
    aA1,A2   vector from atom A1 to atom A2\n\
    dO,H1,H2 axis (dipole) vector of water or similar (O-Hx, Hx=(H1+H2)/2)\n\
    hO,H1[,A1[,A2..]]  vector O-H1, where (water) hydrogen H1 is either\n\
             dangling or single, double or triple bonded to any of acceptors:\n\
             {O on other (water) molecule, A1,A2.. on other molecule}\n\
             H-bond criterion uses distance H1-O (option -r) or H1-A (opt. -R)\n\
             For acceptors on other molecules, use -M,-N,-O,-R\n\
             BUG: should be run twice for both H1 and H2 in water and summed\n\
    iN1,C2,N3,C4,C5  vector perpendicular to a 5-member ring\n\
             (vector product of N1-N2 and C2-Cx, where Cx=(C4+C5)/2)\n\
             e.g., imidazolium ring in standard numbering\n\
    qQ0,Q1,..,QNS  weighed (e.g., charge) density profile\n\
             Q0 is the weight (e.g., charge in e) of site 0, etc.\n\
             BUG: for a mixture, density profiles must be summed over species\n\
    wO,H1,H2 vector perpendicular to a 3-member group (e.g., water H1-O-H2)\n\
             (vector product of H1-H2 and O-Hx, where Hx=(H1+H2)/2)\n\
    * IDs A1,A2,O,... are atom numbers taken from the mol-file\n\
    * lowercase key: geometric center of given group (donor H for key -kh)\n\
      of atoms is calculated and its z-coordinate defines the histogram\n\
      and the z-profile\n\
    * UPPERCASE key: first atom is the z-profile basis (O for key -kh)\n\
      (note that cook* uses either sites or molecule center-of-mass)\n\
  -lLENGTH number of frames analyzed (use after -f, excl. with -t)\n\
  -mM      number of molecules of given kind in a configuration (frame)\n\
  -nNS     number of sites in the molecule measured [1]\n\
  -oOFFSET where given block of M molecules start (in # of sites) [0]\n\
  -M,-N,-O as above for the second acceptor molecule (not water)\n\
  -p       PERIODIC MODE: split z/Lz in [0,1) into GRID bins, periodic b.c.\n\
  -pRECT, -pRECT,GAUSS (or -pRECT+GAUSS):\n\
           As above and smooth periodically by a (2*RECT+1)-point formula\n\
           (rectangle) and GAUSS-times a 3-point formula.\n\
  -PMASK   with -p, sum of: 1=calculate min/max at center and averaged in area\n\
                              farther than 2*(RECT+GAUSS) from the interfaces\n\
                            2=print the convolution window used for -p\n\
  -q       quiet (no output to stdout)\n\
  -rROHW   distance H-O defining an H-bond for dangling hydrogen analysis (-kh)\n\
           normally the minimum on the H-O RDF [%g]\n\
  -RROH2   H-A distance for the second acceptor molecule (not water) [%g]\n\
  -tTO     last frame in NAME.plb analyzed, excl. with -l [default=EOF]\n\
  -vVERBOSE print also zero values of density profiles (good in scripts)\n\
  -zRANGE  z-range (Lz) [%g]\n\
  -Z       periodic b.c. (on input, no tranformation to [0,1])\n\
Examples (BMIM):\n\
  densprof T360 -o0 -n25 -m200 -ka10,14\n\
  densprof T360 -o0 -n25 -m200 -ki10,1,9,21,22\n\
Example (NaCl):\n\
  densprof nacl1000 -o0   -n1 -m500 -kq1\n\
  densprof nacl1000 -o500 -n1 -m500 -kq-1\n\
  mergetab nacl1000+0-q.z:1:3 nacl1000+500-q.z:1=1:3 | plot -:A:C+B\n\
See also:\n\
  plb2rdf plb2diff plbmsd plb2hist plb2cryst plb2dens plb2sqd\n\
  concprof plbqprof\n\
  plbconv plbcut plbfilt plb2asc asc2plb\n",
            COSGRID,ROHW,ROH2,Lz);
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'c': COSGRID=atoi(arg[iarg]+2); break;
        case 'd': DZ=atof(arg[iarg]+2); break;
        case 'f': from=atoi(arg[iarg]+2); break;
        case 'F': F=atoi(arg[iarg]+2); if (F<0) F=-F,smoothF++; break;
        case 'l': to=from+atoi(arg[iarg]+2)-1; break;
        case 'g': GRID=atof(arg[iarg]+2); break;
        case 'h': PHASE=atof(arg[iarg]+2); break;
        case 'k': key=arg[iarg]+2; break;
        case 'm': NM=atoi(arg[iarg]+2); break;
        case 'n': NS=atoi(arg[iarg]+2); break;
        case 'o': OFFSET=atoi(arg[iarg]+2); break;
        case 'M': NMA=atoi(arg[iarg]+2); break;
        case 'N': NSA=atoi(arg[iarg]+2); break;
        case 'O': OFFSETA=atoi(arg[iarg]+2); break;
        case 'P': prtmask=atoi(arg[iarg]+2); break;
        case 'p':
          PER=1;
          ch=strchr(arg[iarg]+1,',');
          /* both -pRECT,GAUSS and -pRECT+GAUSS accepted,
             to be compatible with utility smooth */
          if (!ch) ch=strchr(arg[iarg]+1,'+');
          RECT=atoi(arg[iarg]+2);
          if (ch) GAUSS=atoi(ch+1);
          break;
        case 'q': VERBOSE=0; break;
        case 'r': ROHW=atof(arg[iarg]+2); break;
        case 'R': ROH2=atof(arg[iarg]+2); break;
        case 't': to=atoi(arg[iarg]+2); break;
        case 'v': zerohist++; break;
        case 'z': Lz=atof(arg[iarg]+2); break;
        case 'Z': simpleper=1; break;
        default: Error("unknown option"); }
    else
      strcpy(name,arg[iarg]);

  if (NS<1) Error("NS<1 (no site, see option -n)");
  if (NM<1) Error("NM<1 (no molecule, see option -m)");
  if (from<1) Error("FROM<1");
  if (to<from) Error("empty frame range");
  if (!key) Error("no -k option (consider -kq1 for the simplest z-profile)");
  PHASE+=4;

  SMOOTH=RECT+GAUSS;

  if (COSGRID && SMOOTH)
    fprintf(stderr,"WARNING: smothing ignored for cos(theta)\n");

  if (DZ!=0) {
    if (GRID!=0) Error("both -g and -d specified");
    GRID=1/DZ; }
  if (!GRID) GRID=1;

  if (!(plb=fopen(string("%s.plb",name),"rb"))) Error(name);
  nsites=getindx(key);
  KEY=key[0];

  fread(hdr,4,2,plb);
  ns=hdr[0];
  if (hdr[1]!=-3) Error("only new (L3) format of plb-files is accepted");
  fseek(plb,(long)(ns+1)*(from-1)*sizeof(fvector),SEEK_CUR);
  if (OFFSET+NM*NS>ns)
    Error("selected sites do not exist in the configuration\n\
*** check options -o,-n,-m");
  if (OFFSETA+NMA*NSA>ns)
    Error("acceptor: selected sites do not exist in the configuration\n\
*** check options -o,-n,-m");

  allocarray(R,ns*3);

  switch (tolower(KEY)) {
    case 'q': if (nsites>ns) Error("KEY=Q,q: more weights than sites in .plb"); break;
    case 'i': if (nsites!=5) Error("KEY=I,i: # of sites must be 5"); break;
    case 'a': if (nsites!=2) Error("KEY=A,a: # of sites must be 2"); break;
    case 'd': if (nsites!=3) Error("KEY=D,d: # of sites must be 3"); break;
    case 'w': if (nsites!=3) Error("KEY=W,w: # of sites must be 3"); break;
    case 'h': if (nsites<2) Error("KEY=H,h: # of sites must be >=2"); break;
    default: Error("bad KEY"); }

  if (KEY=='q') {
    out=fopen(string("%s+%d-q.z",name,OFFSET,KEY),"wt");
    fprintf(out,"# [dens.prof.] = 1/AA^3 (number density at given z - Lx,Ly assumed constant)\n\
#   z  [A]    dens.prof. [A-3]  charge dens.prof. in [e.A-3]\n"); }
  else {
    out=fopen(string("%s-%c.lp",name,KEY),"wt");
    fprintf(out,"# [dens.prof.] = 1/AA^3 (number density at given z - Lx,Ly assumed constant)\n\
# Pi = Legendre polynomial\n\
#   z      dens.prof. <P1(cos th)> <P2(cos th)> <P3(cos th)> <P4(cos th)>\n"); }
  fprintf(out,"# args:");
  loop (iarg,1,narg)  fprintf(out," %s",arg[iarg]);
  fprintf(out,"\n");

  RANGE=4*Lz;
  if (F) {
    if (KEY!='q') Error("-F requires -kq");
    allocarrayzero(sf,F+1);
    allocarrayzero(cf,F+1);
    allocarrayzero(ampl,F+1); }

  if (PER) {
    double l,ll,end;

    NH=GRID;
    allocarrayzero(smoothmask,NH);

    if (RECT>=NH/2) Error("rectangular smoothing window RECT too long");
    loopto (iz,NH-RECT,NH+RECT) smoothmask[iz%NH]=1./(1+2*RECT);

    loop (iz,0,GAUSS) {
      l=smoothmask[NH-1];
      end=(smoothmask[NH-2]+l+smoothmask[0])/3;
      loop (i,0,NH-1) {
        ll=l; l=smoothmask[i];
        smoothmask[i]=(smoothmask[i+1]+l+ll)/3; }
      smoothmask[NH-1]=end; } }
  else
    NH=(int)(GRID*RANGE);

  allocarrayzero(hist,NH);
  if (COSGRID) loop (i,0,NH) allocarrayzero(hist[i].coshist,2*COSGRID);

  for (; from<=to; from++) {
    if (fread(L,sizeof(fvector),1,plb)!=1) goto end;
    if (fread(R,sizeof(fvector),ns,plb)!=ns) Error("unexpected EOF");
    if (simpleper)
      loop (n,0,ns) {
        while (R[n][2]<0) R[n][2]+=L[2];
        while (R[n][2]>=L[2]) R[n][2]-=L[2]; }

    if (PER) Aav+=L[0]*L[1]*L[2]/GRID;
    else Aav+=L[0]*L[1];

    if (F) {
      if (PER) Fav+=Sqr(L[0]*L[1]*L[2]/GRID);
      else Fav+=Sqr(L[0]*L[1]*L[2]); }

    loop (n,0,NM) {
      static vector a,b,c,gc;
      double cc;
      int j;

      r=R+OFFSET+n*NS;

#define GC(N) \
  if (isupper(KEY)) { \
    VV(gc,=r[indx[0]]) } \
  else { \
    VO(gc,=0) \
    loop (j,0,N) VV(gc,+=(1./N)*r[indx[j]]) }

      switch (tolower(KEY)) {
        case 'a': /* vector connecting 2 atoms */
          VVV(c,=r[indx[1]],-r[indx[0]])
          cc=sqrt(SQR(c));
          VO(c,/=cc)

          GC(2)
          break;

        case 'd': /* vector of symmetry axis of water */
          VVV(b,=0.5*r[indx[2]],+0.5*r[indx[1]])
          VVV(c,=b,-r[indx[0]])
          cc=sqrt(SQR(c));
          VO(c,/=cc)

          GC(3)
          break;

        case 'h': { /* dangling, H-bonded, double and triple bonded hydrogens */
          int nn;
          int nHB=0;

          /* H1-water O bonds */
          loop (nn,0,NM) if (nn!=n) {
            fvector *rr=R+OFFSET+nn*NS;

            if (DIST(r[indx[1]],rr[indx[0]],L)<ROHW) nHB++; }

          /* H1 bonded to the second molecule */
          loop (nn,0,NMA) {
            fvector *rr=R+OFFSETA+nn*NSA;
            int ii;

            loop (ii,2,nsites)
              if (DIST(r[indx[1]],rr[indx[ii]],L)<ROH2) nHB++; }

          VVV(c,=r[indx[1]],-r[indx[0]]) /* H1-O */
          cc=sqrt(SQR(c));
          VO(c,/=cc)
          if (islower(KEY)) VV(gc,=r[indx[1]])
          else VV(gc,=r[indx[0]])
          gc[2]+=Lz*nHB; }
          break;

        case 'i': /* vector perpendicular to a 5-member (imidazolium) ring */
          VVV(a,=r[indx[2]],-r[indx[0]])
          VVV(c,=0.5*r[indx[4]],+0.5*r[indx[3]])
          VVV(b,=c,-r[indx[1]])
          VECT(c,a,b)
          cc=sqrt(SQR(c));
          VO(c,/=cc)

          GC(5)
          break;

        case 'w': /* vector perpendicular to a plane of 3 atoms (water) */
          VVV(a,=r[indx[2]],-r[indx[1]])
          VVV(c,=0.5*r[indx[2]],+0.5*r[indx[1]])
          VVV(b,=c,-r[indx[0]])
          VECT(c,a,b)
          cc=sqrt(SQR(c));
          VO(c,/=cc)

          GC(3)
          break;

        default:; }

      if (KEY=='q') { /* weighed (e.g., charge) density profile */
        loop (j,0,nsites) {
          if (PER)
            iz=(int)((r[j][2]/L[2]+PHASE)*NH)%NH;
          else
            iz=(int)(r[j][2]*GRID);
          if (F) loopto (k,0,F) {
            sf[k]+=weight[j]*sin(r[j][2]*(2*PI/L[2])*k);
            cf[k]+=weight[j]*cos(r[j][2]*(2*PI/L[2])*k); }
          if (iz>=0 && iz<=NH) {
            hist[iz].n++;
            hist[iz].P1+=weight[j]; } } }
      else {
        if (PER)
          iz=(int)((gc[2]/L[2]+PER)*NH)%NH;
        else
          iz=(int)(gc[2]*GRID);
        if (iz>=0 && iz<=NH) {
          hist[iz].n++;
          hist[iz].P1+=c[2];
          hist[iz].P2+=1.5*Sqr(c[2])-0.5;
          hist[iz].P3+=2.5*Cub(c[2])-1.5*c[2];
          hist[iz].P4+=4.375*Pow4(c[2])-3.75*Sqr(c[2])+0.375;
          if (COSGRID) {
            int ic=(c[2]+1)*COSGRID;

            Max(ic,0) Min(ic,2*COSGRID-1)
            hist[iz].coshist[ic]++; } } } }

    if (F) {
      loopto (k,0,F) {
        ampl[k]+=Sqr(sf[k])+Sqr(cf[k]);
        sf[k]=cf[k]=0; } }

    nf++; }
 end:

  if (F) {
    double w=Fav;

    loopto (k,0,F) {
      ampl[k]=sqrt(ampl[k]/w)*(1+(k>0));
      if (smoothF) ampl[k]*=Sqr(1-(double)k/(1+F));
      fprintf(out,"# %d %12.7f\n",k,ampl[k]); } }

  if (KEY=='q') {
    double m=9e9,M=-9e9,h,x,xx,z;
    int up,down,bad;

    loop (iz,0,NH) if (hist[iz].n || zerohist) {
      double f=0,Lz;

      if (PER) Lz=1; else Lz=L[2]; /* bug: what if change? */

      z=(iz+.5)/GRID;

      if (F) {
        loop (k,0,F) f+=ampl[k]*cos(z*(2*PI/Lz)*k); }

      if (PER)
        fprintf(out,"%8.5f %12.7f %12.7f %g \n",
                z,
                x=smooth(iz,0)/Aav,
                smooth(iz,1)/Aav,
                f);
      else
        fprintf(out,"%8.2f %12.7f %12.7f %g\n",
                z,
                x=hist[iz].n*GRID/Aav,
                hist[iz].P1*GRID/Aav,
                f);
      Min(m,x)
      Max(M,x) }

    /* analyze conv. profile shape - PER only*/
    if (PER && prtmask&1) {
      x=smooth(NH-1,0)/Aav;
      up=-1; down=-1; bad=0;
      h=(m+M)/2;
      loop (iz,0,NH) {
        xx=x;
        x=smooth(iz,0)/Aav;
        if (xx<h && x>=h) { if (up>=0) bad++; up=iz; }
        if (xx>=h && x<h) { if (down>=0) bad++; down=iz; } }

      loop (i,0,2) {
        double sum=0,c;
        int n,j,from,to;

        if (i && down<up || !i && down>=up) iz=(up+down+NH)/2;
        else iz=(up+down)/2;

        c=(smooth(iz,0)+smooth((iz+1)%NH,0))/(2*Aav);
        // (smooth(iz,1)+smooth((iz+1)%NH,1))/(2*Aav); 3rd column

        n=0;
        if (i) { from=up; to=down; }
        else { from=down; to=up; }
        to-=NH;
        while (to<from) to+=NH;
        loopto (j,from+SMOOTH*2,to-SMOOTH*2) {
          n++;
          sum+=smooth(j%NH,0); }
        sum/=(2*n*Aav);

        z=(iz+1)/GRID;
        fprintf(out,"\n%8.5f 0 0\n\
%8.5f %12.7f %12.7f X %s (center, %d points), bad=%d\n",
                z,z,
                c,sum,
                i?"max":"min",n,bad); }
      if (bad) 
        fprintf(stderr,"WARNING: cannot unambiguously determine interface (%d)\n",bad);
      if (prtmask&2) {
        fprintf(out,"\n");
        loop (iz,0,NH)
          fprintf(out,"%8.5f %8.5f X mask\n",
                  iz/GRID,smoothmask[(iz+(NH/2))%NH]); }
    }
  }
  else loop (iz,0,NH) if (hist[iz].n || zerohist) {
    FILE *f;
    double z=(iz+.5)/GRID;
    int sum=0;
    int ic;

    if (PER)
      fprintf(out,"%8.5f %12.7f %12.7f %12.7f %12.7f %12.7f\n",
              z,
              smooth(iz,0)/Aav,
              smooth(iz,1)/Aav,
              smooth(iz,2)/Aav,
              smooth(iz,3)/Aav,
              smooth(iz,4)/Aav);
    else
      fprintf(out,"%8.2f %12.7f %12.7f %12.7f %12.7f %12.7f\n",
              z,
              (double)(hist[iz].n)*GRID/Aav,
              hist[iz].P1/hist[iz].n,
              hist[iz].P2/hist[iz].n,
              hist[iz].P3/hist[iz].n,
              hist[iz].P4/hist[iz].n);

    if (COSGRID) {
      loop (ic,0,2*COSGRID) sum+=hist[iz].coshist[ic];
      f=fopen(string("%s-%c.z=%.2f.cos",name,KEY,z),"wt");

      fprintf(f,"# %s %s %s %s %s z=%.2f\n\
# cos(theta)  prob\n",arg[1],arg[2],arg[3],arg[4],arg[5],z);
      loop (ic,0,2*COSGRID)
        fprintf(f,"%7.2f %g\n",
                (ic+0.5)/COSGRID-1,
                hist[iz].coshist[ic]/(sum/(2.0*COSGRID)));
      fprintf(f,"# %d frames analyzed\n",nf);
      fclose(f); }
  }

  fprintf(out,"# %d frames analyzed\n",nf);
  if (VERBOSE) printf("%d frames analyzed, Aav=%g\n",nf,Aav/nf);
  fclose(out);

  fclose(plb);

  return 0;
}
