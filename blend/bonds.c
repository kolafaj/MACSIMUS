/* cc -Wall -O2 -o bonds bonds.c -lm
*/
#include "../gen/include.h"

#define MAXVAL 8

#define NMAX 16777216
int N=1048576;

struct cfg_s {
  float r[3];
  float q;
  int at;
  int nnbr;
  int nbr[MAXVAL]; 
  char id[8];
} *cfg;

int nval[MAXVAL+1];

int periodic=0;
float L[3];
double Lh[3];

float dist(float *a,float *b)
{
  if (periodic) {
    double sum=0,dr;
    int k;

    loop (k,0,3) {
      dr=a[k]-b[k];
      while (dr>Lh[k]) dr-=L[k];
      while (dr<-Lh[k]) dr+=L[k];
      sum+=Sqr(dr); }
    return sqrt(sum); }
 else
   return sqrt(Sqr(a[0]-b[0])+Sqr(a[1]-b[1])+Sqr(a[2]-b[2]));
}

static char *ats="?NCOSHPF";

int attype(char *c)
{
  char *pos=strchr(ats,*c);
  int t;

  if (!pos) t=0;
  else t=pos-ats;

  return t;
}

#define NAT 8

/* max bond lengths from charmm21; bit less for ?; default=*1.2 */
static float bond[NAT][NAT]={
  /*    ?     N      C     O     S      H     P    F */
  /*?*/{3.1, -1,    -1,   -1,    -1,   -1,   -1,  -1 },                        
  /*N*/{2.1, 1.47,  -1,   -1,    -1,   -1,   -1,  -1 },
  /*C*/{2.1, 1.51,  1.54, -1,    -1,   -1,   -1,  -1 },
  /*O*/{2.2, 1.485, 1.35, 1.415, -1,   -1,   -1,  -1 },
  /*S*/{2.3, 1.74,  1.82, 1.57,  2.04, -1,   -1,  -1 },
  /*H*/{1.5, 1.04,  1.09, 0.966, 1.35, 0.75, -1,  -1 },
  /*P*/{2.3, 1.653, 1.86, 1.58,  2.09, 1.47, 2.21,-1 },
  /*F*/{1.7, 1.35,  1.38, 1.45,  1.7,  0.92, 1.3,  1.41}};

int nbonds;

void add(int i,int j)
{
  if (cfg[i].nnbr>MAXVAL-1 || cfg[j].nnbr>MAXVAL-1)
    fprintf(stderr,"Cannot add bond %d-%d (max valence exceeded).\n",i,j);
  else {
    nbonds++;
    cfg[i].nbr[cfg[i].nnbr++]=j;
    cfg[j].nbr[cfg[j].nnbr++]=i; }
}

int main(int narg,char **arg)
{
  FILE *f;
  int i,j,n=0,ier;
  char line[128],*fn,*dot=NULL,*c,atc[128];
  float d=1.2;
  enum type_e { XYZ,ATM,PDB,PLB } type;
  FILE *mol,*plb=NULL;
  static float hdr[2];
  static int nsin=0,nline=0;

  if (narg<2) {
    fprintf(stderr,"\
Make MACSIMUS-compatible mol/plb files of xyz data in Angstrom. Call by:\n\
  bonds FILE{.xyz,.3dt,.plb,.atm,.pdb} [[-]DIST [VERBOSE]]\n\
Input:\n\
  FILE.xyz: format: x y z, also FILE.3dt\n\
  FILE.plb: plb-file\n\
  FILE.atm: format: 1st line: number of atoms\n\
                    2nd line: empty (=free b.c.) | [box] L (=cube) |\n\
                              [box] Lx Ly Lz (\"box\" is literal keyword)\n\
                    more lines: At x y z [charge...])\n\
            where At = chem. symbol or proton number\n\
  FILE.pdb: PDB format, threshold as above\n\
  DIST (for .xyz,.plb): threshold bond length [default=1.5]\n\
  DIST (for .atm,.pdb): threshold=DIST*(sum of van de Waals radii) [df.=1.2]\n\
  -DIST: forces periodic boundary conditions (.plb,.atm; box must be known)\n\
  VERBOSE: verbose output\n\
Environment:\n\
  MAXN = max number of atoms (not for plb), default=%d\n\
Output:\n\
  FILE.mol: output mol-file\n\
  FILE.plb: output plb-file (not if FILE.plb is input)\n\
See also:\n\
  plb2nbr sh/showatm c/shownbr util/gaussian2plb\n",N);
    exit(0); }

  if (narg>2) d=atof(arg[2]);

  if (getenv("MAXN")) N=atoi(getenv("MAXN"));

  f=fopen(arg[1],"rt");
  if (!f) Error(arg[1]);
  fn=strdup(arg[1]);
  for (c=fn; *c; c++) if (*c=='.') dot=c;
  if (!dot) Error("bonds: no file extension");
  dot++;
  if (!strcmp("3dt",dot) || !strcmp("xyz",dot) ) type=XYZ,d/=bond[attype("C")][attype("C")];
  else if (!strcmp("pdb",dot)) type=PDB;
  else if (!strcmp("atm",dot)) type=ATM;
  else if (!strcmp("plb",dot)) type=PLB,d/=bond[attype("C")][attype("C")];
  else Error("bonds: extension .3dt, .xyz, .pdb, .plb, or .atm expected");

  loop (i,0,NAT) loop (j,0,NAT)
    if (bond[i][j]<0) bond[i][j]=bond[j][i];
  loop (i,0,NAT) loop (j,0,NAT) bond[i][j]*=fabs(d);

  if (type==PLB) {
    float hdr[3];
    int i;
    
    if (2!=fread(hdr,4,2,f)) Error("bonds: plb no header");
    n=hdr[0];
    if (n>=NMAX || n<1) Error("bonds: too many atoms, or this is not a plb-file");
    allocarrayzero(cfg,n);
    if (hdr[1]<0) {
      if (3!=fread(L,4,3,f)) Error("bonds: plb too short");
      Lh[0]=L[0]/2; Lh[1]=L[1]/2; Lh[2]=L[2]/2;
      if (d<0) periodic++; }
    loop (i,0,n) {
      cfg[i].at=attype("C");
      nsin=i;
      if (3!=fread(cfg[i].r,4,3,f)) break; }
    nsin=n;
    goto done; }

  if (!cfg) allocarrayzero(cfg,N);

  nline=0;

  while (nline++,fgets(line,128,f)) if (strlen(line)>5 && line[0]!='#') {
    if (n>=N) Error("bonds: too many atoms - set environment variable MAXN");

    switch (type) {
      case XYZ:
        if (3!=sscanf(line,"%f%f%f",cfg[n].r,cfg[n].r+1,cfg[n].r+2)) goto done;
        cfg[n].at=attype("C");
        n++;
        break;

      case ATM: {
        char *box=strstr(line,"box");
        int i;

        if (nline==1) {
          nsin=atoi(line);
          goto cont; }

        if (type==ATM && nline==2) {
          i=sscanf(line,"%f%f%f",L,L+1,L+2);
          if (i==0 && box) i=sscanf(box+3,"%f%f%f",L,L+1,L+2);
          if (i==1) L[1]=L[2]=L[0];
          if (i==2) Error("bonds: two box sizes found (1 or 3 expected)");
          if (d<0) periodic++; 
          Lh[0]=L[0]/2; Lh[1]=L[1]/2; Lh[2]=L[2]/2;
          goto cont; }

        ier=sscanf(line,"%s%f%f%f%f",atc,cfg[n].r,cfg[n].r+1,cfg[n].r+2,&cfg[n].q);
       
        if (ier<4) goto done;
        if (ier==4) cfg[n].q=0;
        if (isdigit(atc[0])) {
          static char *Mendeleyev[111]={"e",
            "H",                                                                                 "He",
            "Li","Be",                                                  "B" ,"C" ,"N" ,"O" ,"F" ,"Ne",
            "Na","Mg",                                                  "Al","Si","P" ,"S" ,"Cl","Ar",
            "K", "Ca","Sc","Ti","V" ,"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
            "Rb","Sr","Y", "Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I", "Xe",
            "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
                           "Hf","Ta","W" ,"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
            "Fr","Ra","Ac","Th","Pa","U", "Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
                           "Rf","Db","Sg","Bh","Hs","Mt","Uun"};


          i=atoi(atc);
          if (i>110) i=110;
          strcpy(atc,Mendeleyev[i]); }
        cfg[n].at=attype(atc);
        n++;
        break; }

      case PDB: if (!memcmp(line,"ATOM  ",6) || !memcmp(line,"HETATM",6))
        if (line[16]<='A') /* only alt loc A */ {
          cfg[n].at=attype(line+13);
          sscanf(line+12,"%s",cfg[n].id);
          if (sscanf(line+32,"%f%f%f",cfg[n].r,cfg[n].r+1,cfg[n].r+2)!=3) Error(line);
          n++;
          break; }
      }

    cont:; }

 done:

  strcpy(dot,"mol");
  mol=fopen(fn,"wt");
  strcpy(dot,"plb");
  if (type!=PLB) plb=fopen(fn,"w");
  hdr[0]=(float)n;
  if (nsin && n!=nsin) fprintf(stderr,"WARNING NS=%d in header, %d found\n",nsin,n);

  if (plb) {
    hdr[1]=-3;
    fwrite(hdr,4,2,plb);
    fwrite(L,4,3,plb); }

  fprintf(mol,"molecule %s\n\
\n\
parameter_set = dumb\n\
number_of_atoms = %d\n\natoms\n\
! i Crad/atom-id   a-type  charge  chir nb bound_atoms\n",fn,n);

  //  loop (i,0,n)
  //    loop (j,0,i) printf("%d-%d %g %g\n",i,j,dist(r[i],r[j]),bond[at[i]][at[j]]);

  loop (i,0,n)
    loop (j,0,i)
      if (dist(cfg[i].r,cfg[j].r)<bond[cfg[i].at][cfg[j].at]) {
        add(i,j);
        if (narg>3) printf("%d %d %g\n",i,j,dist(cfg[i].r,cfg[j].r)); }
 
  loop (i,0,n) {
    char ati=ats[cfg[i].at];
    char line[64];

    if (!cfg[i].id[0]) {
      sprintf(line,"%c-%03d",ati,i);
      memcpy(cfg[i].id,line,7); }
    fprintf(mol,"%3d %s %c %9.6f 0 %d ",i,cfg[i].id,ati,cfg[i].q,cfg[i].nnbr);
    loop (j,0,cfg[i].nnbr) fprintf(mol," %d",cfg[i].nbr[j]);
    nval[cfg[i].nnbr]++;
    if (plb) fwrite(cfg[i].r,4,3,plb);
    fprintf(mol,"\n"); }

  fclose(mol);
  if (plb) fclose(plb);

  printf("# %d atoms, %d bonds\n",n,nbonds);
  printf("# valence #_of_atoms\n");
  loopto (i,0,MAXVAL) if (nval[i]) printf("# %3d %3d\n",i,nval[i]);
  printf("# box=%g %g %g %s b.c.\n",L[0],L[1],L[2],periodic?"periodic":"free");

  return 0;
}
