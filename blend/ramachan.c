/* cc -Wall -O2 -o ramachan ramachan.c -lm
  update 7/11: output name derived from simname, not sysname
               bug fixed: new plb format
  update 9/01: reads backbone info from a par-file
  update 9/00: peptide bond angle `omega' added
*/
#include "../gen/include.h"

#define VVV(A,B,C) { A[0] B[0] C[0]; A[1] B[1] C[1]; A[2] B[2] C[2]; }
#define SQR(A) (A[0]*A[0]+A[1]*A[1]+A[2]*A[2])
#define SCAL(A,B) (A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define VECT(A,B,C) { \
  A[0]=B[1]*C[2]-B[2]*C[1]; \
  A[1]=B[2]*C[0]-B[0]*C[2]; \
  A[2]=B[0]*C[1]-B[1]*C[0]; }

char name[4][8];

FILE *ble,*ram,*plb,*par;

typedef double vector[3];
typedef float fvector[3];

int fullangle;
int idih;

#include "../gen/dihangle.c"

int ns;

double Dihedral(fvector *r,int i[4])
{
 int h;

 loop (h,0,4) if (i[h]<0 || i[h]>=ns) Error("ramachan: atom index out of range. Possible cause:\n\
  - ble- and plb-files are inconsistent\n\
  - bad option -n (should match N[0] in the simulation def-file)\n\
  - more species (N[1] etc. is not zero)");

 return dihedral(r[i[0]],r[i[1]],r[i[2]],r[i[3]]);
}

struct angle_s {
  struct angle_s *next;
  struct angle_s *omega0,*psi,*omega1; /* for phi-based list only */
  int count; /* for psi,omega: count of assignements to phi */
  enum type_e { NONE,PHI,PSI,OMEGA } type;
  float a;
  int i[4];
} *phi0,*psi0,*omega0;

struct info_s {
  struct info_s *next;
  char s[80];
} *info0;

struct angle_s unknown={NULL,NULL,NULL,NULL,0,NONE,999,{-1,-1,-1,-1}};

/* dihedral angles, atoms in brakets are for preceding/following residue
   omega = [Ca-C]-N-Ca
   phi =      [C]-N-Ca-C
   psi =          N-Ca-C-[N]
   omega =          Ca-C-[N-Ca]
*/

struct bblist {
  struct bblist *next;
  char atom[8]; } *NCaC[3]; /* [0]=N, [1]=Ca, [2]=C */

int is(struct bblist *head,int i)
{
  struct bblist *bb;

  looplist (bb,head)
    if (!strcmp(bb->atom,name[i])) return 1;

  return 0;
}

#define isN(I) is(NCaC[0],I)
#define isCA(I) is(NCaC[1],I)
#define isC(I) is(NCaC[2],I)

char line[256];

void addNCaC(int ii) /********************************************** addNCaC */
{
  char *tok=strtok(line," \n\t");

  while (tok) {
    struct bblist *bb;

    allocone(bb);
    if (strlen(tok)>7) Error("ramachan: too long atom type name in backbone table");
    strcpy(bb->atom,tok);
    bb->next=NCaC[ii];
    NCaC[ii]=bb;
    tok=strtok(NULL," \n\t"); }
}

void addangle(struct angle_s **head,enum type_e type,float a, /**** addangle */
	      int i0,int i1,int i2,int i3)
{
  struct angle_s *x;
  struct info_s *i;

  alloconezero(x);
  alloconezero(i);
  i->next=info0; info0=i;

  if (ram) switch (type) {
    case PHI:
      sprintf(i->s,"%4d   phi = %7.2f  C%d - N%d - Ca%d - C%d\n",idih,a,i0,i1,i2,i3);
      break;
    case PSI:
      sprintf(i->s,"%4d   psi = %7.2f  N%d - Ca%d - C%d - N%d\n",idih,a,i0,i1,i2,i3);
      break;
    case OMEGA:
      sprintf(i->s,"%4d omega = %7.2f  Ca%d - C%d - N%d - Ca%d\n",idih,a,i0,i1,i2,i3);
      break;
    default:
      Error("ramachan: internal"); }

  x->type=type;
  x->a=a;
  x->i[0]=i0; x->i[1]=i1; x->i[2]=i2; x->i[3]=i3;

  x->next=*head; *head=x;
}

void printphi(int verbose,FILE *f) /******************************* printphi */
{
  struct angle_s *phi;
  int n;
  char c1='[',c2=']';

  if (verbose) c1='(',c2=')';

  for (phi=phi0,n=1; phi; phi=phi->next,n++)
    if (verbose==-1) {
      if (fabs(phi->psi->a)>998) fprintf(f," %8.3f",phi->a);
      else fprintf(f," %8.3f %8.3f",phi->a,phi->psi->a); }
    else {
      fprintf(f,"%7.2f %7.2f %7.2f %7.2f %c%d%c%s",
              phi->omega0->a,phi->a,phi->psi->a,phi->omega1->a,c1,n,c2,n>9?"":" ");
      if (verbose)
        fprintf(f,"  %4d %4d %4d %4d %4d %4d %4d",
                phi->omega0->i[0],
                phi->i[0],phi->i[1],phi->i[2],phi->i[3],
                phi->psi->i[3],
                phi->omega1->i[3]);
      if (fputc('\n',f)<0) Error("ramachan: write error"); }

  if (verbose==-1) if (fputc('\n',f)<0) Error("write error");
}

int main(int narg,char **arg) /**************************************** main */
{
  int n[4];
  struct angle_s *phi,*psi,*omega;
  struct info_s *inf;
  double pp;
  char dummy[64];
  char fn[256],*rr,*ch,*dot,*sysname=NULL,*plbname=NULL;
  fvector *r;
  float hdr[2];
  int frame,from=1,to=-1,by=1,optr=0,i,ii,nmol=1,plotit=0,issum=0;
  int varL=0;
  char *parset=NULL;
  FILE *summary=NULL;
  char *point="p";

  if (narg<2) {
    fprintf(stderr,"Ramachandran plot from blend and playback files.\n\
  ramachan [OPTIONS] SYSNAME [SIMNAME.plb] [OPTIONS]\n    \
OPTIONS:\n\
  -a : angles in [0,360) [default=[-180,180]]\n\
  -r : write results to separate files SIMNAME.r#, where #=frame number\n\
       [default=concatenated to SYSNAME.ram]\n\
  -s : write summary phi,psi file (1line=1frame; for columns see SYSNAME.mar)\n\
  -n#: number of molecules (N[0] in the simulation def-file)\n\
       BUG: only identical molecules supported (N[1]=0 required)\n\
  -f#: from frame [default=1]\n\
  -t#: to frame [incl., default=-1=until eof]\n\
  -b#: by frame [default=1=every frame]\n\
  -p[@] : plot the resulting ram-file; @=.pPoO (default=P, see plot)\n\
  -PPARSET: parameter set (relative to BLENDPATH)\n\
            [default=as specified in the ble-file]\n\
FILES:\n\
  SYSNAME.ble: input ble-file\n\
  SIMNAME.plb: input playback file\n\
  SYSNAME.mar: output dihedral angle info extracted from SYSNAME.ble\n\
  SIMNAME.ram: output Ramachandran plot(s) omega,phi,psi,omega:\n\
               first table from SYSNAME.ble with site info\n\
               tables calculated from playback files follow (if not -r)\n\
  SYSNAME.sum: output Ramachandran summary by columns\n\
  SIMNAME.r1,...: tables calculated from playback files (if -r)\n\
  ${BLENDPATH}/PARSET.par: parameter file (only table `backbone' used)\n\
not defined angles/sites are denoted 999/-1\n\
");
   exit(0); }

  loop (i,1,narg)
    if (arg[i][0]=='-') {
      int iarg=atoi(arg[i]+2);

      switch (arg[i][1]) {
        case 'n': nmol=iarg; break;
        case 'f': from=iarg; break;
        case 't': to=iarg; break;
        case 'b': by=iarg; break;
        case 'r': optr++; break;
        case 's': issum++; break;
        case 'a': fullangle++; break;
        case 'p': plotit++; point=arg[i]+2; break;
        case 'P': parset=arg[i]+2; break;
        default: Error("unknown option"); }
    }
    else if (!sysname) sysname=arg[i];
    else if (!plbname) plbname=arg[i];
    else Error("too many file arguments");

  strcpy(fn,sysname);
  strcat(fn,".ble");
  ble=fopen(fn,"rt");
  if (!ble) Error("no ble-file");

  strcpy(fn,sysname);
  strcat(fn,".mar");
  ram=fopen(fn,"wt");
  if (!ram) Error("open inv ram-file for writing");

  strcpy(fn,plbname);
  dot=NULL;
  for (ch=fn; *ch; ch++) if (*ch=='.') dot=ch;
  if (!dot) Error("PLBNAME: no extension (normally .plb expected)");

  if (issum) {
    strcpy(dot,".sum");
    summary=fopen(fn,"wt");
    if (!summary) Error("open summary-file for writing"); }

  if (parset) {
    char *env=getenv("BLENDPATH"),*c;

    if (!env) env="";
    alloczero(c,strlen(env)+strlen(parset));
    sprintf(c,"%s/%s.par",env,parset);
    parset=c; }
  else {
    while ( (rr=fgets(line,256,ble)) ) {
      if ( (ch=strchr(line,'\n')) ) *ch=0;
      if (!memcmp("! parameter file ",line,17)) {
        parset=line+17;
        strcpy(parset+strlen(parset)-4,".par");
        break; } }
    if (!rr) Error("! parameter file info not found in the ble-file"); }

  fprintf(stderr,"parameter file = \"%s\"\n",parset);
  par=fopen(parset,"rb");
  if (!par) Error("cannot open parameter file");
  parset=NULL; /* no longer valid */

  while ( (rr=fgets(line,256,par)) ) {
    if ( (ch=strchr(line,'\n')) ) *ch=0;
    if (!strcmp("backbone",line)) break; }

  if (!rr) Error("no \"backbone\" keyword in the parameter file");

  ii=0;
  while ( (rr=fgets(line,256,par)) ) {
    if (rr[0]!='!') addNCaC(ii++);
    if (ii==3) break; }

  if (ii!=3) Error("wrong backbone info in the parameter file");

  fclose(par);

  /* reading dihedrals */
  while ( (rr=fgets(line,256,ble)) ) {
    if ( (ch=strchr(line,'\n')) ) *ch=0;
    if (!strcmp("dihedrals",line)) break; }

  if (!rr) Error("no \"dihedrals\" keyword");

  while ( (rr=fgets(line,256,ble)) ) {
    if (line[0]=='!') continue;
    idih++;

    if (12!=sscanf(line,"%d%s%d%s%d%s%d%s%s%s%s%lf",
               n,name[0],n+1,name[1],n+2,name[2],n+3,name[3],
               dummy,dummy,dummy,&pp))
      break;

    if (fullangle) pp=fmod(pp+360,360);

    if (isC(0) && isN(1) && isCA(2) && isC(3)) addangle(&phi0,PHI,pp,n[0],n[1],n[2],n[3]);
    if (isC(3) && isN(2) && isCA(1) && isC(0)) addangle(&phi0,PHI,pp,n[3],n[2],n[1],n[0]);

    if (isN(0) && isCA(1) && isC(2) && isN(3)) addangle(&psi0,PSI,pp,n[0],n[1],n[2],n[3]);
    if (isN(3) && isCA(2) && isC(1) && isN(0)) addangle(&psi0,PSI,pp,n[3],n[2],n[1],n[0]);

    if (isCA(0) && isC(1) && isN(2) && isCA(3)) addangle(&omega0,OMEGA,pp,n[0],n[1],n[2],n[3]);
    if (isCA(3) && isC(2) && isN(1) && isCA(0)) addangle(&omega0,OMEGA,pp,n[3],n[2],n[1],n[0]);
  }

  ii=0;
  looplist (inf,info0) {
    char *omega=strstr(inf->s,"omega");
    int x=0;

    if (!omega) x=++ii;

    fprintf(ram,"%3d %s",x,inf->s); }

  fclose(ram);
  fclose(ble);

  /* assigning adjacent omega,psi to the phi-list */
  /* Ca-C-N-Ca-C-N-Ca
     omega:  0 1 2 3
     phi:      0 1 2  3
     psi:        0 1  2 3
     omega:        0  1 2 3    */

  for (phi=phi0; phi; phi=phi->next) {

    for (psi=psi0; psi; psi=psi->next)
      if (phi->i[1]==psi->i[0] && phi->i[2]==psi->i[1] && phi->i[3]==psi->i[2])
        phi->psi=psi, psi->count++;

   for (omega=omega0; omega; omega=omega->next) {
     if (phi->i[0]==omega->i[1] && phi->i[1]==omega->i[2] && phi->i[2]==omega->i[3])
       phi->omega0=omega, omega->count++;
     if (phi->i[2]==omega->i[0] && phi->i[3]==omega->i[1])
       phi->omega1=omega, omega->count++; } }

  for (psi=psi0; psi; psi=psi->next) if (!psi->count) {
    fprintf(stderr,"WARNING psi=%7.2f  N%d - Ca%d - C%d - N%d NOT ASSIGNED TO phi\n",psi->a,psi->i[0],psi->i[1],psi->i[2],psi->i[3]);
    alloconezero(phi);
    phi->psi=psi; psi->count++;
    phi->a=999;
    phi->next=phi0;
    phi0=phi;

    for (omega=omega0; omega; omega=omega->next) {
      if (psi->i[1]==omega->i[0] && psi->i[2]==omega->i[1] && psi->i[3]==omega->i[2])
        phi->omega1=omega, omega->count++; } }

  for (omega=omega0; omega; omega=omega->next) if (!omega->count) {
    /* generally, these omega should be tried to assign to previous psi */
    fprintf(stderr,"WARNING omega=%7.2f  Ca%d - C%d - N%d - Ca%d NOT ASSIGNED TO phi\n",omega->a,omega->i[0],omega->i[1],omega->i[2],omega->i[3]);
    alloconezero(phi);
    phi->omega0=omega; omega->count++;
    phi->a=999;
    phi->next=phi0;
    phi0=phi; }

  /* not-assigned are marked as unknown to be printable (no NULL pointers) */
  for (phi=phi0; phi; phi=phi->next) {
    if (!phi->omega0) phi->omega0=&unknown;
    if (!phi->omega1) phi->omega1=&unknown;
    if (!phi->psi) phi->psi=&unknown; }

  strcpy(fn,plbname);
  strcpy(dot,".ram");
  ram=fopen(fn,"wt");
  if (!ram) Error("ramachan: open ram-file for writing");
  fprintf(ram,"# Ramachandran plot                           omega0 phi psi omega1\n");
  fprintf(ram,"# angle 999 marks undefined value              /  \\ /  \\ /  \\ /  \\\n");
  fprintf(ram,"# omega0  phi     psi    omega1         Ca    C    N   Ca    C    N   Ca \n");

  printphi(1,ram);
  if (optr) fclose(ram);

  if (!plbname) return 1;

  plb=fopen(plbname,"rb");
  if (!plb) Error("no plb-file");
  if (2!=fread(hdr,4,2,plb)) Error("ramachan: plb bad header"); 
  ns=hdr[0];
  if (ns<1 || ns>16777216) Error("ramachan: bad number of sites");
  varL=hdr[1]<0;
  if (nmol<=0) Error("ramachan: bad option -n");
  ns/=nmol;
  if ((int)hdr[0]!=ns*nmol) Error("number of sites in the plb-file is not an integer multiple\n\
of the number of molecules (option -n) specified");

  fprintf(stderr,"%d sites in one configuration, %d sites in one molecule, %d molecules\n",(int)hdr[0],ns,nmol);

  r=malloc(ns*sizeof(fvector));
  if (!r) Error("no heap");

  strcpy(fn,plbname);
 again:
  rr=fn;
  ch=NULL;
  while (*++rr) if (*rr=='.') ch=rr;
  if (!ch) { strcat(fn,"."); goto again; }

  for (frame=1; to<0||frame<=to; frame++) {
    if (varL) {
      if (fread(r,12,1,plb)!=1) break; }
    if (fread(r,12,ns,plb)!=ns) break;
    if (frame<from || (frame-from)%by) continue;

    sprintf(ch,".r%d",frame);
    if (optr)
      ram=fopen(fn,"wt");
    else {
      if (nmol>1)
        fprintf(ram,"#\n# configuration %d (frame %d, molecule %d)\n",
                frame,(frame-1)/nmol+1,(frame-1)%nmol+1);
      else
        fprintf(ram,"#\n# frame %d\n",frame); }

    fprintf(ram,"# omega0  phi     psi    omega1\n");

    for (phi=phi0; phi; phi=phi->next) {
      int OK=1; /* not used */

      if (phi->omega0 == &unknown) OK=0;
      else phi->omega0->a=Dihedral(r,phi->omega0->i);

      if (phi->type == NONE) OK=0;
      else phi->a=Dihedral(r,phi->i);

      if (phi->psi == &unknown) OK=0;
      else phi->psi->a=Dihedral(r,phi->psi->i);

      if (phi->omega1 == &unknown) OK=0;
      else phi->omega1->a=Dihedral(r,phi->omega1->i); }

    printphi(0,ram);
    if (summary) printphi(-1,summary);

    if (optr) {
      fclose(ram);
      fprintf(stderr,"%s done\r",fn); }
    else
      fprintf(stderr,"frame %d done\r",frame); }

  fprintf(stderr,"\n");

  if (!optr) fclose(ram);
  if (summary) fclose(summary);

  *ch=0;
  strcpy(fn,plbname);
  *dot=0;
  if (plotit) {
    char *s=malloc(strlen(fn)+30+optr*frame*(strlen(fn)+9));
    int i;
    char *rg=fullangle?"[0:360][0:360]":"[-180:180][-180:180]";

    if (!s) Error("no heap");

    sprintf(s,"plot %s %s.ram:2:3:",rg,fn);
    if (optr) {
      sprintf(s,"plot %s %s.ram:2:3:o @:2:3:pC",rg,sysname);
      loopto (i,1,frame) sprintf(s+strlen(s)," %s.r%d",fn,i); }
    else
      strcat(s,point[0]?point:"P");

    fprintf(stderr,"%s\n",s);
    if (system(s)) Error("system call failed"); }

  return 0;
}
