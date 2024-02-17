/* cc -Wall -O2 -o plbsites plbsites.c -lm
*/
#include "../gen/include.h"

int ifollow,iframe,nfollow;

double *old,maxd=-1;
char *fmt=" %9.6f";
FILE *out=NULL;

double follow(double x,double L,FILE *out) /************************* follow */
{
  if (nfollow>0) {
    if (ifollow>=nfollow) Error("too many arguments");

  //  fprintf(stderr,"%d x=%g old=%g %g %g old[0]=%g\n",ifollow,x,old[ifollow],L,maxd,old[0]);

    if (iframe) {
      while (x<old[ifollow]-L/2) x+=L;
      while (x>old[ifollow]+L/2) x-=L;
      Max(maxd,fabs(x-old[ifollow])/L) }
    old[ifollow]=x;
    ifollow++; }

  if (out) fprintf(out,fmt,x);

  return x;
}

int main(int narg,char **arg) /**************************************** main */
{
  int i,k,ns,nsmol=1,nmol=1,off=0,n,divL=0;
  float (*flt)[3],hdr[2];
  double *L,(*r)[3];
  FILE *plb;
  int varL=0;
  FILE **outfiles=NULL,*outL=NULL;
  char *fnfmt=NULL;

  if (narg<2) {
    fprintf(stderr,"\
Print trajectory of selected atoms in ASCII. Call by:\n\
  plbsites PLB-FILE [OPTIONS] ATOM[xyz] [ATOM[xyz] ..]\n\
OPTIONS:\n\
  -oOFFSET = number added to each ATOM number [default=0]\n\
  -nNS = number of sites per molecule [default=1]\n\
  -mNM = number of molecules [default=1]:\n\
         pattern of ATOMs is repeated NM times with multiples of NS added\n\
  -fFMT = output format, default=\"%s\"\n\
  -FFMT = every atom will be printed to a separate xyz file\n\
          FMT must contain a format accepting int for the atom number\n\
          [default=stdout, one line per frame]\n\
  -l    = follow periodic boundary conditions\n\
  -L    = as above, r divided by box (good for NPT)\n\
ATOM = in range [0,NS-1]: atom (site) number\n\
       -1: box size\n\
       [-NS,-2]: center of mass (eq. weights) of atoms 0..|NS|-1\n\
NB: ATOM numbers are changed by options -o,-m\n\
If x,y,z is given, only selected coordinates are printed (default=xyz)\n\
Bug: sites are always in order x,y,z even if different order is specified\n\
Examples:\n\
  plbsites xxx.plb -1z 33zx               # 3 columns: Lz x_33 z_33\n\
  plbsites xxx.plb -1z 33z 33x            # 3 columns: Lz z_33 x_33\n\
  plbsites xxx.plb -Fa%%03d.xyz -n5 -m20 1 # a001.xyz, a006.xyz, ..., a096.xyz\n\
See also:\n\
  plbfilt plb2plb plb2asc\n\
",fmt);
    exit(0); }

  loop (i,2,narg) if (arg[i][0]=='-') switch (arg[i][1]) {
    case 'm': nmol=atoi(arg[i]+2); break;
    case 'n': nsmol=atoi(arg[i]+2); break;
    case 'o': off=atoi(arg[i]+2); break;
    case 'L': divL=1;
    case 'l': nfollow=-1; break;
    case 'f': fmt=arg[i]+2; break;
    case 'F': fnfmt=arg[i]+2;
              if (!strchr(arg[i],'%')) Error("plbsites: -F: no % in format");
              break;
    default: if (!isdigit(arg[i][1])) Error("plbsites: bad option"); }

  plb=fopen(arg[1],"rb");
  if (!plb) Error("plbsites: no such file");

  if (2!=fread(hdr,4,2,plb)) Error("plbsites: plb no header");
  ns=hdr[0];
  allocarray(flt,ns+1); flt++; /* flt[-1]=L (float) */
  allocarray(r,ns+1); r++;     /* r[-1]=L (double) */

  if (nfollow==-1) nfollow=(ns+1)*6; /* normally should not exceed */
  if (nfollow) allocarray(old,nfollow);

  if (fnfmt) allocarrayzero(outfiles,ns);

  if (hdr[1]<0) {
    if (hdr[1]<0) varL=1;
    else Error("plbsites: wrong L in plb file"); }
  else
    r[-1][0]=r[-1][1]=r[-1][2]=hdr[1];
  L=r[-1];

  for (iframe=0;;iframe++) {
    int oo=off;

    if (ns+varL!=fread(flt-varL,12,ns+varL,plb)) break;

    if (divL) {
      loop (i,0,ns) loop (k,0,3) r[i][k]=(double)flt[i][k]/flt[-1][k];
      loop (k,0,3) r[-1][k]=1; }
    else
      loop (i,-varL,ns) loop (k,0,3) r[i][k]=(double)flt[i][k];

    ifollow=0;

    loop (n,0,nmol) {
      loop (i,2,narg) if (arg[i][0]!='-' || isdigit(arg[i][1])) {
        int s=atoi(arg[i]);
        int all=1;
        int site;

        /* bug or feature:
           -o-1 0x gives Lx, too
           -o-2 0x gives center of mass of sites 0 and 1 */
        if (s<0) s-=oo;
        if (s>=ns) Error("plbsites: specified site number out of range");
        if (s<-ns) Error("plbsites: specified site numbers out of range for center of mass");

        site=oo+s;
        if (site>=ns) Error("(calculated) site number out of range");
        if (outfiles) {
          if (site<0) {
            if (outL==NULL) outL=fopen(string(fnfmt,site),"wt");
            out=outL; }
          else {
            if (outfiles[site]==NULL) outfiles[site]=fopen(string(fnfmt,site),"wt");
            out=outfiles[site]; } }
        else
          out=stdout;

        if (!out) Error("plbsites: no output file (wrong format, too many files)");

        if (site<=-2) {
          /* center of mass from 0 to site (not incl.) */
          int isite;
          double sum[3]={0,0,0};

          loop (isite,0,abs(site)) {
            if (strchr(arg[i],'x')) sum[0]+=follow(r[isite][0],L[0],NULL),all=0;
            if (strchr(arg[i],'y')) sum[1]+=follow(r[isite][1],L[1],NULL),all=0;
            if (strchr(arg[i],'z')) sum[2]+=follow(r[isite][2],L[2],NULL),all=0;
            if (all) loop (k,0,3) sum[k]+=follow(r[isite][k],L[k],NULL); }

          all=1;
          if (out) {
            if (strchr(arg[i],'x')) fprintf(out,fmt,sum[0]/abs(site)),all=0;
            if (strchr(arg[i],'y')) fprintf(out,fmt,sum[1]/abs(site)),all=0;
            if (strchr(arg[i],'z')) fprintf(out,fmt,sum[2]/abs(site)),all=0;
            if (all) loop (k,0,3) fprintf(out,fmt,sum[k]/abs(site)); } }
        else {
          /* single site or box size */
          if (strchr(arg[i],'x')) follow(r[site][0],L[0],out),all=0;
          if (strchr(arg[i],'y')) follow(r[site][1],L[1],out),all=0;
          if (strchr(arg[i],'z')) follow(r[site][2],L[2],out),all=0;
          if (all) loop (k,0,3) follow(r[site][k],L[k],out); }
      }

      if (outfiles) if (out) fprintf(out,"\n");
      oo+=nsmol; }
    if (!outfiles) printf("\n"); }

  if (nfollow)
    fprintf(stderr,"max jump between frames = %g L%s\n",
            maxd,maxd>0.4?" ********** TOO MUCH!!! **********":"");

  fclose(plb);

  return 0;
}
