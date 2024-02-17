/* \make cfg2asc

*/

#include "ground.h"
#include "sds.h"
#include "varfile.h"

typedef double vector[3];

int main(int narg,char **arg) /**************************************** main */
{
  FILE *asc=NULL;
  vector L={0,0,0};
  char *ascfn="-",*cfgfn=NULL;
  struct a_s {
    int  size;    /* whole struct in bytes */
    int  padd;    /* padded to 8 bytes */
    double logs;  /* log of the Nose variable s */
    vector lambda;/* log(box) */
    vector shape; /* (reserved for box shape) */
    vector rp[1]; /* contiguous array of all r[ns] (and p[ns]) */
  } *a[3];
  int nspec;
  struct spec_s {
    int N;        /* number of molecules */
    int ns;       /* number of sites */
  } *spec=NULL; /* [nspec] */
  struct rec_s {
    int key;
    int intval;
    vector vecval;
  } rec;
  int i,k,cfgkey,N=0,ns=0,nsx,oldsize,iarg;
  int np; /* 1:nonpolar, 2:polar */
  int optionlist[32];
  double t=0,h=0,Entot=0;

  if (narg<2) {
    fprintf(stderr,"\
Convert MACSIMUS cfg-file CFG to its ASCII image ASC:\n\
  cfg2asc CFG ASC\n\
As above with CFG=NAME.cfg and ASC=NAME.cfa:\n\
  cfg2asc NAME[.cfg]\n\
Print info:\n\
  cfg2asc CFGNAME OPTION [OPTION..]\n\
where\n\
  -a# = x,y,z accelerations of site # (if available) [AA/ps2]\n\
  -d = line of selected info in the format of the def-file\n\
  -E = energy [K]\n\
  -h = timestep [ps]\n\
  -L = box sizes [AA]\n\
  -N = print all numbers of molecules in one line\n\
  -n[#] = number of molecules: total (No.N) [of species #]\n\
  -r# = x,y,z coordinates of site # [AA]\n\
  -s[#] = number of sites: total [of 1 molecule of species #]\n\
  -S = number of species (nspec)\n\
  -t = print time [ps]\n\
  -v# = x,y,z velocities of site # (shifted for Verlet) [AA/ps]\n\
Examples (bash):\n\
  for n in {1..100} ; do cfg2asc -n; cfg2asc SIM.$n SIM.$n.asc ; done\n\
See also:\n\
  asc2cfg cfg2atm plb2cfg cfgconv\n\
");
    exit(0); }

  initscroll(0);
  stdin=NULL;

  cfgfn=arg[1];
  if (narg==3) {
    ascfn=arg[2]; }
  else if (narg==2) {
    char *fn=strdup(arg[1]),*dot=strrchr(fn,'.');

    if (dot && !strcmp(dot,".cfg")) *dot=0;
    cfgfn=string("%s.cfg",fn);
    ascfn=string("%s.cfa",fn); }
  else
    if (ascfn[0]!='-') ERROR(("too many arguments"))

  VarOpen(cfgfn,"r");
  if (ascfn[0]!='-')
    if ( !(asc=fopen(ascfn,"wt")) ) ERROR(("cannot write to \"%s\"",ascfn))

  /* if !asc, we process options */
  if (asc) prt("reading \"%s\"",cfgfn);
  VarRead(&cfgkey,sizeof(cfgkey));
  np = cfgkey&2 ? 2 : 1;
  if (asc) fprintf(asc,"%d cfgkey (%s)\n",cfgkey,np==2?"polar":"nonpolar");

  VarRead(optionlist,sizeof(optionlist));
  if (asc) loop (i,0,32) fprintf(asc,"-%c%d\n",tolower(i+'@'),optionlist[i]);

  VarRead(&nspec,sizeof(nspec));
  if (asc) fprintf(asc,"%d nspec\n",nspec);
  if (nspec) {
    allocarray(spec,nspec);
    VarRead(spec,nspec*sizeof(spec[0]));
    if (asc) loop (i,0,nspec) fprintf(asc,"%d %d N ns (of species %d)\n",spec[i].N,spec[i].ns,i); }
  else
    prt("WARNING: missing species info (OK for cook)");

  while (VarFile.size==sizeof(struct rec_s)) {
    VarRead(&rec,sizeof(rec));
    if (asc) fprintf(asc,"%d %d %.16g %.16g %.16g key %s\n",
            rec.key,rec.intval,rec.vecval[0],rec.vecval[1],rec.vecval[2],
            rec.key==1?"N t h En.tot" : rec.key==2?"ns L[3]" : "intval vecval[3]");

    switch (rec.key) {
      case 1: N=rec.intval;
        if (asc) put2(cfgkey,N)
        t=rec.vecval[0];
        h=rec.vecval[1];
        Entot=rec.vecval[2];
        if (asc) put3(t,h,Entot)
        break;

      case 2:
        ns=rec.intval;
        loop (k,0,3) L[k]=rec.vecval[k];
        if (asc) prt("Ns=%d L=%g %g %g loaded",ns,L[0],L[1],L[2]);
        break;
    } }

  if (asc) fprintf(asc,". end of header: cfg, vel, acc follow, VarFile.size=%u\n",VarFile.size);

  if (asc) put(VarFile.size)
  nsx=(VarFile.size-sizeof(*a[0]))/sizeof(vector)+1;
  //  nsx=(VarFile.size-((char*)(a[0]->rp) - (char*)(&a[0]->size)))/sizeof(vector);
  if (asc && np==2) prt("POLAR configuration detected");
  if (ns && nsx/np!=ns) prt("WARNING: ns=%d from record length, %d from key=2: %d used",nsx,ns,nsx);
  ns=nsx/np;
  if (asc) put(ns)
  if (ns<N) ERROR(("ns<N read from file"))

  loop (i,0,3) sdsalloczero(a[i],VarFile.size);
  loop (i,0,3) {
    oldsize=VarFile.size;
    if (VarFile.size) {
      VarRead(a[i],VarFile.size);
      if (a[i]->size!=oldsize)
        ERROR(("a[%d]: size of SDS in file (%d) does not match VarFile.size (%d)",
               i,a[i]->size,oldsize)) } }

  if (asc) {
    fprintf(asc,"! value h*velocity h^2*acceleration[3]; for leap-frog, h*velocity=value(t)-value(t-h)\n");
    fprintf(asc,"%20.16f %.16g %.16g logs\n",
            a[0]->logs, a[1]->logs, a[2]->logs);
    fprintf(asc,"%20.16f %20.16f %20.16f  %.16g %.16g %.16g  %.16g %.16g %.16g  lambda=ln(L)\n",
            a[0]->lambda[0], a[0]->lambda[1], a[0]->lambda[2],
            a[1]->lambda[0], a[1]->lambda[1], a[1]->lambda[2],
            a[2]->lambda[0], a[2]->lambda[1], a[2]->lambda[2]);

    fprintf(asc,"! r[3] h*velocity[3] h^2*acceleration[3] site%s\n",np==1?"":" and Drude");
    loop (i,0,nsx)
      fprintf(asc,"%20.16f %20.16f %20.16f  %.16g %.16g %.16g  %.16g %.16g %.16g %s[%d]\n",
              a[0]->rp[i][0], a[0]->rp[i][1], a[0]->rp[i][2],
              a[1]->rp[i][0], a[1]->rp[i][1], a[1]->rp[i][2],
              a[2]->rp[i][0], a[2]->rp[i][1], a[2]->rp[i][2],
              i<ns?"site":"Drude",
              i%ns); }

  VarClose();

  if (asc) {
    /* ASCII image has been written */
    fclose(asc);
    prt("ascii image \"%s\" written",ascfn); }
  else
    /* dump data according to optionlist */
    loop (iarg,2,narg) {
      int j=0, i=atoi(arg[iarg]+2);
      double q=1;
      char *info="?";

      switch (arg[iarg][1]) {
        case 'd':
          if (spec) loop (j,0,nspec) prt("N[%d]=%d",j,spec[j].N);
          prt("h=%.16g t=%.16g",h,t);
          prt("L[0]=%.16g L[1]=%.16g L[2]=%.16g",L[0],L[1],L[2]);
          break;
        case 'E':
          prt("%.16g",Entot);
          break;
        case 'h':
          prt("%.16g",h);
          break;
        case 'L':
          prt("%.16g %.16g %.16g",L[0],L[1],L[2]);
          break;
        case 'n':
          if (!arg[iarg][2]) prt("%d",N);
          else if (i>=0 && i<nspec) prt("%d",spec[i].N);
          break;
        case 'N':
          if (spec) loop (j,0,nspec) prt_("%d%c",spec[j].N," \n"[j==nspec-1]);
          break;
        case 'S':
          prt("%d",nspec);
          break;
        case 's':
          if (!arg[iarg][2]) prt("%d",ns);
          else if (i>=0 && i<nspec) prt("%d",spec[i].ns);
          break;
        case 't':
          prt("%.16g",t);
          break;
        case 'r':
          info="positions";
          goto aa;
        case 'v':
          j=1;
          info="velocities";
          q=h;
          goto aa;
        case 'a':
          j=2;
          info="accelerations";
          q=h*h;
        aa:
          if (i<0 || i>=ns) ERROR(("site #=%d out of range (ns=%d)",i,ns))
          if (!a[j]) ERROR(("%s not available",info))
          prt("%.16g %.16g %.16g", a[j]->rp[i][0]/q, a[j]->rp[i][1]/q, a[j]->rp[i][2]/q);
          break;
        default:
          ERROR(("%s is unknown option",arg[iarg])) } }

  return 0;
}
