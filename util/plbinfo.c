/* cc -O2 -o plbinfo plbinfo.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* C-style loop (last not included) */
#define loop(_I,_From,_To) for (_I=_From;_I<_To;_I++)

/* Pascal etc. style loop (last included) */
#define loopto(_I,_From,_To) for (_I=_From;_I<=_To;_I++)

typedef float floatvector[3];

void endian(char *a) /********************************************** endian */
{
  char x;

  x=a[0]; a[0]=a[3]; a[3]=x;
  x=a[1]; a[1]=a[2]; a[2]=x;
}

int main(int narg,char **arg) /*************************************** main */
{
  int iarg,iarg0=1,pair=0,ret=0,key=0,nms=1,printfn=0;
  double m=1e6;

  if (narg<2) {
    fprintf(stderr,"\
Print info on MACSIMUS playback files (binary). Call by:\n\
  plbinfo [OPTION] [-]FILE [[-]FILE ...]\n\
OPTION (1st arg, if any):\n\
  +p    report also min site-site distance (time-consuming)\n\
  +mMAX print all frames with coordinates out of range [-MAX,MAX]\n\
        MAX>0: quit with error code=1 if out of range\n\
        MAX<0: as above with |MAX|, do not quit if out of range\n\
        MAX=0: do not check range of coordinates\n\
  +f    print number of frames (useful in scripts)\n\
  +s    print size in B\n\
  +n[#] print number of sites [integer-divided by #]\n\
  +F +S +N as above prepended by filename\n\
  +u    print list of plb utilities\n\
  +t    print technical info on plb-files\n\
FILES:\n\
  - in front of file name: reverse endian of binary file\n\
  - instead of file name = stdin (can use only once); -- = rev.endian\n\
    (useful on 32 bit architectures to check files longer than 2GiB)\n\
Return value:\n\
  number of errors encountered (0 on success)\n\
See also:\n\
  - plbinfo +u for the the list of plb utilities\n\
  - plbinfo +t for technical info on plb-files\n\
Example (bash):\n\
  export NFRAMES=`plbinfo +f water999.plb`\n\
");
    return 0; }

  if (arg[1][0]=='+') {
    switch (arg[1][1]) {
      case 'p': pair++; break;
      case 'F': printfn=1;
      case 'f': key='f'; break;
      case 'S': printfn=1;
      case 's': key='s'; break;
      case 'N': printfn=1;
      case 'n':
        key='n';
        nms=atoi(arg[1]+2);
        if (nms<=0) nms=1;
        break;
      case 'm': m=atof(arg[1]+2); break;
      case 'u':
        fprintf(stderr,"\
MACSIMUS playback files are binary files with standard extension .plb.\n\
Run `plbinfo +t' for file format.\n\
Conversions of plb-files:\n\
  plb2asc asc2plb plb2cfg cfg2plb crd2plb gaussian2plb plbconv plbpak plbbox\n\
  plbfill wat2wat m2m plb2flt3 crd2plb\n\
Filtering:\n\
  plbframe plbcut plbsites plb2plb plbfilt plbstack plbmerge plbsplit plbdist\n\
  plbmolc plbaround wat2wat\n\
Transformations:\n\
  plbshift plbscale plbtran plbrot plbreplicate plbrev plbcenter plb2box.sh\n\
  plbsmooth\n\
Calculations:\n\
  plb2cryst plb2dens plb2diff plb2hist plb2nbr plb2sqd plbmsd plb2tree\n\
  densprof plbdist plb2dist plbionpair plbqprof hbonds\n\
More:\n\
  plbcheck\n\
");
        return 0;
      case 't':
        fprintf(stderr,"\
The MACSIMUS plb-file is a binary file consisting of 4-byte floats:\n\
# header 2 floats (8 bytes)\n\
  o float 1 = number of sites (the maximum is ns=16777216)\n\
  o float 2 =\n\
    * 0e0: free (3D Cartesian, vacuum) boundary conditions; \"old format\"\n\
    * L: (for L>0) periodic b.c., the box is a cube and L is constant\n\
    * -3e0: variable box size; all atom positions are positive\n\
    * -6e0: as above, NOT FULLY IMPLEMENTED, see cook -n4\n\
            half box size subtracted from all coordinates\n\
# frame 1\n\
  o L[]    (12 bytes) for variable box size format only:\n\
                      the actual box sizes (x,y,z)\n\
  o r[0][] (12 bytes) coordinates of the 1st atom (site) in the configuration\n\
    ...\n\
  o r[ns-1][] coordinates of the last atom\n\
* frame 2\n\
  ...\n\
(plbinfo +u for the list of utilities)\n\
");
        return 0;
      default:
        fprintf(stderr,"plbinfo: bad option\n");
        return 0; }
    iarg0=2; }

  loop (iarg,iarg0,narg) {
    FILE *plb;
    float hdr[2];
    floatvector r,r0,r9,L0,L9;
    int j,ns,reverse=0,is;
    char *fn;
    int varL=0;
    floatvector *cfg,L={0,0,0};
    double minrr=9e99;
    long unsigned int n=0,L0n[3],L9n[3], r0n[3],r9n[3];
    int bad=0;

    fn=arg[iarg]+(reverse=arg[iarg][0]=='-');

    if (!key) printf("file %s: ",fn);
    if (!strcmp(arg[iarg],"-")) {
      reverse=0;
      plb=stdin; }
    else if (!strcmp(arg[iarg],"--")) {
      reverse=1;
      plb=stdin; }
    else
      plb=fopen(fn,"rb");
    if (!plb) {
      fprintf(stderr,"plbinfo: cannot open %s (too long?)\n",arg[iarg]);
      ret++;
      continue; }

    if (fread(hdr,sizeof(hdr),1,plb)!=1) {
      fprintf(stderr,"plbinfo: file too short\n");
      ret++;
      goto Close; }

    if (reverse) {
      endian((char*)hdr);
      endian((char*)(hdr+1)); }
    ns=hdr[0];
    if (!key) printf("%d sites, ",ns);

    if (printfn) printf("%s: ",arg[iarg]);

    if (key=='n') {
      printf("%d\n",ns/nms);
      goto Close; }

    if (key=='f' || key=='s') {
      size_t size;
      rewind(plb);
      fseek(plb,0,SEEK_END);
      size=ftell(plb);
      if (key=='s') printf("%ld\n",size);
      else printf("%ld\n",size/(12*(ns+(hdr[1]<0))));
      goto Close; }

    if (hdr[1]>=0) {
      varL=0;
      printf("fixed box size L=%.8g (old format)\n",hdr[1]); }
    else if (hdr[1]==-3) {
      varL=1;
      printf("variable box (new format), hdr[1]=-3\n"); }
    else if (hdr[1]==-6) {
      varL=1;
      printf("variable box with centered configuration, hdr[1]=-6 (NOT FULLY SUPPORTED)\n"); }
    else if (hdr[1]<0) {
      varL=1;
      printf("variable box with centered configuration, hdr[1]=%g (UNSUPPORTED)\n",hdr[1]); }
    else {
      printf("WRONG PARAMETER L=%g\n",hdr[1]);
      ret++;
      continue; }

    if (ns<=0 || ns>16777216) {
      fprintf(stderr,"plbinfo: wrong number of sites (bad endian or not a playback file)\n");
      ret++;
      continue; }

    if (pair) cfg=malloc(ns*sizeof(floatvector));

    loop (j,0,3) r0[j]=L0[j]=9e9,r9[j]=L9[j]=-9e9;

    while (fread(r,sizeof(floatvector),1,plb)==1) {
      is=n%(ns+varL);
      if (reverse) loop (j,0,3) endian((char*)(r+j));
      if (varL && is==0) {
        loop (j,0,3) {
          L[j]=r[j];
	  if (L0[j]>r[j]) L0[j]=r[j],L0n[j]=n;
	  if (L9[j]<r[j]) L9[j]=r[j],L9n[j]=n; } }
      else
        loop (j,0,3) {
          if (cfg) cfg[is-varL][j]=r[j];
          if (r0[j]>r[j]) r0[j]=r[j],r0n[j]=n;
          if (r9[j]<r[j]) r9[j]=r[j],r9n[j]=n; }

      if (m)
        if (fabs(r[0])>=fabs(m) || fabs(r[1])>=fabs(m) || fabs(r[2])>=fabs(m)) bad++;

      if ((is==ns+varL-1) && bad) {
        printf("%4d: coord or L out of range %d times\n",(int)(n/(ns+varL)+1),bad);
        if (m>0) return 1;
        bad=0; }

      if (cfg && is==ns+varL-1) {
        int i,j,k;
        int irr=-1,jrr=-1;
        double x,rr=0;

        loop (i,0,ns) loop (j,0,i) {
          rr=0;
 	  loop (k,0,3) {
	    x=cfg[i][k]-cfg[j][k];
	    if (L[k]) {
	      while (x<-L[k]/2) x+=L[k];
	      while (x>L[k]/2) x-=L[k]; }
	    rr+=x*x; }
	  if (rr<minrr) {
	    irr=i; jrr=j; minrr=rr; } }
        printf("%3g: minr=%.7f i=%-4d j=%-4d (L=%g %g %g)\n",
	       (double)(n+1)/(ns+varL),sqrt(minrr),irr,jrr,L[0],L[1],L[2]);
        irr=-1,jrr=-1;
        minrr=3e33; }
      n++; }

    if (cfg) free(cfg);
    cfg=NULL;

    printf("%.14g configurations (frames), %ld vectors\n",(double)n/(ns+varL),n);
    is=n%(ns+varL);
    if (is) {
      printf("WARNING: last frame not complete: %d vectors\n",is);
      ret++; }
    if (varL) {
      printf("min L = %10.6f %10.6f %10.6f  f=%d %d %d\n",
        L0[0],L0[1],L0[2],
        (int)(L0n[0]/(ns+varL)+1),(int)(L0n[1]/(ns+varL)+1),(int)(L0n[2]/(ns+varL)+1));
      printf("max L = %10.6f %10.6f %10.6f  f=%d %d %d\n",
        L9[0],L9[1],L9[2],
        (int)(L9n[0]/(ns+varL)+1),(int)(L9n[1]/(ns+varL)+1),(int)(L9n[2]/(ns+varL)+1)); }
#define SF(N) (int)(N/(ns+varL)+1),(int)(N%(ns+varL)-varL)
    printf("min r = %10.6f %10.6f %10.6f  f.s=%d.%d %d.%d %d.%d\n",
      r0[0],r0[1],r0[2],
      SF(r0n[0]),SF(r0n[1]),SF(r0n[2]));
    printf("max r = %10.6f %10.6f %10.6f  f.s=%d.%d %d.%d %d.%d\n",
      r9[0],r9[1],r9[2],
      SF(r9n[0]),SF(r9n[1]),SF(r9n[2]));

    printf("range = %10.6f %10.6f %10.6f  (frame=1,2,... site=0,1,...)\n",
            r9[0]-r0[0],r9[1]-r0[1],r9[2]-r0[2]);

  Close:
    fclose(plb); }

  return ret;
}
