/* make plbpak

  12/2017: destination path (option -p) added
  08/2010: options improved, multiple files supported
  08/2008: version with -GRID
  02/2007: updated to variable-L format

  Compress/decompress of playback files (with headers!)
  1. Reduces precision to 1/GRID
  2. Uses relative coding of consecutive atom positions whenever possible.
     The idea is that typical configurations have consecutive atoms
     often close together (because bonded).

  Now a standalone program, but written in a way which allows easy use
  of procedures appendpakplb() and getpakplb() in programs
*/

#include "ground.h"
#include "bitfile.h"
#include <unistd.h>

#define VO(A,O) { A[0] O; A[1] O; A[2] O; }
#define VV(A,B) { A[0] B[0]; A[1] B[1]; A[2] B[2]; }

typedef float floatvector[3];
int verbose=1;

int appendpakplb( /******************************************* appendpakplb */
  FILE *f,        /* output packed playback - file must be opened */
  int ns,         /* # of sites incl. L[] for new plb format */
  int outtype,    /* 0: old fixed-L (deprecated), 1: new variable-L */
  floatvector *r, /* array float[ns][3], if outtype=1 then starts with L[] */
  int grid)       /* grid per AA (grid=10 means resolution 0.1 AA) */
/*
  appends one configuration r[ns] of playback in the packed format to file f
  # of bytes written is returned
*/
{
  int4 minr[3],maxr[3],ir[3],irj,ir0[3],dr[3],bcutrg,bcutrgh,minnbits,nbits;
  int i,j,rel,bcut,bits[3],maxbits,ibcut,bestbcut;

  if (ns<=0) Error("ns<=0");

  VO(minr,=0x7fffffffL)
  VO(maxr,=-0x7fffffffL)

  assignbit(f,'w');

  loop (i,0,ns)
    loop (j,0,3) {
      irj=(int4)(r[i][j]*(double)grid+16777216.5)-16777216L;
      Min(minr[j],irj)
      Max(maxr[j],irj) }

  maxbits=0;
  loop (j,0,3) {
    bits[j]=ceillog2(maxr[j]-minr[j]+1);
    Max(maxbits,bits[j]) }
  if (maxbits<4) {
    static int pass;
    if (!pass) WARNING(("too low resolution or coordinates the same (more warnings suppressed)"))
    pass=1; }

  minnbits=0x7fffffffL;
  nbits=0;
  for (ibcut=maxbits-1; ibcut; ibcut--) {
    if (ibcut!=1)
      /* updated 8/97: */
      if (ibcut<maxbits-3 && nbits==ns*(bits[0]+bits[1]+bits[2])) continue;

    VO(ir0,=0)
    nbits=0;
    bcut=ibcut==1?bestbcut:ibcut;
    bcutrg=1L<<bcut;
    bcutrgh=bcutrg/2;
    if (!bcutrgh) Error("cutoff");

    if (ibcut==1) {
      bitfile.check=1;
      if (outtype)	{
	putbit(0,20);
	putbit(ns,24); }
      else
	putbit(ns,20);
      putbit(grid,12);
      putbit(bcut,8);
      if (bitfile.err) Error("write"); 

      loop (j,0,3) {
	if (minr[j]<-(1L<<23) || minr[j]>(1L<<23)-1)
	  Error("minr out of range");
	putbit(minr[j]&(1L<<24)-1,24); /* 2 complement if negative */
	putbit(bits[j],8); }

      if (bitfile.err) Error("write"); }

    bitfile.check=1; /* =0 later */
    loop (i,0,ns) {
      rel=1;
      loop (j,0,3) {
	ir[j]=(int4)(r[i][j]*(double)grid+16777216.5)-16777216L-minr[j];
	dr[j]=ir[j]-ir0[j]+bcutrgh;
	if (dr[j]<0 || dr[j]>=bcutrg) rel=0; }

      if (ibcut==1) {
	putbit(rel,1); /* 1st bit denotes whether the vector will be relatively
			  to the previous */
	if (rel) loop (j,0,3) putbit(dr[j],bcut); /* relative */
	else loop (j,0,3) putbit(ir[j],bits[j]);  /* absolute */ }
      else
	nbits += rel ? 3*bcut : bits[0]+bits[1]+bits[2];

      VV(ir0,=ir) }

#if 0
    /* verbose: report of optimization of the number of bits */
    if (ibcut==1) prt("=> %d",bcut);
    else prt_("%d:%ld ",bcut,(long)nbits);
#endif

    if (nbits<minnbits) {
      minnbits=nbits; 
      bestbcut=bcut; } }

  flushbit();

  return bitfile.bytecount;
}

int getpakplb( /************************************************* getpakplb */
  FILE *f,         /* input packed playback - file must be opened */
  int *nsp,        /* # of sites returned here */
  int *intype,     /* 0=old type (no L), 1=new type (L present) */
  floatvector **rp)/* pointer to array float[ns][3] returned here */
/*
  reads one configuration of playback in the packed format from file f
  array r to hold them is allocated - don't forget to free it!
  # of (packed) bytes read is returned
*/
{
#define ns (*nsp)
#define r (*rp)
  int4 minr[3],ir[3],bits[3],bcutrg,bcutrgh,grid;
  int i,j,bcut;

  VO(ir,=0)

  assignbit(f,'r');

  ns=getbit(20);
  if (bitfile.err) return 0; /* EOF */

  if (ns==0) {
    *intype=1;
    ns=getbit(24); }
  else 
    *intype=0;

  if (ns<=0) Error("ns<=0");

  alloc(r,(int4)ns*sizeof(floatvector));

  grid=getbit(12);
  bcut=getbit(8);
  loop (j,0,3) {
    minr[j]=getbit(24);
    if (minr[j] & (1L<<23)) minr[j]-=(1L<<24);
    bits[j]=getbit(8); }

  if (verbose) prt_("ns=%d grid=%d  ",ns,grid);

  bcutrg=1L<<bcut;
  bcutrgh=bcutrg/2;
  if (!bcutrgh) Error("");

  if (bitfile.err) Error("read");

  loop (i,0,ns) {
    if (getbit(1))
      loop (j,0,3) ir[j]+=getbit(bcut)-bcutrgh; /* relative position */
    else
      loop (j,0,3) ir[j]=getbit(bits[j]); /* absolute position */
    loop (j,0,3) r[i][j]=(double)(ir[j]+minr[j])/grid; }

  flushbit(); /* next data from byte boundary (flush) */
  return bitfile.bytecount;
}
#undef ns
#undef r

int main(int narg,char **arg) /*************************************** main */
{
  char nm[256];
  FILE *plb,*pakplb;
  char *c,*dot,*path=NULL;
  float hdr[2];
  int ns,i,iarg;
  unsigned frame;
  int grid=32,plbtype=0,plztype=1;
  int erase=0,overwrite=0,ignore=0;
  floatvector *r,L;
  char line[128];
  enum mode_e { AUTO,COMPRESS,DECOMPRESS } mode=AUTO,dir=AUTO;


  if (narg<2) {
    fprintf(stderr,"\n\
Lossy compression of MACSIMUS playback files (12/2017).  Call by:\n\
  plbpak [OPTIONS] FILE.EXT [[OPTIONS] FILE.EXT ..]\n\
FILES:\n\
  FILE.plb = MACSIMUS playback file\n\
  FILE.plz = packed playback file\n\
OPTIONS (apply to the rest of arguments)\n\
  -gGRID   grid (resolution) in points per Angstrom (default=%d/AA}\n\
  -a       auto mode: compress FILE.plb, decompress FILE.plz (default)\n\
  -d       force decompress\n\
  -c       force compress\n\
  -i       ignore ERROR \"unexpected EOF\" (i.e., frame not complete)\n\
  -o       overwrite the target without warning (default = ask)\n\
  -e       erase the original after (de)compression\n\
  -pPATH   destination directory (default=./)\n\
  -q       quiet (default = verbose info about every frame)\n\
  -y       shorthand for -o -e -q -i\n\
Compatibility (see also old version plbpak0):\n\
  old plz-file (w/o L) is converted to an old plb-file with L=0 (cf. plbbox)\n\
  any plb-file is converted to a new plz-file with L information\n\
The maximum coordinate is +-2^23/GRID; if a compression fails with range error,\n\
  use plbinfo to check the maximum coordinates in the plb-file\n\
See also:\n\
  plbinfo plbcheck plb2asc asc2plb plbconv plb2plb plbbox\n\
  showplz\n",grid);
  exit(0); }

  initscroll(0);

  loop (iarg,1,narg) 
    if (arg[iarg][0]=='-')
      switch (arg[iarg][1]) {
        case 'g': grid=atoi(arg[iarg]+2); break;
        case 'a': mode=AUTO; break;
        case 'c': mode=COMPRESS; break;
        case 'd': mode=DECOMPRESS; break;
        case 'i': ignore=1; break;
        case 'p': path=arg[iarg]+2; break;
        case 'q': verbose=0; break;
        case 'e': erase=1; break;
        case 'y': erase=1; verbose=0; ignore=1;
        case 'o': overwrite=1; break; 
        default: ERROR(("%s: unnown option",arg[iarg])) }
    else {
      frame=0;
      strcpy(nm,arg[iarg]);
      /* determine extension */
      dot=NULL;
      for (c=nm; *c; c++) if (*c=='.') dot=c;
      dir=mode;
      if (mode==AUTO) {
        if (!dot) ERROR(("%s: no extension and none of -c -d specified",arg[iarg])) 
        if (!strcmp(dot,".plb")) dir=COMPRESS;
        if (!strcmp(dot,".plz")) dir=DECOMPRESS;
        if (dir==AUTO) ERROR(("%s: unknown extension and none of -c -d specified",arg[iarg]))  }

      if (dir==COMPRESS) {

        /*** compression ***/

        prt("compressing %s, grid=%d/AA...",nm,grid);
        plb=fopen(nm,"rb");
        if (!plb) Error("no such file");
        if (dot) strcpy(dot,".plz"); else strcat(nm,".plz");
        if (!overwrite) {
          pakplb=fopen(nm,"rb");
          if (pakplb) {
            prt("%s exists. Overwrite (y/n)? ",nm);
            gets(line);
            if (toupper(line[0])!='Y') continue;
            fclose(pakplb); } }

        if (path) pakplb=fopen(string("%s/%s",path,nm),"wb");
        else pakplb=fopen(nm,"wb");

        if (fread(hdr,sizeof(hdr),1,plb)!=1) ERROR(("%s: too short",arg[iarg]))
        L[0]=L[1]=L[2]=hdr[1];
        plbtype=hdr[1]<0;
        ns=hdr[0];
        prt("ns=%d (%s-L format)",ns,plbtype?"variable":"fixed");
        alloc(r,(ns+1)*sizeof(floatvector));

        for (;;) {
          frame++;
          i=fread(r+1-plbtype,sizeof(floatvector),ns+plbtype,plb);
          if (plbtype==0) VV(r[0],=L)
          if (i==0) {
            flushbit();
            goto done; }
          if (i!=ns+plbtype) {
            if (ignore) prt("unexpected EOF - file truncated");
            else Error("unexpected EOF");
            goto framedone; }
          appendpakplb(pakplb,ns+1,plztype,r,grid);
          /* old plz file:	appendpakplb(pakplb,ns,plztype,r+1,grid); */

          if (verbose) prt("frame %3d: %d B packed to %ld B",
                           frame,(ns+plztype)*sizeof(floatvector),bitfile.bytecount);
         framedone:;
        } } 

      else {

        /*** decompression ***/

        prt("decompressing %s...",nm);
        pakplb=fopen(nm,"rb");
        if (!pakplb) Error("no such file");
        if (dot) strcpy(dot,".plb"); else strcat(nm,".plb");
        if (!overwrite) {
          plb=fopen(nm,"rb");
          if (plb) {
            prt("%s exists. Overwrite (y/n)? ",nm);
            gets(line);
            if (toupper(line[0])!='Y') continue;
            fclose(plb); } }

        if (path) plb=fopen(string("%s/%s",path,nm),"wb");
        else plb=fopen(nm,"wb");
        
        while (getpakplb(pakplb,&ns,&plztype,&r)) {
          ns-=plztype;
          if (!frame) {
            hdr[0]=ns;
            hdr[1]=-3*plztype; /* old type plz does not contain L */
            if (hdr[1]<0) prt("ns=%d, variable box format",ns);
            else prt(" ns=%d L=%g",ns,hdr[1]);
            fwrite(hdr,sizeof(hdr),1,plb); }
          else
            if (hdr[0]!=ns) Error("ns change");
          frame++;
          if (fwrite(r,sizeof(floatvector),ns+plztype,plb)!=ns+plztype)
            Error("cannot write");
          free(r);
          if (verbose) prt("  frame %3d: %ld B unpacked to %d B",
                           frame,bitfile.bytecount,(ns+plztype)*sizeof(floatvector));
        }
        flushbit(); } 
      
    done:
      fclose(pakplb);
      fclose(plb);
      if (erase) {
        fflush(stdout);
        remove(arg[iarg]);
        if (verbose) prt("%s erased.\n",arg[iarg]); } }

  return 0;
}
