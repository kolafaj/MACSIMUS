#include "ground.h"
#include "bitfile.h"
#include "pakcp.h"

/* 12/2019: cleaned, reverse endian removed */

/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Reading packed CP.  Use as follows:

  FILE *pakcp=fopen(NAME,"rb");
  int NCP;
  float *r=getcpheader(pakcp,&NCP); // r[NCP] gets allocated

  // now r[NCP] is the file header

  while (getpakcp(pakcp)) {
    while (nextcprec(r)) {
      // now r[NCP] is one record (records with CPmark included!)
    }
    endgetpakcp(pakcp); 
  }
*/

static int efread(void *rec,size_t size,size_t n,FILE *f)
{
  int ret=fread(rec,size,n,f);

  if (size!=4) Error("unsupported architecture");

  return ret;
}

static int NCP; /* used for decompressing routines only */
static int verbose=1; /* verbose=1 requests info printed; set by getcpheader */

float *getcpheader(FILE *f,int *NCPp,int v) /******************* getcpheader */
/*
  reads 1st record in f (packed or unpacked)
  returns:
    *NCPp and (static) NCP
    the record float[NCP] is allocated and read
*/
{
  float hdr;
  int4 ncp;
  float *r; /* [NCP] allocated and returned */

  verbose=v;

  if (1!=fread(&hdr,4,1,f)) ERROR(("cp-file too short"))
  if (hdr!=CPmark) ERROR(("bad header of cp-file: CPmark requested"))

  if (1!=fread(&ncp,4,1,f)) ERROR(("cp-file too short"))

  /* legacy: in versions for both endianess, ncp was coded twice */
  ncp=ncp&0x0000ffff;
  NCP=*NCPp=ncp;
  if (verbose) prt("ncp=%d columns found in the cp-file",NCP);
  allocarray(r,ncp);
  r[0]=hdr;
  *(int4*)(r+1)=ncp; /* as int */
  if (fread(r+2,sizeof(float),NCP-2,f)!=NCP-2) ERROR(("cp-file too short"))

  return r;
}

/* append packed convergence profile to (opened) file f */
int4 appendpakcp(
  FILE *f,    /* file must be opened for writing in binary mode */
  int NCP,    /* # of CP items recorded; WARNING: shades static NCP */
  int nbit,   /* resolution is 2^nbit values in the [min,max] range */
  void *data, /* if norec=0 then the head of rec_t* list
                 if norec>0 then pointer to r[norec][NCP] (one contig.array) */
  int4 norec) /* data format */
{
  int i,bcut,ibcut;
  float *minr,*maxr,*r;
  int4 ir,dir,*ir0,*len,*minlen,no;
  int4 *bestbcut; /* int 4 because written - compatibility! */
  rec_t *rec;

  if (!data) return 0;

  allocarray(minr,NCP);
  allocarray(maxr,NCP);

  allocarray(ir0,NCP);
  allocarrayzero(bestbcut,NCP);
  allocarray(len,NCP);
  allocarrayzero(minlen,NCP);

  loop (i,0,NCP) {
    minr[i]=3e33;
    maxr[i]=-3e33;
    minlen[i]=0x7fffffffL; }

  for (no=0, rec=(rec_t*)data;
       norec ? no<norec : rec!=NULL;
       no++, rec = norec ? NULL : rec->next) {
    r = norec ? (float*)data+no*NCP : rec->r;

    loop (i,0,NCP) {
      Min(minr[i],r[i])
      Max(maxr[i],r[i]) } }

  if (verbose) prt_("no=%d   ",no);
  assignbit(f,'w');

  /* to improve! */
  loop (i,0,NCP) if (maxr[i]!=minr[i]) {
    float x=maxr[i];
    int shift=nbit;

    do maxr[i]+=(maxr[i]-minr[i])/(1<<shift);
    while (shift-- && maxr[i]==x);

    if (++shift!=nbit)
      if (verbose) prt("! %d : [%.9g, %.9g->%.9g] (shift=%d)",
                       i,minr[i],x,maxr[i],shift); }

  if (fwrite(minr,sizeof(float),NCP,f)!=NCP) Error("cannot write");
  if (fwrite(maxr,sizeof(float),NCP,f)!=NCP) Error("cannot write");
  if (fwrite(&no,sizeof(int4),1,f)!=1) Error("cannot write");

  /* this is for all writes not through bitput */
  /*.....bitfile.bytecount = 2*sizeof(float)*NCP + sizeof(int4)*(2+NCP);*/

  for (ibcut=nbit-1; ibcut>=1; ibcut--) {
    memset(ir0,0,NCP*sizeof(int4));

    if (ibcut==1) {
      if (verbose) prt_("bitcut: ");
      loop (i,0,NCP) {
        if (maxr[i]==minr[i]) bestbcut[i]=0;
        if (verbose) prt_("%2d ",bestbcut[i]); }
      if (verbose) prt("/%3d",nbit);

      if (fwrite(bestbcut,sizeof(int4),NCP,f)!=NCP) Error("cannot write");

      ir=(int4)nbit;
      if (fwrite(&ir,sizeof(int4),1,f)!=1) Error("cannot write"); }

    bcut=ibcut;
    memset(len,0,NCP*sizeof(int4));

    for (no=0, rec=(rec_t*)data;
         norec ? no<norec : rec!=NULL;
         no++, rec = norec ? NULL : rec->next) {
      r = norec ? (float*)data+no*NCP : rec->r;

      loop (i,0,NCP) if (maxr[i]!=minr[i]) {
        if (ibcut==1) bcut=(int)bestbcut[i];
        ir=(int4)((r[i]-minr[i])/(maxr[i]-minr[i])*(1<<nbit)+0.5);

        /* fix rounding problems ... */
        if (ir<0 || ir>=(1<<nbit)) {
          if (ir<0) ir=0;
          if (ir>=(1<<nbit)) ir=(1<<nbit)-1; }

        dir=ir-ir0[i]+(1<<(bcut-1));
        if (dir>=0 && dir<(1<<bcut)) {
          len[i]+=1+bcut;
          if (ibcut==1) {
            putbit(1,1); putbit(dir,bcut); } }
        else {
          len[i]+=1+nbit;
          if (ibcut==1) {
            putbit(0,1); putbit(ir,nbit); } }
        ir0[i]=ir; } }

    if (ibcut>1) {
      loop (i,0,NCP) {
        if (len[i]<minlen[i]) {
          minlen[i]=len[i];
          bestbcut[i]=bcut; } } }
  }

  flushbit();
  fwrite(&bitfile.checksum,sizeof(bitfile.checksum),1,f);
  if (verbose) prt("bitfile.bytecount=%ld",(long)bitfile.bytecount);

  free(minlen);
  free(len);
  free(bestbcut);
  free(ir0);

  free(maxr);
  free(minr);

  return bitfile.bytecount;
}

static float *minr,*maxr,*markrec,*markrec2,*markrec3;
static int4 *ir,*bestbcut;
static int bcut,nbit;
static int4 no;

int nextcprec(float *r) /***************************************** nextcprec */
/*
  gets next record of packed CP file
  returns 1 on success,
  0 if eof or CPmark encountered
*/
{
  int i;

  if (markrec) {
    copy(r,markrec,NCP*sizeof(float));
    free(markrec);
    markrec=markrec2;
    markrec2=markrec3;
    markrec3=NULL;
    return 1; }

  if (--no<0) return 0;

  loop (i,0,NCP)
    if (maxr[i]-minr[i]) {
      bcut=(int)bestbcut[i];
      if (getbit(1)) ir[i]+=getbit(bcut)-(1<<(bcut-1));
      else ir[i]=getbit(nbit);
      r[i]=minr[i]+(float)ir[i]/(1<<nbit)*(maxr[i]-minr[i]); }
    else
      r[i]=minr[i];

  return 1;
}

int getpakcp(FILE *f) /******************************************** getpakcp */
/*
  begins reading one CP series (`no' records, i.e. between CPmark's)
  returns 0 on EOF
*/
{
  int i;
  int4 ibc;

  assignbit(f,'r');
  if (feof(f)) return 0;

  allocarray(minr,NCP);
  allocarray(maxr,NCP);

  allocarrayzero(bestbcut,NCP);
  allocarrayzero(ir,NCP);

  if (efread(minr,sizeof(float),NCP,f)!=NCP) return 0;

  while (minr[0]<=CPmark) {
    if (markrec) {
      if (markrec2) {
        if (markrec3) Error("to many marked records in a series");
        else markrec3=minr; }
      else markrec2=minr; }
    else markrec=minr;
    allocarray(minr,NCP);
    if (efread(minr,sizeof(float),NCP,f)!=NCP) return 1; }

  if (efread(maxr,sizeof(float),NCP,f)!=NCP) Error("cannot read");

  /*.....loop (i,0,NCP) prt("%2d %g %g",i,minr[i],maxr[i]);*/

  loop (i,0,NCP) {
    if (maxr[i]<minr[i]) Error("max<min");
    if (fabs(maxr[i])>3e33 || fabs(minr[i])>3e33) Error("too big number"); }

  if (efread(&no,sizeof(no),1,f)!=1) Error("cannot read");

  if (efread(bestbcut,sizeof(int4),NCP,f)!=NCP) Error("cannot read");

  if (efread(&ibc,sizeof(int4),1,f)!=1) Error("cannot read");
  nbit=(int)ibc;

  if (verbose) {
    prt_("bitcut: ");
    loop (i,0,NCP) prt_("%2d ",bestbcut[i]);
    prt("/%3d  no=%d",nbit,no); }

  /* WARNING: this is symptomatic ad hoc bug fix ...
     there is rather something wrong in WRITING -- extra 4B */
  if (no==1) efread(&i,sizeof(i),1,f);

  return 1;
}

void endgetpakcp(FILE *f) /************************************* endgetpakcp */
/*
  ends reading one CP series (`no' records, i.e. between CPmark's)
*/
{
  unsigned4 checksum=0;

  flushbit();

  if (bitfile.bytecount) {
    if (fread(&checksum,sizeof(checksum),1,f)!=1) ERROR(("missing checksum in .cpz file"))

    if (checksum!=bitfile.checksum)
      ERROR(("checksum %lu, expected %lu",
             (long unsigned) bitfile.checksum,checksum))
    if (verbose) prt("bitfile.bytecount=%ld",(long)bitfile.bytecount); }

  if (bestbcut) {
    free(bestbcut);
    free(ir);
    free(maxr);
    free(minr); }
}
