/* cc -Wall -O2 -o cp2cp cp2cp.c
 */
#include "../gen/include.h"
#include "../sim/cpmark.h"

#define NCP 256
double d[NCP];
float  r[NCP],rr[NCP];
char name[NCP][4];
int map[NCP];
FILE *in,*out;
int inn,outn;
int nrec=0,copymarks=0;

int getrec(void) /*************************************************** getrec */
{
  int im;

  loop (im,0,9) {
    if (fread(r,4,inn,in)!=inn) return 0;
    if (r[0]>CPmark) break;
    if (copymarks) fwrite(r,4,outn,out); }
  nrec++;

  return 1;
}

int main(int narg,char **arg) /**************************************** main */
{
  int iarg,i,stride=1,block=1,from=0,to=0x7fffffff;
  int is,ib;
  char *key=NULL,*namename="-h";
  char *infn=NULL,*outfn=NULL;

  memset(name,' ',sizeof(name));

  if (narg<3) {
    fprintf(stderr,"\
Manipulate MACSIMUS convergence profile (.cp) files. Call by:\n\
  cp2cp [OPTIONS] IN.cp OUT.cp [OPTIONS]\n\
OPTIONS:\n\
  -b#\tblock by # records of data (make subaverages)\n\
  -c\tcopy all CPmark records to output (default=only header)\n\
  -f#\tread the input file from record # (records are numbered from 0)\n\
  -hNAMENAME..\theaders of empty columns, every 4 chars correspond to one 0\n\
  -n#\tmax. # of records read (must appear after -f)\n\
  -rKEY\trearrange columns, KEY=string of {1,2,..,9,A,..,0}\n\
  \twhere 1=1st column, A=10th column, 0=empty column (filled by 0)\n\
  -s#\tstride by # (take every #-th value, start with #-th)\n\
  -s# -b#  strided data are blocked\n\
  -t#\tread the input to record # (record # is not included; = sum of -f -n)\n\
BUGS:\n\
  does not accept the very old 5-column format, opposite endian\n\
  does not accept packed (.cpz) files\n\
  -c: not checked whether output NCP is long enough (to be fixed)\n\
Examples:\n\
  Omit first 1000 records, delete columns 6+7 (of 12), and block by 10:\n\
    cp2cp in.cp out.cp -f1000 -r1234589ABC -b10\n\
  Replace column 4 by 0 (header ZERO) and swap columns 6 and 7:\n\
    cp2cp in.cp out.cp -r1230576 -hZERO\n\
  Split a.cp (100000 records long) into chunks by 10000 records\n\
    tab 10000 100000 10000 \"cp2cp a.cp b%%d.cp -f%%d -n10000\" | sh\n\
See also:\n\
  showcp cppak start\n");
    exit(0); }
  
  loop (iarg,1,narg)
    if (arg[iarg][0]=='-') {
      int ii=atoi(arg[iarg]+2);

      switch (arg[iarg][1]) {
        case 'f': from=ii; break;
        case 'c': copymarks=1; break;
        case 'n': to=from+ii; break;
        case 't': to=ii; break;
        case 's': stride=ii; break;
        case 'b': block=ii; break;
        case 'r': key=arg[iarg]+2; break;
        case 'h': namename=arg[iarg]; break;
        default: Error("cp2cp: bad option"); } }
    else
      if (infn) {
	if (outfn) Error("cp2cp: more than 2 files");
	outfn=arg[iarg]; }
      else infn=arg[iarg]; 

  if (key) {
    outn=strlen(key);
    if (outn<2 || outn>NCP) Error("cp2cp: too few or too many columns on output");
    loop (i,0,outn) map[i]=key[i]-'0'-1-(key[i]>='A')*7;
    is=2;
    loop (i,0,outn) if (map[i]<0) {
      map[i]=-is;
      is+=4; }
    if (is!=strlen(namename)) Error("cp2cp: -h and -r are not consistent\n\
(number of zeros in -r must match the number of 4-char fields in -h)"); }
  else  
    loop (i,0,NCP) map[i]=i;

  in=fopen(infn,"rb"); if (!in) Error(infn);
  out=fopen(outfn,"wb"); if (!out) Error(outfn);

  /* record 1 */
  if (fread(r,4,2,in)!=2) Error("cp2cp: too short file");
  if (r[0]!=CPmark) Error("cp2cp: not a CP file (CPmark expected)");
  inn=*(int*)(r+1);
  if (!key) outn=inn;
  if (inn<2 || inn>NCP) {
    fprintf(stderr,"NCP=%d on input\n",inn);
    Error("cp2cp: too few or too many columns on input"); }
  put2(inn,outn)
  put2(from,to)

  if (inn>2) {
    if (fread(r+2,4,inn-2,in)!=inn-2) Error("cp2cp: too short file");
    rr[0]=r[0];
    *(int*)(rr+1)=outn;
    loop (i,2,outn)
      if (map[i]>=0) rr[i]=r[map[i]];
      else copy(rr+i,namename-map[i],4);
    fwrite(rr,4,outn,out); }

  loop (ib,0,from)
    if (!getrec()) Error("cp2cp: not enough data (-f option beyond length)");
  for (;;) {
    loop (i,0,outn) d[i]=0;
    loop (ib,0,block) {
      loop (is,0,stride) if (!getrec() || nrec>to) goto done;
      loop (i,0,outn)
	if (map[i]>=0) rr[i]=r[map[i]];
	else rr[i]=0;
      loop (i,0,outn) d[i]+=rr[i]; }
    loop (i,0,outn) r[i]=d[i]/block;
    fwrite(r,4,outn,out); }

 done:;
  fclose(out);
  fclose(in);

  return 0;
}
