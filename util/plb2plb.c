/* cc -o plb2plb -Wall -O2 plb2plb.c
*/

#include "../gen/include.h"

typedef float vector[3];

int ns;
int from,to,by,excl;
void fromtoby(char *arg)
/* get from:to:by, returns 1 if range is empty */
{
  char *c;
  
  excl=arg[0]=='-';
  arg+=excl;
  from=atoi(arg);
  to=from+1;
  by=1;
  if ((c=strchr(arg,':'))) {
    to=atoi(++c);
    if ((c=strchr(c,':'))) by=atoi(++c); }
  if (from<0) Error("plb2plb: FROM<0");
  if (to<0) to+=ns;
  if (by<=0) Error("plb2plb: BY<=0 or missing");
  if (to>ns) to=ns;
}

int main(int narg,char **arg)
{
  FILE *f1,*f2;
  int i,iarg,fromitem;
  static float header[2];
  vector *r;
  int *yes;
  int varL=0;
  vector L;
  int sort=0;
  int newns=0;
  int iframe,fromframe=1,toframe=0x7fffffff;
  int m=0,warn=0;

  if (narg<3) {
    fprintf(stderr,"\
Extract selected sites from a playback file. Call by:\n\
  plb2plb INPUT_FILE OUTPUT_FILE [OPTIONS] RANGE [[-]RANGE ..]\n\
OPTIONS:\n\
  -m     merge RANGEs (subtract -RANGEs) and then write all selected sites in\n\
         the original order (default)\n\
  -o     write sites in the order of RANGEs specified (sites may repeat)\n\
  -fFROM from frame, default=1 (first frame)\n\
  -tTO   to frame (incl.) (default=FROM if -fFROM precedes, otherwise EOF)\n\
  -iN    write info for first N sites (default=0)\n\
Frames are numbered from 1, frame TO is included in the range (FORTRAN-style)\n\
RANGE = \n\
  FROM:TO    range of sites (site TO is not included); too large TO = NS\n\
  FROM:-TO   the same as NS+TO (= TO from end)\n\
  FROM:0     empty range\n\
  FROM:[-]TO:BY sites FROM, FROM+BY, FROM+BY*2, ..., <TO; BY>0 required\n\
  SITE       given site\n\
  -RANGE     exclude given range (with option -m only)\n\
Sites are numbered from 0 to NS-1, TO is not included in the range (C-style)\n\
If no ITEM is given, NS=0 and only box sizes are selected\n\
Example - the following statements are equivalent:\n\
  plb2plb in.plb out.plb 0:6 -1:5 2\n\
  plb2plb in.plb out.plb 0 2 5\n\
See also:\n\
  plbinfo plbsites plbconv plbmerge plbcut plbfilt\n\
  m2m wat2wat plb2asc asc2plb plb2cfg cfg2plb mol2mol\n");
    exit(0); }

  if (!strcmp(arg[1],arg[2])) Error("plb2plb: INPUT_FILE = OUTPUT_FILE");

  f1=fopen(arg[1],"rb");
  if (f1==NULL) Error("plb2plb: no INPUT_FILE");

  f2=fopen(arg[2],"wb");
  if (f2==NULL) Error("plb2plb: cannot open OUTPUT_FILE");

  loop (iarg,3,narg) 
    if (arg[iarg][0]=='-') 
      switch (arg[iarg][1]) {
        case 'u': sort=0; break;
        case 'o': sort=1; break;
        case 'f': fromframe=toframe=atoi(arg[iarg]+2); break;
        case 't': toframe=atoi(arg[iarg]+2); break;
        case 'i': m=atoi(arg[iarg]+2); break;
        default: Error("plb2plb: bad option"); }
    else
      break;
  fromitem=iarg;
  
  if (1!=fread(header,sizeof(header),1,f1)) Error("plb2plb: no header");
  ns=header[0];

  if (header[1]<0) {
    if (header[1]<0) varL=1;
    else Error("plb2plb: wrong header (L<0, L!=-3)"); }

  printf("plb2plb: # of sites=%d  frame size=%dB  L=%f (%s)\n",
	 ns,(int)(ns*sizeof(vector)),header[1],varL?"variable box":"old format");
  if ((header[0]-ns)!=0 || ns<1 || header[0]>16777216.) Error("plb2plb: bad header");

  allocarray(r,ns);
  allocarrayzero(yes,ns);
 
  loop (iarg,fromitem,narg) {
    fromtoby(arg[iarg]);
    for (i=from; i<to; i+=by) {
      if (sort) {
        if (excl) Error("plb2plb: cannot exclude range in the -m mode");
        yes[i]++; }
      else {
        yes[i]=1; 
        if (excl) yes[i]=0; } } }
  
  loop (i,0,ns) {
    newns+=yes[i];
    warn+=yes[i]>1;
    if (ns<m) {
      putchar('.');
      if (yes[i]<10) putchar('0'+yes[i]);
      else putchar('X'); } }
  if (warn) printf("WARNING: %d sites repeat more than once\n",warn);

  header[0]=newns;
  printf("plb2plb: %d site(s) selected\n",newns);
  fwrite(header,sizeof(header),1,f2);

  loopto (iframe,1,toframe) {
    if (varL) {
      if (!fread(L,sizeof(vector),1,f1)) break;
      if (iframe>=fromframe) fwrite(L,sizeof(vector),1,f2); }
  
    if (ns!=fread(r,sizeof(vector),ns,f1)) break;
    
    if (sort) {
      if (iframe>=fromframe) loop (iarg,fromitem,narg) {
        fromtoby(arg[iarg]);
        for (i=from; i<to; i+=by)
          fwrite(r+i,sizeof(vector),1,f2); } }
    else  {
      if (iframe>=fromframe) 
        loop (i,0,ns) if (yes[i]) fwrite(r+i,sizeof(vector),1,f2); } }
  
  fclose(f2);
  fclose(f1);

  return 0;
}
