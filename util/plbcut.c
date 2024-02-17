/* cc -o plbcut -O2 plbcut.c
(c) J.Kolafa 11/1995-
05/2008: quiet mode
06/2005: new format with L[3]
10/2000: arbitrary selection added
         header[1]=L
08/1997: header[1] printed as t0
         bug fixed (FROM did not work correctly)
         -BY added
*/

#include "../gen/include.h"

int fromby=1; /* 1: uses from,by; 0: uses ranges */
int by,from=1;

struct rg_s {
  struct rg_s *next;
  int from,to; } *rg0,*rg;

int towrite(int n)
{
  if (fromby)
    return n>=from && (n-from)%by==0;
  else {
    for (rg=rg0; rg; rg=rg->next) 
      if (n>=rg->from && n<=rg->to) return 1;
    return 0; }
}

int main(int narg,char **arg)
{
  FILE *f1,*f2,*chk=NULL;
  int ny,nn,n,i,ns,nsf,dump=0,site=0,ident=0,incl,omitted=0;
  static float header[2];
  char *ch;
  typedef float vector[3];
  vector *r;
  struct old_s {
    struct old_s *next;
    int n;
    vector r; } *old0=NULL,*old;
  int varL=0,verbose=1;
  vector L;

  if (narg<4) {
    fprintf(stderr,"Extract selected frames from a MACSIMUS playback file.\n\
- Frames (configurations) are numbered from 1\n\
- Output plb format is the same as the input format\n\
Extract with regular stride:\n\
  plbcut {INPUT|-} {OUTPUT|-} BY [FROM]\n\
    BY = extract each BY-th frame\n\
    FROM = first frame extracted, default=1\n\
Extract arbitrary frames, as union of given ranges:\n\
  plbcut {INPUT|-} {OUTPUT|-} 0 FROM[:TO] [FROM[:TO] ...]\n\
    TO = last frame written, incl., default=FROM, very large=until EOF\n\
Check coordinates of SITE, remove repeating frames, and generate plbcut.chk:\n\
  (useful after similation crash and restart, cf. restart.sh)\n\
  plbcut {INPUT|-} {OUTPUT|-} -BY [SITE]\n\
    SITE = site number\n\
Environment:\n\
  PLBCUT (any value): quiet mode (bash: export PLBCUT=1, tcsh: setenv PLBCUT 1)\n\
Example:\n\
  mv simul.plb simul.plb~\n\
  plbcut - simul.plb 0 1:1000 < simul.plb~\n\
See also:\n\
  plbinfo plbfilt plbsites plbframe plb2plb plbmerge plbpak plbtran plbrot\n");
   exit(0); }

  by=atoi(arg[3]);
  if (by && narg>5) Error("too many parameters");

  if (strcmp(arg[1],"-") && !strcmp(arg[1],arg[2])) 
    Error("INPUT_FILE = OUTPUT_FILE");

  f1=strcmp(arg[1],"-")?fopen(arg[1],"rb"):stdin;
  if (!f1) Error("no INPUT_FILE");

  if (getenv("PLBCUT")) verbose=0;

  if (by==0) {
    fromby=0;
    loop (i,4,narg) {
      allocone(rg);
      rg->from=rg->to=atoi(arg[i]);
      if ( (ch=strchr(arg[i],':')) ) rg->to=atoi(ch+1);
      if (rg->to<rg->from)
	fprintf(stderr,"WARNING %d:%d bad or empty range\n",rg->from,rg->to);
      rg->next=rg0; rg0=rg; }
    if (!rg0) fprintf(stderr,"WARNING no range specified\n"); }

  if (by<0) { dump++; by=-by; }

  if (by && narg>4)
    if (dump)
      site=atoi(arg[4]);
    else {
      from=atoi(arg[4]);
      if (from<1) Error("illegal FROM"); }

  fread(header,sizeof(header),1,f1);
  ns=header[0];

  if (header[1]<0)
    if (header[1]<0) varL=1;
    else Error("wrong header (L<0, L!=-3)");

  printf("# of sites=%d  size of 1 frame=%d  L=%f (%s)",
	 ns,(int)(ns*sizeof(vector)),header[1],varL?"variable box":"old format");
  if ((header[0]-ns)!=0 || ns<1 || header[0]>16777216.) Error("bad header");
 
  if (site<0 || site>=ns) Error("illegal SITE");
 
  f2=strcmp(arg[2],"-")?fopen(arg[2],"wb"):stdout;
  if (dump) chk=fopen("plbcut.chk","wt");
  if (!f2) Error("cannot open OUTPUT_FILE");
 
  if (chk) fprintf(chk,"# ns=%d site %d\n\
#frame      x         y         z  [identical to frame(s)]\n",
		   ns,site);

  fwrite(header,sizeof(header),1,f2);
 
  allocarray(r,ns);
 
  n=ny=nn=0;
  for (;;) {
    if (verbose) {
      if (n%50==0) printf("\n%5d:",n);
      if (n%10==0) putchar(' '); }
    n++;
    if (varL) { if (fread(L,sizeof(vector),1,f1)<1) { nsf=0; break; } }
    nsf=fread(r,sizeof(vector),ns,f1);
    if (nsf<ns) break;
    if (towrite(n)) {
      incl=1;
      if (chk) {
	fprintf(chk,"%4d %9.5f %9.5f %9.5f",
		n,r[site][0],r[site][1],r[site][2]);
	for (old=old0; old; old=old->next) 
	  if (!memcmp(old->r,r[site],sizeof(vector))) {
	    ident++; incl=0;
	    fprintf(chk," %d", old->n); }
	fprintf(chk,"\n");
	allocone(old);
	copy(old->r,r[site],sizeof(vector));
	old->n=n;
	old->next=old0; old0=old; }
      if (incl) {
	if (varL) fwrite(L,sizeof(vector),1,f2);
	nsf=fwrite(r,sizeof(vector),ns,f2);
	if (nsf!=ns) Error("write");
	ny++;
	if (verbose) putchar('W'); }
      else {
	if (verbose) putchar('X');
	omitted++; } }
    else {
      if (verbose) putchar('-');
      nn++; } }
 
  if (chk) fclose(chk);
  fclose(f2);
  fclose(f1);
 
  if (nsf) fprintf(stderr,"\nWARNING last cfg only %d sites - ignored\n",nsf);
 
  printf("\n%d skipped(-)  %d written(W)  %d repeating omitted(X)\n",nn,ny,omitted);
 
  if (ident)
    printf("WARNING: %d frames omitted (%d identical pairs)\n\
         watch column(s) 4,5,... in file plbcut.chk!\n",omitted,ident);

  return 0;
}
