/* cc -O2 -o asc2plb asc2plb.c
*/

#include "../gen/include.h"

typedef float floatvector[3];

void endian(char *a) /********************************************** endian */
{
  char x;

  x=a[0]; a[0]=a[3]; a[3]=x;
  x=a[1]; a[1]=a[2]; a[2]=x;
}

static int nextferr;
static char *tok;
double nextf(void) /************************************************** nextf */
{
  double ret=0;

  if (tok) {
    nextferr=0;
    ret=atof(tok);
    tok=strtok(NULL," \t\n\r,"); }
  else
    nextferr=1;

  return ret;
}

int main(int narg,char **arg) /*************************************** main */
{
  FILE *plb,*pla;
  float hdr[2];
  floatvector r;
  int j,n=0,ns=0,reverse=0,iarg,ask=1,atm=-1,atmnsline,readxyz,noheader=0,nv=0;
  char *infn=NULL,*outfn=NULL;
  char line[128];
  char *ext="plb";

  if (narg<2) {
    fprintf(stderr,"\
Convert text coordinate file to binary playback file.  Call by:\n\
  asc2plb ASC-FILE PLB-FILE [OPTIONS]\n\
  asc2plb NAME [OPTIONS]\n\
File arguments:\n\
  ASC-FILE  text image of playback file, usual .EXT=.pla,.vla; -=stdin\n\
            see plbinfo for the format\n\
            with -a: ATM-file, 2 lines before every xyz block:\n\
                     1. NS\n\
                     2. info or box\n\
  PLB-FILE  MACSIMUS playback file, usual .EXT=.plb,.vlb; -=stdout\n\
  NAME      ASC-FILE=NAME.pla, PLB-FILE=NAME.plb\n\
Options:\n\
  -a0       input is ATM-format, 2nd header line=info, box=(0,0,0)\n\
  -a1       .., 2nd header line='box L', box=(L,L,L)\n\
  -a3       .., 2nd header line='box Lx Ly Lz', box=(Lx,Ly,Lz)\n\
  -eEXT     output EXT (2nd form)\n\
  -r        reverse endian of PLB-FILE\n\
  -y        overwite without asking\n\
  -n        no header on input: by lines .xyz (.3dt) -> .3db\n\
See also:\n\
  plb2asc plbconv plbmerge plbcut plbinfo plbfilt\n\
Deprecated old version: plbasc\n\
");
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-' && arg[iarg][1])
      switch (arg[iarg][1]) {
        case 'a': atm=atoi(arg[iarg]+2); break;
        case 'e': ext=arg[iarg]+2; break;
        case 'r': reverse++; break;
        case 'y': ask=0; break;
        case 'n': noheader++; break;
        default: Error("asc2plb: wrong option"); }
    else {
      if (infn) Error("asc2plb: more than 2 file names");
      infn=outfn;
      outfn=arg[iarg]; }

  if (!outfn) Error("asc2plb: no file/name");

  if (!infn) {
    char *name=outfn;

    if (strlen(ext)>3) Error("asc2plb: EXT loger than 3 chars");
    infn=malloc(strlen(name)+5);
    sprintf(infn,"%s.pla",name);
    outfn=malloc(strlen(name)+5);
    sprintf(outfn,"%s.%s",name,ext); }

  if (strcmp("-",infn)) pla=fopen(infn,"rt");
  else pla=stdin;
  if (!pla) Error(infn);
  fprintf(stderr,"%s (ascii) -> %s (binary)\n",infn,outfn);

  if (strcmp("-",outfn)) {
    plb=fopen(outfn,"rb");
    if (plb && ask) {
      char line[256];
      fprintf(stderr,"File %s exists. Overwrite (y/n)? ",outfn);
      gets(line);
      if (tolower(line[0])!='y') return 1;
      fclose(plb); }
    plb=fopen(outfn,"wb"); }
  else
    plb=stdout;

  if (noheader) {
    /* .3dt (.xyz) -> .3db */
    if (atm>=0) Error("asc2plb: cannot combine options -a and -n");
    while (fgets(line,128,pla)) {
      tok=strtok(line," \t\n\r,");
      r[0]=nextf(); r[1]=nextf(); r[2]=nextf();
      if (nextferr) {
        fprintf(stderr,"asc2plb: noheader: line incomplete\n");
        goto end; }
      nv++;
      if (reverse) loop (j,0,3) endian((char*)(r+j));
      fwrite(r,sizeof(floatvector),1,plb); }
    goto end; }

  else if (atm>=0)
    /* .atm format */
    for (;;) {
      if (!fgets(line,128,pla)) goto end;

      /* ns */
      sscanf(line,"%d",&j);
      if (ns && j!=ns) {
        fprintf(stderr,"asc2plb: ATM: ns cannot change between frames\n");
        goto end; }
      if (!ns) {
        ns=j;
        hdr[0]=ns; hdr[1]=-3;
        fprintf(stderr,"ns=%d\n",ns);
        if (reverse) { endian((char*)hdr); endian((char*)(hdr+1)); }
        fwrite(hdr,sizeof(float),2,plb); }

      loop (j,-1,ns) {
        if (!fgets(line,128,pla)) {
          fprintf(stderr,"ATM: missing line %d (-1=info/box), frame incomplete",j);
          goto end; }

        readxyz=1;
        /* atom or keyword box */
        if (j==-1)
          switch (atm) {
            case 0:
              r[0]=r[1]=r[2]=0;
              readxyz=0;
              break;
            case 1:
              tok=strtok(line," \t\n\r,");
              if (tok) tok=strtok(NULL," \t\n\r,");
              if (!tok) {
                fprintf(stderr,"ATM: missing info in 2nd header line");
                goto end; }
              r[0]=r[1]=r[2]=atof(tok);
              readxyz=0;
              break; }

        if (readxyz) {
          tok=strtok(line," \t\n\r,"); /* skip atom/box */
          if (tok) tok=strtok(NULL," \t\n\r,");
          r[0]=nextf(); r[1]=nextf(); r[2]=nextf();
          if (nextferr) {
            fprintf(stderr,"asc2plb: ATM: line %d incomplete\n",j);
            goto end; }
          nv++;
          if (reverse) loop (j,0,3) endian((char*)(r+j)); }

        fwrite(r,sizeof(floatvector),1,plb); }
        n++; }

  else {
    /* PLA format */
    if (!fgets(line,128,pla)) Error("asc2plb: PLA: missing header");

    sscanf(line,"%d%f",&ns,&hdr[1]);
    fprintf(stderr,"ns=%d L=%f\n",ns,hdr[1]);

    hdr[0]=ns;
    if (reverse) { endian((char*)hdr); endian((char*)(hdr+1)); }
    fwrite(hdr,sizeof(float),2,plb);

    while (fgets(line,128,pla)) {
      tok=strtok(line," \t\n\r,");

      if (!tok)  {
        fprintf(stderr,"plb2asc: PLA: line incomplete\n");
        goto end; }
      r[0]=nextf(); r[1]=nextf(); r[2]=nextf();
      if (nextferr) {
        fprintf(stderr,"asc2plb: PLA: line %d incomplete\n",j);
        goto end; }
      nv++;
      if (reverse) loop (j,0,3) endian((char*)(r+j));

      fwrite(r,sizeof(floatvector),1,plb); }
    n=nv/(ns+(hdr[1]==-3)); }

 end:

  fprintf(stderr,"%d vectors read (incl. box), %d configuration(s) written\n", nv, n);

  fclose(plb);
  fclose(pla);

  return 0;
}
