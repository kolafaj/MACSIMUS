/* cc -O2 -o plb2asc plb2asc.c
*/

#include "../gen/include.h"

typedef float floatvector[3];

void endian(char *a) /********************************************** endian */
{
  char x;

  x=a[0]; a[0]=a[3]; a[3]=x;
  x=a[1]; a[1]=a[2]; a[2]=x;
}

int main(int narg,char **arg) /*************************************** main */
{
  FILE *plb,*pla;
  float hdr[2];
  floatvector r,L;
  int j,n=0,ns=0,reverse=0,is=0,varL,iarg,ask=1,NS;
  enum mode_e {PLA,ATM,COM} mode=PLA; /* asc-image, xyz+header, gaussian input */
  char *infn=NULL,*outfn=NULL;
  char *FMT="%9.6f",*atoms=NULL,*ext="pla";
  char *boxopt=" ";
  char *at;

  if (narg<2) {
    fprintf(stderr,"\
Convert playback (and similar) file to text.  Call by:\n\
  plb2asc PLB-FILE ASC-FILE [OPTIONS]\n\
  plb2asc NAME [OPTIONS]\n\
File arguments:\n\
  PLB-FILE  MACSIMUS playback file, usual .EXT=.plb,.vlb; -=stdin\n\
  ASC-FILE  its text image, usual .EXT=.pla,.vla,.atm; -=stdout\n\
  NAME      PLB-FILE=NAME.plb, ASC-FILE=NAME.pla\n\
Options:\n\
  -eEXT     set EXT\n\
  -r        reverse endian of PLB-FILE\n\
  -y        overwite without asking\n\
  -%%FMT\n\
  -f%%FMT    number format (must accept double), default=%s\n\
  -nNS      input is .3db (no header), NS=number of sites (deprecated)\n\
The default output is an ascii image of PLB-FILE. Modificators:\n\
  -aATOMS   header of NS+box (repeats before every frame), EXT=.atm\n\
            atom names need not be separated (AlCl=Al,Cl), repeated if needed\n\
  -a        as above, reads NAME.mol for atoms (1 letter, only second form)\n\
  -g        the same as -a but no header, EXT=.asc\n\
            to add box, use -l' '\n\
            useful for Gaussian input (do not forget to append empty line)\n\
  -l-       omit box info line [default=-l\" \" (plain box info)]\n\
  -l        blank line instead of box info\n\
  -lHDR     use HDR to start line with box info (e.g., -lbox)\n\
Example (print x,y,z block of last 25 atoms in format compatible with PDB)\n\
  plb2asc simul.plb - \"-%%7.3f\" | tail -n 25\n\
See also:\n\
  asc2plb plbconv plbmerge plbcut plbinfo plbfilt view3db\n\
Deprecated old version: plbasc\n\
",FMT);
    exit(0); }

  loop (iarg,1,narg)
    if (arg[iarg][0]=='-' && arg[iarg][1])
      switch (arg[iarg][1]) {
        case '%': FMT=arg[iarg]+1; break;
        case 'f': FMT=arg[iarg]+2; break;
        case 'e': ext=arg[iarg]+2; break;
        case 'r': reverse++; break;
        case 'y': ask=0; break;
        case 'l': boxopt=arg[iarg]+2; break;
        case 'a': atoms=arg[iarg]+2; ext="atm"; mode=ATM; break;
        case 'g': atoms=arg[iarg]+2; ext="asc"; mode=COM; boxopt="-"; break;
        case 'n': ns=atoi(arg[iarg]+2); break;
        default: Error("plb2asc: wrong option"); }
    else {
      if (infn) Error("plb2asc: more than 2 file names");
      infn=outfn;
      outfn=arg[iarg]; }

  if (!outfn) Error("plb2asc: no file/name");
  if (!infn) {
    char *name=outfn;

    if (strlen(ext)>3) Error("plb2asc: EXT loger than 3 chars");
    infn=malloc(strlen(name)+5);
    outfn=malloc(strlen(name)+5);

    if (atoms && atoms[0]==0) {
      FILE *mol;
      char line[256],*tok;
      int n=0,i;

      sprintf(infn,"%s.mol",name);
      mol=fopen(infn,"rt");
      if (!mol) Error("plb2asc: mol-file not found");

      while (fgets(line,256,mol)) {
        if (!memcmp(line,"number_of_atoms",15)) {
          char *c=line;

          while (!isdigit(*c)) c++;
          n=atoi(c); 
          break; } }

      if (n<1 || n>16777216) Error("plb2asc: bad number_of_atoms in mol-file");
      atoms=malloc(n);

      while (fgets(line,256,mol)) 
        if (!memcmp(line,"atoms",5)) break;

      loop (i,0,n) {
        do {
          if (!fgets(line,256,mol)) Error("plb2asc: unexpected EOF in mol-file");
        } while (line[0]=='!');
        tok=strtok(line," \t");
        tok=strtok(NULL," \t");
        tok=strtok(NULL," \t");
        if (!tok) Error("plb2asc: not enough data in mol-file");
        atoms[i]=tok[0]; }
      fclose(mol); }
    sprintf(infn,"%s.plb",name);
    sprintf(outfn,"%s.%s",name,ext); }

  if (strcmp("-",infn)) plb=fopen(infn,"rb");
  else plb=stdin;
  if (!plb) Error(infn);
  fprintf(stderr,"%s (binary) -> %s (ascii)\n",infn,outfn);

  if (strcmp("-",outfn)) {
    pla=fopen(outfn,"rt");
    if (pla && ask) {
      char line[256];

      getsbufsize=256;
      fprintf(stderr,"File %s exists. Overwrite (y/n)? ",outfn);
      gets(line);
      if (tolower(line[0])!='y') return 1;
      fclose(pla); }
    pla=fopen(outfn,"wt"); }
  else
    pla=stdout;

  if (!ns) {
    if (fread(hdr,sizeof(hdr),1,plb)!=1) Error("plb2asc: file too short");
      if (reverse) {
	endian((char*)hdr);
	endian((char*)(hdr+1)); }
      ns=hdr[0];
      fprintf(stderr,"ns=%d L=%f\n",ns,hdr[1]);
      if (hdr[0]<1. || hdr[0]>16777216.) Error("plb2asc: bad ns in hdr"); }

  varL=hdr[1]<0;
  if (!varL) L[0]=L[1]=L[2]=hdr[1];
  NS=ns+varL;

  if (mode==PLA) {
    fprintf(pla,"%d ",ns);
    fprintf(pla,FMT,hdr[1]);
    fprintf(pla,"\n"); }

  while (fread(r,sizeof(floatvector),1,plb)==1) {

    if (is==0) {
      /* header */
      if ( (at=atoms) && mode==ATM ) {
        fprintf(pla,"%d ",ns);
        fprintf(pla,"\n"); }

      /* box info line */
      if (!varL) {
        /* L=r */
        if (!strcmp(boxopt,"")) goto blank; /* empty line */
        if (!strcmp(boxopt,"-")) goto omit; /* no box line */
        if (mode!=COM) fprintf(pla,"%s",boxopt); /* will write L=r */ }
      else {
        if (!strcmp(boxopt,"")) {
          if (mode!=COM) fprintf(pla,"\n"); /* empty line */
          goto omit; }
        else if (!strcmp(boxopt,"-")) {
          goto omit;
          /* no box line */ }
        else {
          /* write L */
          if (mode!=COM) fprintf(pla,"%s",boxopt); } } }
    else {
      if (atoms) {
        if (*at==0) at=atoms;
        if (*at) {
          fprintf(pla,"%c",*at++);
          if (islower(*at)) fprintf(pla,"%c",*at++);
          fprintf(pla," ");
          if (*at && !isalnum(*at)) at++; } }
      else
        fprintf(pla," "); }

    loop (j,0,3) {
      if (reverse) endian((char*)(r+j));
      if (j || !is) fputc(' ',pla);
      fprintf(pla,FMT,r[j]); }
   blank:
    fprintf(pla,"\n");
   omit:
    is=(is+1)%NS;
    n++; }

  fprintf(stderr,"%d vectors read, %g configuration(s) written\n",n,(double)n/NS);

  fclose(plb);
  fclose(pla);

  return 0;
}
