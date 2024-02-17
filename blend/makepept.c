/* cc -O2 -o makepept makepept.c
*/

#include "../gen/include.h"

#define LINE 256
#define LINES 16

char *blendpath=".",*rsddir;
int nerr;
int mergeonly;

char *amino(char code) /********************************************** amino */
{
  char msg[]="? is invalid 1-letter aminoacid code";

  if (isdigit(code)) return "1"; /* count, not amino */

  switch (code) {
    case 'A': return "ala";
    case 'C': return "cys";
    case 'D': return "asp";
    case 'E': return "glu";
    case 'F': return "phe";
    case 'G': return "gly";
    case 'H': return "his";
    case 'I': return "ile";
    case 'K': return "lys";
    case 'L': return "leu";
    case 'M': return "met";
    case 'N': return "asn";
    case 'P': return "pro";
    case 'Q': return "gln";
    case 'R': return "arg";
    case 'S': return "ser";
    case 'T': return "thr";
    case 'V': return "val";
    case 'W': return "trp";
    case 'Y': return "tyr";
    default: msg[0]=code; Error(msg); }
}

struct rsd_s {
  int n;
  char fn[LINE];
  char line[LINES][LINE]; } A,B;

char *zgets(char *line,FILE *f) /************************************* zgets */
{
  if (!fgets(line,LINE,f)) return NULL;

  if (strlen(line)>0 && line[strlen(line)-1]=='\n') line[strlen(line)-1]=0;
  if (strlen(line)>0 && line[strlen(line)-1]=='\r') line[strlen(line)-1]=0;
  if (strlen(line)>0 && line[strlen(line)-1]=='\n') line[strlen(line)-1]=0;

  return line;
}

void readrsd(struct rsd_s *R,char *r) /***************************** readrsd */
{
  FILE *f;
  static int parset=0;
  int comment=1;
  char line0[LINE];

  sprintf(R->fn,"%s/%s/%s.che",blendpath,rsddir,r);

  f=fopen(R->fn,"rt");
  if (!f) {
    strcat(R->fn,": not found (check BLENDPATH)");
    Error(R->fn); }

  zgets(line0,f);
  fprintf(stderr,"%s\n",line0);
  R->n=0;

  while (zgets(R->line[R->n],f)) if (R->line[R->n][0]!='!') {
      if (!memcmp(R->line[R->n],"parameter_set",13))
        R->n-=parset;
      else if (comment) {
        strcat(R->line[R->n],strlen(R->line[R->n])?" ! ":"! ");
        strcat(R->line[R->n],line0);
        comment=0; }
      if (++R->n>=LINES) {
        strcat(R->fn,": too many lines");
        Error(R->fn); } }

  if (!R->n) {
    strcat(R->fn,": empty residue");
    Error(R->fn); }

  parset=1;

  fclose(f);
}

char *tok[9];
int gettok(char *s) /************************************************ gettok */
{
  int i=0;
  char *t=strtok(s," \t");

  t=strtok(NULL," \t");
  if (!t) return 0;
  tok[0]=strtok(t,",");
  for (i=1;;) {
    t=strtok(NULL,",");
    if (!t) break;
    if (i>=9-1) return 0;
    tok[i++]=t; }

  return i;
}

void paste(char *what,char *into) /*********************************** paste */
{
  char *e=strchr(what,'!');
  char *ei=strchr(into,'!');
  char eil[LINE];

  if (ei) {
    strcpy(eil,ei);
    memset(ei,0,LINE-strlen(ei)); }

  if (e) *e++=0;
  if (strlen(into)<strlen(what)) strcpy(into,what);
  else memcpy(into,what,strlen(what));
  if (e && strlen(e)>1) {
    strcat(into," !");
    strcat(into,e); }
  if (ei) strcat(into,eil);
}

void rsd(char *r,char *mpl) /******************************************* rsd */
/*
  copies residuum r mpl-times
  if mpl starts with letter, then once
  ignored if r starts with digit
*/
{
  int a,b,nt,t,i,Bfrom=0,m=1;
  char *s;
  static char last[LINE];

  if (isdigit(*r)) return;

  if (mpl && *mpl && isdigit(*mpl)) m=atoi(mpl);
  if (m<0 || m>1000) Error("negative or too big residue count");

  while (m--) {

    if (mergeonly) {
      readrsd(&A,r);
      loop (a,0,A.n) printf("%s\n",A.line[a]);
      A.n=0;
      continue; }

    if (!memcmp(last,"nter",4) && (!strcmp(r,"gly") || !strcmp(r,"pro"))) {
      fprintf(stderr,"ERROR %s after nter:\n\
  use glyp[n] or prop[n] instead\n",r);
      nerr++; }
    strcpy(last,r);

    if (A.n==0) {
      readrsd(&A,r);
      continue; }
    readrsd(&B,r);

    loop (a,0,A.n) if ( (s=strstr(A.line[a],"!>>>")) ) {
      nt=gettok(s);
      if (!nt) {
        strcat(A.fn," none or too many forward matching groups");
        Error(A.fn); }
      *s=0;
      loop (b,0,B.n)
        loop (t,0,nt) if (strstr(B.line[b],tok[t])) {
        if (a>=b) {
          loop (i,a+1,A.n) if (!strchr(A.line[i],'|'))
            fprintf(stderr,"WARNING %s:\n\
  data after !>>> ignored",A.fn);
          while (b>=0)
            paste(A.line[a--],B.line[b--]);
          A.n=a+1; }
        else {
          strcat(A.fn,": sorry, I cannot handle shorter patch...");
          Error(A.fn); }
        goto cont; }

      fprintf(stderr,"ERROR %s:\n\
  forward matching group not found\n",A.fn);
      nerr++; }

  cont:
    loop (b,0,B.n) if ( (s=strstr(B.line[b],"!<<<")) ) {
      nt=gettok(s);
      if (!nt) {
        strcat(B.fn," none or too many backward matching groups");
        Error(B.fn); }
      *s=0;
      for (a=A.n-1; a>=0; a--) {
        loop (t,0,nt) if (strstr(A.line[a],tok[t])) {
          if (B.n-b>=A.n-a) {
            loop (i,0,b) if (!strchr(B.line[i],'|'))
              fprintf(stderr,"WARNING %s:\n\
  data before !<<< ignored",A.fn);
            while (a<A.n)
              paste(B.line[b++],A.line[a++]);
            Bfrom=b; }
          else {
            strcat(A.fn,": sorry, I cannot handle shorter patch...");
            Error(A.fn); }
          goto finish; } }

      fprintf(stderr,"ERROR %s:\n\
  backward matching group not found\n",B.fn);
      nerr++; }

  finish:
    loop (a,0,A.n) printf("%s\n",A.line[a]);
    A.n=0;
    loop (b,Bfrom,B.n) strcpy(A.line[A.n++],B.line[b]);
    B.n=0;
  }
}

int main(int narg,char **arg) /**************************************** main */
{
  int iarg;
  char *c;

  if (narg<2) {
    fprintf(stderr,"\
Simple MACSIMUS peptide tool.\n\
Print aminoacid table with 1 and 3 letter codes:\n\
  makepept ANYPARM\n\
Make a peptide of given residues in MACSIMUS che-format:\n\
  makepept [-]RSDDIR {rsd [COUNT] | RRR} [rsd [COUNT] | RRR ...]\n\
where\n\
  RSDDIR = subdirectory of BLENDPATH with residues in che-format\n\
  -RSDDIR = as above, only merges .che files (does not solve patches)\n\
  rsd = residue or terminus name in lowercase\n\
  COUNT = optional count (applies to previous residue)\n\
  RRR = chain of one-letter aminoacid codes (uppercase) and count numbers\n\
Example:\n\
  makepept charmm22 ace A3PT hisn ct1 > pept.che\n\
will make peptide Acetyl-ALA-ALA-ALA-PRO-THR-HIS(neutral)-Methyl\n\
Recommended blend commands:\n\
  blend -e40 pept.che     # initialize, optimize, show\n\
  blend -e40 -r2 pept.che # to run again (overwrite pept.mol,pept.plb)\n");
   exit(0); }

  mergeonly=arg[1][0]=='-';
  rsddir=arg[1]+mergeonly;
  blendpath=getenv("BLENDPATH");
 
  if (narg==2) {
    puts("\
   A ALA alanine\n\
   C CYS cysteine\n\
   D ASP aspartic acid\n\
   E GLU glutamic acid\n\
   F PHE phenylalanine\n\
   G GLY glycine\n\
   H HIS histidine\n\
   I ILE isoleucine\n\
   K LYS lysine\n\
   L LEU leucine\n\
   M MET methionine\n\
   N ASN asparagine\n\
   P PRO proline\n\
   Q GLN glutamine\n\
   R ARG arginine\n\
   S SER serine\n\
   T THR threonine\n\
   V VAL valine\n\
   W TRP tryptophan\n\
   Y TYR tyrosine\n\
makepept FORCEFIELD ace ACDEFGHIKLMNPQRSTVWY ct3 > test.che");
    exit(0); }

  fputs("peptide",stdout);
  loop (iarg,2,narg) printf(" %s",arg[iarg]);
  _n

  loop (iarg,2,narg)
    if (isupper(arg[iarg][0])) for (c=arg[iarg]; *c; c++) rsd(amino(*c),c+1);
    else rsd(arg[iarg],arg[iarg+1]);

  loop (iarg,0,A.n) printf("%s\n",A.line[iarg]);

  if (nerr)
    fprintf(stderr,"%d ERROR%s in the output che-file: it needs editing\n",nerr,nerr==1?"":"S");

 return 0;
}
