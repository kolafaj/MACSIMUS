/* cc -O2 -o plbasc plbasc.c 

Conversion between ascii and binary playback files.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Error(X) do { fprintf(stderr,"ERROR %s\n",X); exit(1); } while(0)

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
  FILE *plb,*pla;
  float hdr[2];
  floatvector r;
  int j,n=0,ns=0,reverse=0,s=-1,is=0,prec=0,varL,NS;
  char *fn1,fn2[128],*c,*dot;
  char line[128];

  if (narg<2) {
    fprintf(stderr,"\
Conversion between ascii and binary playback files.  Call by:\n\
  plbasc [-|+]FILE.EXT [NS] [S]\n\
where\n\
  -: reverse endian of binary file\n\
  +: overwite without asking; ascii also: more dec. digits on output\n\
  EXT==plb or EXT==p00:  ascii image of binary FILE.EXT will be created\n\
  other or no EXT:  binary image FILE.plb of ascii FILE.EXT will be created\n\
  NS: # of sites if missing header, NS=0 or missing: # of sites from header\n\
  S: only output site of given # (S=-1 denotes L[3] - new format only)\n\
This is an old deprecated version! Instead use:\n\
  plb2asc asc2plb plb2plb\n");
    exit(0); }
  
  if (narg>2) ns=atoi(arg[2]);
  if (narg>3) s=atoi(arg[3]);

  fn1=arg[1]+( (reverse=arg[1][0]=='-') || (prec=arg[1][0]=='+') );
  strcpy(fn2,fn1);
  for (c=fn2, dot=NULL; *c; c++) if (*c=='.') dot=c;

  if (dot && (!strcmp(dot,".plb") || !strcmp(dot,".p00"))) {

    /* binary->ascii */
    strcpy(dot,".pla");

    printf("%s (binary) ",fn1);
    plb=fopen(fn1,"rb");
    if (!plb) Error("no such file");
    printf("-> %s (ascii)\n",fn2);
    pla=fopen(fn2,"rt");
    if (pla && !prec) {
      char line[256];

      printf("File %s exists. Overwrite (y/n)? ",fn2);
      gets(line);
      if (tolower(line[0])!='y') return 1; 
      fclose(pla); }
    pla=fopen(fn2,"wt");

    if (!ns) {
      if (fread(hdr,sizeof(hdr),1,plb)!=1) Error("file too short");
      if (reverse) {
	endian((char*)hdr);
	endian((char*)(hdr+1)); }
      ns=hdr[0];
      printf("ns=%d L=%f\n",ns,hdr[1]); }

    varL=hdr[1]<0;
    NS=ns+varL;

    if (prec) fprintf(pla,"%d %.8f\n",s<0?ns:1,hdr[1]);
    else fprintf(pla,"%d %.5f\n",s<0?ns:1,hdr[1]);
    
    while (fread(r,sizeof(floatvector),1,plb)==1) {
      if (reverse) loop (j,0,3) endian((char*)(r+j));
      if (s<0 || is-varL==s)
	if (prec) fprintf(pla,"%s%11.8f %11.8f %11.8f\n",is?"":" ",r[0],r[1],r[2]);
	else fprintf(pla,"%s%8.5f %8.5f %8.5f\n",is?"":" ",r[0],r[1],r[2]);
      is=(is+1)%NS;
      n++; }
    }

  else {
    /* ascii->binary */
    if (dot) strcpy(dot,".plb");
    else strcat(fn2,".plb");

    printf("%s (ascii) ",fn1);
    pla=fopen(fn1,"rt");
    if (!pla) Error("no such file");
    printf("-> %s (binary)\n",fn2);
    plb=fopen(fn2,"rb");
    if (plb && !prec) {
      char line[256];
      printf("File %s exists. Overwrite (y/n)? ",fn2);
      gets(line);
      if (tolower(line[0])!='y') return 1; 
      fclose(plb); }

    plb=fopen(fn2,"wb");
    if (!ns) {
      if (!fgets(line,128,pla)) Error("file too short");
      for (c=line; *c; c++) if (*c=='d' || *c=='D') *c='e';
      sscanf(line,"%d%f",&ns,&hdr[1]);
      printf("ns=%d L=%f\n",ns,hdr[1]); }

    varL=hdr[1]<0;
    NS=ns+varL;

    hdr[0]=(s<0?ns:1);
    if (reverse) {
      endian((char*)hdr); endian((char*)(hdr+1)); }
    fwrite(hdr,sizeof(float),2,plb);

    while (fgets(line,128,pla)) {
      char *tok=strtok(line," \t\n\r,");
      if (!tok) break;  
      r[0]=atof(tok); 
      tok=strtok(NULL," \t\n\r,"); if (!tok) break;  
      r[1]=atof(tok);
      tok=strtok(NULL," \t\n\r,"); if (!tok) break;  
      r[2]=atof(tok);
      n++;
      is=(is+1)%NS;
      if (reverse) loop (j,0,3) endian((char*)(r+j));
      fwrite(r,sizeof(floatvector),1,plb); }
    }

  printf("%d vectors read, %g configuration(s) written\n",n,(double)n/NS);

  fclose(plb);
  fclose(pla);

  return 0;
}
