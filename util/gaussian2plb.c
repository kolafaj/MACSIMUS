/* cc -o gaussian2plb gaussian2plb.c -lm
 */
#include "../gen/include.h"
typedef double real;
#include "../gen/vector3d.h"

double ftof(char *c) /************************************************* ftof */
/* get fortran-style number (accepts 9.9D+9 and 9.9d+9) */
{
  char nr[128],*cc;

  while (*c && *c<=' ') c++;

  if (!*c) return 0;
  cc=nr;
  while (*c>' ') {
    if (toupper(*c=='D')) *cc++='e';
    else *cc++=*c;
    c++; }
  *cc=0;

  return atof(nr);
}

double field(char *line,char *key,int n) /**************************** field */
/* get n-th fortran number after key (or BOL id key="") */
{
  char *c=line;

  if (key && key[0]) {
    c=strstr(line,key);
    if (c) line=c+strlen(key); }

  for (;;) {
    if (!*line) return 0;
    while (strchr(" \t",*line)) line++;
    if (!*line) return 0;
    n--;
    if (!n) return ftof(line);
    while (!strchr(" \t",*line)) line++; }
}

int main(int narg,char **arg) /**************************************** main */
{
  FILE *in;
  char *fn,*dot;
  FILE *plb=NULL,*atm=NULL,*nmff=NULL;
  float r[3];
  char line[256];
  int i,n,N=0;
  vector *center;
  vector **dr;
  double *nmf; /* [cm-1] */
  double A=0.5;
  int isfreq=0;

  if (narg<2) {
    fprintf(stderr,"\
Normal mode vibrations from Gaussian output. Call by:\n\
  gaussian2plb NAME[.log] [AMPLITUDE]\n\
where AMPLITUDE is the normal mode vibration amplitude [default=%f AA]\n\
- The Gaussian output file NAME.log is read\n\
- NAME-gauss.atm is created (hint: use showatm|bonds to watch)\n\
- NAME-gauss.nmf with frequencies is created (excl. first 6 singular)\n\
- NAME-gauss.plb with energy minimization trajectory is generated;\n\
  if NAME.mol exists, show is called with coordinate matcher:\n\
    show NAME NAME-opt.plb -mNAME.plb:1:NAME.keep -I$\n\
- Files nm####.plb are generated;\n\
  if NAME.mol exists, show is called with vibration viewer:\n\
    show NAME nm%%04d.plb -Iww%%\n\
Bug: cannot be used for linear molecules\n",A);
    exit(0); }

  if (narg>2) A=atof(arg[2]);

  fn=malloc(strlen(arg[1])+11);
  strcpy(fn,arg[1]);
  if (strlen(fn)>4) 
    if (!strcmp(fn+strlen(fn)-4,".log")) fn[strlen(fn)-4]=0;
  dot=fn+strlen(fn);
  strcpy(dot,".log");
  in=fopen(fn,"rt");
  if (!in) Error(fn);

  /* reading NAtoms: it may appear AFTER configuration !!! */
  while (fgets(line,256,in)) 
    if (strstr(line,"NAtoms=")) {
      int newN=field(line,"NAtoms=",1);

      if (N) {
        if (N!=newN) Error("NAtoms changed"); }
      else {
        N=newN;
        if (plb) Error("plb already open");
        strcpy(dot,"-gauss.plb");
        plb=fopen(fn,"wb");
        r[0]=N; r[1]=0; /* old plb-file */
        put(N)
        fwrite(r,4,2,plb);
        allocarray(center,N); 
        allocarray(dr,3*N-6);
        loop (n,0,3*N-6) allocarray(dr[n],N);
        allocarray(nmf,3*N-6); } }

  if (!N) Error("keyword \"NAtoms=\" is missing");

  rewind(in);

  while (fgets(line,256,in)) 
    if (strstr(line,"Z-Matrix orientation:") 
        || strstr(line,"Standard orientation")) {
      strcpy(dot,"-gauss.atm");                                               
      atm=fopen(fn,"wt");
      fprintf(atm,"%d\n\n",N);

      fgets(line,256,in);
      fgets(line,256,in);
      fgets(line,256,in); 
      fgets(line,256,in); 
      if (line[3]!='-') Error("format");

      loop (i,0,N) {
        fgets(line,256,in);
        /*
    1          1             0        1.919787   -0.003953    0.753920
        */
        center[i][0]=field(line,"",4);
        center[i][1]=field(line,"",5);
        center[i][2]=field(line,"",6);
        fprintf(atm,"%.0f %10.6f %9.6f %9.6f\n",field(line,"",2),center[i][0],center[i][1],center[i][2]);
        VV(r,=center[i])
        fwrite(r,4,3,plb); } }

  fclose(plb);

  strcpy(dot,".mol");
  plb=fopen(fn,"rt");
  if (plb) {
    *dot=0;
    system(string("show %s %s-gauss.plb -m%s.plb:1:%s.keep \'-|\' \'-I$\'",
      fn,fn,fn,fn));
    fclose(plb); }
  
  rewind(in);
  n=0;

  while (fgets(line,256,in))
    if (strstr(line,"Frequencies --")) {
      isfreq++;
      if (n>=3*N-6) Error("to many NMFs");

      /*
                     1                      2                      3
                     A                      A                      A
 Frequencies --   127.9989               147.9907               150.7587
 Red. masses --     1.0709                 1.3622                 1.0323
 Frc consts  --     0.0103                 0.0176                 0.0138
 IR Inten    --   112.2665               160.1418                59.7591
 Atom AN      X      Y      Z        X      Y      Z        X      Y      Z
   1   1     0.00   0.63   0.00    -0.23  -0.01   0.17     0.00   0.52   0.00
   2   1     0.00   0.20   0.00     0.07   0.00  -0.37     0.00  -0.32  -0.01
   3   8     0.00  -0.06   0.00     0.14   0.00   0.00     0.00   0.00   0.00
   4   1    -0.37   0.06   0.38    -0.53   0.00   0.32     0.49   0.25   0.09
   5   1     0.36   0.07  -0.37    -0.52  -0.01   0.33    -0.50   0.25  -0.08
   6   8     0.00   0.00   0.00    -0.06   0.00  -0.03     0.00  -0.04   0.00
NEW VERSION:
     1   6     0.00   0.00  -0.02    -0.12  -0.16   0.00     0.27   0.14   0.00
      */
      nmf[n]  =field(line,"--",1);
      nmf[n+1]=field(line,"--",2);
      nmf[n+2]=field(line,"--",3);
      while (fgets(line,256,in)) if (strstr(line," Atom ")) break;

      loop (i,0,N) {
        fgets(line,256,in);
        dr[n][i][0]=field(line,"",3);
        dr[n][i][1]=field(line,"",4);
        dr[n][i][2]=field(line,"",5);

        dr[n+1][i][0]=field(line,"",6);
        dr[n+1][i][1]=field(line,"",7);
        dr[n+1][i][2]=field(line,"",8);

        dr[n+2][i][0]=field(line,"",9);
        dr[n+2][i][1]=field(line,"",10);
        dr[n+2][i][2]=field(line,"",11); }

      n+=3; }

  put2(n,3*N-6)

  strcpy(dot,"-gauss.nmf");
  nmff=fopen(fn,"wt");
  fprintf(nmff,"# normal mode frequencies from Gaussian\n\
#\n\
#  i      [THz]         1/cm]\n\
");
    
  if (isfreq) loop (n,0,3*N-6) {
    char *fn=string("nm%04d.plb",n);
    int k;

    plb=fopen(fn,"wb");
    r[0]=N; r[1]=0; /* old plb-file */
    fwrite(r,4,2,plb);

    printf("%s: %g\n",fn,nmf[n]);
    fprintf(nmff,"%4d %13.6f %13.4f\n",n,nmf[n]*0.0299792458,nmf[n]);

    loopto (k,0,6) {
      double q=A*cos(k*PI/6);

      loop (i,0,N) {
        VVV(r,=center[i],+q*dr[n][i])
        fwrite(r,4,3,plb); } }
    fclose(plb); }
  fclose(nmff);
  
  *dot=0;
  plb=fopen(string("%s.mol",fn),"rt");
  if (plb) {
    system(string("SHOWNMF=\"%s-gauss.nmf\" show %s nm%%04d.plb -d50 \'-|\' \'-I$ww\'",fn,fn));
    fclose(plb); }

  return 0;
}
