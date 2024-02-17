/* forces by Gaussian, #included from cook/forces.c */
#if PARALLEL
#  error PARALLEL with GAUSSIAN invalid
#endif /*# PARALLEL */

#ifndef FREEBC
#  error only FREEBC allowed with GAUSSIAN
#endif /*# FREEBC */

void forces(ToIntPtr B, ToIntPtr A) /******************************** forces */
/*
   B=forces, A=configuration
*/
{
  FILE *f;
  molecule_t *mn;
  siteinfo_t *si;
  vector *r;
  int i,j,k,n,nn,tabn,ns,sp,aux,err;
  static unsigned ncom0;
  char *found,*ptr;
#define NCOM 8
  static char *com0[NCOM] /* ncom0 */,line[128];
#define NPRED 4 // ASPC(4) or 5th order predictor or ...
  static double *yC;
  static struct pred_s {
    double *y; /* [ny] orbitals */
    double omega; // damping
    double c[NPRED+2]; // c[0] = coeff for t-h, etc.
  } pred[NPRED+2] = {
    {NULL, 1,    {1}},                // [0] k=-1 (last SCF copied)
    {NULL, 2./3, {2,-1}},             // [1] ASPC k=0
    {NULL, 0.6,  {2.5,-2,0.5}},       // [2] ASPC k=1
    {NULL, 4./7, {2.8,-2.8,1.2,-0.2}}, // [3] ASPC k=2 (npred=4, index=3)
    {NULL, 5./9, {3,-24./7,27./14,-4./7,1./14}}, // [3] ASPC k=2 (npred=4, index=3)
    {NULL, 6./11,{22./7,-55./14,55./21,-22./21,5./21,-1./42}}, // [3] ASPC k=2 (npred=4, index=3)
  };
#define NHEADERS 4
  static struct headers_s {
    char *hdr;
    int pos; // from line 0 = file start
    int N;
  } headers[NHEADERS] = {
    {"Coordinates of each shell                  R",0,0},
    {"Constraint Structure                       R",0,0},
    {"Alpha MO coefficients                      R",0,0},
    {"Total SCF Density                          R",0,0}};
  static int npred=0; // @k=2 denotes ASPC(k=2); npred=k+2=4
  static int known=0; // max. "npred=k+2" known at given time 
  static int ny; /* number of orbitals */
  static int nchklines;
#define LINE 88
  static char *fchkfile; /* LINE*nchklines */

  if (com0[0]==NULL) {
    /* four lines of com-input expected, max 128 bytes, in this order:
1) # method and parms for no-measurement forces
2) # method and parms for forces and measurements
3) charge and multiplicity (as before table of atoms x y z)
4) %-command (typically %mem, 1st line)
   %-command (optional)...
5) @-command (optional, predictor) (is not included in NCOM)
   @-command (optional, predictor)...
NB: %Nproc will be printed from environment variable NSLOTS
    */
    f=fopen(Fn("com0"),"rt");
    if (!f) ERROR(("file %s is missing",lastFn))
    loop (ncom0,0,NCOM) {
      if (!fgets(line,128,f)) break;
      if (line[0]=='@') break;
      com0[ncom0]=strdup(line); }
    if (ncom0<4) ERROR(("%s: at least 3 lines expected",lastFn))
    if (com0[0][0]!='#') ERROR(("%s: 1st line must start with #",lastFn))
    if (com0[1][0]!='#') ERROR(("%s: 2nd line must start with #",lastFn))
    if (!isdigit(com0[2][0])) ERROR(("%s: 3rd line must contain two integers",lastFn))
    if (com0[3][0]!='%') ERROR(("%s: 4th line must start with %",lastFn))

    while (line[0]=='@') {
      /* predictor, format:
         @k omega c[1] c[2] .. c[k+2]
         @k omega
         @k=<max. k>
       where omega = dumping parameter,
             k = order; number of constants = k+2
       A series from k=-1 must be given up to the maximum k
       Example, giving ASPC up to k=2 with unoptimized omega
         @-1 1     1
         @0  0.666 2 -1
         @1  0.6   2.5 -2 0.5
         @2  0.571 2.8 -2.8 1.2 -0.2
       This example is equivalent to ASPC, default omega (except rounding)
         @k=2 */
      if (!memcmp(line,"@k=",3)) {
        npred=atoi(line+3)+2;
        if (npred>NPRED+2) ERROR(("%s: @k=%d too large",lastFn,npred-2)) }
      else {
        int k,i=0;
        char *tok;

        k=atoi(line+1)+1;
        if (k<0 || k>NPRED+1) ERROR(("%s: @%d out of range",lastFn,k-1))
        if (k<0 || k>npred-1) ERROR(("%s: @%d exceeds @k=%d",lastFn,k-1,npred-2))
        tok=strtok(line+1," \t");
        tok=strtok(NULL," \t");
        if (!tok || tok[0]<' ') ERROR(("%s: @%d missing omega",lastFn,k-1))
        else pred[k].omega=atof(tok);
        while ( (tok=strtok(NULL," \t")) ) {
          if (i>k) ERROR(("%s: @%d too many coefficients",lastFn,k-1))
          pred[k].c[i++]=atof(tok); } }

      if (!fgets(line,128,f)) break; }

    if (npred) {
      underline("SCF predictor");

      prt("\
Gaussian orbitals predicted via chk-files.\n\
One step of the algorithm (except start):\n\
   yP(t) := c[0]*y(t-h) + c[1]*y(t-2h) + ... \n\
   y(t)  := omega*SCF(yP(t)) + (1-omega)*yP(t)\n\
       t := t+h\n\
y stands for orbital expansion coefficients\n\
SCF() = Gaussian iteration (one or more)");

      /* NPRED=4 assumed in this code
               0   1.00000   1.00000  0.00000  0.00000  0.00000  */
      header(" k    omega      c[1]     c[2]     c[3]     c[4] ");
      loop (k,0,npred)
        prt("%2d  %8.5f  %8.5f %8.5f %8.5f %8.5f",
            k-1,pred[k].omega,
            pred[k].c[0],pred[k].c[1],pred[k].c[2],pred[k].c[3]);
      header(""); }
    else
      prt("no SCF predictor");
    
    fclose(f); }


  /*** making com-file (Gaussian input) */
  f=fopen(Fn("com"),"wt");
  if (!f) ERROR(("cannot create %s",lastFn))

  loop (n,3,ncom0) fputs(com0[n],f);
  found=getenv("NSLOTS");
  if (found) n=atoi(found); else n=1;
  if (n<1 || n>48) ERROR(("NSLOTS=%d is out of range",n))
  fprintf(f,"%%Nproc=%d\n",n);
  fprintf(f,"%s\n\
MACSIMUS-generated configuration t=%.9g\n\
\n\
%s",com0[!!measure],t,com0[2]);

  loop (n,0,No.N) {
    mn=molec+n;
    ns=mn->ns;
    sp=mn->sp;
    si=spec[sp]->si;
    r=rof(mn,A->rp);

    loop (i,0,ns)
      fprintf(f,"%-2s %9.6f %9.6f %9.6f\n",
              sitedef[si[i].st].name,VARG(r[i])); }

  fprintf(f,"\n");
  fclose(f);

  if (option('v')&8) { prt_("CHK PREDbefore t=%g",t); loop (j,0,ny) prt_(" %g",pred[0].y[j]); _n }

  /*** predict SCF, shift the history ***/
  if (npred && known) {
    char *chkfn=Fn("chk");
    char *fchkfn=Fn("fchk");
    /*
      known=0: nothing known (start, 1st pass)
      known=1: pred[0].y=y[t-h] known
      known=2: pred[0].y and pred[1].y known (corresponding to k=0)
      ...
      NB: yP is actually stored in pred[0].y */

    /* predictor order k=known-2 and shift of the history */
    loop (j,0,ny) {
      double yy=0;

      if (j==0)
        loop (k,0,known)
          prt("t=%g known=%d k=%d pred[known-1].c[k]=%g pred[k].y[0]=%g",
              t,known,k,pred[known-1].c[k],pred[k].y[0]);
      
      loop (k,0,known) yy+=pred[known-1].c[k]*pred[k].y[j];
      for (k=min(npred-1,known); k>0; k--) pred[k].y[j]=pred[k-1].y[j];
      pred[0].y[j]=yy; }

    /* re-create the chk-file */
    f=fopen(fchkfn,"wt");

    /* put data of yP=pred[0].y to fchkfile */
    j=0;
    loop (k,0,NHEADERS) {
      for (i=0; i<headers[k].N; i+=5) {
        ptr=fchkfile+LINE*(headers[k].pos+i/5);
        switch (headers[k].N-i) {
          case 1:  sprintf(ptr," %15.8E\n",pred[0].y[j]); j++; break;
          case 2:  sprintf(ptr," %15.8E %15.8E\n",pred[0].y[j],pred[0].y[j+1]); j+=2; break;
          case 3:  sprintf(ptr," %15.8E %15.8E %15.8E\n",pred[0].y[j],pred[0].y[j+1],pred[0].y[j+2]); j+=3; break;
          case 4:  sprintf(ptr," %15.8E %15.8E %15.8E %15.8E\n",pred[0].y[j],pred[0].y[j+1],pred[0].y[j+2],pred[0].y[j+3]); j+=4; break;
          default: sprintf(ptr," %15.8E %15.8E %15.8E %15.8E %15.8E\n",pred[0].y[j],pred[0].y[j+1],pred[0].y[j+2],pred[0].y[j+3],pred[0].y[j+4]); j+=5; break; } } }
    if (j!=ny) ERROR(("internal %d %d",ny,j)) 

    loop (i,0,nchklines) fputs(fchkfile+LINE*i,f);

    fclose(f);

    if (system(string("./unfchk %s %s",fchkfn,chkfn))) ERROR(("./unfchk %s %s FAILED",fchkfn,chkfn))

  } /* predict SCF */

  if (option('v')&8) { prt_("CHK PREDknown=%d t=%g",known,t); loop (j,0,ny) prt_(" %g",pred[0].y[j]); _n }

  /*** Gaussian called ***/
  if (system(string("./gaussian %s",Fn("com"))))
    ERROR(("system call ./gaussian %s failed",Fn("com")))

  if (npred) {
    /*** read fchk file, headers:
Coordinates of each shell                  R   N=          72
Constraint Structure                       R   N=          27
Alpha MO coefficients                      R   N=        3249
Total SCF Density                          R   N=        1653
    */
    char *chkfn=Fn("chk");
    char *fchkfn=Fn("fchk");

    if (system(string("./formchk %s %s",chkfn,fchkfn)))
      ERROR(("./formchk %s %s FAILED",chkfn,fchkfn))

    f=fopen(fchkfn,"rt");
    if (!f) ERROR(("%s missing",fchkfn))

    if (!fchkfile) {
      /* allocate the local copy of the fchk file */
      while (fgets(line,128,f)) nchklines++;
      rewind(f);
      allocarrayzero(fchkfile,LINE*nchklines); }

    /* read in the fchk file */
    i=0; ny=0;
    while (fgets(line,128,f)) {
      if (i>=nchklines)
        ERROR(("%s: number of lines changed while reading",fchkfn))
        /* patch because line 2 contains padded spaces to 90 chars */
        if (strlen(line)>80) {
          if (line[80]=='\n') line[81]=0;
          else strcpy(line+81,"\n"); }
      memcpy(fchkfile+i*LINE,line,LINE);
      loop (k,0,NHEADERS)
        if (!memcmp(line,headers[k].hdr,strlen(headers[k].hdr))) {
          headers[k].pos=i+1; /* where the data start */
          ny+=headers[k].N=atoi(line+50); }
      i++; }
    if (i!=nchklines) ERROR(("%s: number of lines changed while reading",fchkfn))

    /* 1st pass: allocate all arrays with fchk data */
    if (!yC) {
      allocarrayzero(yC,ny);
      loop (k,0,npred) allocarrayzero(pred[k].y,ny); }

    /* get data from fchkfile, yC=SCF(yP); yP=pred[0].y */
    j=0;
    loop (k,0,NHEADERS) {
      for (i=0; i<headers[k].N; i+=5) {
        ptr=fchkfile+LINE*(headers[k].pos+i/5);
        switch (headers[k].N-i) {
          case 0:  err=1; break;
          case 1:  err=sscanf(ptr,"%lf",            yC+j)-1; j++; break;
          case 2:  err=sscanf(ptr,"%lf%lf",         yC+j,yC+j+1)-2; j+=2; break;
          case 3:  err=sscanf(ptr,"%lf%lf%lf",      yC+j,yC+j+1,yC+j+2)-3; j+=3; break;
          case 4:  err=sscanf(ptr,"%lf%lf%lf%lf",   yC+j,yC+j+1,yC+j+2,yC+j+3)-4; j+=4; break;
          default: err=sscanf(ptr,"%lf%lf%lf%lf%lf",yC+j,yC+j+1,yC+j+2,yC+j+3,yC+j+4)-5; j+=5; break; }
        if (err) ERROR(("table %s format or number of data",headers[k].pos-1)) } }
    if (option('v')&8) prt("DEBUG: t=%g chk read, %d=%d+%d+%d+%d",t,ny,headers[0].N,headers[1].N,headers[2].N,headers[3].N);
    if (j!=ny) ERROR(("internal"))

    if (option('v')&8) { prt_("CHK CORRgaussian t=%g",t); loop (j,0,ny) prt_(" %g",yC[j]); _n }
        
    /* corrector = mixing SCF(yP) with yP=pred[0].y, returned as pred[0].y
       (could be merged with the previous step to spare array yC) */
    if (known>0) omega=pred[known-1].omega; else omega=1;
    loop (j,0,ny) pred[0].y[j]=omega*yC[j]+(1-omega)*pred[0].y[j];

    if (option('v')&8) { prt_("CHKomega=%g CORRknown=%d t=%g",omega,known,t); loop (j,0,ny) prt_(" %g",pred[0].y[j]); _n }

    if (known<npred) known++; }


  /*** get forces and energy from the log-file */
  f=fopen(Fn("log"),"rt");
  if (!f) ERROR(("cannot read %s",lastFn))

/*
-------------------------------------------------------------------
Center     Atomic                   Forces (Hartrees/Bohr)
Number     Number              X              Y              Z
-------------------------------------------------------------------
1          8           -.049849321     .000000000    -.028780519
2          1            .046711997     .000000000    -.023346514
3          1            .003137324     .000000000     .052127033
-------------------------------------------------------------------
*/
  while (fgets(line,128,f)) {
    if (strstr(line,"SCF Done:")) {
      char *c=strchr(line,'=');

      if (c) {
        static double Epot0=0;
        double Epot;

        while (*c++==' ');
        Epot=atof(c)*315775.0248020554;
        if (Epot0==0) {
          Epot0=Epot;
          prt("Epot = %.17g K (from \'SCF Done:\') will be used as zero ref",Epot); }
        En.pot=Epot-Epot0; } }

    if ( (found=strstr(line,"Forces (Hartrees/Bohr)")) ) break; }

  if (!found) ERROR(("Forces not found in %s",lastFn))
  loop (n,0,2) if (!fgets(line,128,f)) ERROR(("Not enough lines in %s",lastFn))

  tabn=0;
  loop (n,0,No.N) {
    mn=molec+n;
    ns=mn->ns;
    sp=mn->sp;
    si=spec[sp]->si;
    r=rof(mn,B->rp);

    loop (i,0,ns) {
      tabn++;
      if (!fgets(line,128,f)) ERROR(("forces for mol=%d site=%d missing in %s",n,i,lastFn))
      sscanf(line,"%d%d%lf%lf%lf",&nn,&aux,r[i],r[i]+1,r[i]+2);
      if (nn!=tabn) ERROR(("forces numbering out of sync: read=%g, expected=%d\n\
*** in line %s",nn,tabn,line))
      /* evu: Eh/a0/1[-N] = 596728.313869787 */
      VO(r[i],*=596728.313869787) } }

  fclose(f);
  if (measure) rename(Fn("log"),Fn(string("%d.log",(int)(t/No.dt+0.5))));
  userforces(B,A);

  if (option('v')&128) dumpall(B,A);
}
