static int indx;

static int strcmpx(const char *a,const char *x) /******************* strcmpx */
/* cheap regexp:
     * in x replaces any chars (incl. none)
     ? replaces one char
   return values:
     0 = match
     1 = no match
   SEE /home/jiri/macsimus/c/alg/strcmpx.c
*/
{
  for (;*a;a++,x++) {
    if (*a==*x || *x=='?') continue;
    if (*x=='*') {
      x++;
      if (!*x) return 0;
      while (*x!=*a) {
        a++;
        if (!*a) return 1; } }
    else return 1; }

  while (*x=='*') x++;

  return *x!=0;
}

static char *next(void) /********************************************** next */
{
  char *c=strtok(NULL," \t\n");

  if (!c) ERROR(("missing data in line with number %d in table atoms",indx))

  return c;
}

void readMOL(char *fn) /******************************************** readMOL */
/***
  reads *.mol file in the MEDIT format.
  fn must have extension .mol
***/
{
  char *tok,*c,*slash;
  static char resnm[4];
  int i,j,nnbr,inbr,nbr,ic,testnc=0,ndonor=0,nacceptor=0;
  struct collist_s *cl;

  file=fopen(fn,"rt");
  if (!file) {
    prt("WARNING: %s: no mol-file\n\
*** the configuration will be shown with dots (x) or white spheres\n\
*** remedy hints: molcfg, cook -y# (will call molcfg), bonds",fn);
    ns=1;
    allocarrayzero(site,ns);
    site[0].id="sphere";
    site[0].type="X";
    site[0].molcol='H';
    site[0].r=50; /* diameter = 1 AA */
    nbond=1;
    if (option('i')<=0) nbond=1;
    else nbond=option('i');
    allocarray(bond,nbond);
    bond[0][0]=bond[0][1]=0;

    nc=testnc=1; }

  else while (mygetline()) if ( (tok=strtok(li," =\t\n")) ) {

    if (!strcmp(tok,"parameter_set")) {
      tok=strtok(NULL," =\t\n");
      parameter_set=strdup(tok); }
    else if (!strcmp(tok,"number_of_atoms")) {
      tok=strtok(NULL," =\t\n");
      ns=atoi(tok);
      if (option('v')) prt("! number_of_atoms = %i",ns);
      allocarrayzero(site,ns); /* must :=0 - used later */ }

    else if (!strcmp(tok,"atoms")) {
      if (ns<=0) ERROR(("%s: undefined or bad number_of_atoms",fn))

      if (option('i')<=0) nbond=ns/8+8; /* first guess: will grow on demand */
      else nbond=option('i');
      allocarray(bond,nbond);

      loop (i,0,ns) {
        if (!mygetline()) ERROR(("%s: missing lines",fn))
        tok=li;
        site[i].mark=0;
        if (*tok=='*') { tok++; site[i].mark=4; }
        indx=atoi(strtok(tok," \t"));
        if (i!=indx) ERROR(("%s: atoms are not contiguously numbered, consider command:\n\
***   mol2mol NAME +a NEWNAME",fn))

        tok=next();
        slash=strchr(tok,'/');
        if (slash) {
          memcpy(site[i].Crad,tok,4);
          tok=slash+1;
          site[i].id=strdup(tok); }
        else
          site[i].id=strdup(tok);
        if (strlen(tok)>3 && memcmp(resnm,tok,3)) {
          memcpy(resnm,tok,3); resnm[3]=0;
          if (option('v')) printf("\r%d:%s",i,resnm); }

        /* id: mark if connect to backbone */
        if (abs(option('b'))==1) {
          if (strstr(aminoacids,resnm)) {
            /* is aminoacid */
            for (c=tok+3; (*c>='0' && *c<='9'); c++); /* skip res. number */
            if (*c>='a' && *c<='z') c++;              /* skip chain letter */
            /* mark backbone C or CA or N (1=not draw) */
            site[i].c=strcmp(c,"C") && strcmp(c,"CA") && strcmp(c,"N"); }
          else {
            /* was not aminoacid - mark if -b, not mark if -b-1 */
            site[i].c=option('b')<0; } }

        site[i].r=16; /* means no color mark found */
        for (cl=collisthead; cl; cl=cl->next)
          if (cl->pattern && strstr(tok,cl->pattern))
            site[i].r=cl->colno;

#ifdef SHELL
        /* protein & draw marked */
//        site[i].shell=!memcmp(tok,"HOH",3) || strchr(tok,'w') ? 1 : 4;
        site[i].shell= strstr(tok,shellkey)?1:4;
#endif
        tok=next(); /* type */
        loop (j,0,i) if (!strcmp(tok,site[j].type)) {
          site[i].type=site[j].type;
          goto tdone; }
        site[i].type=strdup(tok);
       tdone:

#ifdef SHELL
        if (strchr(tok,'O')) site[i].shell |= 2; /* water O marked */
#endif
        site[i].molcol=tok[0];
        /* H-bonds */
        if (!strcmpx(tok,donor)) {
          site[i].mark|=1;
          ndonor++; }
        if (strcmp(donor,acceptor))
          if (!strcmpx(tok,acceptor)) {
            site[i].mark|=2;
            nacceptor++; }

        tok=next(); /* charge */
#ifdef CHARGES
        site[i].charge=atof(tok);
#endif
        next();
        nnbr=atoi(next());
        loop (inbr,0,nnbr) {
          nbr=atoi(next());
#ifdef SHELL
          if (site[i].shell==3) wmol[inbr]=nbr-i;
#endif
          loop (ic,0,nc)
            if (bond[ic][0]==nbr && bond[ic][1]==i) {
              testnc++; goto cont; }
          if (nc>=nbond) nbond=enlarge(&bond,nbond,nbond*3/2+4);
          bond[nc][0]=i; bond[nc][1]=nbr;
          nc++;
        cont:; }
        if (nnbr==0 && option('f')) {
          if (nc>=nbond) nbond=enlarge(&bond,nbond,nbond*3/2+4);
          bond[nc][0]=i; bond[nc][1]=i;
          nc++; testnc++; } } } }

  if (file) fclose(file);

#ifdef SHELL
  if (abs(wmol[0]+wmol[1])>3) {
    if (option('l') && option('j')>0) ERROR(("unsupported water"))
    wmol[0]=wmol[1]=0; }
  if (option('v')) prt("\nwater mol: H1-O=%d H2-O=%d",wmol[0],wmol[1]);
#endif

  _n

#ifdef CHARGES
  loop (i,0,ns) {
    if (negA!=NULL || posA!=NULL) {
      int gc=groupcharge(i);

      if (negA) if (gc<=option('n')) site[i].r=negC;
      if (posA) if (gc>=option('p')) site[i].r=posC;

      if (option('e')) { /* the same for all neighbors */
        int j;

        loop (j,0,nc) {
          int b0=bond[j][0],b1=bond[j][1];

          if (i==b0 || i==b1) {
            gc=groupcharge(i==b0?b1:b0);
            if (negA) if (gc<=option('n')) site[i].r=negC;
            if (posA) if (gc>=option('p')) site[i].r=posC; } } } } }
#endif

  if (nc!=testnc) ERROR(("%s: [%d:%d] bad bonds",fn,nc,testnc))

  loop (ic,0,nc) {
    i=bond[ic][0]; nbr=bond[ic][1];
    bond[ic][2] = colortab[
                           site[i].r & site[nbr].r & 16 ? 0 :
                           site[i].r&16 ? site[nbr].r : site[i].r].turbocolor; }

  if (option('v')) {
    prt("%s: %d sites, %d bonds",fn,ns,nc);
    if (ndonor+nacceptor)
      prt("%s: %d donors (%s) and %d acceptors (%s)",
          fn, ndonor,    donor,  nacceptor,    acceptor); }

  /* bond==bond0 if not h-bonds */
  bond0=bond; nc0=nc; nbond0=nbond;
  bondeqbond0=1;
  ncfix=nc;
}

void mol2mol(int newns) /********************************************* newns */
{
  int i,j,n,newnc;
  site_t *oldsite=site;
  bond_t *oldbond=bond;

  allocarrayzero(site,newns*sizeof(site_t));
  loop (i,0,newns) site[i]=oldsite[i%ns];
  free(oldsite);

  allocarray(bond,nbond);

  newnc=0;
  n=(newns+ns-1)/ns;
  loop (i,0,n)
    loop (j,0,nc)
      if (oldbond[j][0]+i*ns<newns && oldbond[j][1]+i*ns<newns) {
        if (newnc>=nbond) nbond=enlarge(&bond,nbond,nbond*3/2+4);
        bond[newnc][0]=oldbond[j][0]+i*ns;
        bond[newnc][1]=oldbond[j][1]+i*ns;
        bond[newnc][2]=oldbond[j][2];
        newnc++; }
  free(oldbond);

  ns=newns;
  nc=newnc;

  /* cf. readMOL above */
  bond0=bond; nc0=nc; nbond0=nbond;
  bondeqbond0=1;
  ncfix=nc;
}

void PDBbackbone(void) /**************************************** PDBbackbone */
{
  int i0=0,ic,j;

  loop (ic,0,nc)
    if (!site[bond[ic][0]].c && !site[bond[ic][1]].c) {
      bond[ic][2]&=255;
      loop (j,0,3) bond[i0][j]=bond[ic][j];
      i0++; }
  nc=i0;
}

void backbone(void) /****************************************** backbone */
/* warning: readGOL must be called after backbone */
{
  int *tree,i=0,i0,j=0,ic,t=0,ii;

#ifdef DEBUG
  out=fopen("show.dbg","wt");
#endif

  if (option('v')) prts("looking for the backbone...");
  alloc(tree,ns*sizeof(int));
  memset(tree,0,ns*sizeof(int));
  loop (i,0,ns) site[i].c=1;

 again:
    i=j;
    tree[i]=++t;
 again2:
#ifdef DEBUG
    put2(i,t)
#endif
    i0=-1;
    loop (ic,0,ns) {
      if (bond[ic][0]==i) j=bond[ic][1];
      else if (bond[ic][1]==i) j=bond[ic][0];
      else continue;
      /* (j) is neighbour of (i) */
      if (!tree[j]) goto again;
#ifdef DEBUG
      else if (tree[j]<t-1) prt("%i %i: cycle",i,j);
#endif
      else if (tree[j]==t-1) i0=j; }

    i=i0;
    if (--t && i0>=0) goto again2;

#ifdef DEBUG
  loop (i,0,ns) prt("%4i : %i",i,tree[i]);
#endif

  i0=ns-1;
  t=tree[i0];
  if (!t) {
    /* no begin-end backbone: try find max from beginning */
    t=0; i0=0;
    loop (i,0,ns) if (tree[i]>t) { t=tree[i]; i0=i; } }

/* mark bonds+sites of main chain */
  while (--t) {
    loop (ic,0,nc) {
      i=bond[ic][0]; j=bond[ic][1];
      if (i>=0 && j>=0)
        if ( (i==i0 && tree[ii=j]==t) || (tree[ii=i]==t && j==i0) ) {
          bond[ic][2]|=256;
          site[i].c=0; site[j].c=0;
          i0=ii;
          break; } } }

#ifdef DEBUG
  loop (ic,0,nc) prt("%4i - %-4i  %i",bond[ic][0],bond[ic][1],bond[ic][2]);
#endif

  /* remove other bonds */

  i0=0;
  loop (ic,0,nc)
    if (bond[ic][2]&256) {
      bond[ic][2]&=255;
      loop (j,0,3) bond[i0][j]=bond[ic][j];
      i0++; }
  nc=i0;

#ifdef DEBUG
  loop (ic,0,nc) prt("%4i - %i",bond[ic][0],bond[ic][1]);
  fclose(out);
#endif

  free(tree);
}

void hbonds(void) /************************************************** hbonds */
/* calculate hydrogen bonds (from bond[ncfix] to bond[nc]) */
{
  int i,j,k,AD,jpair;
  float x,rr,hbdistq;
  static int lastb;

  /* BUG: still problems with nc0 after removing H-bonds */
  nc0=nc=ncfix; /* remove H-bonds */

  if (hbdist==0) return;

  hbdistq=Sqr(hbdist);

  if (option('v')>2)
    for (i=0; i<ns; i+=4) fprintf(stderr,"%d: %d/%d %d/%d %d/%d %d/%d\n",
                                  i/4,
                                  site[i].c,site[i].r,
                                  site[i+1].c,site[i+1].r,
                                  site[i+2].c,site[i+2].r,
                                  site[i+3].c,site[i+3].r);

  if (strcmp(donor,acceptor)) jpair=0,AD=2;
  else jpair=1,AD=1;

  loop (i,0,ns) if (site[i].mark==1 && site[i].c<NOCOL*COLLEN)
    loop (j,(i+1)*jpair,ns) if (site[j].mark==AD && site[j].c<NOCOL*COLLEN) {
      rr=0;
      loop (k,0,3) {
        x=cfg[i][k]-cfg[j][k];
        if (option('l')) if (L[k]) {
          while (x>Lh[k]) x-=L[k];
          while (x<-Lh[k]) x+=L[k]; }
        rr+=Sqr(x);
        if (rr>hbdistq) goto cont; }
      /* already a bond? */
      loop (k,0,nc)
        if (bond[k][0]+bond[k][1]==i+j && abs(bond[k][0]-bond[k][1])==abs(i-j)) goto cont;
      if (nc>=nbond) nbond=enlarge(&bond,nbond,nbond*5/4+4);
      bond[nc][0]=i;
      bond[nc][1]=j;
      bond[nc][2]=hbcolor;
      nc++;
     cont:; }

  if (option('v')>1 && lastb!=nc-ncfix)
    prt("%d hydrogen bond(s) found (|%s-%s|<%g)",nc-ncfix,acceptor,donor,hbdist);

  lastb=nc-ncfix;
}

void readGOL(void) /*********************************************** readGOL */
/***
  reads *.gol file
  fn must have extension .gol
  warning: readGOL must be called after backbone
***/
{
  char *tok,*fn=Fn("gol");
  int i,nsgol,c;
  static int GOLread=0;

  if (GOLread || option(']')==0) return;

  file=fopen(fn,"rt");
  if (!file) {
    fprintf(stderr,"%s: no goal file, using atom info\n",fn);
    loop (i,0,ns) if (site[i].Crad[0]) {
      char col=site[i].Crad[0];
      loop (c,0,NOCOL-1)
        if (colortab[c].col[0]==col) goto colfound1;
      c=NOCOL-1;
    colfound1:
      if (!(site[i].r&16)) c=site[i].r;
      site[i].r=atoi(site[i].Crad+1)*!site[i].c;
      /* note that sites to remove have been marked by backbone() in site[i].c */
      site[i].c=1+c*COLLEN; }
    else {
      /* table in blendmed.c differs in LIGHTGRAY x LIGHTCYAN for Carbon */
                     /*          111111*/
      static         /*0123456789012345*/
      char colors[17]="xxxxxxx  NXCOMSH";
                     /*BBGCRMBGGBGCRMYW*/
                     /*klryearrrlryeaeh*/

      static unsigned char rads[16]={0,0,0,0,0,0,0,0,145,128,145,140,110,103,132,70};

      /* see also blendmed.c, function atomcolor() */

      char ch=site[i].molcol;
      int col;

      loop (col,7,16) if (ch==colors[col]) break;
      if (col==16) col=LIGHTMAGENTA;
      if (ch=='P') col=BROWN; /* = ORANGE, phosphorus default */
      /* TO BE OPTIMIZED: REPEATED BELOW */
      loop (c,0,NOCOL-1)
        if (col==colortab[c].turbocolor) break;
      if (!(site[i].r&16)) c=site[i].r;
      site[i].r=rads[col]*!site[i].c;
      /* note that sites to remove have been marked by backbone() in site[i].c */
      site[i].c=1+c*COLLEN;
    }
  }
  else {
    /* reading the GOL-file */
    mygetline();
    nsgol=atoi(li);
    if (nsgol!=ns) prt("%d site(s) found in %s (expected %d): will be %s",
      nsgol,fn,ns,nsgol<ns?"replicated":"truncated");

    loop (i,0,min(ns,nsgol)) {
      if (!mygetline()) Error("gol file too short");
      tok=strtok(li," \t\n");
      /* GRAY==CYAN for Carbons */
      if (!strcmp(tok,"GRAY")) tok="CYAN";
      if (!strcmp(tok,"ORANGE")) tok="N";
      loop (c,0,NOCOL-1)
        if (colortab[c].col[0]==tok[0]) goto colfound2;
      c=NOCOL-1;
    colfound2:
      if (!(site[i].r&16)) c=site[i].r;
      site[i].r=(100*atof(strtok(NULL," \t\n"))+0.5)*!site[i].c;
      /* note that sites to remove have been marked by backbone() in site[i].c */
      site[i].c=1+c*COLLEN; }
    if (nsgol<ns) loop (i,nsgol,ns) {
      site[i].r=site[i-nsgol].r;
      site[i].c=site[i-nsgol].c; } }

  GOLread=1;
}


/*** marking ***/

static void floodmark(int center) /******************************* floodmark */
{
  int n,i;

  site[center].mark|=4;

  if (curmode>=BOND) circle(xy0[center][0],xy0[center][1],2);

  loop (n,0,nc) {
    if (bond[n][0]==center) if (!(site[i=bond[n][1]].mark&4)) floodmark(i);
    if (bond[n][1]==center) if (!(site[i=bond[n][0]].mark&4)) floodmark(i); }
}

static void markatom(int key,int i) /****************************** markatom */
{
  if (curmode>=BOND) {                    
    setwritemode(!trace);
    setcolor(MODE-2); }
  
  if (key==MIDCLICK)
    floodmark(i);
  else if (key==LEFTCLICK) {
    if (!(site[i].mark&4)) {
      site[i].mark|=4;
      if (curmode>=BOND) circle(xy0[i][0],xy0[i][1],2); } }
  else {
    if ((site[i].mark&4)) {
      site[i].mark-=4;
      if (curmode>=BOND) circle(xy0[i][0],xy0[i][1],2); } }
  if (curmode>=BOND) setwritemode(0);
}

void printmarkedinfo(int i) {
  char line[128];
  FILE *edt=NULL;
  int maxfw;
  static int llasti=0,lllasti=0,newedt=1;
  float lrr,rr,cf;
  int j;

  /* this part probably out of order... */
  if (option('~')) edt=fopen(Fn("edt"),option('~')==1 || !newedt?"at":"wt");
  if (edt && newedt) fprintf(edt,"! %s\n",Fn("edt"));
  newedt=0;

  prt("--- %d: %s %s %f",i,site[i].id,site[i].type,site[i].charge);
  prt("    [%8.5f %8.5f %8.5f ]",cfg[i][0],cfg[i][1],cfg[i][2]);

  maxfw=maxx/xfont.width;
  setfillstyle(0,0); /* black in both modes */
  setwritemode(0);

  if (curmode>=BOND && clickmode==0) {
    bar(0,0,maxx,2*xfont.height);
    //                    setcolor(palcolor(WHITE));
    setcolor(MODE-2);
    sprintf(line,"%d: %s %s",i,site[i].id,site[i].type);
    maxfw=(maxx-outtextxymax(1,1,line,maxfw))/xfont.width; /* space left in char */
    for (j=4; j>=0; j-- ) {
      sprintf(line,"%.*f %.*f %.*f",j,cfg[i][0],j,cfg[i][1],j,cfg[i][2]);
      if (strlen(line)<=maxfw) break; }
    //                    setcolor(palcolor(LIGHTCYAN));
    outtextxy(maxx-strlen(line)*xfont.width,1,line); }
  
  /* bond length */
  lrr=rr=cf=0;
  loop (j,0,3) {
    float x=cfg[i][j]-cfg[lasti][j];
    float y=cfg[llasti][j]-cfg[lasti][j];
    cf+=y*x;
    rr+=x*x;
    lrr+=y*y; }
  prt_("    %d-%d=%.5f",i,lasti,sqrt(rr));
  sprintf(line,"%d-%d=%.2f",i,lasti,sqrt(rr));
  //                  setcolor(palcolor(YELLOW));
  outtextxy(1,1+xfont.height,line);
  
  /* bond angle */
  lrr*=rr;
  if (lrr) {
    prt_("  %d-%d-%d=%.4f", i,lasti,llasti, (180/PI)*acos(cf/sqrt(lrr)));
    sprintf(line,"%d-%d-%d=%.2f", i,lasti,llasti,(180/PI)*acos(cf/sqrt(lrr)));
    outtextxy(maxx-strlen(line)*xfont.width,1+xfont.height,line); }
  
  /* dihedral angle */
  if (lllasti!=llasti)
    prt_("  %d-%d-%d-%d=%.4f", i,lasti,llasti,lllasti,
         dihedral(cfg[i],cfg[lasti],cfg[llasti],cfg[lllasti]));
  _n
    
  if (edt) {
    char *tok=strtok(line," \t");
    tok=strtok(NULL," \t");
    fprintf(edt,option('~')>0?"ra %s\n":"rm %s\n",tok);
    fclose(edt); }
  lllasti=llasti; llasti=lasti; lasti=i;
}
  
static int markcolor;

static int colormark(int i) /************************************* colormark */
/* returns the color of site i, or the mark color if marked */
{
  if (site[i].mark&4) return markcolor;
  else return site[i].c;
}

/*** clusters ***/

static int cluster,clsize;

static void recursecluster(int i) /************************** recursecluster */
{
  int j;

  if (site[i].cluster) return;

  site[i].cluster=cluster;
  clsize++;

  /* extremely inefficient because of storing bonds separately */
  loop (j,0,nc) {
    if (bond[j][0]==i) recursecluster(bond[j][1]);
    if (bond[j][1]==i) recursecluster(bond[j][0]); }
}

static void findclusters(int nth) /**************************** findclusters */
{
  int i,j;
  struct cl_s { 
    int cluster;
    int size; }
  *cl,xch;

  if (nth<0) nth=1;

  fprintf(stderr,"searching for the %d-th largest cluster\n",nth);

  allocarrayzero(cl,ns); /* subscript shifted by -1 */

  loop (i,0,ns) site[i].cluster=0;

  cluster=0;
  
  loop (i,0,ns) {
    cluster++;
    clsize=0;
    recursecluster(i);
    cl[cluster-1].cluster=cluster;
    cl[cluster-1].size=clsize; }

  /* bubblesort enough */
  j=1;
  while (j) {
    j=0;
    loop (i,1,cluster) {
      if (cl[i-1].size<cl[i].size) {
        j=1;
        xch=cl[i-1],cl[i-1]=cl[i],cl[i]=xch; } } }

  loop (i,0,ns) if (site[i].cluster==cl[nth-1].cluster) site[i].mark|=4;
  //  loop (i,0,ns) fprintf(stderr,"%d %d\n",i,site[i].cluster);
  free(cl);
}
