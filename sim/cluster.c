/*
  cluster analysis
  #included from simmeasx.c #ifdef CLUSTERS
*/
#define PERMUTETYPE int
#include "permute.c"

#ifndef BONDHIST
/* # of previously broken or created bonds (from a molecule to a molecule).
   At least 2 needed for reliable definition of `really new' bonds. */
#  define BONDHIST 32
#endif /*# BONDHIST */

             /* mode format maxn maxcluster */
struct cl_s cl={0,   1,     20,  12};

static int cmpint(const void *a,const void *b) /********************* cmpint */
{
  if (*(int*)a<*(int*)b) return -1;
  else if (*(int*)a>*(int*)b) return 1;
  else return 0;
}

static double cluster_cq;
static int isincq(real *r1,real *r2) /******************************* isincq */
{
  double rr;
  vector dr;

  beginCellImage(0,cluster_cq) 
    return 1;
  endCellImage

  return 0;
}

static struct nbrs_s {
  struct nbrs_s *next;
  int n;
#ifdef XSECTION
  int ithis; /* site of molecule nn engaged in the bond (from nbrs[nn]->..) */
  int inbr; /* sites creating bond of molecule this->n */
#endif /*# XSECTION */
}
  **nbrs, /* nbrs[nn] is the head of the list of neighbors of molecule nn */
  **last; /* last timestep */

static void destroylast(void) /********************************* destroylast */
{
  int n;

  if (last) {
    struct nbrs_s *lst,*llst;

    loop (n,0,No.N)
      for (lst=last[n]; lst; lst=llst) {
        llst=lst->next; free(lst); }

    free(last); }
 else
   ERROR(("internal"))
}

static int cluster_k; /* frame (timestep) for bonds, integer 0,1, */

static struct bondhist_s {
  struct bondhist_s *next;
  float t; /* time */
  int k; /* frame index from 0 */
  int x; /* 1: creation, 0: breakage; 2:flag of one create-break-... series*/
  int n; /* neighbor */
} **bondhist; /* bondhist[i] = head of the list of created or broken bonds
                 for molecule i */

static float tstart,tfinish; /* patchy way how to get time */

static int *color; /* [No.N] 
                      all molecules i in a cluster have the same color[i]=1,2,..
                      color[i]=0 = not assigned yet 
                      max color = mark (incl.) */
static FILE *clustcp;
static char (*golcolors)[64];

typedef struct cluster_s {
  struct cluster_s *next;
  char *name;  /* only if it is a named cluster */
  int col;     /* only if it is a named cluster, column # in SIMNAME.cl.cpa */
  int count;   /* how many times found in one config (active only if name) */
  int total;   /* how many times found total */
  int nm;      /* # of molecules */
  int version; /* 1,2,3... for clusters with the same stoichometry */
  int *stoich; /* stoichometry, [nspec] */
  int **nbrs;  /* array[nm] of neighbor lists;
                  each neighbor list starts with its size */
  float (*r)[3]; /* first cluster positions (of site[0]) stored */
} cluster_t;

static cluster_t *clhead,*scl;
static int maxnbrs=8,*rennbrs,nnbr;
static int *ren;

int cmpcluster(const cluster_t *a,const cluster_t *b) /********** cmpcluster */
/* Returns 0 if the clusters are identical (stoichometry+topology), 1 otherwise
   Considers all clusters larger than cl.maxcluster molecules with
   the same stoichometry as identical (returns 0)
   WARNING: does not define ordering - not to be used for sorting */
{
  int n,i;
  int *clnbrs,*cnbrs;

  /* compare cluster size */
  if (a->nm!=b->nm) return 1;

  /* compare stoichiometry */
  if (memcmp(a->stoich,b->stoich,nspec*sizeof(b->stoich[0]))) return 1;

  /* too large clusters considered identical */
  if (a->nm>cl.maxcluster) return 0;

  /* compare topology */

  initpermutes0(b->nm,ren,nspec,b->stoich);

  do { /* all permutations of molecules of the same type in cluster cl */

    /* we try to match permuted cl : c */
    loop (n,0,b->nm) {

      /* mol.# n in cl is matched to mol.# ren[n] in c */
      clnbrs=b->nbrs[n];
      nnbr=clnbrs[0];
      cnbrs=a->nbrs[ren[n]];

      /* compare # of neighbors */
      if (nnbr!=cnbrs[0]) goto nextperm;

      if (maxnbrs<nnbr) {
        /* update auxiliary array of renumbered neighbors if has grown */
        free(rennbrs);
        allocarray(rennbrs,nnbr+1);
        maxnbrs=nnbr; }

      /* compare neighbors, in standard order */
      loopto (i,1,nnbr) rennbrs[i]=ren[clnbrs[i]];
      if (nnbr>1) qsort(rennbrs+1,nnbr,sizeof(rennbrs[0]),cmpint);
      if (memcmp(rennbrs+1,cnbrs+1,nnbr*sizeof(rennbrs[0]))) goto nextperm;

    } /* n */

    /* now c == renumbered cl */
    return 0;

  nextperm:;
  } while (permutes(ren));

  return 1;
}

static
int cmpstoich(const cluster_t *a,const cluster_t *b) /*********** cmpstoich */
/* Cluster ordering according to stoichiometry,
   to be used for sorting for `nice' print only
   WARNING: topologically different clusters are not distinguished */
{
 /* compare cluster size */
 if (a->nm<b->nm) return -1;
 if (a->nm>b->nm) return 1;

 /* compare stoichiometry */
 return memcmp(a->stoich,b->stoich,nspec*sizeof(b->stoich[0]));
}

static
void prtcluster(FILE *removescript,cluster_t *c) /*************** prtcluster */
{
  int i,j,n,k,sp;
  char pref=' ';
  char formula[64],*f;
  double charge=0;

  formula[0]=0;

  /* PATCH, see also PATCH below */
  if (!removescript && cl.mode>1) prt_("%5d ", c->count);
  prt_("%7d",c->total);

  loop (i,0,nspec) if (c->stoich[i]) {
    sprintf(strend(formula),"%c%s",pref,spec[i]->name);
    pref='.';
    // omit stoichiometric coefficient 1 //    if (c->stoich[i]!=1)
    sprintf(strend(formula),"%d",c->stoich[i]);
    charge+=spec[i]->charge*c->stoich[i]; }

  if (c->version) {
    if (c->version<27) sprintf(strend(formula),"%c",'a'-1+c->version);
    else sprintf(strend(formula),"-%d",c->version); }

  prt_("%-12s ",formula);
  charge/=electron;
  if (fabs(charge)<1e-6) prt_(" 0");
  else prt_("%+2.0f",charge);

  if (c->name) prt_(" %d:%s ",c->col,c->name);
  else prt_("    ");
  if (c->nm<=cl.maxcluster) loop (i,0,c->nm) loopto (j,1,c->nbrs[i][0])
    if ( (n=c->nbrs[i][j]) >=i)
      prt_(" %d-%d", i,n);
  _n

  if (removescript) if (c->r) {
    FILE *mol,*plb,*gol=NULL;
    static float hdr[2];

    f=strend(formula);
    strcpy(f,".mol");
    mol=fopen(formula+1,"wt");
    if (removescript) fprintf(removescript,"rm%s",formula);
    if (golcolors) {
      strcpy(f,".gol");
      gol=fopen(formula+1,"wt");
      fprintf(gol,"!\n!\n%d\n",c->nm);
      if (removescript) fprintf(removescript,"%s",formula); }
    strcpy(f,".plb");
    plb=fopen(formula+1,"wt");
    if (removescript) fprintf(removescript,"%s\n",formula);

    hdr[0]=c->nm;
    fwrite(hdr,2,4,plb);
    fwrite(c->r,c->nm,12,plb);
    fclose(plb);

    fprintf(mol,"cluster");
    if (c->name) fprintf(mol," %s",c->name);
    *f=0;
    fprintf(mol,"%s\n\n\
number_of_atoms=%d\n\n\
atoms\n",formula,c->nm);

    sp=k=0;
    loop (i,0,c->nm) {
      char *nm;

      while (k>=c->stoich[sp]) sp++,k=0;
      k++;

      if (gol) fprintf(gol,"%s\n",golcolors[sp]);
      nm=spec[sp]->name;
      fprintf(mol,"%2d %d-%s %s 0 0",i,i,nm,nm);

      loopto (j,0,c->nbrs[i][0])
        fprintf(mol," %d",c->nbrs[i][j]);
      fprintf(mol,"\n"); }

    fclose(mol);
    if (gol) fclose(gol); }
}

#ifdef XSECTION
static void addnbr(int n1,int i1,int n2,int i2) /******************** addnbr */
/*
   neighbor n2, site s2 is added to the list of neighbors of n1
   multiple entries possible
*/
{
  struct nbrs_s *nbr;

  allocone(nbr);
  nbr->n=n2;
  nbr->ithis=i1;
  nbr->inbr=i2;
  nbr->next=nbrs[n1];
  nbrs[n1]=nbr;
}
#else /*# XSECTION */
static void addnbr(int n1,int n2) /********************************** addnbr */
/*
   neighbor n2 is added to the list of neighbors of n1
   multiple entries possible
*/
{
  struct nbrs_s *nbr;

  allocone(nbr);
  nbr->n=n2;
  nbr->next=nbrs[n1];
  nbrs[n1]=nbr;
}
#endif /*#!XSECTION */

static int mark; /* number of (so far) found clusters */

/* recursive */
static void findcluster(int n) /******************************** findcluster */
{
  struct nbrs_s *nbr;

  color[n]=mark;
  looplist (nbr,nbrs[n])
    if (!color[nbr->n]) {
#ifdef XSECTION
      /* always? re-use for mode&2 ? */
      vector *r1=rof(molec+n,xs.A->rp);
      int ns=molec[n].ns,i;
      vector *r2=rof(molec+nbr->n,xs.A->rp);
      int k;

      loop (k,0,3) {
        while (r1[nbr->ithis][k]<r2[nbr->inbr][k]-box.Lh[k]) {
          loop (i,0,ns) r2[i][k]-=box.L[k]; }
        while (r1[nbr->ithis][k]>r2[nbr->inbr][k]+box.Lh[k]) {
          loop (i,0,ns) r2[i][k]+=box.L[k]; } }
#endif /*# XSECTION */
      findcluster(nbr->n); }
}

static int getsitetype(char *nm) /****************************** getsitetype */
{
  int i;

  loop (i,0,nsites) if (!strcmp(nm,sitedef[i].name)) return i;
  ERROR(("%s is unknown site type",nm))

  return -1;
}

static int spectype(char *c) /************************************* spectype */
{
  int i;

  loop (i,0,nspec) if (!strcmp(c,spec[i]->name)) return i;
  ERROR(("no species (molecule) %s",c))

  return 0;
}

/*****************************************************************************/
/*               chemical format - modified from blendmed.c                  */
/*****************************************************************************/

/* see also code in blend/blendmed.c */

/* 'x' no longer means crossing bonds: use '+' */
const char CHESEP[]=" -/\\=|+\t\n";
static struct coor_s { int st,n,x0,x1,y; } *coor;
static int ns;
static char *clustername;
static int siten(int x,int y) /*************************************** siten */
{
  int i;

  loop (i,0,ns)
    if (y==coor[i].y)
      if (coor[i].x0<=x && x<=coor[i].x1) return coor[i].n;

  ERROR(("%s: bad bond in rel.line=%d, col.=%d",clustername,y,x))

  return 0;
}

static char **screen;
static int si;

static void findnbr(int x,int y,int c,int dx,int dy) /************** findnbr */
{
  int *nbr,from,to;

  if (screen[y][x]!=c) return;

  do {
    x+=dx; y+=dy;
  } while (screen[y][x]==c || screen[y][x]=='+'); /* 'x' removed */
  to=siten(x,y);
  from=coor[si].n;
  if (from>=ns || to>=ns) ERROR(("internal"))
  allocarray(nbr,scl->nbrs[from][0]+2);
  copy(nbr,scl->nbrs[from],(scl->nbrs[from][0]+1)*sizeof(nbr[0]));
  nbr[++nbr[0]]=to;
  free(scl->nbrs[from]);
  scl->nbrs[from]=nbr;
}

cluster_t *readCHE(FILE *f) /*************************************** readCHE */
{
  int nx,ny,i,is,x,st;
  char *c;
  int *stsum;
  cluster_t *cl;
  char line[256];

  struct line_s {
    struct line_s *next;
    char *l; } *head=NULL,*l;

  alloconezero(cl);
  scl=cl;
  allocarrayzero(cl->stoich,nspec);

  rallocarrayzero(stsum,nspec);

  /***
    read and store to a list
    # of sites calculated
    # of lines & columns calculated
  ***/
  nx=ny=ns=0;
  while (fgets(line,256,f) && line[0]!='\n')
    if (!strchr("!#",line[0])) {
      rallocone(l);
      l->l=strdup(line);
      if (!l->l) ERROR(("strdup"))
      l->next=head;
      head=l;
      Max(nx,strlen(line));

      c=strtok(line,CHESEP);
      while (c) {
        if (isalpha(c[0])) ns++;
        c=strtok(NULL,CHESEP); }

      ny++; }

  if (option('v')&4) prt("cluster %d sites  %d lines  %d columns",ns,ny,nx);
  cl->nm=ns;

  nx+=2;

  ralloc(screen,(ny+2)*sizeof(*screen));

  ralloc(screen[0],nx); memset(screen[0],' ',nx);
  ralloc(coor,sizeof(struct coor_s)*ns);

  allocarrayzero(cl->nbrs,ns);

  is=0;
  for (i=ny,l=head; i>=1; i--,l=l->next) {
    if (!l) ERROR(("internal"))
    ralloc(screen[i],nx);
    memset(screen[i],' ',nx);
    if (strlen(l->l)>=nx) ERROR(("internal"))
    strcpy(screen[i]+1,l->l);

    c=strtok(l->l,CHESEP);
    while (c) {
      if (isalpha(c[0])) {
        if (is>=ns) ERROR(("internal"))
        coor[is].x0=(int)(c-l->l)+1;
        coor[is].x1=(int)(c-l->l)+strlen(c);
        coor[is].y=i;
        st=spectype(c);
        coor[is].st=st;
        coor[is].n=cl->stoich[st]++;
        is++; }
      c=strtok(NULL,CHESEP); } }

  if (l || is!=ns) ERROR(("internal"))

  loop (is,1,nspec) stsum[is]=stsum[is-1]+cl->stoich[is-1];
  loop (is,0,ns) coor[is].n+=stsum[coor[is].st];

  ralloc(screen[ny+1],nx); memset(screen[ny+1],' ',nx);

  /*** finding bonds ***/
  loop (si,0,ns) alloconezero(cl->nbrs[si]);

  loop (si,0,ns) {

    /* shorter sides */
    findnbr(coor[si].x0-1,coor[si].y,'-',-1,0);
    findnbr(coor[si].x1+1,coor[si].y,'-',1,0);

    /* longer sides + corners */
    loopto (x,coor[si].x0,coor[si].x1) {
      findnbr(x-1,coor[si].y-1,'\\',-1,-1);
      findnbr(x,coor[si].y-1,'|',0,-1);
      findnbr(x+1,coor[si].y-1,'/',1,-1);

      findnbr(x+1,coor[si].y+1,'\\',1,1);
      findnbr(x,coor[si].y+1,'|',0,1);
      findnbr(x-1,coor[si].y+1,'/',-1,1); } }

  loop (si,0,ns)
    if (cl->nbrs[si][0]>1)
      qsort(cl->nbrs[si]+1,cl->nbrs[si][0],sizeof(cl->nbrs[0][0]),cmpint);

  release(stsum);
  return cl;
} /* readCHE */

static double **dist; /* [nsites][nsites] table of distances^2 defining bonds */
static int namedclusters;

void readclusterdef(char *fn) /****************************** readclusterdef */
{
  char line[256];
  char nm1[16],nm2[16];
  double d;
  int nd=0,i;
  FILE *f,*oldf=NULL;

  if ((cl.mode&1)==0) {
    prt("NOTE: CLUSTER calculation is off because (cl.mode&1)==0");
    return; }

  f=fopen(fn,"rt");
  if (!f) ERROR(("readclusterdef: %s not found",fn))

  if (!rennbrs)
    allocarray(rennbrs,maxnbrs+1); /* permanent, growing on demand */

  for (;;) {
    char *tok;

    while (!fgets(line,256,f)) {
      if (!oldf) goto prtinfo;
      f=oldf; oldf=NULL; }

    tok=strtok(line," \t\n=");

    if (!line[0] || strchr("!#\n",line[0]))
      continue;

    else if (line[0]=='$' && line[1]=='i') {
      char *fn=line+2;

      if (oldf)
        WARNING(("nested $i in reading cluster definitions: returns only once"))
      if (!*fn) fn=strtok(NULL," \t\n=");
      if (!fn) ERROR(("$i missing file name"))
      prt("opening include file \"%s\"",fn);
      oldf=f;
      f=fopen(fn,"rt");
      if (!f) {
        ERROR(("no include file \"%s\"",fn))
        f=oldf; } }

    else if (!strcmp("clusters",tok)) cl.mode|=2;

    else if (!strcmp("configurations",tok)) cl.mode|=4;

    else if (!strcmp("bonddynamics",tok)) cl.mode|=8;

    else if (!strcmp("colors",tok)) {
      int i;

      allocarray(golcolors,nspec);
      loop (i,0,nspec) strcpy(golcolors[i],"WHITE 1");
      while (fgets(line,256,f)) {
        char *tok=strtok(line," \n\t");
        int sp;

        if (!tok) break;
        sp=spectype(tok);
        tok=strtok(NULL,"\n");
        if (tok) strcpy(golcolors[sp],tok); } }

    else if (!strcmp("bonds",tok)) {
      alloc2Darrayzero(dist,nsites,nsites);

      while (fgets(line,256,f)) {
        int sp1,sp2;

        i=sscanf(line,"%s%s%lf",nm1,nm2,&d);
        if (i<=0) break;
        if (i!=3) ERROR(("%sATOM ATOM DIST expected",line))
        sp1=getsitetype(nm1);
        sp2=getsitetype(nm2);

        dist[sp1][sp2]=dist[sp2][sp1]=d*d;
        nd++; } }

    else if (!strcmp("maxcluster",tok)) {
      cl.maxcluster=atoi(strtok(NULL," \t\n="));
      if (cl.maxcluster<1 || cl.maxcluster>2048)
        ERROR(("cl.maxcluster out of bounds 1..2048")) }

    else if (!strcmp("mincluster",tok)) {
      cl.mincluster=atoi(strtok(NULL," \t\n="));
      if (cl.mincluster<1 || cl.mincluster>1024)
        ERROR(("cl.mincluster out of bounds 1..1024")) }

    else if (!strcmp("maxn",tok)) {
      cl.maxn=atoi(strtok(NULL," \t\n="));
      if (cl.maxn<1 || cl.maxn>65536)
        ERROR(("cl.maxn out of bounds 1..65536")) }

    else if (!strcmp("format",tok)) {
      cl.format=atoi(strtok(NULL," \t\n="));
      if (cl.format==0)
        WARNING(("cl.format=0: no output"))
      if (cl.format<0 || cl.format>63)
        ERROR(("cl.format=%d out of bounds 0..63",cl.format)) }

    else if (!strcmp("cluster",tok)) {
      cluster_t *cl,*c;
      char *name=strdup(strtok(NULL," \t\n="));
      int inm;

      clustername=name; /* static, for err msg */
      cl=readCHE(f);
      namedclusters++;

      allocarray(ren,cl->nm);
      cl->name=name;
      inm=namedclusters;
      looplist (c,clhead) {
        inm--;
        if (!cmpcluster(c,cl))
          ERROR(("named clusters #%d=%s and #%d=%s are equivalent\n\
*** (check whether cl.maxcluster is large enough)",
                 inm,c->name,namedclusters,cl->name)) }
      free(ren);

      cl->next=clhead;
      clhead=cl; } }

 prtinfo:
  prt("cluster analysis data: %d distances, cl.maxcluster=%d, %d named clusters",
      nd,cl.maxcluster,namedclusters);
  /* warnings and opening clustcp postponed to analyzeclusters() */

  prt("Lines with cluster data contain keyword CLUSTERCOUNT:\n\
SEE THE MANUAL FOR DETAILS");

  fclose(f);
}

void analyzeclusters(int frame) /*************************** analyzeclusters */
{
  int i,i1,i2,n,n1,n2,sp1,sp2,ns1,ns2,imark,nbonds=0;
  molecule_t *mn1,*mn2,*mn;
  siteinfo_t *si1,*si2;
  vector *r1,*r2;
  vector lastL;
  int maxinframe=0;
  cluster_t *c;
  static int pass;
  int *clustercount=NULL; /* initialized to suppress compile warning */

  if ((cl.mode&1)==0) return;

  if (cl.maxn) allocarrayzero(clustercount,cl.maxn);

  // probably forgoten..removed
  //  if (!iscube()) ERROR(("analyzeclusters: cube expected"))

  VV(lastL,=box.L)

  if (!dist) {
    if (!pass++) prt("no cluster analysis");
    namedclusters=0;
    return; }
  if (cl.maxcluster<0 || cl.maxcluster>32) WARNING(("suspicious cl.maxcluster"))

  rallocarrayzero(color,No.N);
  allocarrayzero(nbrs,No.N); /* to be renamed to last and freed later */

  loop (n1,0,No.N) {
    mn1=molec+n1;
    sp1=mn1->sp;
    r1=rof(mn1,cfg[0]->rp);
    ns1=mn1->ns;
    si1=spec[sp1]->si;

    loop (n2,0,n1) {
      mn2=molec+n2;
      sp2=mn2->sp;
      ns2=mn2->ns;
      r2=rof(mn2,cfg[0]->rp);
      si2=spec[sp2]->si;

      loop (i1,0,ns1) loop (i2,0,ns2)
        if ( (cluster_cq=dist[si1[i1].st][si2[i2].st]) )
          if (isincq(r1[i1],r2[i2])) {
#ifdef XSECTION
            addnbr(n1,i1,n2,i2); addnbr(n2,i2,n1,i1);
#else /*# XSECTION */
            addnbr(n1,n2); addnbr(n2,n1);
#endif /*#!XSECTION */
            nbonds++;
            goto nextmols; }

      nextmols:; } }

  /* all neighbors nbrs[]-> known now */

#ifdef XSECTION
  /* cross-section interface (bad programming style) */
  if (!xs.A) sdsalloc(xs.A,cfg[0]->size);
  sdscopy(xs.A,cfg[0]);
  /* .. findcluster will periodically adjust molecules */
#endif /*# XSECTION */

  mark=0;
  loop (n,0,No.N)
    if (!color[n]) {
      mark++;
      findcluster(n); }

#ifdef XSECTION
  if (!xs.color) allocarray(xs.color,No.N);
  copyarray(xs.color,color,No.N);

#  ifndef FREEBC
  /* find periodic clusters */
  loopto (imark,1,mark) { // BUG found in V 3.5b: was loop (imark,0,mark)
    int ncl=0,periodic=0;

    loop (n,0,No.N) if (imark==color[n]) {
      struct nbrs_s *nbr;

      ncl++;
      looplist (nbr,nbrs[n]) {
        vector *r1=rof(molec+n,xs.A->rp);
        vector *r2=rof(molec+nbr->n,xs.A->rp);
        int k;

        loop (k,0,3) {
          if (r1[nbr->ithis][k]<r2[nbr->inbr][k]-box.Lh[k]) periodic++;
          if (r1[nbr->ithis][k]>r2[nbr->inbr][k]+box.Lh[k]) periodic++; }
      } /* neighbors of n */
    } /* n of color imark */ 
    if (periodic)
      ERROR(("%d-cluster %d is periodic (%d)",ncl,imark,periodic))
  } /* imark */
#  endif /*# FREEBC */
#endif /*# XSECTION */

  if (cl.mode&4) { /* keyword configurations */
    /* Good for monoatomic molecules only!
       Export the whole configuration (plb, mol; if colors are given in the
       cli-file, then also gol) with connection information.
       Sites #0 in molecules are connected instead of the sites really
       connected, original bonds ignored.
       Subclusters are periodically shifted - may fail in case of a
       periodic cluster.  */
    FILE *mol=fopen(Fn(string("%d.mol",frame)),"wt");
    FILE *plb=fopen(Fn(string("%d.plb",frame)),"wb");
    FILE *gol=NULL;
    float r[DIM],(*r0)[DIM];
    int is,*renum; /* [No.N] : molecule -> site #0 renumbering */

    if (golcolors) gol=fopen(Fn(string("%d.gol",frame)),"wb");
    alloc(r0,mark*4*DIM);
    loop (i,0,mark) r0[i][0]=-9e9;
    allocarray(renum,No.N);

    is=0;
    loop (n,0,No.N) {
      renum[n]=is;
      is+=molec[n].ns; }

    /* new format, keeping periodic b.c - why not to use free b.c.? */
    r[0]=No.s;
    r[1]=-3;
    fwrite(r,4,2,plb);
    VV(r,=box.L)
    fwrite(r,4,DIM,plb); /* DIM=3 only ... */

    if (gol) fprintf(gol,"!sphere#1\n\
! %s\n\
%d\n",lastFn,No.s);

    fprintf(mol,"frame %d\n\n\
number_of_atoms = %d\n\n\
atoms\n",frame,No.s);

    is=0;
    loop (n,0,No.N) {
      int ns=molec[n].ns,i;
      struct nbrs_s *nbr;
      int col=color[n];

      if (col<0 || col>mark) ERROR(("mark=%d color[%d]=%d",mark,n,col))
      col--;

      mn=molec+n;

      if (gol) fprintf(gol,"%s\n",golcolors[mn->sp]);

      nnbr=0;
      looplist (nbr,nbrs[n]) nnbr++;

      loop (i,0,ns) {
        fprintf(mol,"%4d %2s-%dc%*s %-2s 0 0 %d",
                is,
                sitedef[spec[mn->sp]->si[i].st].name,
                col,4-(int)strlen(string("%d",col)),"",
                sitedef[spec[mn->sp]->si[i].st].name,
                nnbr);

        VV(r,=rof(mn,cfg[0]->rp)[i])

        if (!i)
          looplist (nbr,nbrs[n]) fprintf(mol," %d",renum[nbr->n]);

        if (r0[col][0]<-8e8) {
          copy(r0[col],r,4*DIM); }
#ifndef FREEBC
        else {
          int k;

          loop (k,0,DIM) {
            while (r[k]-r0[col][k]> box.Lh[k]) r[k]-=box.L[k];
            while (r[k]-r0[col][k]<-box.Lh[k]) r[k]+=box.L[k]; }
        }
#endif /*# FREEBC */
        fwrite(r,4,DIM,plb);
        is++;
        nnbr=0; /* sites >0 not connected */
        fprintf(mol,"\n"); } }

    fclose(mol);
    fclose(plb);
    if (gol) fclose(gol);
    free(r0);
    free(renum); }

  loopto (imark,1,mark) {
    cluster_t *b;
    int inm=0;

    alloconezero(b);
    loop (n,0,No.N) b->nm+=color[n]==imark;

    allocarrayzero(b->stoich,nspec);
    allocarrayzero(b->nbrs,b->nm);
    allocarray(ren,No.N);

    loop (n,0,No.N) if (color[n]==imark) {
      ren[n]=inm++;
      b->stoich[molec[n].sp]++; }

    if (cl.mode&2) allocarray(b->r,inm);

    inm=0;
    loop (n,0,No.N) if (color[n]==imark) {
      struct nbrs_s *nbr;
      vector *r=rof(molec+n,cfg[0]->rp);

      nnbr=0;
      looplist (nbr,nbrs[n]) nnbr++;

      allocarray(b->nbrs[inm],nnbr+1);
      b->nbrs[inm][0]=nnbr;
      nnbr=0;
      looplist (nbr,nbrs[n]) b->nbrs[inm][++nnbr]=ren[nbr->n];
      if (nnbr>1) qsort(b->nbrs[inm]+1,nnbr,sizeof(b->nbrs[0][0]),cmpint);

      if (b->r) {
        VV(b->r[inm],=r[0])
#ifndef FREEBC
        if (inm) {
          int k;

          loop (k,0,DIM) {
            while (b->r[inm][k]-b->r[0][k]<-box.Lh[k]) b->r[inm][k]+=box.L[k];
            while (b->r[inm][k]-b->r[0][k]> box.Lh[k]) b->r[inm][k]-=box.L[k]; } }
#endif /*# FREEBC */
      }

      inm++; }

    if (cl.maxn) {
      Max(maxinframe,inm)
      if (inm<cl.maxn) clustercount[inm]++;
      else clustercount[0]++; }

    /* .. now b is a `standard' image of the cluster, unique but permutation */
    /* array ren[] will re-used */

    /* finding whether b matches any of the clusters already detected */
    looplist (c,clhead)
      if (!cmpcluster(c,b)) {
        c->count++;
        c->total++;
        loop (i,0,b->nm) free(b->nbrs[i]);
        free(b->nbrs);
        free(b->stoich);
        if (b->r) {
          if (!c->r) {
            allocarray(c->r,c->nm);
            loop (i,0,c->nm) VV(c->r[ren[i]],=b->r[i]) }
          free(b->r); }
        free(b);
        goto match; }

    b->count=b->total=1;
    b->next=clhead;
    clhead=b;
  match:;
    free(ren); }

  if (cl.maxn) {
    int i;

    /* NB: clustercount[0] = all maxn and bigger */
    if (cl.format&1) {
      if (cl.format&16) prt_("%.3f ",t);
      loop (i,1,cl.maxn) prt_("%d ",clustercount[i]);
      prt(" %d %d CLUSTERCOUNT1 %d",clustercount[0],maxinframe,mark); }
    if (cl.format&2) {
      if (cl.format&16) prt_("%.3f ",t);
      prt_("%d  ",maxinframe);
      loop (i,1,cl.maxn) prt_("%d ",clustercount[i]);
      prt(" %d CLUSTERCOUNT2 %d",clustercount[0],mark); }
    if (cl.format&4) {
      int k;
      if (cl.format&16) prt_("%.3f ",t);
      loop (k,0,clustercount[0]) prt_("%d+ ",cl.maxn);
      for (i=cl.maxn-1; i>0; i--) loop (k,0,clustercount[i]) prt_("%d ",i);
      prt(" CLUSTERCOUNT4 %d",mark); }
    if (cl.format&8) {
      if (cl.format&16) prt_("%.3f ",t);
      if (clustercount[0]) prt_("%d+*%d ",cl.maxn,clustercount[0]);
      for (i=cl.maxn-1; i>0; i--) if (clustercount[i]) prt_("%d*%d ",i,clustercount[i]);
      prt(" CLUSTERCOUNT8 %d",mark); }
    free(clustercount); }

  /* list all clusters (print molecule numbers) */
  if (cl.format&32) {
    loopto (imark,1,mark) {
      int clsize=0;
      
      loop (n,0,No.N) if (imark==color[n]) clsize++;
      if (clsize>=cl.mincluster) {
        if (cl.format&16) prt_("%.3f ",t);
        prt_("%d %d ",imark,clsize);
        loop (n,0,No.N) if (imark==color[n]) prt_(" %d",n);
        prt(" CLUSTERCOUNT32 %d",mark); } } }
  
  {
    /* sort clusters with respect to stoichiometry */
    cluster_t cl00,*cl0=&cl00,*aux;
    int changed;

    cl0->next=clhead;

    if (cl0->next) do {
      changed=0;
      for (c=cl0; c->next->next; c=c->next)
        if (cmpstoich(c->next,c->next->next)>0) {
          aux=c->next->next;
          c->next->next=aux->next;
          aux->next=c->next;
          c->next=aux;
          changed++; }
      } while (changed);
    clhead=cl0->next;
  }

  /* mark clusters with the same stoichiometry */
  if (clhead) for (c=clhead; c->next; c=c->next)
    if (cmpstoich(c,c->next)==0) {
      cluster_t *cc,*zero=NULL;
      int maxversion=0;

      for (cc=c; cc && cmpstoich(cc,c)==0; cc=cc->next) {
        Max(maxversion,cc->version);
        if (!cc->version) zero=cc; }
      if (zero) zero->version=maxversion+1; }

  /* number named clusters
     (=columns in SIMNAME.cl.cpa, used as info in .sfd only) */
  if (clhead) {
    int n=0;

    for (c=clhead; c; c=c->next) if (c->name) c->col=++n; }

  if (pass==0) {
    int l,i;

    pass++;
    /* open clustcp if necessary */
    if (namedclusters) clustcp=fopen(Fn("cl.cpa"),"wt");

    if (clustcp) loop (l,0,2) {
      fprintf(clustcp,"#");
      i=1;
      for (c=clhead; c; c=c->next)
        if (c->name) {
          if (l) fprintf(clustcp," %s",c->name);
          else fprintf(clustcp," %*d",(int)strlen(c->name),i++); }
      if (l) fprintf(clustcp," REST  bonds\n");
      else fprintf(clustcp," %4d  bonds\n",i); } }

  {
    int rest=0;

    for (c=clhead; c; c=c->next) {
      if (option('v')&4) prtcluster(NULL,c);
      if (clustcp) {
        if (c->name) fprintf(clustcp," %*d",(int)strlen(c->name),c->count);
        else rest+=c->count; }
      c->count=0; }
    if (clustcp) {
      fprintf(clustcp," %4d",rest);
      fprintf(clustcp," %6d\n",nbonds); }
  }

  release(color);

  /*** break bond analysis ***/

  if (cl.mode&8) {
    if (last) {
      int n;
      struct bondhist_s *bh;
      
      if (!bondhist) {
        allocarrayzero(bondhist,No.N);
        tstart=t; }
      tfinish=t;
      
      loop (n,0,No.N) {
        struct nbrs_s *nbr,*lst;
        
        for (nbr=nbrs[n]; nbr; nbr=nbr->next) {
          for (lst=last[n]; lst; lst=lst->next)
            if (nbr->n==lst->n) goto contnbr;
          
          /* created bond to a neighbor found */
          //          if (option('v')&4) prt("created bond: k=%d %d-%d",k,n,nbr->n);
          // ??? - to be checked here
          if (option('v')&4) prt("created bond: k=%d %d-%d",cluster_k,n,nbr->n);
          
          allocone(bh);
          bh->t=t;
          bh->k=cluster_k;
          bh->n=nbr->n;
          bh->x=1;
          bh->next=bondhist[n];
          bondhist[n]=bh;
          
        contnbr:; }
        
        for (lst=last[n]; lst; lst=lst->next) {
          for (nbr=nbrs[n]; nbr; nbr=nbr->next)
            if (lst->n==nbr->n) goto contlst;
          
        /* broken bond to a neighbor found */
          //          if (option('v')&4) prt("broken bond: k=%d %d-%d",k,n,lst->n);
          // ??? - to be checked here
          if (option('v')&4) prt("broken bond: k=%d %d-%d",cluster_k,n,lst->n);
          
          allocone(bh);
          bh->t=t;
          bh->k=cluster_k;
          bh->n=lst->n;
          bh->x=0;
          bh->next=bondhist[n];
          bondhist[n]=bh;

        contlst:; } }
      
      destroylast();
      cluster_k++; } 
    last=nbrs; }
  else {
    last=nbrs;
    destroylast(); }

  VV(box.L,=lastL)
  VV(box.Lh,=0.5*box.L)
}

static void bonddynamics(void)
{
  int n,i,j;

  typedef int (*cp_t)[BONDHIST+1];
  /* created[k][0]: created
     created[k][1]: created with create-break-... removed
     created[k][2]: created with create-(other)-break .. removed ... */
  cp_t created,broken;

  if (!dist) return;
  if (!bondhist) return;

  destroylast();

  prt("create/break bond convergence profiles from t=%.3f to %.3f (k=%d)",
      tstart,tfinish,cluster_k);

  allocarrayzero(created,cluster_k);
  allocarrayzero(broken,cluster_k);

  loopto (i,0,BONDHIST)
    loop (n,0,No.N) {
      struct bondhist_s *bh,*c=NULL;

      for (bh=bondhist[n]; bh; bh=bh->next) bh->x &= 1;

      for (bh=bondhist[n]; bh; bh=bh->next) if (!(bh->x&2)) {
        int x=bh->x,nflip=1;

        bh->x |= 2;
        /* this loop empty (no check) for i=0
 	  for i>0, finding create-break... or break-create... sequences
 	  separated by at most i-1 other bond changes */
        for (c=bh->next,j=i; c && j; c=c->next,j--) if (c->n == bh->n) {
          x ^= 1;
          if (c->x != x)
            ERROR(("INTERNAL (n=%d i=%d c->t=%f c->n=%d)",n,i,c->t,c->n))
          c->x |= 2;
 	 j=i+1; /* .. looks strange, but not thet there's j-- immediately */
          nflip++; }

        if (nflip&1) {
          /* odd number of flips, e.g. create-break-create */
          j=0;
          for (c=bh; c; c=c->next)
            if (c->n == bh->n)
              if (j++==nflip/2) break;

          if (!c) ERROR(("internal"))

          if ( (option('v')&4) && nflip>1)
            prt("i=%d %d-%d %s %d flips center=%d %f",
                i,n,c->n,bh->x&1 ? "created" : "broken",nflip,c->k,c->t);

          if (c->k<0 || c->k>=cluster_k)
            ERROR(("internal: c->k=%d not in range [0,%d]",c->k,cluster_k))
          (bh->x&1 ? created : broken)[c->k][i]++; } } }

  loop (n,0,2) {
    char *info=n?"created":"broken";
    FILE *f=fopen(Fn(n?"created.cpa":"broken.cpa"),"wt");
    cp_t cp=n?created:broken;

    if (!f) ERROR(("open %s",lastFn))
    fprintf(f,"# t 2*%s[0] ... 2*%s[%d]\n",info,info,BONDHIST);
    loop (i,0,cluster_k) {
      fprintf(f,"%7.3f",tstart+i*(tfinish-tstart)/(cluster_k-1));
      loopto (j,0,BONDHIST) {
        fprintf(f," %d",cp[i][j]);
        if (j%BONDHIST) StaSet(0,2,2,0);
        else StaSet(0,lag.err,2,lag.n);
        StaAdd(string("%s[%d]",info,j),cp[i][j]*0.5); }
      fprintf(f,"\n"); }
  fclose(f);
  free(cp); }
}

void printclusters(void) /************************************ printclusters */
{
  if (cl.mode&8) bonddynamics();

  if (clhead) {
    cluster_t *c;
    FILE *removescript=NULL;

    /* PATCH, to be checked !!! */
    if (cl.mode>1) removescript=fopen(Fn("rm.sh"),"wt");

    header(" total  formula   charge [name] bonds (0 = 1st molecule from chem.formula..)");
    looplist (c,clhead) prtcluster(removescript,c);
    if (clustcp) fclose(clustcp);
    if (removescript) fclose(removescript);
    header(""); }
}
