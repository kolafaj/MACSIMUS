#include "ground.h"
#include "options.h"
#include "vector.h"
#include "intrapot.h"
#include "sitesite.h"
#include "blendpar.h"
#include "blendmed.h"
#include "rndgen.h"
#include "blendedt.h"

int ble_file;

static int match(char *s1, char *s2) /******************************** match */
/*
  returns 1 if s1 matches s2, 0 otherwise
  s2 may contain a wildcard `?' in arbitrary position (exactly 1 char)
  or `*' at the end (any # of characters)
  as soon as I find regexp somewhere, I'll use it...
*/
{
  char *c1=s1,*c2=s2;

  for (;;c1++,c2++) {
    if (*c1==0 && *c2==0) return 1;
    if (*c2=='*') return 1;
    if (*c1==0 || *c2==0) return 0;
    if (*c2!='?' && *c1!=*c2) return 0; }
}

static char Line[128];
static site_t *site,*fsite;

int findsite(species_t *spec,char *id,int chk) /****************** findsite */
{
  int i;
  int n=-1;
  loop (i,0,spec->ns) {

    if (match(spec->site[i].id,id)) {
      if (chk) if (n>=0) WARNING(("%s: %s not unique",Line,id))
      n=i; } }

  /* if no match found, analyzed for site # */
  if (n<0) {
    char *c;

    for (c=id; *c; c++) if (!isdigit(*c)) goto skip;
    n=atoi(id);
    if (n>=spec->ns) ERROR(("site %s >= ns=%d",id,spec->ns));
   skip:; }

  if (chk)
    if (n<0) ERROR(("%s does not match any site",id))
      else fsite=&spec->site[n]; /* ignored if called globally */

  return n;
}

static void prtsite(char *s1,site_t *s,int n,char *s2)
{
  prt_("%s %s [%d] %s %s", s1,s->id,n,atom[s->type].name,s2);
  if (ble_file) fprintf(stderr,"%s %s [%d] %s %s", s1,s->id,n,atom[s->type].name,s2);
}

static int maxns;

static char *wateratom(int i)
{
  static char waterinfo[32];

  if (water0->ns!=3)
    strcpy(waterinfo,"ERROR only ns=3 water supported");
  else
    sprintf(waterinfo,"%s %g",atom[water0->type[i]].name,water0->charge[i]);

  return waterinfo;
}

void addH(species_t *spec,int i,char *itype,int nH,char *Htype,double Hcharge)
/* adds nH hydrogens (type=Htype) to atom i (new type=itype) */
{
  site_t *sh;
  site_t *si=spec->site+i;
  char *c;
  int iH;

  loopto (iH,1,nH)
    if (spec->ns>=maxns) ERROR(("addH"))
    else {
      sh=spec->site+spec->ns;
      alloczero(sh->id,strlen(si->id)+3);
      strcpy(sh->id,si->id);
      if (strlen(sh->id)>4) {
        c=strend(sh->id)-2;
        if (*c!='C') c--;
        if (*c!='C') c[3]='H',c[4]=iH+'0';
        else c[5]=c[3],c[4]=c[2],c[3]=c[1],c[2]='H',c[1]=iH+'0',c[0]='~'; }

      sh->type=atomn(Htype);
      sh->charge=Hcharge;
      si->charge-=Hcharge;
      si->type=atomn("CT");
      sh->nnbr=1;
      sh->nbr[0]=i;
      if (si->nnbr>=MAXVAL) {
        ERROR(("extended C too many bonds"))
        return; }
      si->nbr[si->nnbr++]=spec->ns;
      sh->keep=WANTED;
      sh->r[0]=UNDEFATOM;
      spec->ns++; }

  si->type=atomn(itype);
}

void movesite(species_t *spec,int from,int to) /****************** movesite */
{
  int ns=spec->ns,i,n;
  site_t *s,aux;

  if (from<0 || from>=ns || to<0 || to>=ns) ERROR((""))
  else if (from<to) {
    aux=site[from];
    loop (i,from,to) site[i]=site[i+1];
    site[to]=aux;
    loop (i,0,ns) {
      s=site+i;
      loop (n,0,s->nnbr)
        if (s->nbr[n]==from) s->nbr[n]=to;
        else if (s->nbr[n]>from && s->nbr[n]<=to) s->nbr[n]--; } }
  else if (from>to) {
    aux=site[from];
    for (i=from; i>to; i--) site[i]=site[i-1];
    site[to]=aux;
    loop (i,0,ns) {
      s=site+i;
      loop (n,0,s->nnbr)
        if (s->nbr[n]==from) s->nbr[n]=to;
        else if (s->nbr[n]<from && s->nbr[n]>=to) s->nbr[n]++; } }

  /*.....checkneighbors(spec,spec->ns);*/
}

void moledit(species_t *spec) /************************************* moledit */
{
  char *line,*l,*oldline;
#if MAXVAL<5
#define MAXNBRID 5
#else
#define MAXNBRID MAXVAL
#endif
  char *id,*wid,*wdist,*wabsy=NULL,*wabsz, *nbrid[MAXNBRID];
  char *type,*charge,*key;
  int i,n,nnbr,addwater=4,water=0;
  site_t *s;
  const char *sep=" \t";
  FILE *edit;
  int probe=0;
  int ff=1; /* force field match */
  int parsed;

  maxns=spec->ns;

  site=spec->site;
  alloc(line,option('l'));
  strcpy(spec->ext,".edt");
  edit=fopen(spec->fn,"rt");
  if (!edit) {
    edit=stdin;
    fprintf(stderr,"\n! no %s : enter editing commands now (`?' for help)\n",
	    spec->fn); }
  else
    prt("! editing molecule using %s",spec->fn);

  *(spec->ext)=0;

  for (;;) {
    switch (++addwater) {
      case 1:
	if (probe) {
	  sprintf(line,"pf W%d-O %s %s %s %s",
		  water,wateratom(2),wid,wdist,wabsy);
	  spec->probe.ns=2; /* DIRTY: 1 will be added when "pf" parsed */
	  probe=0; }
	else
	  if (!strcmp(wid,"="))
	    sprintf(line,"af W%d-O %s %s %s %s %s",
		    water,wateratom(2),wid,wdist,wabsy,wabsz);
	  else {
	    if (wabsy && wabsy[0]) water=atoi(wabsy);
	    sprintf(line,"af W%d-O %s %s %s",water,wateratom(2),wid,wdist); }
        break;
      case 2: sprintf(line,"aa W%d-H1 %s W%d-O",
		      water,wateratom(0),water); break;
      case 3: sprintf(line,"aa W%d-H2 %s W%d-O",
		      water,wateratom(1),water); break;
      case 4: /* end of including water */
	free(line); line=oldline; oldline=NULL;
	/*** ugly patch: reordering TIP3P water from OHH to HHO ***/
	{ site_t *se=site+spec->ns,sx=se[-3];
	  se[-3]=se[-2],se[-2]=se[-1],se[-1]=sx;
	  if (se[-1].nnbr!=2) ERROR(("internal"))
				else se[-1].nbr[0]--,se[-1].nbr[1]--;
	  if (se[-2].nnbr!=1) ERROR(("internal")) else se[-2].nbr[0]+=2;
	  if (se[-3].nnbr!=1) ERROR(("internal")) else se[-3].nbr[0]+=2;
	  prt("! water reordered to H1=%d H2=%d O=%d",
	      spec->ns-3,spec->ns-2,spec->ns-1); }
      default:
	line[0]=0;
	if (!fgets(line,128,edit)) strcpy(line,"."); /* EOF=end mark */
	prt_("! %s",line); }
    if ( (id=strchr(line,'\n')) ) *id=0;
    strcpy(Line,line);
    for (l=line; *l==' ' || *l=='\t'; l++);
    key=strtok(l,sep);

    /*** main edt parser ***/
    parsed=1;
    if (!key || !*key || key[0]=='!') /* no action */;

    else if (key[0]=='?') fputs(
				"\n\
add atom ID, bond to ID1, ID2..., position calc. later (=mark as missing):\n\
                aa ID TYPE CHARGE [ID1 [ID2 [...]]]\n\
add free atom ID, place DIST from ID1 out of its nbrs, or at (X,Y,Z):\n\
                af ID TYPE CHARGE [ID1 [[-]DIST]]\n\
                af ID TYPE CHARGE = X Y Z\n\
add TIP3 water, place DIST from ID1 out of its nbrs, number from #:\n\
                aw [ID1 [DIST [#]]]\n\
                aw = X Y Z\n\
probe by free atom or TIP3 water:\n\
                pf ID TYPE CHARGE [GRID [SHELL [NforPLB]]] \n\
                pw [GRID [SHELL [NforPLB]]]\n\
add aliphatic hydrogens:\n\
                aah\n\
move hydrogens to their atoms (only reorders mol-file):\n\
                mh\n\
remove atom ID: ra ID\n\
remove molecule:rm ID_of_any_atom\n\
add bond:       ab ID1 ID2\n\
remove bond:    rb ID1 ID2\n\
change ID:      i ID NEWID\n\
change charge:  q ID CHARGE\n\
change type:    t ID TYPE\n\
print info:     p ID\n\
reaction:       react REA-FILE\n\
force field (editing ignored if no match, empty matches everything):\n\
                ff PARSET\n\
end of editing: end  OR  .\n\
wildcards are  ?  for 1 char,  *  at the end of ID\n"
				,stderr);

    else if (!strcmp(key,"end") || !strcmp(key,".")) {
      spec->wrmol++;
      checkneighbors(spec,spec->ns);
      free(line);
      return; }

    else if (!strcmp(key,"ff") || !strcmp(key,"parameter_set")) {
      /* Force Field */
      char *c=strtok(NULL," \t=\n");
      if (!c || *c==0) c="NONE";
      ff=!strcmp(c,"NONE") || !strcmp(PARfn,c);
      prt("! force field set to %s - editing %s",c,ff?"on":"off"); }
    else parsed=0;

    if (parsed) ;
    else if (!ff) /* all commands ignored because force field does not match */;
    else if (!strcmp(key,"aa") || !strcmp(key,"af")
	                       || (probe=!strcmp(key,"pf"))
	     ) { /* Add [Free] Atom */

    int free=key[1]=='f';

    if (spec->ns>=maxns) { /* must reallocate! */
      site_t *old=site;
      maxns += (spec->ns+23)/8;
      alloc(spec->site,maxns*sizeof(site_t));
      copy(spec->site,old,spec->ns*sizeof(site_t));
      site=spec->site;
      free(old); }

    if (!(id=strtok(NULL,sep))) {
      ERROR(("%s: nothing to add",Line)) continue; }
    if (!(type=strtok(NULL,sep))) {
      ERROR(("%s: no type",Line)) continue; }
    if (!(charge=strtok(NULL,sep))) {
      ERROR(("%s: no charge",Line)) continue; }

    loop (i,0,MAXNBRID) nbrid[i]=NULL;
    i=0;
    while ( (nbrid[i]=strtok(NULL,sep)) ) {

#if 0 /* giving up parameter test... */
      if (++i>=(free?MAXVAL:2+probe)) {
        ERROR(("%s: too many %s",Line,free?"parameters":"atoms to bond"))
          continue; }
#else /* new simplified test */
      if (++i>=MAXNBRID) {
        ERROR(("%s: too many parameters",Line))
        continue; }
#endif
    }
    nnbr=min(i,MAXVAL);

    if ((i=findsite(spec,id,0))>=0) {
      ERROR(("%s already present (site [%d])",id,i))
      continue; }

    s=&site[spec->ns];
    memset(s,0,sizeof(site_t));

    alloczero(s->id,strlen(id)+1); strcpy(s->id,id);
    s->type=atomn(type);
    s->charge=atof(charge);
    s->nnbr=free?0:nnbr;
    s->keep=WANTED;
    s->r[0]=UNDEFATOM;

    if (free) { /* add free (non-bonded) atom */
      float dist=3,rr;
      int absr=nbrid[0] && !strcmp(nbrid[0],"=");

      if (probe) {
        if (nbrid[0]) spec->probe.d=atof(nbrid[0]);
        else spec->probe.d=1;
        if (nbrid[1]) spec->probe.shell=atof(nbrid[1]);
        else spec->probe.shell=5;
        if (nbrid[2]) spec->probe.show=atoi(nbrid[2]);
        else spec->probe.show=0;
        spec->probe.i=spec->ns;
        if (spec->probe.ns && !spec->opt_m) {
          ERROR(("molecular probe requires nonzero option -m"))
          spec->opt_m=-100; }
        spec->probe.ns++; /* DIRTY: command pw sends here 2 */
        if (spec->opt_y) {
          WARNING(("option -y (center molecule) turned off"))
          spec->opt_y=0; }
        if (spec->opt_w) {
          WARNING(("option -w (write) turned off"))
          spec->opt_w=0; }
        probe=0; }
      else if (nbrid[0] && (absr || (i=findsite(spec,nbrid[0],1))>=0)) {

        if (!absr) if (nbrid[1] && nbrid[1][0]) dist=atof(nbrid[1]);

        if (absr) {
          int j,err=1;
          loop (j,0,3)
            if (nbrid[j+1]) s->r[j]=atof(nbrid[j+1]);
            else {
              if (err) WARNING(("%s: missing %c coordinate",Line,j+'x'))
              err=0;
              s->r[j]=0; } }

        else if (dist>0 && fsite->nnbr) {
          float R[3];
          int nnbr=0,j;

          VO(R,=0)
          loop (j,0,fsite->nnbr)
            if (fsite->keep!=WANTED) {
              VV(R,+=site[fsite->nbr[j]].r)
              nnbr++; }
          if (!nnbr || fsite->keep==WANTED) {
            ERROR(("selected site(s) have unknown positions"))
            continue; }
          VVV(R,=fsite->r,-1.0/nnbr*R)
          rr=sqrt(SQR(R));

          if (rr<0.2) {
            /* flat case */
            double r[3],rx[3],dmax=6;
          newr:
/*.....            rr=dist/sqrt(rndball(r)); VO(r,*=rr)*/
            rndsphere(r); VO(r,*=dist)
            loop (j,0,fsite->nnbr)
              if (fsite->keep!=WANTED) {
                VVVV(rx,=fsite->r,+r,-site[fsite->nbr[j]].r)
                if (SQR(rx)<dmax) { dmax*=0.9999; goto newr; } }
            VV(R,=r) }
          else
            dist/=rr;
          VVV(s->r,=fsite->r,+dist*R) }

        else {
          double R[3];
          dist/=sqrt(rndball(R));
          VVV(s->r,=fsite->r,+dist*R)
          prts_("! randomly on sphere\n");
          if (ble_file) fputs("! randomly on sphere\n",stderr); }

        s->keep=FREE; }
      }

    else { /* add (bonded) atom */
      loop (n,0,nnbr) {
        i=s->nbr[n]=findsite(spec,nbrid[n],1);
        if (site[i].nnbr>=MAXVAL) {
          ERROR(("%s [%d]: too many bonds",site[i].id,i))
          break; }
        else
          site[i].nbr[site[i].nnbr++]=spec->ns; } }

    n=spec->ns++;
    prtsite("!",s,n,"added\n"); }

  else if (!strcmp("aw",key) || (probe=!strcmp("pw",key)) ) { /* Add Water */
    char *c;
    water++; addwater=0;
    oldline=line; alloc(line,option('l'));
    if (!(wid=strtok(NULL,sep))) { wid=wdist=""; continue; }
    if (!(wdist=strtok(NULL,sep))) { wdist="2"; continue; }
    if (!(wabsy=strtok(NULL,sep))) { wabsy=""; continue; }
    if (!(wabsz=strtok(NULL,sep))) { wabsz=""; continue; }
    if ( (c=strtok(NULL,sep)) ) water=atoi(c); }

  else if (!strcmp(key,"ra")) { /* Remove Atom */
    int ndel=0,m,j;

    if (!(id=strtok(NULL,sep))) {
      ERROR(("%s: nothing to remove",Line))
      continue; }

    while ((n=findsite(spec,id,0))>=0) {
      ndel++;

      loop (i,0,spec->ns) {
        s=&site[i];
        loop (m,0,s->nnbr) {
          if (s->nbr[m]>n) s->nbr[m]--;
          else if (s->nbr[m]==n) {
            s->nnbr--;
            loop (j,m,s->nnbr)
              s->nbr[j]=s->nbr[j+1];
            m--; } } }

      spec->ns--;
      prtsite("!",&site[n],n,"removed\n");
      loop (i,n,spec->ns) site[i]=site[i+1]; }

    if (ndel>1) {
      prt_("! %d atoms removed, ns=%d\n",ndel,spec->ns);
      if (ble_file) fprintf(stderr,"! %d atoms removed, ns=%d\n",ndel,spec->ns); }
    if (ndel==0) {
      prts_("!? no match\n");
      if (ble_file) fputs("!? no match\n",stderr); } }

  else if (!strcmp(key,"rm")) { /* Remove Molecule (cluster) */
    int ndel=0,m;

    if (!(id=strtok(NULL,sep))) {
      ERROR(("%s: nothing to remove",Line))
      continue; }

    site=spec->site;
    cluster();

    while ((n=findsite(spec,id,0))>=0) {
      int cl=site[n].clust,ncl,icl;

      ncl=0;
      loop (i,0,spec->ns) if (cl==site[i].clust) ncl++;

      ndel+=ncl;
      loop (i,0,spec->ns) {
        s=&site[i];
        if (s->clust==cl) icl=i;
        loop (m,0,s->nnbr)
          if (s->nbr[m]>n) s->nbr[m]-=ncl; }

      icl++;
      prtsite("!",&site[n],n,"molecule removed\n");
      loop (i,icl,spec->ns) site[i-ncl]=site[i];
      spec->ns-=ncl; }

    if (ndel>1) {
      prt_("! %d atoms removed, ns=%d\n",ndel,spec->ns);
      if (ble_file) fprintf(stderr,"! %d atoms removed, ns=%d\n",ndel,spec->ns); }
    if (ndel==0) {
      prts_("!? no match\n");
      if (ble_file) fputs("!? no match\n",stderr); } }

  else if (!strcmp(key,"rb")) { /* Remove Bond */
    int n1=findsite(spec,strtok(NULL,sep),1);
    int n2=findsite(spec,strtok(NULL,sep),1);

    if (n1<0 || n2<0) continue;

    s=&site[n1];
    loop (i,0,s->nnbr) if (s->nbr[i]==n2) {
      s->nnbr--;
      while (i<s->nnbr) { s->nbr[i]=s->nbr[i+1]; i++; }
      goto OK1; }

    WARNING(("not connected"))
    continue;

    OK1:
    prtsite("!", s,n1,"--");

    s=fsite;
    loop (i,0,s->nnbr) if (s->nbr[i]==n1) {
      s->nnbr--;
      while (i<s->nnbr) { s->nbr[i]=s->nbr[i+1]; i++; }
      goto OK2; }

    ERROR(("internal: bad bond table"))
    OK2:
    prtsite("", s,n2,"disconnected\n"); }

  else if (!strcmp(key,"ab")) { /* Add Bond */
    int n1=findsite(spec,strtok(NULL,sep),1);
    int n2=findsite(spec,strtok(NULL,sep),1);

    if (n1<0 || n2<0) continue;

    /*
      can connect() from blendedt be used here?
      -- consider (re)numbering of sites!
    */
    s=&site[n1];
    if (s->nnbr>=MAXVAL) {
      ERROR(("%s: too many atoms to bond",s->id))
      continue; }
    s->nbr[s->nnbr++]=n2;
    prtsite("!", s,n1,"--");

    s=fsite;
    if (s->nnbr>=MAXVAL) {
      ERROR(("%s: too many atoms to bond",s->id))
      continue; }
    s->nbr[s->nnbr++]=n1;
    prtsite("", s,n2,"connected\n"); }

  else if (!strcmp(key,"i")) { /* change Id */
    n=findsite(spec,strtok(NULL,sep),1);
    if (n<0) continue;
    id=strtok(NULL,sep);
    if (!id || !*id)
      ERROR(("no new ID"))
    else {
      strcpy(fsite->id,id);
      prtsite("!",fsite,n,"\n"); } }

  else if (!strcmp(key,"q")) { /* change charge */
    n=findsite(spec,strtok(NULL,sep),1);
    if (n<0) continue;
    charge=strtok(NULL,sep);
    if (!charge || !*charge)
      ERROR(("no new CHARGE"))
    else {
      fsite->charge=atof(charge);
      prtsite("!",fsite,n,"\n"); } }

  else if (!strcmp(key,"t")) { /* change Type */
    n=findsite(spec,strtok(NULL,sep),1);
    if (n<0) continue;
    type=strtok(NULL,sep);
    if (!type || !*type)
      ERROR(("no new TYPE"))
    else {
      fsite->type=atomn(type);
      prtsite("!",fsite,n,"\n"); } }

  else if (!strcmp(key,"react")) { /* Print */
    id=strtok(NULL,sep);
    checkneighbors(spec,spec->ns);
    react(spec,id);
    _n }

  else if (!strcmp(key,"p")) { /* Print */
    id=strtok(NULL,sep);

    nnbr=0;
    loop (i,0,spec->ns) if (match((s=&site[i])->id,id)) {
      prt_("%4d: %8s %4s  q=%7.4f  chir=%d  -> ",
           i,s->id,atom[s->type].name,s->charge,s->chir);
      loop (n,0,s->nnbr) prt_(" %d",s->nbr[n]);
      if (ble_file) {
        fprintf(stderr,"%4d: %8s %4s  q=%7.4f  chir=%d  -> ",
                i,s->id,atom[s->type].name,s->charge,s->chir);
        loop (n,0,s->nnbr) fprintf(stderr," %d",s->nbr[n]); }
      nnbr++;
      _n if (ble_file) fputc('\n',stderr); }
    if (!nnbr) {
      prts_("sorry, no such atom\n");
      if (ble_file) fputs("sorry, no such atom\n",stderr); } }

  else if (!strcmp(key,"aah")) { /* Add Aliphatic Hydrogens */

    /* this is CHARMM: to be read one day from a file */
#define NHtab 6
    static struct {
      char *ext;  /* extended atom type */
      char *full; /* full atom type */
      char *h;    /* H atom type */
      float charge; /* H charge */
      int n;      /* # of H */ } Htab[]={
      {"CH1E","CT", "HA",0.1,1},
      {"CH2E","CT", "HA",0.1,2},
      {"CH3E","CT", "HA",0.1,3},
      {"C6RE","C6R","HA",0.1,1},
      {"C5RE","C5R","HA",0.1,1},
      {"SH1E","ST", "H",0.2,1}};

    int j,t;

    if (spec->ns>=maxns) { /* must reallocate! */
      site_t *old=spec->site;
      maxns = spec->ns*3; /* PESSIMISTIC: IMPROVE !!! */
      alloc(spec->site,maxns*sizeof(site_t));
      copy(spec->site,old,spec->ns*sizeof(site_t));
      site=spec->site;
      free(old); }

    loop (i,0,spec->ns) { /* WARNING: spec->ns not constant */
      s=site+i;
      t=s->type;
      loop (j,0,NHtab) if (t==atomn(Htab[j].ext))
        addH(spec,i,Htab[j].full,Htab[j].n,Htab[j].h,Htab[j].charge);
      }
    } /* aah */

  else if (!strcmp(key,"mh")) { /* Move Hydrogens */

    again:

    loop (i,0,spec->ns) {
      s=site+i;
      if (atom[s->type].mass<1.5 && s->nnbr==1) {
        /* H found */
        int k;

        n=s->nbr[0];
        if (n<i) { movesite(spec,i,n); goto again; }
        else loop (k,i+1,n) if (site[k].nnbr!=1 || site[k].nbr[0]!=n) {
          movesite(spec,i,n); goto again; } }

      }
    }
  else
    ERROR(("%s unknown command",key))
  }

/* no return point here! */
}

void readjet(species_t *spec,bond_t **b0ptr,fixdih_t **fd0ptr) /**** readjet */
/***
    adds artificial constraints from file fn.jet
    (fn must not have ext, but space for it)
***/
{
  FILE *f;
  bond_t *b;
  char *name[5],*tok;
  int indx[4];
  char *line;
  int ii,itok,i,j,ncon=0;
  site_t *site=spec->site;
  fixdih_t *fd;

  strcpy(spec->ext,".jet");

  if ( !(f=fopen(spec->fn,"rt"))) {
    fprintf(stderr,"\
No %s file with additional bonds or dihedrals to fix: add them now!\n\
Define bond: (ID1,ID2=ID's of atoms, K=force const [kcal/mol], BOND_LENGTH [AA])\n\
  ID1 ID2 K BOND_LENGTH\n\
Define dihedral to fix: (chain ID1/ID2\\ID3/ID4, angle PHI in degree)\n\
  ID1 ID2 ID3 ID4 PHI\n\
End data by . or end or Ctrl-D\n",
            spec->fn);
    f=stdin; }
  alloc(line,option('l'));

  while (fgets(line,option('l'),f)) {
    tok=strtok(line," \n\r\t");
    if (!strcmp(tok,".") || !strcmp(tok,"end")) break;

    itok=0;
    for (itok=0; tok && itok<5 && tok[0]!='!'; itok++,tok=strtok(NULL," \n\r\t"))
      name[itok]=tok;
    if (!itok) break;
    if (itok<4)  {
      ERROR(("%s: bad format of line:\n%s",spec->fn,line))
      goto ignore; }
//      ierr=sscanf(line,"%s%s"irealfmt irealfmt,name[0],name[1],&K,&r);
    if (itok<4) {
      ERROR(("%s: bad format of line:\n%s",spec->fn,line))
      goto ignore; }

    ii=itok==4?2:4;

    loop (j,0,ii) {
      indx[j]=-1;
      loop (i,0,spec->ns)
        if (!strcmp(site[i].id,name[j])) indx[j]=i; }

    loop (j,0,ii) if (indx[j]<0) {
      prts_(line);
      ERROR(("bad atom id in %s",spec->fn))
      goto ignore;}

    ncon++;

    if (ii==2) {
      ralloc(b,sizeof(bond_t));
      loop (j,0,2) b->indx[j]=indx[j];
      b->parm.K=atof(name[2]); b->parm.Ki2=-2*b->parm.K;
      b->parm.length=atof(name[3]);
      b->next=*b0ptr;
      *b0ptr=b; }
    else {
      ralloc(fd,sizeof(fixdih_t));
      loop (j,0,4) fd->indx[j]=indx[j];
      fd->phi=PI/180*atof(name[4]);
      fd->next=*fd0ptr;
      *fd0ptr=fd; }
    
   ignore:; }

  free(line);
  fclose(f);
  prt("! %s: %d constraints added",spec->fn,ncon);
  *(spec->ext)=0;
}

void react(species_t *spec,char *REAname) /*************************** react */
/***
    performs a special form of `chemical reaction'
    the file REAname should contain lines as e.g.
    HP + O <1.2 *C OAC += H-.6 OT+.35 C+.05 OA+.2 > -12
    denotes the following reaction:

            O              OT--H
           /              /
    HP + -C        ->   -C
           \              \
            OAC            OA

    the reaction takes place if the distance of HP and O is less than the
      distance given by the `<' command (here 1.2A, default = previous line,
      if nothing given = 1A)
    *C means that C is the central atom of the group
      (if no * is given, 1st atom after `+' is central)
    += means that partial charges (numbers after atom types) will be added
      to already present partial charges
    = instead of += means assignment of partial charges
    > energy of the reaction follows (negative=exothermic)

    EXTRA lines:
    if the charge redistribution extends over one group, `continuation'
    lines may be added.  They start with `\', e.g.
    \C5RE + N5R H C5R += C5RE+.05 N5R+.1 H+0.05 C5R+0.05
    Here, `+' does not denote a reaction but an atom that has been marked
    in previous reaction step (=line without `\' + any number of lines
    with `\'); it must be bonded to the next atom (here C5RE--N5R) which
    itself must not be marked.  Marking of other atoms in the group is
    irrelevant. Note that any such marking is canceled by a new reaction
    (=line without `\').  The above example changes charges in the group:

    H--N5R---C5R
        /
    C5RE

    where C5RE has been afected by a reaction but N5R not
***/
{
  FILE *f;
  int nnbr,center;
  char *A;
  char *B[MAXVAL];
  char *newA,*newB[MAXVAL];
  double newchargeA,newchargeB[MAXVAL],dist=1,zero_energy=0;
  char *line,*linea,*tok,*line0,*eq,*sg;
  int i,inbr,j,add,extra;
  site_t *site=spec->site;

  if (!REAname) return;

  if ( !(f=fopen(REAname,"rt"))) {
    ERROR(("No reaction file %s\n",REAname))
    REAname=NULL; return; }

  alloc(linea,option('l')); line=linea;
  alloc(line0,option('l'));

  while (fgets(line,option('l'),f)) if (line[0]!='!' && line[0]!='\n') {
    strcpy(line0,line);
    if (!strcmp(line,".\n") || !strcmp(line,"end\n")) break;

    if ( (extra=line[0]=='\\') ) line++;

    eq=strchr(line,'=');
    if (!eq) {
      ERROR(("%s: `=' or `+=' expected in line:\n%s",spec->fn,line0))
      continue; }
    if ( (add=eq[-1]=='+') ) eq[-1]=0;
    *eq++=0;

    tok=strtok(line," \t+");
    A=tok;

    nnbr=center=0;
    while ( (tok=strtok(NULL," \t+"))) {
      if (tok[0]=='*') {
        center=nnbr; tok++; }
      if (tok[0]=='<') {
        if (tok[1]) dist=sqr(atof(tok+1));
        else dist=sqr(atof(strtok(NULL," \t"))); }
      else {
        if (nnbr>=MAXVAL)
          ERROR(("%s: too many neighbors in line:\n%s",spec->fn,line0))
          B[nnbr++]=tok; } }

    tok=strtok(eq," \t");
    newA=tok;
    if ( (sg=strpbrk(tok,".+-")) ) {
      newchargeA=atof(sg); *sg=0; }
    else
      newchargeA=atof( (tok=strtok(NULL," \t")) );

    loop (i,0,nnbr) {
      tok=strtok(NULL," \t\n");
      if (!tok) ERROR(("%s: missing data on r.h.s. in line:\n%s",spec->fn,line0))
        newB[i]=tok;
      if ( (sg=strpbrk(tok,".+-")) ) {
        newchargeB[i]=atof(sg); *sg=0; }
      else
        newchargeB[i]=atof( (tok=strtok(NULL," \t\n")) ); }

    /* scanning the reaction energy (> zero_energy) */
    if (tok) if ( (tok=strtok(NULL," \t\n")) ) if (tok[0]=='>') {
      if (tok[1]) zero_energy=atof(tok+1);
      else zero_energy=atof(strtok(NULL," \t")); }
    /* ... end of line scan */

    loop (i,0,spec->ns) {
      if (!extra) site[i].rea=0;

      if (match(atom[site[i].type].name,A) && (!extra || site[i].rea)) {
        double minrr=9e9,rr,*r=site[i].r;
        int ii=-1,ic,k,cnbr,l;
        int matchi[MAXVAL];

        if (extra) loop (inbr,0,site[i].nnbr) {
          j=site[i].nbr[inbr];
          if (match(atom[site[j].type].name,B[0]) && !site[j].rea) ii=j; }

        else loop (j,0,spec->ns)
          if (i!=j && match(atom[site[j].type].name,B[0])) {
            rr=SQRD(r,site[j].r);
            if (rr<dist && rr<minrr) {
              minrr=rr;
              if (ii>=0 && option('v')&8)
                prt("! %d:%dx%d not unique for %s",i,j,ii,line0);
              ii=j; } }

        if (ii<0) continue;

        if (center==0) ic=ii;
        else if (site[ii].nnbr!=1) continue;
        else ic=site[ii].nbr[0];

        loop (j,0,nnbr) matchi[j]=-1;
        if (!match(atom[site[ic].type].name,B[center])) continue;
        matchi[center]=ic;
        if (center) matchi[0]=ii;

        loop (j,0,nnbr) if (matchi[j]<0) {
          cnbr=-1;
          loop (k,0,site[ic].nnbr) {
            cnbr=site[ic].nbr[k];
            loop (l,0,nnbr) if (matchi[l]==cnbr) goto cont0; /* already matched */
            if (match(atom[site[site[ic].nbr[k]].type].name,B[j])) goto matched;
           cont0:; }
          goto cont;
         matched:
          if (cnbr<0) ERROR((""))
          matchi[j]=cnbr; }

        site[i].type=atomn(newA);
        if (!add) site[i].charge=0;
        site[i].charge+=newchargeA;

        loop (j,0,nnbr) {
          site[matchi[j]].type=atomn(newB[j]);
          if (!add) site[matchi[j]].charge=0;
          site[matchi[j]].charge+=newchargeB[j];
          site[matchi[j]].rea=1; }

        if (!extra)
          if (site[ii].nnbr>=MAXVAL || site[i].nnbr>=MAXVAL)
            ERROR(("%s, line\n%s site %d too many neighbors",spec->fn,line0,i))
          else {
            site[ii].nbr[site[ii].nnbr++]=i;
            site[i].nbr[site[i].nnbr++]=ii;
            prt_("\n! %s: %3d %-4s %s - %3d %-4s %s bonded, zero_energy=%g",
                 REAname,
                 i,atom[site[i].type].name,site[i].id,
                 ii,atom[site[ii].type].name,site[ii].id,zero_energy);
            spec->zero_energy+=zero_energy; }

        spec->wrmol++;
       cont:; }
      } /* all sites in reaction line */
    line=linea; } /* fget(line) */

  free(line0);
  free(linea);
  fclose(f);
  *(spec->ext)=0;
  checkneighbors(spec,spec->ns);
}

#define T(I) site[I].type

void autobonds(species_t *spec,int rg)
/* it somehow repeats the code in buildforcefield... */
{
  site_t *site=spec->site;
  int ns=spec->ns,i,j,n,nb=0;
  double q=fabs(rg/100.),m=0,d;
  bond_t *bond;

  if (rg==0) return;

  if (q<=1) WARNING(("create bonds for distances <= 100%% bond lengths"))

  for (bond=bond0; bond; bond=bond->next)
    if (bond->parm.K>0)
      Max(m,bond->parm.length)
  m=Sqr(m*q);

  /* erase all bonds (to re-create them) */
  if (rg<0) loop (i,0,ns) site[i].nnbr=0;

  loop (i,0,ns) if (site[i].r[0]<UNDEFATOM)
    loop (j,0,i) if (site[j].r[0]<UNDEFATOM) {

      /* test whether already bonded */
      loop (n,0,site[i].nnbr) if (site[i].nbr[n]==j) goto alreadybonded;

      d=Sqr(site[i].r[0]-site[j].r[0])
       +Sqr(site[i].r[1]-site[j].r[1])
       +Sqr(site[i].r[2]-site[j].r[2]);

      if (d<m) {
        for (bond=bond0; bond; bond=bond->next) if (bond->parm.K>0)
          if ( (bond->indx[0]==T(i) && bond->indx[1]==T(j))
            || (bond->indx[1]==T(i) && bond->indx[0]==T(j)) ) {
            if (d<Sqr(bond->parm.length*q)) {
              if (site[i].nnbr>=MAXVAL || site[j].nnbr>=MAXVAL)
                ERROR(("autobonds: site %d or %d too many bonds\n\
(bad configuration or too long bond limit (option -d))",i,j))
              else {
		nb++;
		if (rg>0 || option('v')&8) {
		  prtsite("! autobond:",&site[i],i,"-");
		  prtsite("",&site[j],j,"\n"); }
		site[i].nbr[site[i].nnbr++]=j;
		site[j].nbr[site[j].nnbr++]=i; } }
	    break; } }
    alreadybonded:; }
  prt("! autobond: %d bond%s %s for dist<%g*r0",
      nb,"s"+(nb==1),rg<0?"created":"added",q);
}

