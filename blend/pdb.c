/*
  pdb -> mol converter (c) J.Kolafa 1994-96, 2000, 2001, 2009, 2024 (bug fix)
  This software is obsolete and will not be developed
  VERSION is #defined in pdbbasic.h
*/

#include "ground.h"
#include "options.h"
#include "forscanf.h"

#include "pdbbasic.h"
#include "pdbpdb.h"
#include "pdbinfo.h"
#include "pdbmol.h"
#include "pdbconv.h"
#include "pdbrep.h"

#define x badoption
int optionlist[32] =
/* ` a b c d e f g h i j k l m n o p   q r s t u v w x y z { | } ~   */
  {x,0,2,0,0,0,0,0,0,0,x,x,4,1,0,0,0,100,0,0,0,0,1,x,0,x,1,x,0,x,x,x};
/* @ A B C D E F G H I J K L M N O P   Q R S T U V W X Y Z [ \ ] ^ _ */
#undef x


/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

int main(int narg, char **arg) /************************************** main */
{
  int i,j,n=0;
  site_t *s,*pat;
  struct id_s *id;
  ssbond_t *ss;
  connect_t *cc;
  reconnect_t *re;
  residue_t *res,*resnext,**lastptr,*r;
  int breakchain=1;
  int naltloc,merge,no;
  char *cter="CTER",*nter="NTER";
  char *dir;
  int plbs=0,toframe=-1,byframe=1;

  initscroll(0);

  /*** argument parsing ***/
  loop (i,1,narg)

    if (arg[i][0]=='-') switch (arg[i][1]) {
      case 'p': {
	char *c=strchr(arg[i],':');

	option('p')=atoi(arg[i]+2);
	if (c) {
	  plbs++; c++;
	  toframe=atoi(c);
	  c=strchr(c,':');
	  if (c) byframe=atoi(c+1); 
	  if (byframe<=0) ERROR(("illegal byframe=%d",byframe)) } 
	break; }
      case 'h': nter=arg[i]+2; option('h')++; break;
      case 'e': cter=arg[i]+2; option('e')++; break;
      case 'o': {
        omitrsd_t *omitrsd;
        if (strlen(arg[i]+2)<1 || strlen(arg[i]+2)>3) {
          ERROR(("bad -oRSD : 1-3 char expected"))
          break; }
        alloczero(omitrsd,sizeof(omitrsd_t));
        strcpy(omitrsd->resnm,arg[i]+2); option('o')++;
        omitrsd->next=omitrsd0;
        omitrsd0=omitrsd;
        break; }
      case 'f': parameter_set=arg[i]+2; break;
      case 'r': rsddir=arg[i]+2; break;
      default: getoption(arg[i]) /* see options.h */ }
  
    else {
      /* no -option ==> file name */
      pdbname=molname;
      molname=arg[i]; }
  
  initscroll((option('s')+(option('s')==1)*30)*1024);
  
  if (option('v') || !molname)
    prts("*** pdb->mol converter V"VERSION" *** (c) J.Kolafa 1994- ***");
  
  if (!molname) help();
  
  if (rsddir && !parameter_set) parameter_set=rsddir;
  
  if (parameter_set) {
    dir=getrsddir();
  
    if (rsddir && strcmp(dir,rsddir))
      prt("option -r%s overrides `rsddir %s' in parameter file %s",rsddir,dir,fn);
    if (!rsddir) rsddir=strdup(dir); }
  else {
    ERROR(("unknown parameter set, use -fPARSET and/or -rRSDDIR"))
    return 1; }
  
  if (option('p')) {
    if (option('b')*option('t'))
      ERROR(("cannot both bin and text input if -p"))
    if (!pdbname || !strcmp(molname,pdbname))
      if (!plbs)
	ERROR(("must give separate input and output .pdb names if -p")) }

  if (plbs) {
    if (option('b') && option('b')!=2) {
      prt("option -b2 forced because range of frames specified in option -p");
      option('b')=2; }
    if (option('t') && option('t')!=2) {
      prt("option -t2 forced because range of frames specified in option -p"); 
      option('t')=2; } }
  
  if (option('v'))
    prt("force field (parameter set) = %s, residues = %s/*%s",
        parameter_set,rsddir,rsdext);
  
  {
    char *swapname=molname;
    if (pdbname) molname=pdbname;
    readPDB();
    molname=swapname; }
  
  readSEL();
  
  if (option('n')>0) {
    /* select neutral terminal residues (NEW in V1.3i) */
    if (!option('e')) cter="CTERH";
    if (!option('h')) nter="NTERN"; }
  else if (option('n')<0) {
    /* select terminal residues with counterions */
    if (!option('e')) cter="CTERNA";
    if (!option('h')) nter="NTERCL"; }
  /* ... nter will be modified for PRO,GLY */
  
  if (option('a')<0 || option('a')>15) {
    if (option('v')) {
      WARNING(("invalid option -a%d",option('a')))
      prt("-a0 set"); }
    option('a')=0; }
  
  if (option('v')) {
    prt_(
      option('n')>1 ? "use neutral residues" : 
      option('n')<0 ? "use charged residues neutralized by counterions Na+ Cl-" :
      "use charged residues");
    prt(" (incl. termini)");
    if (option('x')) prt("wildcard A for rotating groups enabled");
    prt("head (N-ter) patch=%s",nter);
    prt("end  (C-ter) patch=%s",cter);
  
    prt("col. 17 is %s for HETATM, %s for ATOM",
      option('a')&4?"part of id":"alt.loc.",
      option('a')&8?"part of id":"alt.loc."); }
  
  strcpy(optioninfo,"!");
  loop (i,0,5) {
    char c="acgnq"[i]; /* five interesting options */
    sprintf(optioninfo+strlen(optioninfo)," -%c%d",c,option(c)); }
  sprintf(optioninfo+strlen(optioninfo)," %s %s\n",nter,cter);
  
  /* solve alternate locations */
  naltloc=altlocations();
  
  /* find CYS-CYS bonds from S-S distances */
  if (option('c')) findCYSCYS();
  
  /* read mol-files */
  looplist (res,reshead) {
  
    if (!strcmp(res->resnm,"CYS")) {
      /* use CYSS if CYS creates the S-S bridge */
      looplist (ss,sshead) loop (j,0,2)
        if (ss->chain[j]==res->chain && ss->resno[j]==res->resno) {
          res->site=readRSD("CYSS",&res->moltype,res->resno,res->chain);
          if (res->moltype!=AMINOACID) ERROR(("CYSS is not AMINOACID"))
          /* remember SG sites to connect */
          looplist (s,res->site) looplist (id,s->ids)
            if (!strcmp(id->id,"SG")) ss->site[j]=s;
          goto CYSSread;
          }
      if (option('v')) { prtres(res); prt(" is not SS bonded"); }
      }
    res->site=readRSD(res->resnm,&res->moltype,res->resno,res->chain);
    CYSSread:; }
  
  /* CTER and NTER patches */
  for (res=reshead,lastptr=&reshead; res; lastptr=&res->next,res=resnext) {
    resnext=res->next;
  
    if (res->moltype==NTER) {
      /* N-terminus specified in the PDB file */
      if (resnext->moltype!=AMINOACID) 
        ERROR(("residue %4s %c %d %c: aminoacid does not follow NTER",
               res->resnm,res->chain,res->resno,res->resins))
      merge=option('z')==0;
      res->site=patch(res->site,resnext->site,merge);
      if (merge) *lastptr=resnext;
      if (option('v'))
        prt("pdb-specified patch %s applied to following %s",res->resnm,resnext->resnm);
      res->ter=-1; /* chain head marked */ }
  
    else if (breakchain) {
      if (res->moltype==CTER)
        ERROR(("residue %4s %c %d %c: C-terminus first in chain",
               res->resnm,res->chain,res->resno,res->resins))
  
      if (res->moltype==AMINOACID) {
        /* N-terminus (auto or by -h option) */
        char locnter[8];
        char *actualnter=nter;
  
        if (!strcmp(res->resnm,"PRO") || !strcmp(res->resnm,"GLY")) {
          /* PRO and GLY as 1st residue require special patches */
          if (!option('h')) {
            sprintf(locnter,"%sP%s",res->resnm,nter+4); /* e.g., PROPN */
            actualnter=locnter;
            if (option('v'))
              prt("%s as N-terminus: %s replaced by %s",
                  res->resnm,nter,actualnter); }
          else if (option('v'))
            prt("Note: -h%s specified for %s as N-terminus",
                actualnter,res->resnm); }

        merge=option('z')!=2;
        pat=readRSD(actualnter,&test_nterp,-1,0);
        pat=patch(pat,res->site,merge);
        if (!merge) {
          r=newpatchres(res,actualnter);
          *lastptr=r;
          r->next=res;
          r->site=pat;
          r->atom=NULL;
          r->moltype=NTER;
          r->ter=-1; /* chain head marked */ } } }
  
    breakchain = !resnext ||
      resnext && ( option('g') && resnext->resno-res->resno!=1
                   || resnext->chain != res->chain
                   || resnext->moltype != res->moltype 
                      && !(res->moltype==NTER && resnext->moltype==AMINOACID)
                      && !(res->moltype==AMINOACID && resnext->moltype==CTER) );
  
    if (res->moltype==CTER) {
      breakchain=1;
      WARNING(("unexpected CTER"))
      res->ter=1; }
      
    else if (resnext && resnext->moltype==CTER) {
      /* C-terminus specified in the PDB file */
      if (res->moltype!=AMINOACID) 
        ERROR(("residue %4s %c %d %c: aminoacid does not precede CTERP",
               res->resnm,res->chain,res->resno,res->resins))
  
      merge=option('z')==0;

      resnext->site=patch(resnext->site,res->site,merge);

#if 0
  { /* DEBUGGING DUMP */
    residue_t *x; looplist (x,reshead) {
      site_t *s;
      atom_t *a;
    prt("========= %s:%d ==========",x->resnm,x->moltype);
    for (a=x->atom,s=x->site; a||s; ) {
      if (a) { prt_("%3d:%-3d %-3s",a->line,a->atomno,a->id); a=a->next; }
      else prt_("           ");
      if (s) {
	prt_("  %3d %-3s [%s]",s->n,s->ids->id,s->type);
	if (s->patch[0]) prt_("   patch=%s",s->patch);
	s=s->next; }
      _n
    } } }
#endif

      if (option('v'))
        prt("pdb-specified patch %s applied to previous %s",resnext->resnm,res->resnm);
      if (merge) res->next=resnext->next;
      resnext=resnext->next;
      breakchain=1; }
  
    else if (breakchain) {
      if (res->moltype==NTER)
        ERROR(("residue %4s %c %d %c: N-terminus last in chain",
               res->resnm,res->chain,res->resno,res->resins))
  
      if (res->moltype==AMINOACID) {
        /* C-terminus (auto or by -e option) */
        pat=readRSD(cter,&test_cterp,-1,0);
        merge=option('z')!=2;
	pat=patch(pat,res->site,merge);
        if (!merge) {
	  r=newpatchres(res,cter);
	  res->next=r;
	  r->next=resnext;
	  r->site=pat;
	  r->atom=NULL;
	  r->ter=1; 
	  r->moltype=CTER; } } }
    }
  
  
  /***
    looks for match of atom-id in the PDB ATOM-line
    with a site in last read residue
  ***/
  looplist (res,reshead) {
    site_t *s;
    atom_t *a;
    char *aid;
    int l,j,f,ff;
  
    looplist (a,res->atom) {
      looplist (s,res->site) looplist (id,s->ids) {
        /* was patch: exact match required, but H:HN accepted */
	if (s->patch[0]) {
	  if (!strcmp(s->patch,a->id)) goto match;
	  if (!strcmp(s->patch,"HN")
	      && (!strcmp(a->id,"H")) || !strcmp(a->id,"1H") || !strcmp(a->id,"H1") )
	      goto match; }
  
	/* NEW in V1.3m (because of hydrogens): 
	   1HX in PDB matches HX1 in RES.rsd  */
	if (a->id[1]=='H' && strchr("123",a->id[0])) {
	  l=strlen(id->id);
	  if (id->id[l-1]==a->id[0]
	      && strlen(a->id)==l
	      && !memcmp(a->id+1,id->id,l-1)) goto match; }

        /* try 1XXX in PDB (a->id) to match to XXX in *.rsd (s->id) */
        ff=a->id[0]=='1' && id->id[0]>='A';
  
        loopto (f,0,ff) {
          aid=a->id+f;
  
          /* exact match */
          if (!strcmp(id->id,aid)) goto match;

          /* XX1 in PDB (a->id)[w/o leading 1] matches XX in RES.rsd (id->id) */
          l=strlen(aid)-1;
          /* NEW: test aid[l-1] added */
          if (aid[l-1]>='A' 
              && aid[l]=='1'
              && id->id[l]==0
              && !memcmp(id->id,aid,l)) goto match;
  
          /* XX in PDB (a->id)[w/o leading 1] matches XX1 in RES.rsd (id->id) */
          l=strlen(id->id)-1;
          if (id->id[l-1]>='A'
              && id->id[l]=='1'
              && aid[l]==0
              && !memcmp(id->id,aid,l)) goto match;
          }
  
        /* unique element s in the res file may match a->id without numbers */
        if (option('u')) if (s->unique)
          if (!strcmp(removedigits(id->id),removedigits(a->id))) {
            if (option('v')) {
              prtres(res);
              prt(": %s matches %s because %s is unique",
                  id->id,a->id,removedigits(id->id)); }
            goto match; } 
  
        /* (added 6/95) enable A to match any atom for rotating groups */
        if (option('x'))
          if (a->id[0]=='A' && !strcmp(a->id+1,id->id+1)) goto match; 
  
        /* other matching conditions ? */
  
        /*** no match ***/
        continue; /* s->next, id->next */
  
        match:
  
          /* matches! */
          s->r=a->r; /* reference 3D coordinates of atom */
          s->line=a->line; /* pdb line number */
  
          /* find sites to connect remembered from the CONECT statement */
          looplist (cc,connecthead) loop (j,0,2)
            if (cc->atomno[j]==a->atomno) cc->site[j]=s;
          goto done;
        }
  
      /* some non-matches ignored */
      if (sitebyoption(a->id,res,option('i'))) {
        omit_t *omit;
        if (option('v')) {
          prtres(res);
          prt(": atom %d:%s not matched in *.rsd file - ignored",
              a->atomno,a->id); }
        alloc(omit,sizeof(omit_t));
        omit->atomno=a->atomno;
        omit->next=omithead; omithead=omit;
        goto done; }
  
      WARNING(("residue %4s %c %d %c: atom %d:%s not matched in *.rsd file\n\
*** this may be caused by:\n\
  * unknown atom names in the RSD file - check the matching atoms in *.rsd\n\
  * atoms in a patch (like OXT=OC2) which are included in its resuidue\n\
  * CYS-CYS bond problem (check option -%c)\n\
(note: a few non-matched atoms can be filled in the blend stage)",
             res->resnm,res->chain,res->resno,res->resins,a->atomno,a->id,
             option('x') || a->id[0]!='A' ? 'i' : 'x'))
      done:;
      }
  
    /*
      calculate counterion positions
    */
    if (option('n')<0) {
      vector rc,rr;
      float pushout = ( option('n')==-1 ? 100 : -option('n') ) * 0.01;
      int i;
      site_t *sion,*s;
  # define NION 5
      static struct {
        char resnm[4];
        char center[4];  /* --(reflex)--(center)....(ion) */
        char reflex[4];
        char ion[4]; } table[NION+2]={
      {"ARG","CZ","NE","CL"},    
      {"LYS","NZ","CE","CL"},    
      {"HIS","CE1","CG","CL"},    
      {"ASP","CG","CB","NA"},
      {"GLU","CD","CG","NA"},
      {"",   "N", "CA","CL"}, /* NTER */
      {"",   "C", "CA","NA"}, /* CTER */ };
  
      loop (i,0,NION) 
        if (!strcmp(res->resnm,table[i].resnm)) goto fillion;
  
      if (!option('h') && res->ter==-1) i=NION;       /* NTER with counterion */
      else if (!option('e') && res->ter==1) i=NION+1; /* CTER with counterion */
      else goto nofill;
  
      fillion:
        VO(rc,=UNDEFATOM)
        VO(rr,=UNDEFATOM)
        sion=NULL;
  
        looplist (s,res->site) looplist (id,s->ids) {
          if (!strcmp(id->id,table[i].ion)) sion=s;
          if (!strcmp(id->id,table[i].center)) VV(rc,=s->r)
          if (!strcmp(id->id,table[i].reflex)) VV(rr,=s->r) }
      
        if (rc[0]==UNDEFATOM || rr[0]==UNDEFATOM) 
          ERROR(("residue %4s %c %d %c\n\
  *** atoms %s or %s not found or undefined (cannot add counterion)",
                 res->resnm,res->chain,res->resno,res->resins,
                 table[i].center,table[i].reflex))
        else if (!sion) ERROR(("%s missing in %s%s.rsd",
                               table[i].ion,res->resnm,table[i].ion))
        else {
          float *r;
          alloc(r,sizeof(vector));
          VVV(r,=(pushout+1)*rc,-pushout*rr)
          sion->r=r; }
  
      nofill:;
      } /* option('n')<0 */
    } /* res */
  
  /* renumber sites globally */
  /* this allows applying connect across residues */
  n=0;
  looplist (res,reshead) n+=renumber(res->site,n);
  
  if (option('v')) prt("%d sites",n);
  
  /* connect aminoacid residues */
  no=0;
  for (res=reshead; res; res=resnext) {
    resnext=res->next;
    if (res->ter<=0
        && resnext
        && ((res->moltype==AMINOACID && resnext->moltype==AMINOACID)
          ||(res->moltype==AMINOACID && resnext->moltype==CTER)
          ||(res->moltype==NTER && resnext->moltype==AMINOACID)) ) {
      site_t *Cter=findid(res,"C",NULL),*Nter=findid(resnext,"N","NT");
  
      if (Cter && Nter) {
        connect(Cter,Nter); /* existence checked in findid() */
        no++; }
      } }
  
  if (option('v')) prt("%d peptide C-N bonds",no);
  
  /* connect SG of CYS using data in SSBOND command */
  no=0;
  if (sshead) if (option('v')>0) prt("SSBOND bonds:");
  looplist (ss,sshead) {
    int c=0,j;
    loop (j,0,2) if (!ss->site[j]) {
      c++;
      ERROR(("SSBOND: residue CYS %c %d or atom SG to connect not found",
             ss->chain[j],ss->resno[j])) }
    if (!c) {
      connect(ss->site[0],ss->site[1]);
      no++; } }
  if (no) if (option('v')) prt("total %d SSBOND bonds",no);
  
  /*
    deletes (some of) sites not found in *.pdb
  */
  if (option('d')) {
    no=0;
    looplist (res,reshead)
      looplist (s,res->site) looplist (id,s->ids)
        if (sitebyoption(id->id,res,option('d')) && s->r[0]==UNDEFATOM)
          no+=delete(s);
    if (no) prt("%d atoms deleted",no);
    }
  
  /* connect sites using CONECT data */
  no=0;
  if (connecthead) if (option('v')>0) prt("CONECT bonds:");
  looplist (cc,connecthead) {
    int c=0,j;
    loop (j,0,2) if (!cc->site[j]) {
      omit_t *omit;
  
      c++;
      looplist (omit,omithead)
        if (omit->atomno==cc->atomno[j]) goto site_omitted;
      ERROR(("CONECT: atom number %d to connect not found",cc->atomno[j]))
      site_omitted:; }
    if (!c) {
      connect(cc->site[0],cc->site[1]);
      no++; } }
  if (no) if (option('v')) prt("total %d CONECT bonds",no);
  
  /* reconnect sites using patch data */
  no=0;
  if (reconnecthead) if (option('v')>0) prt("reconnect bonds:");
  looplist (re,reconnecthead) {
    connect(re->site[0],re->site[1]);
    no++; }
  if (no) if (option('v')>0) prt("total %d reconnected bonds",no);
  
  /*
    doublecheck of SS distances: prints warning if there are
    unconnected SG atoms of CYS with distances less than 4A
    tested on final coordinates (with alt loc solved)
  */
  
  {
  residue_t *r1,*r2;
  site_t *s1,*s2;
  int i;
  float dist;
  vector dr;
  
  /* WARNING: no aliases allowed for SG */
  looplist (r1,reshead) if (!strcmp(r1->resnm,"CYS")) {
    looplist (s1,r1->site) if (!strcmp(s1->ids->id,"SG")) {
      for (r2=reshead; r2!=r1; r2=r2->next) if (!strcmp(r2->resnm,"CYS"))
        looplist (s2,r2->site) if (!strcmp(s2->ids->id,"SG")) {
          VVV(dr,=s1->r,-s2->r)
          dist=sqrt(SQR(dr));
          if (dist<=4) {
            loop (i,0,s1->nnbr) if (s1->nbr[i]==s2->n) goto SSbond;
            if (option('v')) {
              prt_("WARNING "); prtres(r1); prtres(r2);
              prt("S-S distance=%.3f and no bond\n        (check option -c!)",
                  dist); }
            SSbond:; } }
      goto SSfound; }
    ERROR(("residue %4s %c %d %c: atom SG not found",
           r1->resnm,r1->chain,r1->resno,r1->resins))
    SSfound:; }
  }
  
  if (option('v')) {
    prt_("%d alternate locations",naltloc);
    if (naltloc) {
      i=option('a')&3;
      prts(i ?
           i==2 ? ": max. occupancy used" : ": average used"
           : ": location A used (check option -a!)"); }
    else _n }
  
  if (abs(option('d'))==2) removesites(freeH,"option -d2 or -d-2");
  
  replacepattern();
  
  checkneighbors();

  if (!option('p')) {
    if (option('m')) writeMOL();
    if (option('t')) write3D("wt",option('t'));
    /* 'b' must be the last because of possible endian change in place */
    if (option('b')) write3D("wb",option('b')); }
  else if (!plbs) {
    if (option('t')) read3D("rt",option('t'),option('p'));
    if (option('b')) read3D("rb",option('b'),option('p'));
    pastePDB(pdbname,-1); }
  else {
    int iframe;

    retonerror=1;
    for (iframe=option('p'); toframe<0 || iframe<=toframe; iframe+=byframe) {
      if (option('t')) if (read3D("rt",option('t'),iframe)) break;
      if (option('b')) if (read3D("rb",option('b'),iframe)) break;
      pastePDB(pdbname,iframe); } }
  
  if (option('v')) {
    /* residue statistics */
    prtcount();
    prtomitted(); }
  
  if (option('s')) {
    fprintf(stderr,"quit by \';\' or scroll (\'$?\' for help)\n");
    getdata enddata }
  
  return 0;
} /* main */
