#include "ground.h"
#include "options.h"
#include "forscanf.h"

#include "pdbbasic.h"
#include "pdbpdb.h"
#include "pdbinfo.h"

struct A_s A;

FILE *pdb;
int pdbline; /* line of pdb-file:
                to be advanced if \n appears in forscanf format */

static struct {
  omitrsd_t *rsd;
  char chain,resins;
  int resno;
} lastomit;

static int omit(void) /************************************************ omit */
{
  omitrsd_t *omitrsd;

  looplist (omitrsd,omitrsd0)
    if (!strcmp(omitrsd->resnm,A.resnm)) {
      omitrsd->l++;
      if (lastomit.rsd!=omitrsd
          || lastomit.chain!=A.chain
          || lastomit.resins!=A.resins
          || lastomit.resno!=A.resno) {
        omitrsd->n++;
        lastomit.rsd=omitrsd;
        lastomit.chain=A.chain;
        lastomit.resins=A.resins;
        lastomit.resno=A.resno; }
      return 1; }

  return 0;
}

static void tilde(char *id) /***************************************** tilde */
/*
  If id contains space(s) inside (I consider names containing spaces
  inside as strange, but the Fortranists from PDB don't mind) they are
  replaced by tildes (~).
  id should not contain trailing spaces - OK for id obtained by forscanf.
*/ 
{
  char *c;

  for (c=id; *c; c++) if (*c==' ') *c='~';
}

static void readATOM(int hetatm) /********************************* readATOM */
/*
  read one ATOM or HETATM line
  pdb file pointer must be just after the keyword
  it is at eol after reading

  fortran format for ATOM:
   (a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,1x,i3,2x,a4,i4) 

  fortran format for HETATM:
   (a6,i5,1x,a4,a1,a3,1x,a1,i5,   3x,3f8.3,2f6.2,1x,i3,2x,a4,i4) 

  NOTE: see forscanf.c for function forscanf, analogous to fscanf with
  fortran-like fixed-length fields
*/
{
  forscanf(pdb,
           "%5i %4s%1a%3s %1a",
           &A.atomno,
                A.id,&A.altloc,A.resnm,
                          &A.chain);

  if ( (hetatm && option('a')&4) || (!hetatm && option('a')&8) )
    if (A.altloc>' ') {
      static char app[2];
      app[0]=A.altloc;
      strcat(A.id,app);
      A.altloc=' '; }

  tilde(A.id);
  tilde(A.resnm);

#if 0
  if (hetatm) {
    forscanf(pdb,"%5i",&A.resno);
    A.resins=0; }
  else
#endif

    forscanf(pdb,"%4i%1a",&A.resno,&A.resins);

  if (forscanf(pdb,"   %8f%8f%8f%6f%6f %3i  %4s%4i\n",
                       &A.r[0],&A.r[1],&A.r[2],&A.occup,&A.bvl,
                                       &A.footnote,
                                            A.ydent,&A.recno)<0)
    ERROR(("unexpected EOF while reading ATOM/HETATM"))
  A.hetatm=hetatm;
  A.line=pdbline;
  pdbline++;
}

static atom_t *newatom(void) /************************************** newatom */
/* to allocate ATOMs read from pdb in longer chunks */
{
  static atom_t *buf;
  static int ibuf=ATOMBUFLEN;

  if (ibuf>=ATOMBUFLEN) {
    allocarray(buf,ATOMBUFLEN);
    ibuf=0; }

  return buf+ibuf++;
}

residue_t *appendres(residue_t *res) /**************************** appendres */
/*
  new residue is allocated and put after res (end of list normally);
  if res==NULL, new list reshead-> is created
  info is copied from PDB-line A and 1st atom is included
  new res is returned
*/
{
  residue_t *r;

  allocone(r);
  if (res) { r->next=res->next; res->next=r; }
  else { reshead=r; r->next=NULL; }

  strcpy(r->resnm,A.resnm); /* residue type in 3-letter code */
  r->chain=A.chain;         /* chain identity letter */
  r->resno=A.resno;         /* residue number in sequence */
  r->resins=A.resins;       /* residue insertions (?) */

  /* the following info will be filled later */
  r->moltype=UNKNOWN;
  r->atom=NULL;
  r->site=NULL;
  r->ter=0;

  count(r->resnm); /* residue statistics */

  if (!includeA(r)) DISASTER((""))

  return r;
}

residue_t *newpatchres(residue_t *res,char *ter)
/*
  new residue, to hold a patch
*/
{
  int l=strlen(ter);
  residue_t *r;

  if (l>3) l=3;
  allocone(r);
  *r=*res; /* name copied from rsd */

  if (option('z')>=1) {
    /* name from patch: 3 letters in UPPERCASE */
    copy(r->resnm,"XXX",3);
    copy(r->resnm,ter,l);
    loop (l,0,3) r->resnm[l]=toupper(r->resnm[l]); }

  return r;
}

int includeA(residue_t *r) /*************************************** includeA */
/*
  if atom A is a member of residue r, atom a is allocated and included
    to r->a list and 1 is returned
  if atom A is not a member of residue r, 0 is returned
*/
{
  atom_t *a,*newa;

  if (r)
    /* if (A.hetatm==(r->moltype!=AMINOACID)) -oops, moltype not known here! */
    if (r->chain==A.chain
        && r->resno==A.resno
        && r->resins==A.resins
        && !strcmp(r->resnm,A.resnm)) {
      newa=newatom();
      newa->atomno=A.atomno; /* atom number in sequence */
      strcpy(newa->id,A.id); /* atom type in BDB convention */
      newa->altloc=A.altloc; /* alternative location */
      newa->occup=A.occup;   /* occupancy */
      newa->line=A.line;     /* pdb line # */
      VV(newa->r,=A.r)       /* cartesian coordinates */
        /* incl. to list head: problems with alt.loc.
           a->next=r->atom;
           r->atom=a; */
        /* incl. to the list end (alt loc calculations require the same order of
           atoms as in PDB ) */
      looplist (a,r->atom)
        if (a->next==NULL) { a->next=newa; goto included; }
      r->atom=newa;
    included:
      newa->next=NULL;
    return 1; }
  
  return 0;
}


void readPDB(void) /************************************************ readPDB */
/*
  Reads the PDB file into program structures
  Proteins are stored in the units by residues which are recognized by
  chain letter, number, and residue change.
*/
{
  residue_t *res=NULL,*r;
  char keyword[7];  /* ATOM, HETATM etc. */

  pdb=openfile("pdb","rt",1);
  pdbline=0;

  while (forscanf(pdb,"%6s",keyword)==1) {

    if (!strcmp(keyword,"ATOM")) { /* ///////////////////////////////// ATOM */
  
      readATOM(0);
      if (omit()) continue;

      if (!strcmp(A.resnm,"OXY")) {
        /* exception: OXY merged to the previous residue */
        if (!res || A.chain!=res->chain) {
          ERROR(("%s %c %d : misplaced OXY terminus",A.id,A.chain,A.resno))
          continue; }
        
        strcpy(A.resnm,res->resnm);
        A.resno=res->resno;
        strcpy(A.id,"OC2"); /* my name of C-ter oxygen is OC2 */
        res->ter=1; /* terminate chain ! */
        if (!includeA(res)) ERROR(("pdb internal error")) }

      else {

        /* exception: OXT (C-term oxygen) is translated to OC2 */
        if (!strcmp(A.id,"OXT")) strcpy(A.id,"OC2");

        if (includeA(res)) goto Aincluded;
        looplist (r,reshead) if (includeA(r)) goto Aincluded;
        res=appendres(res);
      Aincluded:; } }    

    else if (!strcmp(keyword,"HETATM")) { /* //////////////////////// HETATM */
  
      readATOM(1);
      if (omit()) continue;
  
      if (includeA(res)) goto Aincluded;
      looplist (r,reshead) if (includeA(r)) goto HETAincluded;
      res=appendres(res);
      HETAincluded:; }
  
    else if (!strcmp(keyword,"SSBOND")) { /* //////////////////////// SSBOND */
      int i,ssnr;
      ssbond_t *ss;
      char CYSname[4];
    
      alloc(ss,sizeof(ssbond_t));
      /* format:
         SSBOND   1 CYS      3    CYS     40
         SSBOND   2 CYS A    6    CYS A   11
         IIII SSS A IIII    SSS A IIII
      */
  
      forscanf(pdb,"%4i ",&ssnr);
      loop (i,0,2) {
        forscanf(pdb,"%3s %1a %4i    ",CYSname,&ss->chain[i],&ss->resno[i]);
        if (strcmp(CYSname,"CYS"))
          ERROR(("SSBOND connecting %s (should be CYS)",CYSname)) }
      includeSS(ss);
  
      forscanf(pdb,"\n"); /* clear to eol */
      pdbline++; }
  
    else if (!strcmp(keyword,"CONECT")) { /* //////////////////////// CONECT */
      int j,indx[5];
      connect_t *cc;
  
      loop (j,0,5) {
        indx[j]=0;
        forscanf(pdb,"%5d",&indx[j]); }
  
      forscanf(pdb,"\n");
      pdbline++;
  
      if (!indx[0] || !indx[1]) {
        ERROR(("bad CONECT format"))
        continue; }
  
      loop (j,1,5) if (indx[j]) {
        alloc(cc,sizeof(connect_t));
        cc->next=connecthead; connecthead=cc;
        cc->atomno[0]=indx[0]; cc->atomno[1]=indx[j];
        cc->site[0]=cc->site[1]=NULL; /* determined later */ } }
  
    else if (!strcmp(keyword,"TER")) { /* ////////////////////////////// TER */
      if (!res) ERROR(("misplaced TER: no residue to terminate"))
      else res->ter=1;
      forscanf(pdb,"\n");
      pdbline++; }
  
    else if (!strcmp(keyword,"END")) /* //////////////////////////////// END */
      /* the same as EOF */
      break;
  
    else if (!strcmp(keyword,"ENDMDL")) /* ////////////////////////// ENDMDL */
      /* the same as EOF */
      break;
  
    else if (!strcmp(keyword,"HEADER")) { /* //////////////////////// HEADER */
      forscanf(pdb,"%80s\n",pdbheader);
      pdbline++; }
  
    else if (!strcmp(keyword,"COMPND")) { /* //////////////////////// COMPND */
      forscanf(pdb,"%80s\n",pdbcompnd);
      pdbline++; }
    else {
      /* unknown keyword: skip to eol */
      forscanf(pdb,"\n");
      pdbline++; }
  
    } /* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> */
    
  if (res) res->ter=1;
  
  if (!option('p')) fclose(pdb);
} /* readPDB */
  
void pastePDB(char *pdbname,int iframe) /************************** pastePDB */
{
  FILE *out;
  char line[LINELEN],r[LINELEN];
  residue_t *res;
  site_t *s;

  if (iframe<0) {
    out=openfile("pdb","wt",0);
    fprintf(out,
	    "REMARK %s.pdb = recalculated 3D configuration pasted to %s.pdb\n",
	    molname,pdbname); }

  else {
    char ext[16];
#ifdef FAT
    sprintf(ext,"%03d",iframe);
#else
    sprintf(ext,"%04d.pdb",iframe);
#endif
    out=openfile(ext,"wt",0); }
    
  rewind(pdb);
  pdbline=0;

  while (fgets(line,LINELEN,pdb)) {
    looplist (res,reshead)
      looplist (s,res->site)
        if (s->line==pdbline) {
	  sprintf(r,"%8.3f%8.3f%8.3f",s->r[0],s->r[1],s->r[2]);
	  if (strlen(r)>24) ERROR(("line %d: too large coord",pdbline))
          if (strlen(r)<24) ERROR(("line %d: bad format",pdbline))
          memcpy(line+30,r,24); }
    fputs(line,out);
    pdbline++; }

  fclose(out);
}
