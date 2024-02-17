#include "ground.h"
#include "options.h"

#include "pdbbasic.h"
#include "pdbrep.h"
#include "pdbconv.h"

static float charge(char *atom) /************************************ charge */
{
  char *c=strpbrk(atom,"+-");
  float f;

  if (c) {
    f=atof(c);
    *c=0; }
  else
    f=10e9;

  return f;
}

static site_t *s[3];
static char *from[3],*to[3];
static float fromq[3],toq[3];

static int match(int i) /********************************************* match */
{
  return
    !strcmp(s[i]->type,from[i])
    && (fromq[i]>9e9 || fabs(fromq[i]-s[i]->charge)<0.001);
}

int markedbyDEL(site_t *s) /************************************ markedbyDEL */
{
  return !strcmp(s->type,"DEL");
}

void replacepattern(void) /********************************** replacepattern */
/*
  replaces patterns of atom types by 3 according to data in file *.rep
  example of file *.rep:
! this is comment
CH1E OT     CH1E    = CH1E OE CH1E
CH2E OT-0.4 CH1E+.1   CH2E OE-0.3 CH1E+.1

  Numbers after atom types mean optional partial charges, positive must have `+'
  If not given than they are irrelevant for atomtype match or not changed

  ! Only bonded triplets supported, no wildcards !
  ! EXTREMELY inefficient implementation !
*/
{
  residue_t *res;
  char r[82],rr[82],*swap;
  int i,j,k,no;

  FILE *pattern=openfile("rep","rt",0);

  if (!pattern && pdbname && strcmp(molname,pdbname)) {
    if (option('v'))
      prt("no MOLFILE.rep pattern replacement file, trying PDBFILE.rep");
    swap=molname; molname=pdbname;
    pattern=openfile("rep","rt",0);
    molname=swap; }

  if (!pattern) {
    if (option('v')) prt("no pattern replacement file");
    return; }

  while (fgets(r,82,pattern)) if (r[0]!='!') {

    strcpy(rr,r);
#define ERR { ERROR(("bad line in *.rep file:\n%s",rr)) goto nextline; }

    if ( !(from[0]=strtok(r," \t")) ) goto nextline; /* empty line */
    if ( !(from[1]=strtok(NULL," \t")) ) ERR
    if ( !(from[2]=strtok(NULL," \t:=")) ) ERR
    loop (i,0,3) if ( !(to[i]=strtok(NULL," \t\n")) ) ERR

    prt_("replacing pattern %s %s %s -> %s %s %s",
         from[0],from[1],from[2],
         to[0],to[1],to[2]);
    no=0;

    loop (i,0,3) {
      fromq[i]=charge(from[i]);
      toq[i]=charge(to[i]);
      if (strlen(to[i])>=STRLEN) {
        ERROR(("%s too long in *.rep file",to[i])) goto nextline; } }

    looplist (res,reshead)
      looplist (s[1],res->site)
        if (match(1))
          loop (i,0,s[1]->nnbr)  {
            s[0]=findsite(s[1]->nbr[i]);
            if (!s[0]) {
              ERROR(("site number %d not found",s[1]->nbr[i]))
              continue; }
            loop (j,0,s[1]->nnbr) if (i!=j) {
              if (match(0)) {
                s[2]=findsite(s[1]->nbr[j]);
                if (!s[2]) {
                  ERROR(("site number %d not found",s[1]->nbr[j]))
                  continue; }
                if (match(2)) {
                  if (option('v'))
                    prt_("\n %d %s - %d %s - %d %s  replaced",
                         s[0]->n,s[0]->ids->id,
                         s[1]->n,s[1]->ids->id,
                         s[2]->n,s[2]->ids->id);
                  no++;
                  loop (k,0,3) {
                    strcpy(s[k]->type,to[k]);
                    if (toq[k]<9e9) s[k]->charge=toq[k]; } } } } }

    prt(" : total %d replaced",no);

    nextline: ; }

  removesites(markedbyDEL,"marked by DEL in *.rep file");
}

#undef ERR


select_t *select0;

void readSEL(void) /********************************************** readSEL */
/*
  Selecting different charged residues according to info in a file

  Example:

  option -n  ! if this is given then the runtime value of option -n is checked

  !RSD No chain repl
   ARG 11   -   argn ! means that ARG11 will be replaced by argn

  The format is free,
  `-' in the chain field means no chain
  `*' in the chain field means any chain
  -1 as No means any residue

  RSD MUST be upper case, replacement refers to .rsd files and for UNIX
  is case sensitive; lowercase recommended!
*/

{
  char r[82],rr[82];
  char *resnm,*resnotok,*chaintok,*resfn;
  char chain;
  int resno;
  int nsel=0;
  select_t *select;
  FILE *sel=openfile("sel","rt",0);

  if (!sel) {
    if (option('v')) prt("no residue selection file");
    return; }

  while (fgets(r,82,sel)) if (r[0]!='!') {

    strcpy(rr,r);
#define ERR { ERROR(("bad line in *.sel file:\n%s",rr)) goto nextline; }

    if ( !(resnm=strtok(r," \t")) ) goto nextline; /* empty line */

    if ( !(resnotok=strtok(NULL," \t")) ) ERR

    if (!strcmp(resnm,"option")) {
      int n;

      if (memcmp(resnotok,"-n",2)) ERR
      n=atoi(resnotok+2);
      if (resnotok[2]=='-' && resnotok[3]==0) n=-1;
      if (option('n')!=n)
        WARNING(("option -n%d in *.sel differs from actual: -n%d",n,option('n')))
      goto nextline; }

    if ( !(chaintok=strtok(NULL," \t")) ) ERR
    if ( !(resfn=strtok(NULL," \t\n")) ) ERR

    if (!strchr("-+0123456789*",resnotok[0])) ERR
    if (chaintok[1]) ERR
    chain=chaintok[0];
    if (chain=='-' || chain=='~') chain=' ';
    if (resnotok[0]=='*') resno=-1; else resno=atoi(resnotok);

    alloc(select,sizeof(select_t)+strlen(resfn));
    select->next=select0;
    select0=select;
    strcpy(select->resnm,resnm);
    select->chain=chain;
    select->resno=resno;
    strcpy(select->resfn,resfn);
    nsel++;

    nextline:; }

  prt("%d selections in *.sel file",nsel);
  fclose(sel);
}
