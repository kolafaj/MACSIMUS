/*
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% residue statistics %%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include "ground.h"
#include "options.h"
#include "pdbinfo.h"

omitrsd_t *omitrsd0;

static struct stat_s {
  struct stat_s *next;
  char resnm[4];
  int n;
} *stathead=NULL;

void count(char *resnm) /******************************************** count */
{
  struct stat_s *news,*s,*lasts=NULL;
  int cmp;

  looplist (s,stathead) {
    cmp=strcmp(resnm,s->resnm);
    if (cmp==0) { s->n++; return; }
    else if (cmp<0) break;
    lasts=s; }

  alloc(news,sizeof(struct stat_s));
  if (lasts) lasts->next=news;
  else stathead=news;
  news->next=s;
  memcpy(news->resnm,resnm,4); news->resnm[3]=0;
  news->n=1;
}

void prtcount(void) /********************************************* prtcount */
{
  struct stat_s *s;
  int i=0,n=0;

  prts_("residue statistics:");
  looplist (s,stathead) {
    if (i++%5==0) _n
    n+=s->n;
    prt_("  %3s %-4i",s->resnm,s->n); }
  prt("\ntotal %d residues",n);
}

void prtomitted(void) /***************************************** prtomitted */
{
  omitrsd_t *omitrsd;

  if (omitrsd0) {
    prts("omitted residues:");
    looplist (omitrsd,omitrsd0)
      prt("  %s %d (total %d lines)",omitrsd->resnm,omitrsd->n,omitrsd->l); }
}

void help(void) /***************************************************** help */
{
  prts_("\
  pdb [OPTIONS] PDBNAME [MOLNAME]\n\
FILE ARGUMENTS:\n\
  PDBNAME.pdb  source PDB-file (default MOLNAME=PDBNAME - cannot if -p)\n\
  MOLNAME.pdb  output PDB-file (if option -p)\n\
  MOLNAME.rep  optional pattern replacement file\n\
               (example of .rep line=\"C OT C  C OE C\")\n\
  PDBNAME.rep  try this if MOLNAME.rep not found\n\
  MOLNAME.sel  optional residue selection file\n\
               (example of .sel line=\"ASP 11 A asph\")\n\
  MOLNAME.mol  output mol-file\n\
\n\
  output (if not -p) or input (if -p) configuration files, see also -b,-t\n\
  MOLNAME.3db  bin format float[][3], deprecated\n\
  MOLNAME.plb  playback format (header float[2]+bin)\n\
  MOLNAME.3dt  text 3 column format (NB: elsewhere denoted .xyz)\n\
  MOLNAME.plt  text playback format (NB: elsewhere denoted .pla)\n\
FORCE FIELD OPTIONS:\n\
  -rRESDIR  directory of residue files *.rsd, relative to BLENDPATH\n\
  -fPARSET  parameter_set = PARSET.par (in BLENDPATH, used by blend)\n\
            if not given, RESDIR.par is assumed\n\
OPTIONS\n\
  -a#   alternate locations: default 0=loc.A  1=average  2=max occupancy\n\
       4=no alt.l. and wider id (col.17) for HETATM, 8=for ATOM; #'s can sum\n\
  -b#   binary config.: #=1 MOLNAME.3db, #=2 MOLNAME.plb (default), #=3 both\n\
        #=0 none  #<0 reverse endian\n\
  -c#   make CYS-CYS bond if |SG-SG|<# AA\n\
  -d#   don\'t add atoms missing in PDBNAME.pdb\n\
  -i#   ignore additional atoms in PDBNAME.pdb\n\
        -d,-i apply for:\n\
         #<0:all residues  #>0:molecules only only\n\
         |#|=1:H  |#|=2:heavy atoms  |#|=3:all atoms\n\
  -hRSD head: use RSD.rsd as N-term patch (default=according to -n)\n\
  -eRSD end: use RSD.rsd as C-term patch (default=according to -n)\n\
  -g    gap in residue numbering terminates chain\n\
  -l#   max valence\n\
  -m0   don\'t write MOLNAME.mol file\n\
  -n    use neutral termini+residues (ASP GLU protonated, ARG LYS HIS w/o H+)\n\
  -n-1  add counterions Na+,Cl-\n\
  -n-%  distance of counterion in % of bond\n\
  -oRSD omit residue (usually HOH or WAT)\n\
  -p#[:TO:BY] make MOLNAME.pdb [MOLNAME.####.pdb] by pasting #-th config[s]\n\
              from MOLNAME.plb into PDBNAME.pdb\n\
  -q%   charge multiplication factor\n\
  -s#   scrolling (#kB buffer)\n\
  -t#   text config: #=1 MOLNAME.3dt, #=2 MOLNAME.plt, #=3 both (default=none)\n\
  -u    enable PDB:rsd atom match if rsd name unique even if not the same #\n\
  -v#   0=no info/warnings, -1=no (re)connect report, 1=default, 2=verbose\n\
  -x    enable A in atom name to match any atom for rotating groups\n\
  -z    use patch names for residues: 0=never 2=always 1=only from pdb (df.)\n\
  -\\0   don\'t remember residues - always re-read file\n");
  exit(0);
}
