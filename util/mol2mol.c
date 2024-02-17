/* cc -O2 -Wall -o mol2mol mol2mol.c
*/

#include "../gen/include.h"

typedef float vector[3];

#define MAXVAL 8
#define MAXMOL 100

struct mol_t {
  vector r;
  int nnbr;
  int nbr[MAXVAL];
  char *charge;
  float R;
  int n;
  char *id;
  char *type;
  int chir;
  char color;
  char incl;
} *mol[MAXMOL];
int ns[MAXMOL];

char *parameter_set="UNKNOWN";
char molhdr[10240]="mol2mol";

char line[256];
char *getfield(char *key)
{
  char *c=strchr(line,'\n');

  if (c) *c=0;

  if ( (c=strstr(line,key)) ) {
    c=strchr(c+strlen(key),'=');
    if (!c) return NULL;
    c++;
    while (*c && *c<=' ') c++;
    return c; }

  return NULL;
}

#define MOL 1
#define GOL 2
#define PLB 4

int main(int narg,char **arg)
{
  int iarg,mode=31,n,frame=-1;
  static float header[2];
  static vector maxL;
  vector L;
  char *outfn;
  int nm=0,im=0,i,ii,rns,inbr,nnbr,nbroken=0;
  FILE *in,*out;
  char *lastname="UNDEFINED";

  if (narg<3) {
    fprintf(stderr,"\
Extract selected sites from mol,gol, and plb files\n\
  mol2mol [PROCESS-OPTIONS] \\\n\
    [-fFRAME] INPUT_NAME [OPTIONS]  [[-fFRAME] [INPUT_NAME [OPTIONS] ..]\n\
    OUTPUT_NAME\n\
PROCESS-OPTIONS:\n\
  -m process only mol files        [if none of -m -g -p is specified,\n\
  -g process only gol files        all files are processed]\n\
  -p process only plb files\n\
INPUT_NAME:\n\
  refers to INPUT_NAME.mol, INPUT_NAME.gol, INPUT_NAME.plb\n\
OUTPUT_NAME:\n\
  refers to OUTPUT_NAME.mol, OUTPUT_NAME.gol, OUTPUT_NAME.plb\n\
FRAME selection (applies to all following INPUT_NAME  unless changed):\n\
  -fFRAME     the frame of plb-file, -f1=first, -f-1=last (default)\n\
OPTIONS (apply to the previous INPUT_NAME, are executed from left):\n\
  +a          include all atoms\n\
  +NUMBER     include atom (site) of given number\n\
  +FROM,TO    include atoms in given range [FROM,TO); also +FROM:TO\n\
  +FROM,TO,BY include atoms FROM, FROM+TO, ..; TO is the first not included\n\
  +iSTRING    include atoms with STRING in their ID (as substring)\n\
  +tTYPE      include atoms with given TYPE\n\
  +@FILE      include atoms listed in FILE, one number per line\n\
  +@          include atoms listed in INPUT_NAME.mark (from blend, hotkey=@m)\n\
  -a etc.     - instead of + means exclude\n\
  +xDX        shift all x coordinates by DX and adjust box size\n\
  +yDY        shift all y coordinates by DY   (DX,DY,DZ should not be negative)\n\
  +zDR        shift all z coordinates by DZ\n\
  -x -y -z    as above, do not adjust box size (max of overlapping boxes)\n\
Sites are numbered from 0 to NS-1\n\
Sites are in the original order irrespective of ITEM ordering (use plbmerge)\n\
Example:\n\
  * from cfg1.mol and cfg1.plb, select all ALA and GLY\n\
  * merge with cfg2.mol,cfg2.plb, shifted by 10 A in z, with excluded atoms\n\
    listed in cfg2.mark (obtained by blend -g cfg2, hot key @m)\n\
  mol2mol -m -p  cfg1 +iALA +iGLY  cfg2 +a -@ +z10\n\
See also:\n\
  plb2plb plbinfo plbsites plbconv plbmerge plbcut plbfilt\n\
  plb2asc asc2plb plb2cfg cfg2plb  molren\n");
    exit(0); }

  outfn=arg[narg-1];
  if (strchr("-+",outfn[0])) Error("mol2mol: bad OUTPUT_NAME (starts with - or +)");

  loop (iarg,1,narg-1) {
    strcat(molhdr," ");
    strcat(molhdr,arg[iarg]);

    if (strchr("+-",arg[iarg][0])) {
      int i,from=atoi(arg[iarg]+2),to,by=1;
      float d=atof(arg[iarg]+2);
      char *sep;
      int incl=arg[iarg][0]=='+';

      switch (arg[iarg][1]) {
        case 'm': if (mode==31) mode=0; mode|=MOL; break;
        case 'g': if (mode==31) mode=0; mode|=GOL; break;
        case 'p': if (mode==31) mode=0; mode|=PLB; break;

        case 'f': frame=from; break;

        case 'a': loop (i,0,ns[im]) mol[im][i].incl=incl;
                  break;
        case 'i': loop (i,0,ns[im])
                    if (strstr(mol[im][i].id,arg[iarg]+2))
                      mol[im][i].incl=incl;
                  break;
        case 't': loop (i,0,ns[im])
                    if (!strcmp(mol[im][i].type,arg[iarg]+2))
                      mol[im][i].incl=incl;
                  break;
        case '@': {
          char *fn=arg[iarg][2]?arg[iarg]+2:string("%s.mark",lastname);
          FILE *in=fopen(fn,"rt");
          if (in) {
            char *c;

            fprintf(stderr,"%s %s\n",incl?"including":"excluding",fn);
            while (fgets(line,256,in)) if (line[0]!='!') {
              c=line; while (*c && *c<=' ') c++;
              i=atoi(c);
              if (i>=0 && i<ns[im]) mol[im][i].incl=incl; }
            fclose(in); }
          else
            Error("mol2mol: no mark file (option -@)"); }
        case 'x': loop (i,0,ns[im]) mol[im][i].r[0]+=d;
                  if (incl) maxL[0]+=d;
                  break;
        case 'y': loop (i,0,ns[im]) mol[im][i].r[1]+=d;
                  if (incl) maxL[1]+=d;
                  break;
        case 'z': loop (i,0,ns[im]) mol[im][i].r[2]+=d;
                  if (incl) maxL[2]+=d;
                  break;
        default: if (isdigit(arg[iarg][1])) {
                    sep=strpbrk(arg[iarg],",:");
                    if (sep) {
                      to=atoi(sep+1);
                      sep=strpbrk(sep+1,",:");
                      if (sep) by=atoi(sep+1); }
                    else
                      to=from+1;
                    for (i=from; i<to; i+=by)
                      if (i>=0 && i<ns[im]) mol[im][i].incl=incl; }
                  else
                    Error("mol2mol: bad option"); } }

    else { /* arg[iarg] not beginning by + - */
      char *c,*fn;

      im=nm;
      if (im>=MAXMOL) Error("mol2mol: too many molecules (max MAXMOL)");
      ns[im]=-1;
      lastname=arg[iarg];
      fprintf(stderr,"=== processing %s ===\n",lastname);

      if (mode&MOL) {
        if ( (in=fopen(fn=string("%s.mol",arg[iarg]),"rt")) ) {
          int *ren=NULL,err=0;

          fprintf(stderr,"reading %s.mol\n",arg[iarg]);
          while (fgets(line,256,in)) if (line[0]!='!') {
              if ( (c=getfield("parameter_set")) ) parameter_set=strdup(c);
              if ( (c=getfield("number_of_atoms")) ) ns[im]=atoi(c);
              if (!memcmp(line,"atoms",5)) {
                if (ns[im]<=0) Error("mol2mol: undefined or zero number of sites");
                if (!mol[im]) allocarrayzero(mol[im],ns[im]);
                allocarray(ren,ns[im]);

                loop (i,0,ns[im]) {
                  char id[128],type[32],charge[32];

                  do if (!fgets(line,256,in)) Error("mol2mol: unexpected EOF");
                  while (line[0]=='!');

                  sscanf(line,"%d%s%s%s%d%d%d%d%d%d%d%d%d%d",
                         &ii,
                         id,
                         type,
                         charge,
                         &mol[im][i].chir,
                         &mol[im][i].nnbr,
                         mol[im][i].nbr,
                         mol[im][i].nbr+1,mol[im][i].nbr+2,mol[im][i].nbr+3,mol[im][i].nbr+4,mol[im][i].nbr+5,mol[im][i].nbr+6,mol[im][i].nbr+7);
                  mol[im][i].id=strdup(id);
                  mol[im][i].type=strdup(type);
                  mol[im][i].charge=strdup(charge);
                  if (ii!=i) err++;
                  ren[i]=ii; } } }
          if (err) {
            int j,k;

            fprintf(stderr,"%d numbering errors detected, will try to renumber\n",err);
            loop (i,0,ns[im])
              loop (j,0,i) if (ren[i]==ren[j]) Error("mol2mol: numbering not unique, bailing out");
            loop (i,0,ns[im]) {
              loop (j,0,mol[im][i].nnbr) {
                loop (k,0,ns[im]) if (ren[k]==mol[im][i].nbr[j]) {
                  mol[im][i].nbr[j]=k; goto found; }
                fprintf(stderr,"%d(%d) neighbor %d not found - omitting\n",
                        i,ren[i],mol[im][i].nbr[j]);
                mol[im][i].nnbr--;
                for (k=mol[im][i].nnbr; k>0; k--) mol[im][i].nbr[k-1]=mol[im][i].nbr[k];
              found:; } }
            fprintf(stderr,"successfully renumbered\n"); }
          free(ren);

          fclose(in); }
        else {
          fprintf(stderr,"cannot open %s (=> no output)\n",fn);
          mode &= 31-MOL; } }

      if (mode&GOL) {
        if ( (in=fopen(fn=string("%s.gol",arg[iarg]),"rt")) ) {
          fprintf(stderr,"reading %s.gol\n",arg[iarg]);
          fgets(line,256,in);
          fgets(line,256,in);
          fgets(line,256,in);
          rns=atoi(line);
          if (ns[im]<0) ns[im]=rns;
          if (rns<ns[im]) Error("mol2mol: not enough sites in the gol file");
          if (ns[im]<=0) Error("mol2mol: undefined or zero number of sites");
          if (!mol[im]) allocarrayzero(mol[im],ns[im]);
          loop (i,0,ns[im]) {
            char longcolor[128];
            if (!fgets(line,256,in)) Error("mol2mol: unexpected EOF");
            sscanf(line,"%s%f",longcolor,&mol[im][i].R);
            mol[im][i].color=longcolor[0]; }
          fclose(in); }
        else {
          fprintf(stderr,"cannot open %s (=> no output)\n",fn);
          mode &= 31-GOL; } }

      if (mode&PLB) {
        if ((in=fopen(fn=string("%s.plb",arg[iarg]),"rb")) ) {
          if (frame==0) Error("mol2mol: -f0 is invalid (use -f1 for the first frame)");

          if (frame<0) frame=0x7fffffff;
          fprintf(stderr,"reading %s.plb, frame=%s\n",
                  arg[iarg],
                  frame==0x7fffffff?"last":string("%d",frame));

          fread(header,sizeof(header),1,in);
          rns=header[0];
          if (ns[im]<0) ns[im]=rns;
          if (rns<ns[im]) Error("mol2mol: not enough sites in the plb file");
          if (ns[im]<=0) Error("mol2mol: undefined or zero number of sites");
          if (!mol[im]) allocarrayzero(mol[im],ns[im]);

          loop (i,0,frame) {
            if (header[1]<0) {
              if (1!=fread(L,sizeof(L),1,in)) goto ex; }
            else L[0]=L[1]=L[2]=header[1];

            loop (i,0,ns[im])
              if (1!=fread(mol[im][i].r,sizeof(vector),1,in)) goto ex; }
        ex:
          Max(maxL[0],L[0])
          Max(maxL[1],L[1])
          Max(maxL[2],L[2])

            fclose(in); }
        else {
          fprintf(stderr,"cannot open %s (=> no output)\n",fn);
          mode &= 31-PLB; } }

      nm++; } }

  n=0;
  loop (im,0,nm)
    loop (i,0,ns[im])
      if (mol[im][i].incl) mol[im][i].n=n++;

  if (!(mode&7)) Error("mol2mol: nothing to write because of missing files");
  fprintf(stderr,"=== writing output files %s.* (ns=%d) ===\n",outfn,n);
  if (!n) Error("mol2mol: no atom selected - nothing to do (consider +a)");

  if (mode&MOL && (out=fopen(string("%s.mol",outfn),"wt")) ) {
    fprintf(out,"%s\n\
\n\
parameter_set = %s\n\
number_of_atoms = %d\n\
\n\
atoms\n\
! i   atom-id  a-type charge chir nbonds bound_atoms\n\
",molhdr,parameter_set,n);
    ii=0;
    loop (im,0,nm)
      loop (i,0,ns[im])
      if (mol[im][i].incl) {
        fprintf(out,"%3d %10s %4s %9s %d",
                ii, mol[im][i].id, mol[im][i].type,mol[im][i].charge,
                mol[im][i].chir);

        nnbr=0;
        loop (inbr,0,mol[im][i].nnbr)
          nnbr+=mol[im][mol[im][i].nbr[inbr]].incl;
        fprintf(out," %d",nnbr);
        nbroken+=mol[im][i].nnbr-nnbr;

        nnbr=0;
        loop (inbr,0,mol[im][i].nnbr)
          if (mol[im][mol[im][i].nbr[inbr]].incl)
            fprintf(out," %d",mol[im][mol[im][i].nbr[inbr]].n);

        fprintf(out,"\n");
        ii++; }

    fclose(out);
    if (ii!=n) Error("mol2mol: internal: ii!=n"); }

  if (mode&GOL && (out=fopen(string("%s.gol",outfn),"wt")) ) {
    fprintf(out,"!sphere#1\n! %s\n%d\n",molhdr,n);
    loop (im,0,nm)
      loop (i,0,ns[im])
        if (mol[im][i].incl) fprintf(out,"%c %.3f\n",mol[im][i].color,mol[im][i].R);
    fclose(out); }

  if (mode&PLB && (out=fopen(string("%s.plb",outfn),"wb")) ) {
    header[0]=n;
    header[1]=-3;
    fwrite(header,sizeof(header),1,out);
    fwrite(maxL,sizeof(vector),1,out);
    loop (im,0,nm)
      loop (i,0,ns[im]) if (mol[im][i].incl)
        fwrite(mol[im][i].r,sizeof(vector),1,out);
    fclose(out); }

  if (nbroken) fprintf(stderr,"WARNING: %d broken bonds removed\n",nbroken/2);

  return 0;
}
