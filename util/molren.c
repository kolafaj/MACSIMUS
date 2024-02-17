/* cc -O2 -Wall -o molren molren.c
*/

#include "../gen/include.h"

typedef float vector[3];

/* check the code if MAXVAL changes! */
#define MAXVAL 8

struct mol_t {
  vector r;
  int nnbr;
  int nbr[MAXVAL];
  char *charge;
  float R;
  char *id;
  char *type;
  int chir;
  char color;
} *mol;
int ns;

char *parameter_set="UNKNOWN";
char molhdr[10240]="molren";

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

int main(int narg,char **arg)
{
  int frame,i,ii,inbr;
  static float header[2];
  vector L;
  char *infn,*outfn,*c;
  FILE *in,*out;
  int *ren,*ner,*swp;

  if (narg<4) {
    fprintf(stderr,"\
Renumber mol,plb,gol files. Call by:\n\
  molren {-d|-i} INPUT_NAME OUTPUT_NAME n0 n1 ...\n\
  refers to INPUT_NAME.mol, INPUT_NAME.gol, INPUT_NAME.plb\n\
OUTPUT_NAME:\n\
  refers to OUTPUT_NAME.mol, OUTPUT_NAME.gol, OUTPUT_NAME.plb\n\
OPTIONS:\n\
  -d	direct renumbering: site 0 of INPUT_NAME becomes n0 in OUTPUT_NAME\n\
  -i	inverted renumbering: site 0 of OUTPUT_NAME is site n0 of INPUT_NAME\n\
See also:\n\
  mol2mol wat2wat\n");
    exit(0); }

  infn=arg[2];
  outfn=arg[3];

  /* renumbering tables */
  allocarray(ren,narg-4);
  allocarray(ner,narg-4);

  ns=narg-4;
  fprintf(stderr,"molren: ns=%d in the renumbering table\n",ns);
  loop (i,0,ns) {
    int is=atoi(arg[i+4]),j;

    if (is<0 || is>=ns) Error("molren: atom number out of range");
    loop (j,0,i) if (is==ren[j]) Error("molren: atom number repeats");
    ren[i]=is;
    ner[i]=-1; }

  loop (i,0,ns) ner[ren[i]]=i;
  loop (i,0,ns) if (ner[i]<0) Error("molren: numbering");
  
  if (!strcmp(arg[1],"-d")) swp=ren,ren=ner,ner=swp;
  else if (!strcmp(arg[1],"-i")) ;
  else Error("molren: bad option");

  loop (i,0,ns) printf("%d %d\n",ren[i],ner[i]);

  /* MOL */
  
  if ( (in=fopen(string("%s.mol",infn),"rt")) ) {
    fprintf(stderr,"molren: reading %s.mol .. ",infn);
    while (fgets(line,256,in)) if (line[0]!='!') {
        if ( (c=getfield("parameter_set")) ) parameter_set=strdup(c);
        if ( (c=getfield("number_of_atoms")) ) ii=atoi(c);
        if (!memcmp(line,"atoms",5)) {
          if (ns!=ii) Error("molren: ns in mol-file does not match the renumbering table");
          if (ns<=0) Error("molren: undefined or zero number of sites");
          if (!mol) allocarrayzero(mol,ns);
          
          loop (i,0,ns) {
            char id[128],type[32],charge[32];

            do if (!fgets(line,256,in)) Error("molren: mol unexpected EOF");
            while (line[0]=='!');
            
            sscanf(line,"%d%s%s%s%d%d%d%d%d%d%d%d%d%d",
                   &ii,
                   id,
                   type,
                   charge,
                   &mol[i].chir,
                   &mol[i].nnbr,
                   mol[i].nbr,
                   mol[i].nbr+1,mol[i].nbr+2,mol[i].nbr+3,mol[i].nbr+4,mol[i].nbr+5,mol[i].nbr+6,mol[i].nbr+7);
            mol[i].id=strdup(id);
            mol[i].type=strdup(type);
            mol[i].charge=strdup(charge);
            if (ii!=i) Error("molren: incorrect numbering of mol-file. Try mol2mol to fix first."); }
        } }

    fclose(in);
    
    out=fopen(string("%s.mol",outfn),"wt");
    if (!out) Error("molren: cannot write mol-file");
    fprintf(out,"%s\n\
\n\
parameter_set = %s\n\
number_of_atoms = %d\n\
\n\
atoms\n\
! i   atom-id  a-type charge chir nbonds bound_atoms\n\
",molhdr,parameter_set,ns);
    
    loop (i,0,ns) {
      fprintf(out,"%3d %10s %4s %9s %d",
              i,
              mol[ren[i]].id,
              mol[ren[i]].type,
              mol[ren[i]].charge,
              mol[ren[i]].chir);

      fprintf(out," %d",mol[ren[i]].nnbr);

      loop (inbr,0,mol[ren[i]].nnbr)
        fprintf(out," %d",ner[mol[ren[i]].nbr[inbr]]);

      fprintf(out,"\n"); }
    fprintf(stderr,"%s.mol written\n",outfn);
    fclose(out);
  }

  /* GOL */
  
  if ( (in=fopen(string("%s.gol",infn),"rt")) ) {
    fprintf(stderr,"reading %s.gol\n",infn);
    loop (i,0,3)
      if (!fgets(line,256,in)) Error("molren: gol-file too short");
    ii=atoi(line);
    if (ns!=ii) Error("molren: ns in gol-file does not match the renumbering table");
    if (!mol) allocarrayzero(mol,ns);
    loop (i,0,ns) {
      char longcolor[128];
      
      if (!fgets(line,256,in)) Error("mol2mol: unexpected EOF");
      sscanf(line,"%s%f",longcolor,&mol[i].R);
      mol[i].color=longcolor[0]; }
    fclose(in);
    
    out=fopen(string("%s.gol",outfn),"wt");
    if (!out) Error("molren: cannot write gol-file");
    fprintf(out,"!sphere#1\n! %s\n%d\n",molhdr,ns);
    loop (i,0,ns) fprintf(out,"%c %.3f\n",mol[ren[i]].color,mol[ren[i]].R);
    fclose(out);
  }

  /* PLB */
  
  if ( (in=fopen(string("%s.plb",infn),"rb")) ) {
    fprintf(stderr,"reading %s.plb ... ",infn);

    if (1!=fread(header,sizeof(header),1,in))
      Error("molren: plb-file too short");
    ii=header[0];
    if (ns!=ii) Error("molren: ns in plb-file does not match the renumbering table");
    if (!mol) allocarrayzero(mol,ns);
    
    out=fopen(string("%s.plb",outfn),"wb");
    if (!out) Error("molren: cannot write plb-file");
    fwrite(header,sizeof(header),1,out);

    for (frame=0;;frame++) {
      if (header[1]<0) {
        if (1!=fread(L,sizeof(L),1,in)) goto ex; }
      else L[0]=L[1]=L[2]=header[1];
      
      loop (i,0,ns)
          if (1!=fread(mol[i].r,sizeof(vector),1,in)) goto ex; 
      
      if (header[1]<0) fwrite(L,sizeof(vector),1,out);
      loop (i,0,ns)
        fwrite(mol[ren[i]].r,sizeof(vector),1,out); }
    
  ex:
    fprintf(stderr,"%s.plb %d frames\n",outfn,frame);    
    fclose(out); }
  
  return 0;
}
