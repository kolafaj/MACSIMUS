/* cc -O2 -o shownear shownear.c -lm 
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  char line[128],*ch;
  FILE *gol,*plb;
  int from,to,ns,i,j,frame,incl,varL;
  double dist;
  char *colornear,*colormark=NULL;
  float (*r)[3];
  float hdr[2];
  struct gol_s {
    char col[8];
    float R; } *g;
  
  if (narg<7) {
    fprintf(stderr,"Mark neighbors. Call by:\n\
  shownear FILE.plb[:FRAME] FILE.gol FROM TO [-]DIST COLORNEAR [COLORMARK] \\\n\
    > OUTFILE.gol\n\
Marks the atoms closer than |DIST| from atoms [FROM..TO) by COLORNEAR,\n\
optionally the atoms [FROM..TO) by COLORMARK\n\
  -DIST means that atoms [FROM..TO) are not included in the tests\n\
  TO<=0 means to the end\n\
BUG: only free b.c., irrespective of L\n");
    exit(0); }

  strcpy(line,arg[1]);
  if ( (ch=strchr(line,':')) ) *ch=0,frame=atoi(ch+1);
  else frame=1;
  if (frame<1) Error("frame<1");
  
  if (!(plb=fopen(line,"rb"))) Error(line);
  if (!(gol=fopen(arg[2],"rt"))) Error(arg[2]);
  from=atoi(arg[3]);
  to=atoi(arg[4]);
  dist=atof(arg[5]); 
  incl=dist>0;
  dist*=dist;
  colornear=arg[6];
  if (narg>7) colormark=arg[7];
  
  fread(hdr,4,2,plb);
  ns=hdr[0];
  if (to<=0) to=ns;
  if (from>=to) Error("from>=to");
  varL=hdr[1]<0;
  
  fgets(line,128,gol); fputs(line,stdout);
  fgets(line,128,gol); fputs(line,stdout);
  fgets(line,128,gol); fputs(line,stdout);
  if (ns!=atoi(line)) Error("ns[plb] != ns[gol]");
  
  alloc(r,sizeof(r[0])*(ns+varL));
  loop (i,0,frame) if (fread(r,sizeof(r[0]),ns+varL,plb)!=ns+varL) Error("read plb");
  
  alloc(g,sizeof(g[0])*ns);
  loop (i,0,ns) {
    fgets(line,128,gol);
    sscanf(line,"%s%f",g[i].col,&g[i].R); }
  
  fclose(gol);
  fclose(plb);

  if (colormark) loop (i,from,to) strcpy(g[i].col,colormark);

  r+=varL;

  loop (i,0,ns) if (incl || i<from || i>=to)
    loop (j,from,to) if (i!=j) {
      double rr=Sqr(r[i][0]-r[j][0])+Sqr(r[i][1]-r[j][1])+Sqr(r[i][2]-r[j][2]);
      if (rr<dist) {
        strcpy(g[i].col,colornear); break; } }
  
  loop (i,0,ns) printf("%8s %5.3f\n",g[i].col,g[i].R);

  return 0;
}
