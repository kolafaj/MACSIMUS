/* cc -Wall -O2 -o shiftzprof shiftzprof.c -lm
 */
#include "../gen/include.h"

int main(int narg,char **arg)
{
  int iarg;
  double MAX;

  if (narg<3) {
    fprintf(stderr,"\
Shift/extend MACSIMUS z-density profiles periodically. Call by:\n\
  shiftzprof MAX DENSPROF.z [DENSPROF.z..]\n\
where\n\
  0<MAX<1   z will be in range [(MAX-1)*Lz,MAX*Lz); i.e., shifted\n\
  -1<MAX<0  z will be in range [MAX*Lz,Lz); i.e., padded left\n\
The output extensions are .Z\n\
Example:\n\
  shiftzprof -.5 SIM.cm.kgm-3.z; plot SIM.cm.kgm-3.Z\n\
");
    exit(0); }

  MAX=atof(arg[1]);
  if (MAX<-1 || MAX>1) Error("MAX out of range [-1,1]");

  loop (iarg,2,narg) {
    FILE *in=fopen(arg[iarg],"rt"),*out;
    char *fn=strdup(arg[iarg]),*last=strlast(fn);
    char line[1024];
    struct line_s {
      struct line_s *next;
      double z;
      char rest[1];
    } *head=NULL,*tail,*n;
    double Lz=0,dz=0;
    int nl=0;

    if (*last!='z') Error("wrong extension: z expected");
    *last='Z';
    out=fopen(fn,"wt");

    fprintf(stderr,"writing to %s\n",fn);

    fprintf(out,"# shifted periodically by %g Lz\n",MAX);

    while (fgets(line,1024,in)) {
      char *c=line;
      double z;

      if (*c=='#') {
        fputs(line,out);
        if ( (c=strstr(line,"<Lz>=")) ) Lz=atof(c+5);
        if ( (c=strstr(line,"dz=")) ) dz=atof(c+3); }
      else {
        while (*c==' ') c++;
        if (!*c) Error("line too short");
        z=atof(c);
        while (*c>' ') c++;
        if (!*c) Error("line too short");
        while (*c==' ') c++;
        if (!*c) Error("line too short");
        n=malloc(sizeof(*n)+strlen(c));
        strcpy(n->rest,c);
        n->z=z;
        nl++;
        n->next=NULL;
        if (head) {
          tail->next=n;
          tail=n; }
        else
          head=tail=n; } }

    fprintf(stderr,"%d histogram lines read Lz=%g dz=%g\n",nl,Lz,dz);
    if (Lz<=0) Error("keyword <Lz>= missing or wrong");

    //    looplist (n,head) printf("DEBUG %g %s",n->z,n->rest);
    
    if (MAX>=0) {
      looplist (n,head) if (n->z>=MAX*Lz) fprintf(out,"%g %s",n->z-Lz,n->rest);
      looplist (n,head) if (n->z<MAX*Lz) fprintf(out,"%g %s",n->z,n->rest); }
    else {
      looplist (n,head) if (n->z>=(MAX+1)*Lz) fprintf(out,"%g %s",n->z-Lz,n->rest);
      looplist (n,head) fprintf(out,"%g %s",n->z,n->rest); }

    fclose(in);
    fclose(out);
  }

  return 0;
}
