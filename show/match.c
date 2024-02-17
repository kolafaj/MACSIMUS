/* make show
 */
fvector centeR;

char Matchfile[]="match.plb:1", *matchfile=Matchfile;

double cfgdist(fvector *r,fvector *R,int ns,int *fix) /************* cfgdist */
/* 
  returns the (averaged) distance of configurations r and R of ns sites
  if fix!=NULL, then only of selected sites
*/
{
  float sum=0;
  fvector cr,dr;
  int i,nfix=0;

  loop (i,0,ns) if (!fix || fix[i]) {
    cr[0]=r[i][0]-center[0]+centeR[0];
    cr[1]=r[i][1]-center[1]+centeR[1];
    cr[2]=r[i][2]-center[2]+centeR[2];
    dr[0]=rot[0][0]*cr[0]+rot[0][1]*cr[1]+rot[0][2]*cr[2]-R[i][0];
    dr[1]=rot[1][0]*cr[0]+rot[1][1]*cr[1]+rot[1][2]*cr[2]-R[i][1];
    dr[2]=rot[2][0]*cr[0]+rot[2][1]*cr[1]+rot[2][2]*cr[2]-R[i][2];
    sum+=Sqr(dr[0])+Sqr(dr[1])+Sqr(dr[2]); 
    nfix++; }

  return sqrt(sum/nfix);
}

/* cheap random numbers */
double rnd31(void) /************************************************** rnd31 */
{
  static double seed=1111;

  seed=fmod(seed*16807e0,2147483647e0);

  return seed;
}

double gauss(void) /************************************************** gauss */
{
  return (rnd31()-rnd31()+rnd31()-rnd31())/2147483647e0;
}

void match(fvector *r,int ns,double eps) /**************************** match */
/* 
  matches plb-file matchfile (of given frame) to *r with precision eps
  example of matchfile: "to-match.plb:2:fixed.fix" (see -m of show)
  if mrev, the endian is reversed 
*/
{
  FILE *plb;
  int frame=1,iframe,i;
  static fvector *R;
  static int *fix;
  char *fixname=NULL;
  float hdr[2];
  static double dr,da;
  double drot0[3][3];
  float center0[3];
  double dist0,dist;
  int it=30; /* minimum number of iterations */
  int matchvarL=0;

  if (eps<1e-5) it=100;

  if (!R) {
    char *c=strchr(matchfile,':');

    if (!c) 
      c=strend(matchfile);
    else {
      *c++=0;
      frame=atoi(c);
      c=strchr(c,':');
      if (c) fixname=c+1; }      
      
    if (frame<1) frame=1;

    fprintf(stderr,"reading match file %s, frame %d\n",matchfile,frame);

    if (!(plb=fopen(matchfile,"rb"))) {
      fprintf(stderr,"no such file\n");
      return; }

    if (fread(hdr,4,2,plb)!=2) {
      fprintf(stderr,"too short file\n");                                       
      return; } 
      
    if (hdr[1]==-3) matchvarL=1;

    if (hdr[0]>16777216. || ns!=(int)hdr[0])
      ERROR(("current ns=%d, match file %s: ns=%g\n",ns,matchfile,hdr[0]))

    alloc(R,(ns+matchvarL)*sizeof(fvector));

    loop (iframe,0,frame) {
    
      if (matchvarL) 
        if (fread(R,sizeof(fvector),1,plb)!=1) ERROR(("%s too short",matchfile))
      if (fread(R,sizeof(fvector),ns,plb)!=ns) ERROR(("%s too short",matchfile))

      if (iframe==frame-1) {
	int j;
	fvector maxr,minr;

	loop (j,0,3) maxr[j]=-(minr[j]=9e9);
	
	loop (i,0,ns) if (R[i][0]<99999.0)
	  loop (j,0,3) {
            Min(minr[j],R[i][j])
	    Max(maxr[j],R[i][j]) }

	loop (i,0,3) centeR[i]=(maxr[i]+minr[i])/2; } }
    fclose(plb); 
    
    if (fixname) {
      char line[256];
      int i;
      
      fprintf(stderr,"reading list of sites to match %s\n",fixname);

      if (!(plb=fopen(fixname,"rt"))) {
        fprintf(stderr,"\ano such file\n");
        return; }
      
      allocarrayzero(fix,ns);
      while (fgets(line,256,plb)) if (line[0]!='!') 
        if (sscanf(line,"%d",&i)) if (i>=0 && i<ns) fix[i]++; }
    } /* !R */

  if (!OPTION_C)
    loop (i,0,3) center[i]-=centeR[i],centeR[i]=0;   

  copy(drot0,drot,sizeof(drot));
  copy(center0,center,sizeof(center));
  dist=dist0=cfgdist(r,R,ns,fix);
  printf("av.dist = %9.7f --> ",dist);

  dr=0.1; da=0.01;

  do { 
    copy(drot,drot0,sizeof(drot));
    loop (i,0,3) rotate(i,da*gauss());
    dist=cfgdist(r,R,ns,fix);
    if (dist<dist0) {
      dist0=dist;
      copy(drot0,drot,sizeof(drot));
      da*=1.1; }
    else {
      copy(drot,drot0,sizeof(drot));
      rotate(0,0); /* rot:=drot */
      da*=0.98; }

    loop (i,0,3) center[i]=center0[i]+dr*gauss();
    dist=cfgdist(r,R,ns,fix);
    if (dist<dist0) {
      dist0=dist;
      copy(center0,center,sizeof(center));
      dr*=1.1; }
    else {
      copy(center,center0,sizeof(center));
      dr*=0.98; }
    } while (fabs(dr)+fabs(da*10)>eps || --it>0);

  copy(drot,drot0,sizeof(drot));
  copy(center,center0,sizeof(center));
  printf("%9.7f da=%8.2e dr=%8.2e\n",dist0,da,dr);
}
