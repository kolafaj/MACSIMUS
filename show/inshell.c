double distq=9;
int wmol[4];

// inline
float ndq(int j,float x)
/* squared nearest image distance */
{
  /* return Sqr(Lh-fabs(Lh-fabs(x))); */

  if (x>Lh[j]) do x-=L[j]; while (x>Lh[j]);
  else while (x<-Lh[j]) x+=L[j];

  return x*x;
}

void inshell(void)
{
  int i,j,k;
  fvector r;
  float rr;
  float sh;
  int ninshell=0;

  /* site[].shell: 1=water, 2=water Oxygen, 4=to draw */

  if (distq==0) {
    loop (i,0,ns) site[i].shell|=4;
    return; }

  loop (i,0,ns) {
    if (site[i].shell&1) {
      if (site[i].shell&2) {
	site[i].shell=3;
	loop (k,0,3) r[k]=cfg[i][k];

	if (option('l')) {
	  loop (j,0,ns) if (!(site[j].shell&1)) {
	    if ((rr=ndq(0,r[0]-cfg[j][0]))>distq) continue;
	    if ((rr+=ndq(1,r[1]-cfg[j][1]))>distq) continue;
	    if ((rr+ndq(2,r[2]-cfg[j][2]))<distq) {
	      ninshell++;
	      site[i].shell=7;
	      if (option('j')>0) loop (k,0,3) {
		sh=floor((cfg[j][k]-r[k]+Lh[k])/L[k])*L[k];
		cfg[i][k]+=sh;
		cfg[wmol[0]+i][k]+=sh;
		cfg[wmol[1]+i][k]+=sh; }
	      goto nexti; } } }

	else loop (j,0,ns) if (!(site[j].shell&1)) {
	  if ((rr=Sqr(r[0]-cfg[j][0]))>distq) continue;
	  if ((rr+=Sqr(r[1]-cfg[j][1]))>distq) continue;
	  if ((rr+Sqr(r[2]-cfg[j][2]))<distq) {
	    ninshell++;
	    site[i].shell=7;
	    goto nexti; } } } 
      else 
        site[i].shell=1; }
  nexti:; }

  /* marking also n.n. */
  loop (i,0,nc) loop (j,0,2)
    if (site[bond[i][j]].shell&4) site[bond[i][!j]].shell|=4;

  prt("%d waters in %f-shell",ninshell,sqrt(distq));
}
