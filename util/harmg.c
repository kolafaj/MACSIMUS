/* make harmg
see harmonic.h for more info
*/


#include "ground.h"
#define SDS
#include "alloc.h"
#include "varfile.h"
#define HARMONICS 3
#include "harmonic.h"

#include "vector.h" /* FLOAT #defined in simopt.h ==> vector = float[3] */

/* angle H-O-H = 107, dist O-H = 0.98 (compromise between ST2 and TIP4P) */
const vector O ={0, 0    ,0    };
const vector H1={0, 0.788,0.583};
const vector H2={0,-0.788,0.583};
#define INITw VV(w[2],=O) VV(w[0],=H1) VV(w[1],=H2)

int NM=100;
int OOFRAME=0;
int merged=0;
int pint=1;

struct harmonics_sds *harmonics;

#ifdef G2DGRID
int g2d[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];

void rgbcolor(char rgb[3],double x)
{
#define QPPM (1.0)
#define NCOL 6
  const double RGB[NCOL+1][3]={

#if NCOL==7
    { 1,1,1 },
    { 1,0,0 },
    { 1,1,0 },
    { 0,1,0 },
    { 0,1,1 },
    { 0,0,1 },
    { 1,0,1 },
    { 0,0,0 } 
#elif NCOL==6
    { 1,1,0 },
    { 0,1,0 },
    { 0,1,1 },
    { 0,0,1 },
    { 1,0,1 },
    { 1,0,0 },
    { 1,1,0 }
#else
#error bad NCOL
#endif
  };

  int i,j;

#if NCOL==7
  if (x<0) x=0;
  if (x>=1) x=.9999999;
#else
  x=fmod(x,1);
#endif
  i=x*NCOL;
  x=x*NCOL-i;
  /* rounding to have limited number of colors for GIF */
  x=((int)(x*32))/32.;
  loop (j,0,3) rgb[j]=(RGB[i][j]*(1-x)+RGB[i+1][j]*x)*255;
}

#endif

#if NHARM>1
int anglehist
  /*     r     phi   theta1 theta2 alpha1 alpha2  */
    [R_NHIST][NHARM][NHARM][NHARM][NHARM][NHARM];

#define cosd(X) cos((X)*(PI/180))
#define sind(X) sin((X)*(PI/180))
void rotd(vector r,int axis,double angle) /******************** rotd */
{
  double sa=sind(angle), ca=cosd(angle);
  double o[3][3];
  double newr[3];
  int i=(axis+1)%3,j=(i+1)%3,k=(j+1)%3,a;

  o[i][i]= ca; o[i][j]=-sa; o[i][k]=0;
  o[j][i]= sa; o[j][j]= ca; o[j][k]=0;
  o[k][i]=  0; o[k][j]=  0; o[k][k]=1;

  loop (a,0,3) newr[a]=SCAL(o[a],r);
  VV(r,=newr)
}
#endif

int main(int narg, char **arg)
{
  double q,qq,g,gc,gc2,gcc,ga,gaa,r,qxyz,R=0,RQ;
  double Kirkwood=0;
  int ir,rel;
  char *fmt="%8.4f";
  char fn0[128],fn[128],*c;
#define IF(C) if (strpbrk(arg[2],C))

  initscroll(0);

  if (narg<3) {
    prt("*** full and 2D corr.f. analysis for water *** (c) J.Kolafa ***\n\
Usage:\n\
  %s SIMNAME[.har] KEYS [OPTIONS]\n\
KEYS:\n\
  h = g(r) and harmonics, SIMNAME.hg\n\
  k = running Kirkwood factor, SIMNAME.rkf\n\
  g = harmonics divided by g(r), SIMNAME.hgr\n\
  a = list of max g(FULL) in shells, SIMNAME.ah\n\
  r = raw playbacks, XYr-####.plb (X={st}, Y={cf})\n\
  s = smart playbacks (grouped histograms), xxs-####.plb\n\
  2 = 2D distribution functions as xyz tables, SIMNAMEsH.xyz (s=O or H)\n\
  p = 2D distribution functions as color pictures, SIMNAME[r]sH.ppm\n\
  w = write anglehist after particle interchange resolved (unless -x)\n\
OPTIONS:\n\
 -f@ = format of g-functions [default=%%8.4f]\n\
 -k# = for KEY k, will add horizontal line at # to mark the limit\n\
 -o  = O-O frame [default=one molecule centered, other symmetrically repeated]\n\
 -m  = one merged `molecule' and corresponding gol-file [default=trajectory]\n\
 -n# = for KEY asr, # of maximum correlation function values [default=100]\n\
 -r# = for KEY sr, interpolated output at given r instead of all distances\n\
 -x  = turn off particle interchange preprocessing (some cfgs twice)",arg[0]);
  exit(1); }

  if (arg[2][0]=='-') ERROR(("KEY expected"))

  loop (ir,3,narg) if (arg[ir][0]=='-') switch (arg[ir][1]) {
    case 'f': fmt=arg[ir]+2; break;
    case 'm': merged++; break;
    case 'o': OOFRAME++; break;
    case 'k': Kirkwood=atof(arg[ir]+2); break;
    case 'r': R=atof(arg[ir]+2); break;
    case 'n': NM=atoi(arg[ir]+2); break;
    case 'x': pint=0; break;
    default: ERROR(("unknown option")); }
  else ERROR(("option expected"))

  strcpy(fn0,arg[1]);
  if ( (c=strstr(fn0,".har")) )
  if (strend(fn0)-c==4) *c=0;

  fprintf(stderr,"system=%s  KEY=%s  FRAME=%s  NM=%d  FORMAT=%s\n",
          fn0,arg[2],OOFRAME?"O-O":"center",NM,fmt);

  strcpy(fn,fn0); strcat(fn,".har");
  VarOpen(fn,"r");
  harmonics=VarGetSds();

  put(harmonics->nmeas)

  if (!harmonics->nmeas || !harmonics->N) ERROR(("nothing to do"))

  qxyz = harmonics->V / (harmonics->N*(harmonics->N-1)
                        *((double)harmonics->nmeas)*(8*PI/Sqr(G2DGRID)));

  loop (rel,0,3) {

    strcpy(fn,fn0); 

    switch (rel) {
      case 0: IF("h") strcat(fn,".hg"); else goto skip; break;
      case 1: IF("g") strcat(fn,".hgr"); else goto skip; break;
      case 2: IF("k") strcat(fn,".rkf"); else goto skip; }
    
    out=fopen(fn,"wt");

    prt("# %d %s molecules  %d measurements\n\
# cutoff=%g grid=%g nhist=%d\n\
# T=%g rho=%g\n\
# Note: glLm=<g YlLm> (angle aver.), Y000=1, Y100=sqrt(3) cos(theta1), ...",
        harmonics->N,harmonics->info,harmonics->nmeas,
        harmonics->cutoff,harmonics->grid,harmonics->nhist,
        harmonics->T,harmonics->rho);

    if (rel<2) {
      static char hdr[320]="\
# 1=A     2=B     3=C     4=D     5=E     6=F     7=G    8=H    #\n\
#                -g010            g020    g11-1  <cos    <cos^2 #\n\
# r       g000    g100    g110    g200    g111    alpha>  alpha>#";
      if (rel) strcat(hdr,"\n\
#                /g000   /g000   /g000   /g000   /g000   /g000  #");
      header(hdr);

      q = harmonics->V / (2*PI*harmonics->N*(harmonics->N-1)
                          *Sqr((double)harmonics->nmeas)/Cub(harmonics->grid));

      loop (ir,0,harmonics->nhist) {
        r=(ir+0.5)/harmonics->grid;
        qq=q/(ir*(ir+1.0)+1.0/3);
        g   = harmonics->hist[ir].n *qq;
        gc  = harmonics->hist[ir].c *qq;
        gc2 = harmonics->hist[ir].c2*qq;
        gcc = harmonics->hist[ir].cc*qq;
        ga  = harmonics->hist[ir].a *qq;
        gaa = harmonics->hist[ir].aa*qq;
        if (harmonics->hist[ir].n) {
          prt_("%6.3f ",r);
          prt_(fmt,g);
          qq=1; if (rel) qq=1/g;
          prt_(fmt,sqrt(0.75)*gc*qq);
          prt_(fmt,3*gcc*qq); 
          prt_(fmt,sqrt(1.25)*(1.5*gc2-g)*qq); 
          prt_(fmt,-1.5*(ga-gcc)*qq);
          prt_(fmt,ga*qq);
          prt(fmt,gaa*qq); }
      } }
    else {
      double
        rkf=1,r0=0,
        qkf=2./((harmonics->N-1)*(double)harmonics->nmeas);

      r=0;
      header("# r     running Kirkwood factor");

      loop (ir,0,harmonics->nhist) {
        r=ir/harmonics->grid;
        if (harmonics->hist[ir].a) {
          prt("%5.2f %7.4f",r,rkf);
          if (!r0) r0=r; }
        rkf+=qkf*harmonics->hist[ir].a; }
      if (Kirkwood) prt("\n%5.2f %7.4f\n%5.2f %7.4f",r0,Kirkwood,r,Kirkwood);
    }

    header("");

    fclose(out);
  skip:; }

#if NHARM>1

  VarRead(anglehist,sizeof(anglehist));
  fprintf(stderr,"anglehist read\n");

  strcpy(fn,fn0); strcat(fn,".ah");
  IF ("a") out=fopen(fn,"wt");
  else out=stdout;

  IF ("asrw") {
    int ir,iphi,xphi,itheta1,itheta2,ialpha1,ialpha2,i,j,zero;
    int ip,dia1,dia2,ia1,ia2,smart,totheta2;
    int irfrom,irto;
    double all=NHARM*NHARM*NHARM*NHARM*NHARM;
    double r;
    double maxomega,omegat1,sumh,sumhsmart,qharm;
    FILE *plb,*gol;
    float hdr[2];
    struct {
      double g,omega,phi,theta1,theta2,alpha1,alpha2;
      int h,dalpha1,dalpha2;
    } *G,g;
    struct {
      double omega,theta1,theta2;
      int itheta1,itheta2; 
    } Minomega,Maxomega;
    char *fnx;
    /* cos(theta) weights (NHARM==12):
       ihist weight    weight/max
       0 11 0.034074  0.1316525 
       1 10 0.099900  0.38598559 
       2 9  0.158919  0.61401441
       3 8  0.207107  0.80019915   
       4 7  0.241181  0.93185165
       5 6  0.258819  1
  */

    /* divisors of NHARM/2 allowed because of symmetries */
    static int alphastride[NHARM]={4,2,1,1,1,1}; /* and symmetrically */

    alloc(G,NM*sizeof(G[0]));

    loop (i,0,NHARM/2) alphastride[NHARM-1-i]=alphastride[i];
    IF ("a") {
      prt("alphastride, as function of thetai (used for \"smart angles\"):");
      loop (i,0,NHARM) prt_("%2d ",alphastride[i]);
      _n }

    maxomega=sqr(sin(PI/NHARM));

    if (pint) {
      loop (ir,0,R_NHIST) 
        loop (iphi,0,NHARM)
          loop (itheta1,0,NHARM) loop (itheta2,0,NHARM)
            loop (ialpha1,0,NHARM) loop (ialpha2,0,NHARM)
              if (itheta1+itheta2>=NHARM) {
                /* merging histogram bins equivalent due to particle interchange */
                int *unnorm=&anglehist[ir][iphi][itheta1][itheta2][ialpha1][ialpha2];
                int *norm=&anglehist[ir][iphi][NHARM-1-itheta2][NHARM-1-itheta1][ialpha2][ialpha1];
                if (unnorm==norm) ERROR(("internal"))
                *norm+=*unnorm;
		*unnorm=0; }
      fprintf(stderr,"particle interchange normalized\n"); }
    else 
      fprintf(stderr,"WARNING: particle interchange not normalized\n");

    IF ("w") {
      FILE *h;
      int ir;
#define WS(T,X) { T x; x=X; fwrite(&x,sizeof(T),1,h); }

      strcpy(fn,fn0); strcat(fn,".fcf");
      h=fopen(fn,"wb");
      WS(float,harmonics->T);
      WS(float,harmonics->rho);
      WS(float,harmonics->V/harmonics->nmeas);
      WS(int,harmonics->N);
      WS(int,harmonics->nmeas);
      WS(int,R_NHIST);
      WS(int,NHARM);
      loopto (ir,0,R_NHIST) WS(float,R_INV_INDEX(ir))
	if (fwrite(anglehist,sizeof(anglehist),1,h)!=1)
	  ERROR(("writing %s",fn))
#undef WS
        fclose(h); }

    irfrom=0;
    irto=R_NHIST;
    if (R)
      loop (ir,irfrom,irto-1) {
        double r0=sqrt( ( sqr(R_INV_INDEX(ir+1)) + sqr(R_INV_INDEX(ir)) )/2 );
	double r1=sqrt( ( sqr(R_INV_INDEX(ir+2)) + sqr(R_INV_INDEX(ir+1)) )/2 );

	if (r0<R && R<r1) {
	  RQ=(r1-R)/(r1-r0);
	  irfrom=ir; irto=ir+1;
	  break; } }

    loop (ir,irfrom,irto) {
      r=sqrt( ( sqr(R_INV_INDEX(ir+1)) + sqr(R_INV_INDEX(ir)) )/2 );
      if (R) r=R;
      else fprintf(stderr,"%2.0f%%\r",(double)ir/(0.01*R_NHIST));
      qharm = (harmonics->V 
	   / (harmonics->N*(harmonics->N-1)*Sqr((double)harmonics->nmeas)) )
	/( (cub(R_INV_INDEX(ir+1))-cub(R_INV_INDEX(ir)))*(2*PI/3)
	   *Cub(PI/NHARM)) * (4*Cub(PI));

      g.dalpha2=g.dalpha1=1;
      sumh=sumhsmart=0;
      zero=0;
    
      loop (smart,0,2) {
	int isplb=0;
	double unity;

	IF (smart?"s":"r") isplb=1;
      
	Minomega.omega=9e9;
	Maxomega.omega=0;
	memset(G,0,NM*sizeof(G[0]));
      
	loop (itheta1,0,NHARM) {
	  g.theta1=itheta1*(180./NHARM)+(90./NHARM);
	  if (smart) g.dalpha1=alphastride[itheta1];
	  omegat1=(cosd(g.theta1-90./NHARM)-cosd(g.theta1+90./NHARM))*g.dalpha1;
        
	  if (pint) totheta2=NHARM-itheta1; else totheta2=NHARM;

	  loop (itheta2,0,totheta2) {
	    g.theta2=itheta2*(180./NHARM)+(90./NHARM);
	    if (smart) g.dalpha2=alphastride[itheta2];
	    g.omega=omegat1*(cosd(g.theta2-90./NHARM)-cosd(g.theta2+90./NHARM))
	      *g.dalpha2;
	    if (pint && itheta2+itheta1<NHARM-1) g.omega*=2;
          
	    if (g.omega>Maxomega.omega+1e-9) {
	      Maxomega.omega=g.omega;
	      Maxomega.itheta1=itheta1;
	      Maxomega.itheta2=itheta2; 
	      Maxomega.theta1=g.theta1;
	      Maxomega.theta2=g.theta2; }
          
	    if (g.omega<Minomega.omega-1e-9) {
	      Minomega.omega=g.omega;
	      Minomega.itheta1=itheta1;
	      Minomega.itheta2=itheta2; 
	      Minomega.theta1=g.theta1;
	      Minomega.theta2=g.theta2; }
          
	    for (ialpha1=0; ialpha1<NHARM; ialpha1+=g.dalpha1) {
	      g.alpha1=ialpha1*(180./NHARM)+(90./NHARM)*g.dalpha1;
            
	      for (ialpha2=0; ialpha2<NHARM; ialpha2+=g.dalpha2) {
		g.alpha2=ialpha2*(180./NHARM)+(90./NHARM)*g.dalpha2;
              
		loop (iphi,0,NHARM*(1+smart)) {

		  if (smart) {
		    xphi=0;
		    g.h=0;
		    loop (dia1,0,g.dalpha1)
		      loop (dia2,0,g.dalpha2) {
		      ia1=ialpha1+dia1;
		      ia2=ialpha2+dia2;
		      ip=iphi + (itheta1<NHARM/2?-1:1)*dia1
			+ (itheta2<NHARM/2?-1:1)*dia2;
		      xphi+=ip;
		      ip=(ip+NHARM*2)%(NHARM*2);
		      if (ip<NHARM)
			if (R) 
			  g.h+=anglehist[ir][ip][itheta1][itheta2][ia1][ia2]*RQ
			    +anglehist[ir+1][ip][itheta1][itheta2][ia1][ia2]*(1-RQ);
			else
			  g.h+=anglehist[ir][ip][itheta1][itheta2][ia1][ia2];
		      /* ???: ip>=NHARM turned off not to include twice... */
                    }
		    sumhsmart+=g.h;
		    g.phi=fmod(xphi/(double)(g.dalpha1*g.dalpha2)*(180./NHARM)+(90./NHARM)+720.,360.);
                  } 
		  else {
		    if (R) 
		      g.h+=anglehist[ir][iphi][itheta1][itheta2][ialpha1][ialpha2]*RQ
			+anglehist[ir+1][iphi][itheta1][itheta2][ialpha1][ialpha2]*(1-RQ);
		    else
		      g.h=anglehist[ir][iphi][itheta1][itheta2][ialpha1][ialpha2];
		    sumh+=g.h;
		    if (!g.h) zero++;
		    g.phi=iphi*(180./NHARM)+(90./NHARM);
                  }
                
		  g.g=g.h/g.omega*qharm;
		  if (g.g>G[NM-1].g) {
		    loop (i,0,NM) if (g.g>G[i].g) {
		      for (j=NM-1; j>i; j--) G[j]=G[j-1];
		      G[i]=g;
		      /* particle exchange REMOVED, see exchange.old */
		      break; } } 
		}
              
	      } } } }
      
	fnx=string("%c%c%c-%.0f.plb",
		   fn[0],
		   strchr(fn,'c')?'c':'f',
		   "rs"[smart],r*1000); /* WAS A BUG HERE (+.5) */
	IF ("a") prt("\nr=%.3f file=%s %s angles\n\
__i________g_________h____omega____phi_theta1_theta2_alpha1_alpha2/a1_a2",
		     r,fnx,smart?"smart":"normal");
	gol=NULL;
	if (isplb) {
	  plb=fopen(fnx,"wb");
	  hdr[0]=6; hdr[1]=0;
	  if (merged) {
	    static int first=1;
	    int ns=NM*4*(pint+1);
	    hdr[0]=6*ns;
	    strcpy(strend(fnx)-3,"gol");
	    gol=fopen(fnx,"wt");
	    fprintf(gol,"!sphere#1\n! %s\n%.0f\n",fnx,hdr[0]); 
	    if (first) {
	      first=0;
	      fprintf(stderr,"to use show, make a `merged' mol file:\n\
  molcfg -%d w2 %d\n\
and `ln -s' it to the name of plb and gol file selected\n\
don't use `molcfg -%d w2 NAME' directly, gol-file would be overwritten\n"
		      ,ns,ns,ns);
	    } }
	  fwrite(hdr,4,2,plb); }
	
	loop (i,0,NM) if (G[i].h) {
	  vector w[3]; /* H1,H2,O */
	  int xmirror,ymirror;
	  int ig;

	  if (i==0) unity=G[i].g;        
	  if (gol) loop (ig,0,4*(1+pint))
		     fprintf(gol,"\
YELLOW  %.3f\n\
YELLOW  %.3f\n\
BLUE    %.3f\n\
WHITE   %.3f\n\
WHITE   %.3f\n\
RED     %.3f\n",
			     G[i].g/unity,G[i].g/unity,G[i].g/unity*1.9,
			     G[i].g/unity,G[i].g/unity,G[i].g/unity*1.9);

	  IF ("a") prt("%3d %11.3f %6d %8.5f %6.1f %6.1f %6.1f %6.1f %6.1f %2d %2d",
		       i+1,G[i].g,G[i].h,G[i].omega/maxomega,
		       G[i].phi, G[i].theta1, G[i].theta2, G[i].alpha1, G[i].alpha2,
		       G[i].dalpha1, G[i].dalpha2);

	  if (OOFRAME) {
	    INITw
	    rotd(w[0],2,G[i].alpha1); rotd(w[0],0,G[i].theta1); 
	    rotd(w[1],2,G[i].alpha1); rotd(w[1],0,G[i].theta1);
	    loop (j,0,3) w[j][2]-=r/2;
	    if (isplb) fwrite(w,4,9,plb);             
	    INITw
	    rotd(w[0],2,G[i].alpha2); rotd(w[0],0,G[i].theta2); rotd(w[0],2,G[i].phi);
	    rotd(w[1],2,G[i].alpha2); rotd(w[1],0,G[i].theta2); rotd(w[1],2,G[i].phi);
	    loop (j,0,3) w[j][2]+=r/2;
	    if (isplb) fwrite(w,4,9,plb); }
	  else { /* !OOFRAME */
	    loop (xmirror,0,2) loop (ymirror,0,2) {
	      INITw              
	      if (isplb) fwrite(w,4,9,plb);             
	      loop (j,0,3) {
		rotd(w[j],2,G[i].alpha2); 
		rotd(w[j],0,G[i].theta2);
		rotd(w[j],2,G[i].phi);
		w[j][2]+=r;
		rotd(w[j],0,-G[i].theta1);
		rotd(w[j],2,-G[i].alpha1); }
	      if (xmirror) loop (j,0,3) w[j][0]*=-1;
	      if (ymirror) loop (j,0,3) w[j][1]*=-1;
	      if (isplb) fwrite(w,4,9,plb); 
          
	      if (pint) {
		/* particle interchange */
		INITw
		if (isplb) fwrite(w,4,9,plb);             
		loop (j,0,3) {
		  rotd(w[j],2,G[i].alpha1);
		  rotd(w[j],0,-G[i].theta1);
		  rotd(w[j],2,-G[i].phi);
		  w[j][2]-=r;
		  rotd(w[j],0,G[i].theta2);
		  rotd(w[j],2,-G[i].alpha2); }
		if (xmirror) loop (j,0,3) w[j][0]*=-1;
		if (ymirror) loop (j,0,3) w[j][1]*=-1;
		if (isplb) fwrite(w,4,9,plb); } }
          } /* !OOFRAME */
        }
	if (isplb) fclose(plb); 
	if (gol) fclose(gol);
      
	IF ("a") {
	  if (ir==R_NHIST-1) {
	    prt("");
	    prt("max omega: %7.5f %7.5f for %d %d (theta1=%g theta2=%g)",
		Maxomega.omega, Maxomega.omega/maxomega,
		Maxomega.itheta1,Maxomega.itheta2,
		Maxomega.theta1,Maxomega.theta2);
	    prt("min omega: %7.5f %7.5f for %d %d (theta1=%g theta2=%g)",
		Minomega.omega, Minomega.omega/maxomega,
		Minomega.itheta1,Minomega.itheta2,
		Minomega.theta1,Minomega.theta2); }
	  put3(sumh,sumhsmart,zero)
          put3(all,sumh/all,zero/all) }
      } /* smart */
      if (!R) if (sumh!=sumhsmart)
	ERROR(("sumh=%.0f != sumhsmart=%.0f",sumh,sumhsmart))
      } /* ir */
  }
#endif

  free(harmonics);

#ifdef G2DGRID
  {
    int i,ir,j,m,method,ng2d;
    FILE *ppm,*flt;
    unsigned char rgb[3];
    double dm;
    
    out=stdout;

    if (VarFile.size) {
      
      VarRead(&ng2d,sizeof(ng2d));
      put(ng2d)
      qxyz/=ng2d;

      loop (i,0,2) {
	VarRead(g2d,sizeof(g2d));
	m=0;
	dm=0;
    
	loop (ir,0,G2DGRID*G2DMAXR) 
	  loop (j,0,G2DGRID*G2DMAXR) {
	  Max(m,g2d[ir][j])
	    Max(dm,g2d[ir][j]/((ir+0.5)*(j+0.5))) }
        
	loop (method,0,2) {
	  strcpy(fn,fn0);
	  if (method) strcat(fn,"r");
	  strcat(fn,i?"HH.ppm":"OH.ppm");
	  IF ("p") ppm=fopen(fn,"wb");
	  strcpy(fn,fn0);
	  strcat(fn,i?"HH.xyz":"OH.xyz");
	  if (!method) IF ("2") flt=fopen(fn,"wb");
	  IF ("p") fprintf(ppm,"P6\n# CREATOR: %s %s\n%d %d\n255\n",arg[0],arg[1],
			   G2DGRID*G2DMAXR,G2DGRID*G2DMAXR);  
    
	  loop (ir,0,G2DGRID*G2DMAXR) 
	    loop (j,0,G2DGRID*G2DMAXR) {

	    double gray;
	    float fg=(double)g2d[ir][j]*qxyz;
        
	    if (method) {
	      gray=g2d[ir][j]/((ir+0.5)*(j+0.5))/dm;
#if NCOL!=7
	      gray=fg*QPPM*Sqr(G2DGRID)/((ir+0.5)*(j+0.5));
#endif
            }
	    else {
	      gray=g2d[ir][j]/(double)m;
	      IF ("2") fwrite(&fg,4,1,flt); 
#if NCOL!=7
	      gray=fg*QPPM;
#endif
            }

	    if (ir>G2DGRID*5/4 && ir<G2DGRID*7/4) gray=j/(G2DGRID*G2DMAXR-1.0);
          
	    rgbcolor(rgb,gray);
	    if ((ir+1)%G2DGRID<2 && (j+1)%G2DGRID<2)
	      /*.....     if (rgb[0]==255)*/
	      rgb[0]=rgb[1]=rgb[2]=0;
	    IF ("p") fwrite(rgb,3,1,ppm); }
	  if (!method) IF ("2") fclose(flt);
	  IF ("p") fclose(ppm); }
      } }
    else
      fprintf(stderr,"g2d missing\n");
  }
#endif

  VarClose();

  return 0;
}
