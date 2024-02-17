/* 
  #include file for simmeas
  harmonic analysis of linear molecules
*/

struct harmonics_sds *harmonics;
#include "orienta.c"
#define xDEBUG

#if NHARM>1
int anglehist
  /*     r     phi   theta1 theta2 */
    [R_NHIST][NHARM][NHARM][NHARM];
#endif
#ifdef G2DGRID
int ng2d;
int g2doh[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
int g2dhh[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
#endif

void measureharmonics(void)
{
  int ns=spec[0]->ns,n1,n2;
#if CO2VERSION==2
  int O=1,H1=0,H2=2; /* CO2: in fact, H1,H2 are O1,O2 and O stands for C */
#else
  int O1=1,O2=2,H1=0,H2=3; /* O:charges, H:LJ */
#endif
  vector *R1,*R2,r1,r2,dr;
  double rr,x;
  vector *axis;
  double sumcosa=0;
  double sumcc=0;
#ifdef G2DGRID
  vector r2h[2]; /* corresponds to R1 */
#endif

#if NHARM>1
  vector *hh;
  double fromr=R_INV_INDEX(0)+1e-12,tor=R_INV_INDEX(R_NHIST)-1e-12;
  alloc(hh,No.N*sizeof(hh[0]));
#endif

#ifdef G2DGRID
  ng2d++;
#endif

  alloc(axis,No.N*sizeof(axis[0]));

  /* normalized axes of molecules */
  loop (n1,0,No.N) {
    R1=a[0]->rp+n1*ns;
#if NHARM>1
    VVV(hh[n1],=R1[H1],-R1[H2])
#endif
    VVV(axis[n1],=R1[H1],-R1[H2]) 
    rr=sqrt(SQR(axis[n1]));
    VO(axis[n1],/=rr) }

  harmonics->V+=PROD(box.L);
  harmonics->nmeas++;

  loop (n1,0,No.N) {
    R1=a[0]->rp+n1*ns;
#if CO2VERSION==2
    VV(r2,=R1[O])
#else
    VVV(r2,=0.5*R1[O1],+0.5*R1[O2])
#endif
#ifdef G2DGRID
    VV(r2h[0],=R1[H1])
    VV(r2h[1],=R1[H2])
#endif
    loop (n2,0,n1) {
      R2=a[0]->rp+n2*ns;
#if CO2VERSION==2
      VV(r1,=R2[O])
#else
      VVV(r1,=0.5*R2[O1],+0.5*R2[O2])
#endif
/* WARNING: using intermac.h here:
   dr is vector from r2 to r1 ==> making r1<-->r2 swap! */
#ifdef G2DGRID
#ifdef FREEBC
#error FREEBC not supported
#endif
      {
	int i,ir;
	vector r1h[2]; /* corresponds to R2 */

	rr=0;
	loop (i,0,3) {
	  dr[i]=r1[i]-r2[i];
	  r1h[0][i]=R2[H1][i];
	  r1h[1][i]=R2[H2][i];
	  if (dr[i]>box.Lh[i]) dr[i]-=box.L[i],r1[i]-=box.L[i],r1h[0][i]-=box.L[i],r1h[1][i]-=box.L[i];
	  if (dr[i]<-box.Lh[i]) dr[i]+=box.L[i],r1[i]+=box.L[i],r1h[0][i]+=box.L[i],r1h[1][i]+=box.L[i];
	  rr+=Sqr(dr[i]); }

	x=sqrt(rr);
	ir=(int)(x*G2DGRID);
#define ADDH(H,V1,V2) { int ii=(int)(sqrt(SQRD(V1,V2))*G2DGRID); \
  if (ii<G2DGRID*G2DMAXR) H[ir][ii]++; }
	if (ir<G2DGRID*G2DMAXR) loop (i,0,2) {
	  ADDH(g2doh,r1,r2h[i])
	  ADDH(g2doh,r2,r1h[i])
          ADDH(g2dhh,r2h[0],r1h[i])
          ADDH(g2dhh,r2h[1],r1h[i]) }
#undef ADDH
      }
#else
      NI(0) NI(1) NI(2) rr=SQR(dr);
      x=sqrt(rr);
#endif
      { 
	double cosa=SCAL(axis[n1],axis[n2]),cc=Sqr(cosa);
	int i=(int)(x*harmonics->grid);
      
	sumcosa+=cosa;
	sumcc+=cc;

	if (i<harmonics->nhist) {
	  double cost1=SCAL(axis[n1],dr)/x;
	  double cost2=SCAL(axis[n2],dr)/x;
	  harmonics->hist[i].n++;
	  harmonics->hist[i].c +=cost1-cost2;
	  harmonics->hist[i].c2+=Sqr(cost1)+Sqr(cost2);
	  harmonics->hist[i].cc+=cost1*cost2;
	  harmonics->hist[i].a +=cosa;
	  harmonics->hist[i].aa+=cc; 

#if NHARM>1
	  if (x<tor && x>fromr) {
	    double phi,theta1,theta2;
	    int iphi;

	    i=R_INDEX(x);
	    if (i<0 || i>=R_NHIST) ERROR(("internal")) /* temporary */

            theta1=acos(cost1);
            theta2=acos(cost2);
            ORIENTEDANGLE(phi,axis[n1],axis[n2],dr,rr,x)

            /* mirror symmetry */
            iphi=(int)((phi+2*PI)*(NHARM/PI))%(NHARM*2);
            if (iphi>=NHARM) {
#ifdef DEBUG
              prt_("mirror ");
	      phi=-phi;
#endif
              iphi=NHARM*2-1-iphi;
            }
	    /* NOTE: particle exchange symmetry is more complex and is
	       postponed to harmg.c */

	    anglehist[i] [iphi]
#define PI_NHARM(X) [(int)((X+3*PI)*(NHARM/PI))%NHARM]
	      PI_NHARM(theta1)
	      PI_NHARM(theta2) ++;
	    }
#endif
        } }
    } }

  StaSet(0,lag.err,2,lag.n);
  StaAdd("<cosangle(axis1,axis2)>",sumcosa/(No.N*(No.N-1)/2));
  StaAdd("<cosangle(axis1,axis2)^2>",sumcc/(No.N*(No.N-1)/2));

  free(axis);
#if NHARM>1
  free(hh);
#endif
}

void initharmonics(double gr,double cutoff,double T,double rho)
{
  int nh=(int)(cutoff*abs(gr)+1.0000001);

  if (nspec!=1) ERROR(("only nspec=1 supported"))

  sdsalloczero(harmonics,sizeof(struct harmonics_sds)+(nh-1)*sizeof(struct harmitem_s));
  harmonics->grid=gr;
  harmonics->nhist=nh;
  harmonics->cutoff=cutoff;
  harmonics->cq=Sqr(cutoff);
  harmonics->T=T;
  harmonics->rho=rho;
  harmonics->N=No.N;
  /* init static arrays ... this added after the Water II paper */
  memset(anglehist,0,sizeof(anglehist));
  ng2d=0;
  memset(g2doh,0,sizeof(g2doh));
  memset(g2dhh,0,sizeof(g2dhh));
}

void loadharmonics(double gr,double cutoff)
{
  VarOpen(Fn("har"),"r");
  VarReadSds(harmonics);
#if NHARM>1
  VarRead(anglehist,sizeof(anglehist));
#endif
#ifdef G2DGRID
  if (VarFile.size) {
    VarRead(&ng2d,sizeof(ng2d));
    VarRead(g2doh,sizeof(g2doh));
    VarRead(g2dhh,sizeof(g2dhh)); }
  else
    WARNING(("initializing g2d"))
#endif
  VarClose();
  prt("%s loaded",lastFn);
  if (harmonics->grid!=gr || harmonics->cutoff!=cutoff)
    ERROR(("harmonics grid=%g cutoff=%g, loaded grid=%g cutoff=%g",
	   gr,cutoff,harmonics->grid,harmonics->cutoff))
}

void saveharmonics(void)
{
  VarOpen(Fn("har"),"w");
  VarPutSds(harmonics);
#if NHARM>1
  VarPut(anglehist,sizeof(anglehist));
#endif
#ifdef G2DGRID
  VarPut(&ng2d,sizeof(ng2d));
  put(ng2d)
  VarPut(g2doh,sizeof(g2doh));
  VarPut(g2dhh,sizeof(g2dhh));
#endif
  VarClose();
  prt("%s saved",lastFn);
}
