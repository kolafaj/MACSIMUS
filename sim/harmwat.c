/* 
  #include file for simmeas
  harmonic analysis of water models and dipolar dumbells
  WARNING: <cosa^2> not included in measurements -- see harmlin.c
  error with init=1 fixed, though

  ANGLE CONVENTIONS

  phi is oriented angle of plane OO-axis2 vs. plane OO-axis1 as if
  rotated around vector 1-->2 (i.e., the angle of the projection of
  OO-axis2 vs. the projection of OO-axis1 when watched from the side
  if molecule 2)

  thetai is the angle between axisi and vector 1-->2

  alphai is oriented angle between the water molecule plane and the
  OO-axisi plane

  mirror symetry is used (phi=-phi, alpha1=-alpha1, alpha2=-alpha2) to
  normalize phi into [0,180).  Should not be used in chiral
  environment.  

  WARNING : WATERPLUS adds also OM : OO 2D rdf (for acetone, this is CC : OC)
*/

struct harmonics_sds *harmonics;
#include "orienta.c"
#define xDEBUG

/* 3 TIP3P: order of sites = HHO */
/* 4 TIP4P: order of sites = HOMH */
/* 5 ST2: order of sites=OLHLH - deprecated */
/* 5 TIP5P: LHOLH */

#if NHARM>1
int anglehist
  /*     r     phi   theta1 theta2 alpha1 alpha2  */
    [R_NHIST][NHARM][NHARM][NHARM][NHARM][NHARM];
#endif
#ifdef G2DGRID
int ng2d;
int g2doh[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
int g2dhh[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
#ifdef WATERPLUS
int g2doc[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
int g2dcc[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
#endif
#endif

void measureharmonics(void)
{
  int ns=spec[0]->ns,n1,n2;
  int O=0,H1=1,H2=2;
  vector *R1,*R2,r1,r2,dr;
  double rr,x;
  vector *axis;
  double sumcosa=0;
#ifdef G2DGRID
  vector r2h[3]; /* corresponds to R1; for normal water only [2] used */
#endif

#if NHARM>1
  vector *hh;
  double fromr=R_INV_INDEX(0)+1e-12,tor=R_INV_INDEX(R_NHIST)-1e-12;
  alloc(hh,No.N*sizeof(hh[0]));
#endif

#ifdef G2DGRID
  ng2d++;
#endif

  switch (ns) {
    case 2: O=0,H1=1,H2=1; strcpy(harmonics->info,"dumbell"); break;
    case 3: O=2,H1=0,H2=1; strcpy(harmonics->info,"TIP3P"); break;
    case 4: O=1,H1=0,H2=3; strcpy(harmonics->info,"TIP4P"); break;
      /* ...this works for acetone too, 1=C=center */
      //    case 5: O=0,H1=2,H2=4; strcpy(harmonics->info,"ST2"); break;
    case 5: O=2,H1=1,H2=4; strcpy(harmonics->info,"TIP5P"); break;
    default: ERROR(("ns=%d: not known water model",ns)) }
  
  alloc(axis,No.N*sizeof(axis[0]));

  /* normalized axes of molecules */
  loop (n1,0,No.N) {
    R1=a[0]->rp+n1*ns;
    VVVV(axis[n1],=R1[H1],+R1[H2],-2*R1[O]) 
      rr=sqrt(SQR(axis[n1]));
#if NHARM>1
    VVV(hh[n1],=R1[H1],-R1[H2])
#endif
    VO(axis[n1],/=rr) }

  harmonics->V+=PROD(box.L);
  harmonics->nmeas++;

  loop (n1,0,No.N) {
    R1=a[0]->rp+n1*ns;
    VV(r2,=R1[O])
#ifdef G2DGRID
    VV(r2h[0],=R1[H1])
    VV(r2h[1],=R1[H2])
#ifdef WATERPLUS
    VV(r2h[2],=R1[2]) /* this is O while O=1 is actually center=C */
#endif
#endif
    loop (n2,0,n1) {
      R2=a[0]->rp+n2*ns;
      VV(r1,=R2[O])
	/* WARNING: using intermac.h here:
	   dr is vector from r2 to r1 ==> making r1<-->r2 swap! */
#ifdef G2DGRID
#ifdef FREEBC
#error FREEBC not supported
#endif
      {
	int i,ir;
	vector r1h[3]; /* corresponds to R2 */

	rr=0;
	loop (i,0,3) {
	  dr[i]=r1[i]-r2[i];
	  r1h[0][i]=R2[H1][i];
	  r1h[1][i]=R2[H2][i];
#ifdef WATERPLUS
	  r1h[2][i]=R2[2][i]; /* O */
#endif
	  if (dr[i]>box.Lh[i]) dr[i]-=box.L[i],r1[i]-=box.L[i],r1h[0][i]-=box.L[i],r1h[1][i]-=box.L[i];
	  if (dr[i]<-box.Lh[i]) dr[i]+=box.L[i],r1[i]+=box.L[i],r1h[0][i]+=box.L[i],r1h[1][i]+=box.L[i];
	  rr+=Sqr(dr[i]); }

	x=sqrt(rr);
	ir=(int)(x*G2DGRID);
#define ADDH(H,V1,V2) { int ii=(int)(sqrt(SQRD(V1,V2))*G2DGRID); \
  if (ii<G2DGRID*G2DMAXR) H[ir][ii]++; }
    if (ir<G2DGRID*G2DMAXR) {
      loop (i,0,2) {
	ADDH(g2doh,r1,r2h[i])
        ADDH(g2doh,r2,r1h[i])
	ADDH(g2dhh,r2h[0],r1h[i])
	ADDH(g2dhh,r2h[1],r1h[i]) }
#ifdef WATERPLUS
      ADDH(g2doc,r1,r2h[2])
      ADDH(g2doc,r2,r1h[2])
      ADDH(g2dcc,r2h[2],r1h[2])
#endif
      }
#undef ADDH
      }
#else
      NI(0) NI(1) NI(2) rr=SQR(dr);
      x=sqrt(rr);
#endif
      { 
	double cosa=SCAL(axis[n1],axis[n2]);
	int i=(int)(x*harmonics->grid);
      
	sumcosa+=cosa;

	if (i<harmonics->nhist) {
	  double cost1=SCAL(axis[n1],dr)/x;
	  double cost2=SCAL(axis[n2],dr)/x;
	  harmonics->hist[i].n++;
	  harmonics->hist[i].c +=cost1-cost2;
	  harmonics->hist[i].c2+=Sqr(cost1)+Sqr(cost2);
	  harmonics->hist[i].cc+=cost1*cost2;
	  harmonics->hist[i].a +=cosa;
	  harmonics->hist[i].aa+=Sqr(cosa); 

#if NHARM>1
	  if (x<tor && x>fromr) {
	    double phi,alpha1,alpha2,theta1,theta2;
	    int iphi;

	    i=R_INDEX(x);
	    if (i<0 || i>=R_NHIST) ERROR(("internal")) /* temporary */

            theta1=acos(cost1);
	    theta2=acos(cost2);
	    ORIENTEDANGLE(phi,axis[n1],axis[n2],dr,rr,x)
            ORIENTEDANGLE(alpha1,dr,hh[n1],axis[n1],1,1)
	    ORIENTEDANGLE(alpha2,dr,hh[n2],axis[n2],1,1)

	      /* mirror symmetry 
		 (should NOT be used in chiral environment, as e.g. water
		 around an optical active solute) */
	      iphi=(int)((phi+2*PI)*(NHARM/PI))%(NHARM*2);
	    if (iphi>=NHARM) {
#ifdef DEBUG
	      prt_("mirror ");
	      phi=-phi;
#endif
	      iphi=NHARM*2-1-iphi;
	      alpha1=-alpha1; alpha2=-alpha2;
            }
	    /* NOTE: particle exchange symmetry is more complex and is
	       postponed to harmg.c */

	    anglehist[i] [iphi]
#define PI_NHARM(X) [(int)((X+3*PI)*(NHARM/PI))%NHARM]
	      PI_NHARM(theta1)
	      PI_NHARM(theta2)
	      PI_NHARM(alpha1)
	      PI_NHARM(alpha2) ++;
#ifdef DEBUG
	    if (alpha1<-PI/2) alpha1+=PI;
	    if (alpha2<-PI/2) alpha2+=PI;
	    prt("r=%.4f ph=%.3f t1=%.3f t2=%.3f a1=%.3f a2=%.3f",
		x,phi*180/PI,theta1*180/PI,theta2*180/PI,alpha1*180/PI,alpha2*180/PI);
#endif
	  }
#endif
        } }
    } }

  StaSet(0,lag.err,2,lag.n);
  StaAdd("<cosangle(axis1,axis2)>",sumcosa/(No.N*(No.N-1)/2));

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
#ifdef WATERPLUS
  memset(g2doc,0,sizeof(g2doc));
  memset(g2dcc,0,sizeof(g2dcc));
#endif
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
    VarRead(g2dhh,sizeof(g2dhh)); 
#ifdef WATERPLUS
    VarRead(g2doc,sizeof(g2doc));
    VarRead(g2dcc,sizeof(g2dcc));
#endif
  }
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
  backup("har");
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
#ifdef WATERPLUS
  VarPut(g2doc,sizeof(g2doc));
  VarPut(g2dcc,sizeof(g2dcc));
#endif

#endif
  VarClose();
  prt("%s saved",lastFn);
}
