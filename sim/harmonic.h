/*
HARMONICS = 3 : water models
HARMONICS = 2 : linear molecules
HARMONICS = 1 : MeOH (partly implemented)

again, WATERPLUS turns on acetone instead of water (affects R_NHIST...)

  ax1
   /
  / theta1
 +------------------+- - - 
            dr       \ theta2
                      \
                      ax2
theta1 = angle(ax1,dr)
theta2 = angle(ax2,dr)
alpha = angle(ax1,ax2)
phi = angle(plane ax1-dr,plane ax2-dr)
cos(alpha) = sin(theta1) sin(theta2) cos(phi) + cos(theta1) cos(theta2)

g(r,ax1,ax2) = SUM glLm(r) YlLm   [l>=0, L>=0]
               lLm

YlLm = Yl(m1) YL(-m2), m=m1-m2, phi=phi1-phi2

glLm(r) = <g YlLm>  [angle average]

Y000 = 1
Y100 = sqrt(3) cos(theta1)
Y010 = -sqrt(3) cos(theta2) [use g100=-g010]
Y110 = 3 cos(theta1) cos(theta2)
Y200 = sqrt(1.25) (3cos^2(theta1)-1)
Y020 = sqrt(1.25) (3cos^2(theta2)-1) [use g200=g020]
Y111 = -1.5 sin(theta1) sin(theta2) cos(phi) [Im part irrelevant]
Y11-1=Y111 [use g11-1=g111]
g001=g00-1=g101=g011...=0 for isotropic fluid
*/

#if HARMONICS==3 /* water */
#define NHARM 12
#define G2DMAXR 8 /* must be integer */
#elif HARMONICS==2 /* CO2 */
#define NHARM 24
#define G2DMAXR 10 /* must be integer */
#elif HARMONICS==1 /* MeOH - partly */
#define NHARM 1
#define G2DMAXR 10 /* must be integer */
#else
#error bad HARMONICS
#endif

/*
  NHARM==1: only harmonic expansion upto cos^2 (P2)
  NHARM>1: +angle analysis
*/

#define G2DGRID 50
/*
  defined: 2D site-site rdf
*/
  
struct harmitem_s {
  int n;     /* 1 */
  double c;  /* cos(theta1) - cos(theta2) */
  double c2; /* cos^2(theta1) + cos^2(theta2) */
  double cc; /* cos(theta1) cos(theta2) */
  double a;  /* cos(alpha) */
  double aa; /* cos^2(alpha) */
  };

#if NHARM>1
#if HARMONICS==3
#define R_NHIST 32

#ifdef WATERPLUS 
/* this gives range of r 3.22..9.5 */
#define R0_4 1.34
#else
/* this gives range of r 2.36421 .. 7.51286 */
#define R0_4 1.24
#endif

#define R_INDEX(X) ((int)((sqrt(sqrt(X))-R0_4)*77.0))
#define R_INV_INDEX(I) Pow4((I)/77.0+R0_4)
extern int anglehist
  /*     r     phi   theta1 theta2 alpha1 alpha2  */
    [R_NHIST][NHARM][NHARM][NHARM][NHARM][NHARM];
#define HARMSIZE (R_NHIST*NHARM*NHARM*NHARM*NHARM*NHARM)
#elif HARMONICS==2
#define R_NHIST 80
/* this gives range of r 2.36421..9.48 */
#define R_INDEX(X) ((int)((sqrt(sqrt(X))-1.24)*160.0))
#define R_INV_INDEX(I) Pow4((I)/160.0+1.24)
extern int anglehist
  /*     r     phi   theta1 theta2  */
    [R_NHIST][NHARM][NHARM][NHARM];
#define HARMSIZE (R_NHIST*NHARM*NHARM*NHARM)
#endif
#if HARMSIZE*4 > 60*1024*1024
#error too large HARMONICS histogram
#endif
#endif

#if CO2VERSION==1
/* DLJD in reduced units */
#undef R_INDEX
#undef R_INV_INDEX
#define R_INDEX(X) ((int)((sqrt(sqrt(X))-0.95)*160.0))
#define R_INV_INDEX(I) Pow4((I)/160.0+0.95)
#endif

extern struct harmonics_sds {
  int size;
  int nmeas,nhist,N;
  char info[8];
  double V,grid,cutoff,cq,T,rho;
  struct harmitem_s hist[1/*variable*/]; } *harmonics;

#ifdef G2DGRID
extern int ng2d; /* # of measurements */
extern int g2doh[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
extern int g2dhh[G2DGRID*G2DMAXR][G2DGRID*G2DMAXR];
#endif
