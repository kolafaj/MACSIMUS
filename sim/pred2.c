/*
  This module implements the ASPC and PSPC methods with Verlet+SHAKE.
  ASPC (option -p2XX):
  J. Kolafa:
    "Time-Reversible Always Stable Predictor-Corrector Method for Molecular Dynamics of Polarizable Molecules",
    J. Comput. Chem. 25, 335-342 (2004)
  PSPC (Partially Stable Predictor-Corrector):
  J. Kolafa:
    "Numerical integration of equations of motion with a self-consistent field given by an implicit equation
    Mol. Simulat. 18 193-212 (1996)

  It is #included by constrd.c, function Shake()
  Options:
    -pKPC: P=0 = no action (Car-Parrinello-like method, see shakecp.c)
           P=1 = no action (no prediction=last value)
           P=2 = selects ASPC (2nd order), K=k (predictor length)
           P=3 = 3rd-order predictors
           P>3 = ERROR (see gear2.c for higher-order methods, which are hardly useful)
*/

  switch (option('p')/10%10) {
    case 0: /* Car-Parrinello-like */
    case 1: /* no prediction */
      break;
    case 3: /* PSPC (start by ASPC) */
    case 2: { /* ASPC */
      real
        *a0=(real*)polarrof(molec,cfg[0]->rp),
        *a1=(real*)polarrof(molec,cfg[1]->rp),
        pred;

      int j,nreal=DIM*No.s;

      int k=option('p')/100%10;
      int k0=option('p')/10%10;
      static int kk;
      static double *AA; /* [k] coeff. for ASPC: B[i] of the paper */
      static real **aa; /* history for ASPC;
                           note that aa[0]=a0 (for time t) and aa[1]=a1 (for time t-h) 
                           are stored in SIMNAME.cfg -- the older history is not */

      if (kk==0) {
        /* allocate all arrays */
        int i;

        if (k>8) WARNING(("ASPC: too long predictor (k=%d)",k))
        if (k<0 || k>=50) ERROR(("ASPC: k=%d out of bounds",k))

        allocarray(AA,k+2);
        AA[0]=2; AA[1]=-1; /* k=0 predictor */
        allocarray(aa,k+2);
        loop (i,2,k+2) alloczero(aa[i],nreal*sizeof(real)); }

      aa[0]=a0; aa[1]=a1; /* must be here - a0,a1 may change */

      /* the predictor + history shift; note that kk=0,1,...,k */

      loop (j,0,nreal) {
        int ik;

        pred = 0;
        loop (ik,0,kk+2) pred += AA[ik]*aa[ik][j];
        for (ik=k; ik>=0; ik--) aa[ik+1][j]=aa[ik][j];
        a0[j]=pred; }

      /* coefficients of ASPC (kk>0) */
      if (kk<k) {
        int i,PSPC;
        double o;

        if (kk==0) _n
        kk++;
        PSPC=(k==kk) && (k0>2);
        /* omega guaranteeing stability for next step(s) */
        o=(kk+2)/(double)(2*kk+3);

        if (option('^')==0) scf.omega=o; /* default = stability */

        prt("@@@@@@@ initializing %cSPC for k=%d (%s) @@@@@@@",
            PSPC?'P':'A',
            kk,k==kk?"final":"1 step");
        if (PSPC) {
          switch (kk) {
            case 1: AA[0]=3,AA[1]=-3,AA[2]=1; break;
            case 2: AA[0]=2.4,AA[1]=-1.2,AA[2]=-0.8,AA[3]=0.6; break;
            case 3: AA[0]=1.94,AA[1]=-0.46/3,AA[2]=-1.18,AA[3]=0.06,AA[4]=1./3; break;
            default: ERROR(("unknown PSPC (3rd order) predictor (k=%d)",kk)); }
          prt_("PSPC(3,%d):",kk);
          loop (i,0,kk+2) prt_(" %g",AA[i]);
          _n }
        else {
          prt("scf.omega=%f (stability guaranteed for scf.omega<=%f)",scf.omega,o);

          o=(double)(4*kk+6)/(double)(kk+3);
          loop (i,0,kk+2) {
            AA[i]=o*(i+1);
            prt("B[%d]=%f",i+1, AA[i]);
            o*= (double)(i-kk-1)/(double)(kk+i+4); } }
        if (k==kk) _n } }
      break;

    default:
      ERROR(("illegal option -p%d (with Verlet+SHAKE)",option('p')))
  }
