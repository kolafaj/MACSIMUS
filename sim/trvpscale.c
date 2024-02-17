/* 
  #included from constrd.c, Shake():
  - TRVP (Time-Reversible Velocity Predictor):
    dx.doi.org/10.1021/ct200108g
    J. Chem. Theory Comput. 2011, 7, 3596-3607
  - constraint (bond) scaling predictor for T_NPT
  - follow given V(t) code (tau.rho<0)
*/

  if (thermostat==T_NOSE || thermostat>=T_NPT || tau.rho<0 || Eext.isB) {
    static int nose_k=-1,lastk=-1;

    if (nose_k>=0 && No.pred!=nose_k)
      ERROR(("TRVP predictor length unexpectedly changed (%d->%d)",nose_k,No.pred))
    if (nose_k<0 && No.pred==0)
      WARNING(("SHAKE+Nose or V(t) without a predictor (option -p0) is imprecise"))
    nose_k=No.pred;

    if (!Vpred) {
      allocarray(b,No.pred+1); /* TRVP predictor coefficients */
      allocarray(Vpred,No.eq); /* predicted velocities (actually, v*h as ususal) 
                                  Vpred[0] -> logs (denoted as \xi in the manual)
                                  Vpred[1,2,3] -> lambda 
                                  vpred -> rp (real velocities) */
      vpred=(vector*)(Vpred+RPOFFSET); /* atom velocities (->rp) start here */
      allocarray(Vhist,No.pred+1); /* the velocity history */
      loopto (j,1,No.pred) allocarray(Vhist[j],No.eq);
      k=0; /* no history known now */ }

    Vhist[0]=&cfg[1]->logs; /* pointer passed - array [No.eq] not allocated */

    /* generate TRVP predictor coefficients */
    if (k!=lastk) {
      b[0]=(double)(2*k+1)/(k+1);
      loop (j,0,k) b[j+1]=-b[j]*(double)(k-j)/(k+j+2);
      prt("TRVP for Verlet+SHAKE with Nose-Hoover, k=%d (max %d):",k,No.pred);
      loopto (j,0,k) prt("  b[%d]=%.15g",j,b[j]); }
    lastk=k;

    /* @4 the velocity predictor for all degrees of freedom: 
       xi=logs, lambda[3], v[No.s][3]) */
    loop (i,0,No.eq) {
      Vpred[i]=0;
      loopto (j,0,k) Vpred[i]+=b[j]*Vhist[j][i]; }

    if (thermostat==T_NOSE)
      Vfactor=Vpred[0];
    else if (thermostat>=T_NPT) {
      /* predicting the change of background scaling (constraints do not change) */
      /* lambda[0] = Vpred[1], lambda[1] = Vpred[2], lambda[2] = Vpred[3] */
      if (rescale&RESCALE_PT) {
        /* anisotropic scaling (used in SHAKE) @53 */
        loop (j,0,3) if (rescale&(1<<j)) {
          if (k>1) scale[j]=exp(2*Vhist[0][j+1]-Vhist[1][j+1]);
          else scale[j]=exp(Vhist[0][j+1]); }
        else
          Vpred[j+1]=0,scale[j]=1;
        // @43
        Vfactor=Vpred[0]+(Vpred[1]+Vpred[2]+Vpred[3])/No.f; }
      else {
        Vfactor=Vpred[0]+No.NPT*Vpred[1]; // @41; No.NPT: see siminit.c
        /* @51 isotropic scaling (used in SHAKE) */
        if (k>1) scale[0]=scale[1]=scale[2]=exp(2*Vhist[0][1]-Vhist[1][1]);
        else scale[0]=scale[1]=scale[2]=exp(Vhist[0][1]); } }
    else if (tau.rho==-3) {
      /* 3 box sizes from file */
      static int pass=1;
      if (pass) {
        WARNING(("tau.rho=-3: unfinished and/or not tested"))
        pass--; }
      int k;

      if (!Vpred) allocarray(Vpred,4);
      loop (k,0,3) {
        // @5b
        scale[k]=Lfromfile(k,t+h)/Lfromfile(k,t);
        Vfactor=0;
        Vpred[0]=0; /* probably not needed */
        /* @4b ? check that Vpred is not used beyond [4] */
        Vpred[k+1]=log(Lfromfile(k,t+h/2)/Lfromfile(k,t-h/2)); } }
    else if (tau.rho==-1) {
      /* cube cell size from file */
      Vpred[0]=Vpred[1]=Vpred[2]=Vpred[3]=0; /* probably not needed */
      scale[0]=scale[1]=scale[2]=Lfromfile(0,t+h)/Lfromfile(0,t); // @5b
      // @4b:
      Vfactor=log(Lfromfile(0,t+h/2)/Lfromfile(0,t-h/2)); }
    
    // DO NOT USE (unfinished attempt of CM-based version)
    //    if (rescale&RESCALE_CM) scale[0]=scale[1]=scale[2]=1;
    
    /* increase k up to No.pred */
    k=min(No.pred,k+1);

    /* history shift */
    for (j=k; j; j--)
      loop (i,0,No.eq) Vhist[j][i]=Vhist[j-1][i]; }
#if VERLET==4
  else
    ERROR(("implementation limitation\n\
*** version VERLET=4 (predicted velocities) requires a Nose-Hoover thermostat"))
#endif /*# VERLET==4 */
