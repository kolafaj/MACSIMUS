/*
   Charge density profile of point charges, both POLAR and nonpolar.
   #included from slabmeas.c #ifndef GAUSSIANCHARGES for geom3=2
*/
#ifdef POLAR
  if (si[i].chargepol==0)
    dpr[nspec+sp]->hist.dhist[nh]+=si[i].charge;
  else if (slabmode&8) {
    /* Drude charge density profile for short dipoles:
       using triangle (instead of rectangle) histogram bin
       algorithm changed in V3.4n */
    double zx=rc*Grid+nbins-0.5;
    int nm=(int)zx;
    double dz=zx-nm;

    double zxD=zx+rpol[i][2]*Grid;
    int nmD=(int)zxD;
    double dzD=zxD-nmD;
    
    dpr[nspec+sp]->hist.dhist[nm%nbins]+=si[i].charge*(1-dz);
    dpr[nspec+sp]->hist.dhist[(nm+1)%nbins]+=si[i].charge*dz;
    
    dpr[nspec+sp]->hist.dhist[nmD%nbins]+=si[i].chargepol*(1-dzD);
    dpr[nspec+sp]->hist.dhist[(nmD+1)%nbins]+=si[i].chargepol*dzD; }
  else {
    /* charge density profile: Drude oscillator directly as 2 charges
       not good for very short dipoles */
    int ndrude=(int)((rc+rpol[i][2])*Grid+nbins)%nbins;

    dpr[nspec+sp]->hist.dhist[nh]+=si[i].charge;
    dpr[nspec+sp]->hist.dhist[ndrude]+=si[i].chargepol; }
#else /*# POLAR */
  dpr[nspec+sp]->hist.dhist[nh]+=si[i].charge;
#endif /*#!POLAR */
