/*
   Charge density profile of Gaussian charges, both POLAR and nonpolar.
   The true Gaussian distribution is used to calculate the histogram.
   New in V3.4o, bug in V3.4n, charges treated as point in older versions.
   #included from slabmeas.c #ifdef GAUSSIANCHARGES for geom3==2
*/
  if (si[i].charge) {
    int range=si[i].esig*5.9*Grid+1,k;
    int iz=rc*Grid+nbins;

    loopto (k,iz-range,iz+range) {
      double t1=((k-nbins+1)/Grid-rc)/si[i].esig;
      double t0=((k-nbins  )/Grid-rc)/si[i].esig;

      // factor /2 moved to printdpr; dpr[nspec+sp]->hist.dhist[k%nbins]+=si[i].charge/2*(erf(t1)-erf(t0)); 
      dpr[nspec+sp]->hist.dhist[k%nbins]+=si[i].charge*(erf(t1)-erf(t0)); } }

#ifdef POLAR
/* the Gaussian sigma is the same for the atom charge and Drude charge */
  if (si[i].chargepol) {
    int range=si[i].esig*5.9*Grid+1,k;
    double zpol=rc+rpol[i][2];
    int iz=zpol*Grid+nbins;

    loopto (k,iz-range,iz+range) {
      double t1=((k-nbins+1)/Grid-zpol)/si[i].esig;
      double t0=((k-nbins  )/Grid-zpol)/si[i].esig;

      // factor /2 moved to printdpr; dpr[nspec+sp]->hist.dhist[k%nbins]+=si[i].chargepol/2*(erf(t1)-erf(t0)); 
      dpr[nspec+sp]->hist.dhist[k%nbins]+=si[i].chargepol*(erf(t1)-erf(t0)); } }
#endif

