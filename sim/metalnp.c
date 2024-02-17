/*
  metal - not pair
*/

struct eldens_s eldens;

#ifdef SS_MEASURE_rep
#  error "SS_MEASURE_rep not supported in nonpolar version"
#endif

#ifdef NIBC
#  error NIBC+METAL not supported
#endif

double LJM(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2 GLOB)
/* LJ site-site interaction, measurements incl.                     **** LJM */
{
  double U=0,rr,f,x,y;
  real fdr;
  rdf_t *ssrdf;
  vector dr;

  int i=(vector*)r1-eldens.rref;
  int j=(vector*)r2-eldens.rref;
  real rho1=eldens.rho[i];
  real rho2=eldens.rho[j];
  real rho12=eldens.rhorho[i][j];

  beginCellImage(0,box.cq)
    if ( (ssrdf=ss->rdf) ) {
      int i=(int)(sqrt(rr)*ssrdf->grid);
      if (i<ssrdf->nhist) (ssrdf->hist[i])++; }
    SS_MEASURE
      ENVIR+=rr*f;
      MADDFORCES(f1,f2,f,dr)
  endCellImage

  MARK_PATCH
  return U;
} /* LJM */

void LJ(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2) /******* LJ */
/* LJ site-site interaction, no measurements */
{
  double rr,f,x;
  real fdr;
  vector dr;

  int i=(vector*)r1-eldens.rref;
  int j=(vector*)r2-eldens.rref;
  real rho1=eldens.rho[i];
  real rho2=eldens.rho[j];
  real rho12=eldens.rhorho[i][j];

  beginCellImage(,box.cq)
    SS_NOMEASURE(=)
    ADDFORCES(f1,f2,f,dr)
  endCellImage

  MARK_PATCH
}

double TIP4P(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
{
  ERROR((""))
  return 0;
}
double ST2(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
{
  ERROR((""))
  return 0;
}

double TIP3P(vector f1[],vector f2[],molecule_t *m1,molecule_t *m2,vector *rp)
{
  ERROR((""))
  return 0;
}

void LJQQ(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2,double qq)
{
  ERROR((""))
}

double LJQQM(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2,double qq GLOB)
{
  ERROR((""))
  return 0;
}

void XQQ(vector r1,vector r2,vector f1,vector f2,double qq)
{
  ERROR((""))
}

double XQQM(vector r1,vector r2,vector f1,vector f2,double qq GLOB)
{
  ERROR((""))
  return 0;
}

void LJQQ14(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2,double qq)
{
  ERROR((""))
}

double LJQQ14M(sitesite_t *ss,vector r1,vector r2,vector f1,vector f2,double qq GLOB)
{
  ERROR((""))
  return 0;
}
