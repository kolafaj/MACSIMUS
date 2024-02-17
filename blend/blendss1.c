/* version allowing cutoff */

#  ifdef POLAR
#    error POLAR not supported with cutoff version
#  endif /*# POLAR */

static struct {
  double C1,C2,C1q,C2q,alpha,beta,betaq,s0,s1,s2,f;
  vector r0,r0minusC2,r0plusC2;     /* to optimize out-of-cutoff tests */
  int is_not;
} cutoff;

#  ifdef OMITKEPT
static enum keep_e keepmask;
#  endif /*# OMITKEPT */

static double sspot( /************************************************ sspot */
                    /* vector r0: see structure cutoff */
                    vector r1,
                    vector f0,vector f1,
                    site_t *s0,site_t *s1,
                    int onefour)
/***
    site-site potential (Lennard-Jones + charge-charge)
    cutoff.r0 is vector of site s0, the force to site s0 will be summed to f0
    r1 is vector of site s1, the force to site s1 will be summed to f1
    energy is returned
    NOTES:
     r0 has moved from the parameter list to the static structure cutoff
     incompatible with sim/intrapot.c - always works as measure==1
***/
{
  double U=0,rr,f,r,x,y,z,s,qq;
  vector dr;
  sitesite_t *ss;
#ifdef METAL
  static double rho1,rho2,rho12; /* for METAL only, TO BE EXTENDED */
#endif
  
  VVV(dr,=cutoff.r0,-r1)

#  ifdef DEBUG
  if ( !cutoff.is_not && SQR(dr)/3>cutoff.C2q) ERROR(("debug: cutoff cube"))
#  endif /*# DEBUG */

  if ( (rr=SQR(dr))<cutoff.C2q ) {

    ss=&sstab[s0->type][s1->type][onefour];
    SS_MEASURE /* #defined in sim/XXX/sitesite.h */
    ULJ+=U;

    /* electrostatic potential */
    if ( (qq=s0->charge*s1->charge)!=0 ) {
      r=sqrt(rr);
      qq*=obfuscate/r;
      if (onefour) qq*=factor14;
      f+=qq/rr;
      U+=qq;
      Uel+=qq; }

    if (checkpairs) {
      int i,j;
      char infos[PAIRINFOLEN],*info=NULL;

      if (signpairs) {
        U*=signpairs;
        loop (i,0,maxpairs)
          if (U>maxpair[i].U) {
            for (j=maxpairs-1;j>i;j--) maxpair[j]=maxpair[j-1];
            maxpair[i].U=U;
            info=maxpair[i].info;
            break; }
        U*=signpairs; }

      if (!info) if (verbose2) info=infos;

      if (info) {
        // in the new version don't know whether fixed...
        //        info[0] = isfix ? 'f' : onefour ? '*' : ' ';
        info[0] = onefour ? '*' : ' ';
        sprintf(info+1,
                "%4i %-4s %-9s - %4i %-4s %-9s  r=%6.3f",
                (int)(s0-site),atom[s0->type].name,s0->id,
                                 (int)(s1-site),atom[s1->type].name,s1->id,
                                                  sqrt(rr));

        if (verbose2) prt("!%s  U=%g LJ=%g el=%g",info,U,U-qq,qq); }

    } /* checkpairs */

    /* multiply by the switch function */
    if (rr>cutoff.C1q) {
      y=rr-cutoff.alpha,
        z=Sqr(y)-cutoff.betaq;
      s=((cutoff.s2*z+cutoff.s1)*z+cutoff.s0)*y+0.5;
      f=f*s-U*cutoff.f*Sqr(z);
      U*=s; }

    /* forces summed to f0,f1 */
    VVO(f0,+=dr,*=f) VV(f1,-=dr)
  } /* in cutoff sphere */

  return U;
}
