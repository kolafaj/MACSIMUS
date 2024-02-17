/* #include code for Shake - fixing sites without measurements
   to be used for not ANCHOR version with -k0
   see initfix() and fixforces()
 */

#ifndef ANCHOR

#if SHAKE!=1
#  error SHAKE=1 required for (this version of) ANCHOR
#endif /*# SHAKE!=1 */

if (option('k')==0 && fixa) {
  fixsites_t *fix;
  vector *rpf=fixa->rp;
  
  looplist (fix,fixsites0)
    loop (i,fix->from,fix->to) {
      VV(rp[i],=rpf[i])
      VO(rp1[i],=0)
      VO(rp2[i],=0) } }
#endif
