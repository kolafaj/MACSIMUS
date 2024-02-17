/* #include code for Shake - anchoring CM of groups (in the Shake loop)
   (there are r(t+h) in p[*])
   note: group->cf cleared in anchorm.c and printed in anchora.c
 */

#if SHAKE!=1
#  error SHAKE=1 required for (this version of) ANCHOR
#endif /*# SHAKE!=1 */

if (mn->anchor & ANCHOR_g) {
  struct anchorgroup_s *group;

  looplist (group,mn->group) {
    vector cm={0,0,0};
    double rr;
    double m=0;
    int j;

    loop (j,0,group->ns) {
      i=group->site[j];
      VV(cm,+=si[i].mass*p[i])
      m+=si[i].mass; }
  VO(cm,/=m)
  VV(cm,-=group->r0)
  normalizer(cm);
  /* NOTE: called repeatedly => -=
     (- because we want to have a constraint force) */
  VV(group->cf,-=(-m/(h*h))*cm)
  rr=SQR(cm);
  if (rr>anchor.rr) {
    if (it>1) goto skipanchorg;
    notinplace++;
    rr=sqrt(anchor.rr/rr);
    VO(cm,*=rr) }
  loop (j,0,group->ns) {
    i=group->site[j];
    VV(p[i],-=cm) }
  skipanchorg:; 
  }
}
