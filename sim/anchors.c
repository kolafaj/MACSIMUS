/* #include code for Shake - anchoring site (in the Shake loop)
   (there are r(t+h) in p[*])
 */

#if SHAKE!=1
#  error SHAKE=1 required for (this version of) ANCHOR
#endif /*# SHAKE!=1 */

if (mn->anchor & ANCHOR_s) {
  /* anchor given site */
  vector cm;
  double rr;

  VVV(cm,=p[mn->iaxis[0]],-mn->r0)
  normalizer(cm);
  rr=SQR(cm);
  if (rr>anchor.rr) {
    if (it>1) goto skipanchor;
    notinplace++;
    rr=sqrt(anchor.rr/rr);
    /* whole molecule moved if site not in place */
    VO(cm,*=rr) 
    loop (i,0,mn->ns) VV(p[i],-=cm) }
  else {
    /* NOTES: called repeatedly ==> -=
       (- because we want to have a constraint force) */
    VV(cf,-=si[mn->iaxis[0]].mass/(h*h)*cm)
    VV(p[mn->iaxis[0]],-=cm) }
 skipanchor:;
}
