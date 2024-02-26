/* #include code for Shake
   (called after SHAKE iterations have been finished)

   constrain the center-of-mass of all not-marked molecules
*/

if (anchor.xyz) {
  vector cm={0,0,0};
  double M=0,rr;
  
  loop (n,FROM,No.N) {
    mn=molec+n;
    if (!mn->anchor) {
      p=rof(mn,rp2);
      si=spec[mn->sp]->si;
      ns=mn->ns;
      M+=spec[mn->sp]->mass;
      loop (i,0,ns) VV(cm,+=si[i].mass*p[i]) } }
  VO(cm,/=M)
  VV(cm,-=anchor.r0)

  //  normalizer(cm);
  /* NOTES: called only once (=> no +=)
     -m here because we want to have a constraint force */
  if (anchor.xyz==7) {
    VV(cf,=(-M/(h*h))*cm)
    rr=SQR(cm); }
  else {
    /* for selected coordinates only 
       NB: other parts of the code are inefficiently calculated for whole vectors) */
    int k,b=1;

    rr=0;
    loop (k,0,3) {
      if (anchor.xyz&b) {
        rr+=Sqr(cm[k]);
        cf[k]=(-M/(h*h))*cm[k]; }
      else
        cf[k]=cm[k]=0;
      b<<=1; } }

  if (rr>anchor.rr) {
    notinplace++;
    rr=sqrt(anchor.rr/rr);
    VO(cm,*=rr) }
  
  loop (n,FROM,No.N) {
    mn=molec+n;
    if (!mn->anchor) {
      p=rof(mn,rp2);
      loop (i,0,mn->ns) VV(p[i],-=cm) } }

  printanchor(FORCEUNIT,cf,anchor.xyz);
}
