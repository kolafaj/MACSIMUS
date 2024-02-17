/* BJERRUM only - to be #included to main.c */

if (bj.mode&1)
  bjerrum_read_oo(); /* read in the 1st pass only */

if (bj.mode&128)
  bjerrum_topology(cfg[0]); /* DEPRECATED */

if (bj.mode&2) {

  depend_r(cfg[0],0);
  
  if (bj.mode&4) {
    /* Calculate defects at end of job from the averaged configuration.
       NPT-compatible: The averaging takes into account the periodic
       boundary conditions; i.e., it is performed with the configuration
       scaled to box $1^3$ and the averaged box is rescaled back to the
       averaged box sizes. */
    int i,k;
    static vector Lsum;

    if (icyc==0) {
      if (!bj.cfg) {
      /* allocate if needed */
        sdsalloc(bj.cfg,cfg[0]->size);
        sdscopy(bj.cfg,cfg[0]); }
      /* set 1st config, init box size averaging */
      VV(Lsum,=box.L)
#ifdef POLAR
      loop (i,No.s,2*No.s) VV(bj.cfg->rp[i],=cfg[0]->rp[i])
#endif
      loop (i,0,No.s) VVV(bj.cfg->rp[i],=cfg[0]->rp[i],/box.L) }
    else {
      /* icyc>0: sum up and follow periodic b.c. (scaled to 1^3) */
      VV(Lsum,+=box.L)
#ifdef POLAR
      loop (i,No.s,2*No.s) VV(bj.cfg->rp[i],+=cfg[0]->rp[i])
#endif
      loop (i,0,No.s) {
        loop (k,0,3) {
          double x=cfg[0]->rp[i][k]/box.L[k];
          double xav=bj.cfg->rp[i][k]/icyc;

          if (x<xav-0.5) x+=1.0;
          if (x>xav+0.5) x-=1.0;
          bj.cfg->rp[i][k]+=x; } } }
    if (sig || icyc==no-1) {
      /* average and process */
#ifdef POLAR
      loop (i,No.s,2*No.s) VO(bj.cfg->rp[i],/=(icyc+1))
#endif
      loop (i,0,No.s) VVO(bj.cfg->rp[i],*=Lsum,/Sqr(icyc+1))
# if 0
      prt("%d -3 mainbjerrum.c",No.s);
      prt("%g %g %g mainbjerrum.c",Lsum[0]/(icyc+1),Lsum[1]/(icyc+1),Lsum[2]/(icyc+1));
      loop (i,0,No.s) prt("%g %g %g mainbjerrum.c",VARG(bj.cfg->rp[i]));
#endif
      bjerrum_gauss(bj.cfg,bj.from,bj.to,bj.q,bj.eps); } }
  else
    /* every cycle */
    bjerrum_gauss(cfg[0],bj.from,bj.to,bj.q,bj.eps);
}
