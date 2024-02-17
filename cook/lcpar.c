/* run slabbs (of the linked-cell list method) in parallel */

#ifdef BARRIER
static
void *pll_linkcell_bar(void *arg) /**************** pll_linkcell_bar:threads */
{
  int dix;
  /* NB: (pthread_t*)arg=No.thread+ith (ith=thread number)*/

  loop (dix,0,par1_lc.maxdix) {
    /* dix passed as a multiple of No.th added to arg */
    pll_linkcell((void*)((pthread_t*)arg+dix*No.th));
    if (!par1_lc.forceserial) pthread_barrier_wait(&par1_lc.barrier); }

  return arg;
}
#endif /*# BARRIER */

static double linkcell(linklist_t ****list) /********************** linkcell */
{
  double U=0;
  int ith,dix;
  pll_global_t *pll_global; /* array [No.cell[0]] */

  par1_lc.forceserial=measure && No.measureserial;
  par1_lc.maxdix=2+(int)((box.cutoff-1e-5)/(box.L[0]/No.cell[0]));

  if (No.bonded) {
    max14=0;
    if (box.excrlimit==0) excrrlimit=Sqr(box.cutoff);
    else excrrlimit=Sqr(box.excrlimit); }

  if (No.th<1 || No.th>No.cell[0]) ERROR(("invalid number of No.th"))

  par1_lc.by=No.th;
  par1_lc.list=list;
  allocarrayzero(pll_global,No.th);
  par1_lc.glob0=pll_global;

#ifdef BARRIER
  if (par1_lc.forceserial) {
    /* serial loop */
    loop (ith,0,No.th)
      pll_linkcell_bar(No.thread+ith); }
  else {
    if (pthread_barrier_init(&par1_lc.barrier,NULL,No.th)) ERROR(("cannot create barrier"))
    PARALLELIZE(pll_linkcell_bar,No.th)
                                                             if (pthread_barrier_destroy(&par1_lc.barrier)) ERROR(("cannot destroy barrier")) }
#else /*# BARRIER */
  /* dix=0: inter-slab interactions, dix=1: neighboring slabs, etc. */
  loop (dix,0,par1_lc.maxdix) {
    if (par1_lc.forceserial) {
      /* serial loop */
      loop (ith,0,No.th) {
        par1_lc.dix=dix;
        pll_linkcell(No.thread+ith); } }
    else {
      par1_lc.dix=dix;
      PARALLELIZE(pll_linkcell,No.th) }
  } /* end of dix loop */
#endif /*#!BARRIER */

  /* sum up all instances of En.el, En.vir ... */
  if (measure) loop (ith,0,No.th) {
#if PRESSURETENSOR&PT_VIR
    loop (dix,0,PT_DIM) En.Pvir[dix]+=pll_global[ith].Pvir[dix];
#endif /*# PRESSURETENSOR&PT_VIR */
    U += pll_global[ith].U;
    En.el += pll_global[ith].Enel;
    En.vir += pll_global[ith].Envir; }

  if (partimes.rspace) loop (ith,0,No.th)
    partimes.rspace[ith]+=pll_global[ith].t/(double)CLOCKS_PER_SEC;

  loop (ith,0,No.th) Max(max14,pll_global[ith].max14)

  free(pll_global);

  if (No.c) update14();

  return U;
} /* linkcell */
