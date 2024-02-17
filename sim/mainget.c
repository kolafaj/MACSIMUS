/* to be directly #included from main.c
  main module reading the def- and get-files
*/
    do {
      /* In the first invocation (pass=0 and 1), this loop is execuded twice:
         - read SYSNAME.def until ';', close SYSNAME.def
         - open SIMNAME.get and read it until ';'
         In any following invocation (pass>1):
         - and read SIMNAME.get until the next ';'
      */

      int *N; /* for pass=0 only */

      if (pass==0) {
        allocarray(N,nspec);
        allocarray(group,nspec);
        loop (i,0,nspec) {
          /* numbers of molecules (because spec[i]->N might have been set in the ble-file) */
          N[i]=spec[i]->N;
          /* groups of species */
          group[i]=i; } }

      if (pass==1) {
        loop (i,0,nspec) spec[i]->N=N[i];
        free(N);

        // initfix, checkfixed moved in V 3.5i because the # of molecules may change after loadcfg
        //        if (option('k')>0) initfix(Fn("fix"));
        initNo(); /* now we know N[] to call it,
                     but will be recalculated if N,s from cfg applies */
        //        checkfixed();

        /* close .def file, open .get file or use stdin */
        fclose(in);
        if (option('s')) in=stdin;
        else in=fopen(Fn("get"),"rt");
        if (in==NULL) ERROR(("open %s",lastFn)) }

      if (pass>=1) {
        _n putline('%',79);
        if (option('s')) prt("Ctrl-D to quit / ? for help / enter data as VAR=EXPR / ; to end data:"); }

      oldt=t;
#ifdef SLAB
      if (!wall.LJ && nsites) {
        allocarray(wall.LJ,nsites);
        loop (i,0,nsites) wall.LJ[i].neps=wall.LJ[i].sig=-3e33; }
#endif /*# SLAB */
      getdata
        /* auxiliary variables for calculations */
        static double
          aux=0, a=0, b=0, c=0, x=0, y=0, z=0;
        static int
          i=0,j=0,k=0,n=0;

        if (pass==0) {
          getvec(N,,nspec)
#ifdef ECC
          get(el.ecc) get(el.epsf)
#endif /*# ECC */
          get(equalize.mol) get(equalize.cfg) get(equalize.sp) }

        get(AllocSizeLim)

        get(el.centroid) get(el.bg) get(el.epsq)
#ifndef FREEBC
        get(el.epsk) get(el.epsr)
        get(el.alpha) get(el.kappa)
        get(el.epsinf) get(el.test) get(el.diag) get(el.diff)
        get(el.sf)
        get(el.rplus) get(el.rshift)
        get(el.sat)
        get(virial)
        get(tau.sat)
#endif /*# FREEBC */

        get(t)

        getvec(group,,nspec)
        getvec(center.K,,DIM)
        getvec(center.r0,,DIM)
        get(el.Perr)

#ifdef SLAB
        get(slab.grid) get(slab.max)
        get(slab.sp) get(slab.mode) get(slab.prt)
        get(slab.T) get(slab.Tz0) get(slab.Tz1)
        get(slab.sym)
        getkey(slab.geom,slabgeomkey)
        get(slab.out)
        getvec(slab.n,,NCENTER)
        getvec(slab.ns,,NCENTER)
        getvec(slab.z,,NCENTER)
        getvec(slab.Kz,,NCENTER)
        getvec(slab.z0,,NCENTER)
        getvec(slab.z1,,NCENTER)
        get(slab.K) get(slab.range)
        get(slab.ext.center) get(slab.ext.span) get(slab.ext.zero)
        get(wall.n) get(wall.g)
        get(wall.rho) getvec(wall.z,,2)
        if (pass>0) {
          getvec(wall.LJ,.sig,nsites)
          getvec(wall.LJ,.neps,nsites) }
#  if SLAB & 1
        get(tau.L) getvec(slab.Lx,,4) getvec(slab.Ly,,4)
#  endif /*# SLAB & 1 */
#  if SLAB & 2
        get(cleave.sigma)
        get(cleave.K)
        get(cleave.n)
        get(cleave.i)
        get(cleave.z0)
        get(cleave.z1)
        get(cleave.init)
#  endif /*# SLAB & 2 */
#  if SLAB & 4
        get(slab.wall.sig) get(slab.wall.epsrho)
#  endif /*# SLAB & 4 */
#endif /*# SLAB */

        getvec(center.cmK,,DIM)
        get(center.cmn)
        getvec(el.E,,DIM) get(el.f) getvec(el.phase,,DIM)
        getvec(el.B,,DIM)
        get(el.m.sp) get(el.m.plus) get(el.m.minus) get(el.m.m)
        get(dV)
        getkey(init,initkey)
        get(mirror) get(unsplit)
        getkey(sort,sortkey)
        getkey(pins,pinskey) get(initrho) getkey(MC,yesno)
        get(E) get(rho) getvec(L,,DIM)
#if 0
        getvec(box.Lx,,4) getvec(box.Ly,,4) getvec(box.Lz,,4)
        getvec(box.rho,,4)
#endif /*# 0 */
        get(nshift) getvec(shift,,DIM) getvec(vshift,,DIM)
        get(tau.E) get(tau.rho)
        getkey(rescale,rescalekey)
        get(tau.T) get(T) get(initvel)
        getkey(thermostat,thermostatkey)
        get(T_tr_in)
        get(tau.P) get(P) get(maxscale) get(bulkmodulus)
        get(tau.sig) get(tau.i) get(tau.j)
        getvec(load.n,,DIM)
        get(load.N)
        getvec(load.L,,DIM)
        getkey(load.tr,loadtrkey)
        getkey(load.zero,loadzerokey)
        get(nm.method) get(nm.dr) get(nm.eps) get(nm.key) get(nm.zero)
        get(nm.frames) get(nm.modes) get(nm.ampl) get(nm.mem)
        get(omega) get(omegac)
        get(eps) get(epsc)
        get(no) get(noint) get(h)
        get(nomax) get(stop) get(Tstop)
        get(reread.from) get (reread.to) get (reread.by)
        get(dt.prt) get(dt.cfg) get(dt.plb)
        get(lag.err) get(lag.n)
        get(lag.J) get(lag.M)
        get(lag.CM) get(lag.LM) get(lag.AM)
        /* now always autoset:        getvec(box.center,,DIM) */

        get(gear.init)
        getvec(gear.C,,MAXGEARORDER)

        get(removemol.n)

#if PRESSURETENSOR&PT_OFF
        get(lag.Pt)
#endif /*# PRESSURETENSOR&PT_OFF */
        get(Emax) get(drmax)
        get(rdf.grid) get(rdf.cutoff) get(rdf.onefour)
        get(el.grid) get(el.minqq) getkey(el.corr,rescalekey)
        get(cutoff) get(LJcutoff) get(corr) get(poteps)
        get(CPnbit)
#ifdef POLAR
        get(scf.eps) get(scf.epsq) get(scf.omega) get(scf.maxit)
        get(scf.epsx) get(scf.omegax) get(scf.test)
        get(scf.E) get(scf.Estride)
        get(tau.dip)
        get(scf.domega) get(scf.margin)
#endif /*#!POLAR */
#ifdef SHEAR
        get(shear)
#endif /*# SHEAR */
#ifdef LINKCELL
        getvec(No.cell,,DIM)
        get(No.percell)
        get(box.over14) get(box.max14)
        get(box.rmin)
#endif /*# LINKCELL */
#ifdef PERSUM
        get(No.molspan)
#endif /*# PERSUM */
        /* structure wall defined in simdef */
#ifdef LOG
        get(No.first)
#endif /*# LOG */
        get(No.rotatefrom)
        get(drift) get(conserved)
        getkey(quit,yesno)
#ifdef WATERPLUS
        {
#  include "watercut.h"
          get(cut.from) get(cut.to)
        }
#endif /*# WATERPLUS */
        if (quit==-4) {
          prt("WARNING: quit=\"shell\"=-4 is deprecated, use key=\"shell\"=8 instead");
          quit=0;
          key=8; }

        getkey(key,keykey)
        if (key) {
          prt("executing key=%d:",key);
          switch (key) {
            case 1: case 2: case 3:
            case -1: case -2: case -3:
              sortmolecules(key);
              break;
            case 4:
              if (system(string("showcp -p99 %s",simils.simname))) fprintf(stderr,"showcp failed\n");
              break;
            case 5:
              if (system(string("show %s",simils.simname))) fprintf(stderr,"show failed\n");
              break;
            case 6:
              if (system(string("rdfg %s -g -p",simils.simname))) fprintf(stderr,"rdfg -g failed\n");
              break;
            case 7:
              if (system(string("rdfg %s -c -p",simils.simname,simils.simname))) fprintf(stderr,"rdfg -c failed\n");
              break;
            case 8:
              /* the shell (DOS,EMX removed...) */
              fprintf(stderr,"\n *** Type `exit' to return to %s! ***\n",arg[0]);
              i=system(getenv("SHELL")?getenv("SHELL"):"bash");
              prt("shell return code=%d",i);
              break;
            case 9: goto TheEnd;
          }

          key=0; }

        get(cache)
#ifdef DIHHIST
        get(dih.res) get(dih.grid) get(dih.cp) get(dih.dcp)
#endif /*# DIHHIST */

#ifdef WIDOM
        get(widom.sp) get(widom.spreal) get(widom.n)
#  ifdef SLAB
        get(widom.z0) get(widom.z1) get(widom.dz) get(widom.mode)
#  endif /*# SLAB */
#endif /*# WIDOM */

        get(diff.mode)

#ifdef BJERRUM
          get(bj.mode) get(bj.from) get(bj.to) get(bj.q) get(bj.eps)
#endif /*# BJERRUM */

#ifdef CLUSTERS
        get(cl.mode) get(cl.format)
        get(cl.maxn) get(cl.maxcluster) get(cl.mincluster)
#endif /*# CLUSTERS */

#ifdef XSECTION
        getkey(xs.mode,xsmodekey)
        get(xs.grid) get(xs.RvdW) get(xs.Rscale)
        get(xs.freq) get(xs.maxs) get(xs.sizelimit)
#  ifdef CLUSTERS
        get(xs.mincluster)
#  endif /*# CLUSTERS */
#endif /*# XSECTION */

#ifdef RGYR
        getvec(rg.end,,2) get(rg.cp)
        get(lag.v) get(lag.nv) get(lag.dim)
#endif /*# RGYR */

        /* auxiliary variables */
        get(aux)
        get(a) get(b) get(c)
        get(x) get(y) get(z)
        get(i) get(j) get(k)
        get(n)

        checkdata
      enddata
      newt=t;
      _n
      onefourinrdf=rdf.onefour;

      if (quit>0) {
        prt("WARNING: quit>0 is deprecated, use key=\"quit\" or key=9 instead");
        goto TheEnd; }
    } while (pass++==0);
