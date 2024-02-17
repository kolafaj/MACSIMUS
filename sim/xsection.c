/* cross-section measurements
   #ifdef XSECTION only:
     is #included from simmeasx.c TWICE:
   XSECTION_M==1: original molecule-based version
   XSECTION_M==1 && CLUSTERS: cluster version instead of molecule version
   XSECTION_M==0: whole configuration version
*/

#if XSECTION_M==2
double clXsection(ToIntPtr A) /*********************************** mXsection */
#elif XSECTION_M==1
double mXsection(ToIntPtr A) /************************************ mXsection */
#else /*# XSECTION_M */
double cXsection(ToIntPtr A) /************************************ cXsection */
#endif /*#!XSECTION_M */
/*
  Calculates cross sections based on spheres with radii given by the
  van der Waals radii + probe radius xs.RvdW, multiplied by xs.Rscale;
  the default xs.Rscale=2^(-1/6).

  xs.grid: number of sample points per Angstrom
  xs.freq: see the manual
  xs.mode:
    0 = 32-point Gaussian integration (probably by A.H. Stroud, The
        Approximate Calculation of Multiple Integrals), see gen/sphint.c
        (Actually only 16 of them are used because of symmetry.)
        The points are are first rotated by a random rotational matrix
    1: one random direction
    2: (x+y),(x-y),z separately
    3: x,y,z separately
    4: average over 3 orthogonal axes (rotated by a random rot matrix)
    5,6,7: x,y,z
  NOTES:
    modes 0,1,4 give the angle-averaged cross-section
    modes 2,3 give cross-sections in specified directions: suitable for
          recognizing orientational order, e.g. use mode=2 to see whether
          shear stress orders elongated molecules
  return value: sum of all cross sections; for modes 2,3 the 1st of three
    directions only
  xs.sizelimit: max size of the grid array
  xs.maxsize: max size actually encountered
  xs.Rscale: factor to get the collision diameter from the van der Waals
             radii sum
  xs.RvdW: probe van der Waals radius
*/
{
  static int    modes[8]=  {16,1,3,3,3,   1,1,1};
  static double weights[8]={0, 1,1,1,1./3,1,1,1};
  static int warned=0;
  int overflow;

  molecule_t *mn;
  siteinfo_t *si;
  vector *r;

  int imode,ns,sp,n,i,dimx,dimy,fromx,tox,fromy,toy,ix;
  int freq=xs.freq,mode=xs.mode&7;
  real rr,dy,minx,maxx,miny,maxy,w;

  unsigned summesh;
  unsigned char *mesh;

  vector ray,rvec,nx,ny;
  typedef real vector2[2];
  vector2 *proj;
#if XSECTION_M==2
  int imark,clsize,nssum,maxnssum=0,maxclsize,maxmark;
  double sum,maxsum=0; /* initialized to suppress compile warning */
#elif XSECTION_M==1
  vector *sum;
  int maxs=xs.maxs;
  int *sumn;
  double Sum=0;
#else
  int maxs=xs.maxs;
  vector sum={0,0,0};
#endif
  int is,ifreq;

  static matrix sheardir={
    {0.70710678118654757, 0.70710678118654757, 0},
    {0.70710678118654757,-0.70710678118654757, 0},
    {0,                   0,                   1} };
  static matrix xyzdir={
    {1,0,0},
    {0,1,0},
    {0,0,1}};
  matrix R;

  if (freq<0) freq=-freq; else freq=1;

#if XSECTION_M==2
  prt("\n_____t_______#/tot__sites__molecs__Xsec[AA2]____KEY_____");
       /*    0.000    5/30      36   12      159.185 XSECCLUSTERS */

  loopto (imark,1,mark) { // BUG found in V 3.5b: was loop (imark,0,mark)
    clsize=nssum=0;
    sum=0;
    loop (n,0,No.N) if (imark==xs.color[n]) {
      clsize++;
      nssum+=molec[n].ns; }
    if (clsize>=xs.mincluster) {
      allocarray(proj,nssum);

#elif XSECTION_M==1
  if (maxs<0) maxs=nspec;
  Min(maxs,nspec)
  allocarrayzero(sum,maxs);
  allocarrayzero(sumn,maxs);

  loop (n,0,No.N) {
    mn=molec+n;
    ns=mn->ns;
    if (ns<2) continue; /* ignore monoatomics */
    sp=mn->sp;
    if (sp>=maxs) continue; /* ignore species */
    si=spec[sp]->si;
    r=rof(mn,A->rp);
    allocarray(proj,ns);
#else
    if (maxs<0) maxs=No.s;
    Min(maxs,No.s)
    allocarray(proj,No.s);
#endif

    loop (ifreq,0,freq) {

      w=weights[mode]; /* overwritten for mode=0 */

      loop (imode,0,modes[mode]) {

        if (mode==0 || mode==4) rndmat(R);

        switch (mode) {
          /* only each second vector of sphint is taken because of symmetry */
          case 0: mpl(ray,R,sphint[imode*2].n); w=sphint[imode*2].w*2; break;
          case 1: rndsphere(ray); break;
          case 2: VV(ray,=sheardir[imode]) break;
          case 3: VV(ray,=xyzdir[imode]) break;
          case 4: VV(ray,=R[imode]); break;
          case 5: VV(ray,=xyzdir[0]); break;
          case 6: VV(ray,=xyzdir[1]); break;
          case 7: VV(ray,=xyzdir[2]); break; }
        /* now, ray is the direction of the projection */

        do rndball(rvec); while (fabs(SCAL(rvec,ray))<0.0625);
        VECT(ny,ray,rvec)
        rr=sqrt(SQR(ny));
        VO(ny,/=rr)
        VECT(nx,ray,ny)
        /* now, nx and ny are axes in the projection plane, randomly rotated */

        /* PASS 1: determinimg min/max sizes of the mesh */
        maxx=maxy=-3e33;
        minx=miny=3e33;

#if XSECTION_M==2
        is=0;
        loop (n,0,No.N) if (imark==xs.color[n]) {
          mn=molec+n;
          ns=mn->ns;
          sp=mn->sp;
          si=spec[sp]->si;
          r=rof(mn,A->rp);
          loop (i,0,ns) {

#elif XSECTION_M==1
        loop (i,0,ns) {
          is=i;
#else /*# XSECTION_M */
        is=0;
        loop (n,0,No.N) {
          mn=molec+n;
          ns=mn->ns;
          sp=mn->sp;
          si=spec[sp]->si;
          r=rof(mn,A->rp);
          loop (i,0,ns) {
#endif /*#!XSECTION_M */
            rr=(sitedef[si[i].st].LJ[0].RvdW+xs.RvdW)*xs.Rscale;

            proj[is][0]=SCAL(nx,r[i]);
            proj[is][1]=SCAL(ny,r[i]);

            Max(maxx,proj[is][0]+rr)
            Max(maxy,proj[is][1]+rr)
            Min(minx,proj[is][0]-rr)
            Min(miny,proj[is][1]-rr)
#if XSECTION_M==2
              is++; } }
        if (is!=nssum) ERROR(("internal: is!=nsum %d %d",is,nssum))
#elif XSECTION_M==1
        }
#else
            is++; } }
        if (is!=No.s) ERROR(("internal: is!=No.s"))
#endif

        dimx=xs.grid*(maxx-minx)+1.0001;
        dimy=xs.grid*(maxy-miny)+1.0001;

        /* to prevent interference errors */
        minx-=rnd()/xs.grid;
        miny-=rnd()/xs.grid;

        overflow=dimx*dimy>xs.sizelimit;
        Max(xs.maxsize,dimx*dimy)
        if (overflow && warned<=4) {
          WARNING(("Xsec: grid array %dx%d too large - truncated, Xsec underestimated%s",dimx,dimy,
                   warned==4?"\n(more occurrences not reported)":""))
          /* not clear what's the best strategy if grid not enough */
          /* this one tries centering, good for some atoms trying to escape */

          do {
            double M=fmax(fmax(maxx,maxy),fmax(-minx,-miny));

            if (maxx/M>0.9) maxx*=.99;
            if (maxy/M>0.9) maxy*=.99;
            if (-minx/M>0.9) minx*=.99;
            if (-miny/M>0.9) miny*=.99;
            dimx=xs.grid*(maxx-minx);
            dimy=xs.grid*(maxy-miny);
          } while (dimx*dimy>xs.sizelimit);

/*.....prt("%g %g  %g %g",minx,maxx,miny,maxy);*/
          warned++; }

        /* CREATE MESH
           - indices of the mesh array are (ix,iy), 0<=ix<dimx, 0<=iy<dimy
           - they represent grid points ((ix+1)/grid, (iy+1)/grid)
           - the map function is mesh[ix*dimy+iy] */
        alloczero(mesh,dimx*dimy); /* 2D array, my own map function */

        /* PASS 2: projecting particles on the mesh */
#if XSECTION_M==2
        is=0;
        loop (n,0,No.N) if (imark==xs.color[n]) {
          mn=molec+n;
          ns=mn->ns;
          sp=mn->sp;
          si=spec[sp]->si;
          r=rof(mn,A->rp);
          loop (i,0,ns) {
#elif XSECTION_M==1
        loop (i,0,ns) {
          is=i;
#else
        is=0;
        loop (n,0,No.N) {
          mn=molec+n;
          ns=mn->ns;
          sp=mn->sp;
          si=spec[sp]->si;
          r=rof(mn,A->rp);
          loop (i,0,ns) {
#endif
            rr=(sitedef[si[i].st].LJ[0].RvdW+xs.RvdW)*xs.Rscale;
            proj[is][0]-=minx;
            proj[is][1]-=miny;
            fromx=(int)(xs.grid*(proj[is][0]-rr));
            tox  =(int)(xs.grid*(proj[is][0]+rr));
            if (fromx<0 || tox>dimx) {
              fromx=tox=0;
              if (!overflow) ERROR(("Xsec: xrange=%d %d",fromx,tox)) }

            rr*=rr*Sqr(xs.grid);
            loop (ix,fromx,tox) {
              dy=rr-Sqr(ix+1-proj[is][0]*xs.grid);
              if (dy>=0) {
                dy=sqrt(dy);
                fromy=(int)(xs.grid*proj[is][1]-dy);
                toy  =(int)(xs.grid*proj[is][1]+dy);
                if (fromy<0 || toy>dimy) {
                  fromy=toy=0;
                  if (!overflow) ERROR(("Xsec: yrange=%d %d",fromy,toy)) }
                if (toy-fromy) memset(mesh+ix*dimy+fromy,1,toy-fromy); } }

#if XSECTION_M==2
              is++; } }
        if (is!=nssum) ERROR(("internal: is!=nssum"))
#elif XSECTION_M==1
              }
#else
            is++; } }
        if (is!=No.s) ERROR(("internal: is!=No.s"))
#endif

      summesh=0;
      loop (i,0,dimx*dimy) summesh+=mesh[i];

#if 0
      /* DEBUG: simple print of screen */
      {
        /* DEBUG */
        static int count;
        int i,j;
        loop (j,0,dimy) {
          _n
          loop (i,0,dimx)
          if (mesh[i*dimy+j]) prtc('%'); else prtc('.'); }
        _n
        if (count++>2000) exit(0);
      }
#endif /*# 0 */

#if XSECTION_M==2
      if (mode==2 || mode==3) ERROR(("xs.mode=2,3 not supported for clusters"))
      else sum+=summesh*w;
#elif XSECTION_M==1
      if (mode==2 || mode==3) sum[sp][imode]+=summesh*w;
      else sum[sp][0]+=summesh*w;
#else
      if (mode==2 || mode==3) sum[imode]+=summesh*w;
      else sum[0]+=summesh*w;
#endif

      free(mesh); } /* imode */
    } /* ifreq */

#if XSECTION_M==2
    prt("%9.3f %4d/%-4d %5d %4d %12.3f XSECCLUSTERS",
         t,imark,mark,nssum,clsize,sum/(freq*Sqr(xs.grid)));
    if (nssum>maxnssum) {
      maxnssum=nssum;
      maxsum=sum/(freq*Sqr(xs.grid));
      maxmark=imark;
      maxclsize=clsize; }
    free(proj);
  } /* clsize>=xs.mincluster */
  } /* imark */
  if (maxnssum)
    prt("%9.3f %4d/%-4d %5d %4d %12.3f MAXCLUSTER (by # of sites)",
        t,maxmark,mark,maxnssum,maxclsize,maxsum);
  _n
  StaAdd("Xsection[max.cluster]",maxsum);
  return maxsum;

#elif XSECTION_M==1
    sumn[sp]++;
    free(proj); } /* n */

  loop (sp,0,maxs) if (sumn[sp]) {
    Sum+=sum[sp][0];
    if (mode==2 || mode==3) loop (i,0,3) {
      char *s=string("%s Xsection[%d]",
                     i*4+(mode==2 ? "x+y\0x-y\0z" : "x\0HAy\0HAz"),sp);
      rr=sum[sp][i]/(freq*sumn[sp]*Sqr(xs.grid));
      StaAdd(s,rr);
      En.xsec[i]=rr;
      if (0) prt_("%s=%g  %c",s,rr,"\0\n"[i==2]); }
    else {
      char *s=string("Xsection[%d]",sp);
      rr=sum[sp][0]/(freq*sumn[sp]*Sqr(xs.grid));
      StaAdd(s,rr);
      if (0) prt("%s=%g",s,rr); } }

  free(sumn);
  free(sum);

  return Sum/Sqr(xs.grid);

#else
  if (mode==2 || mode==3)
    loop (i,0,3) {
      char *s=string("%s Xsection[cfg]",
               i*4+(mode==2 ? "x+y\0x-y\0z" : "x\0HAy\0HAz"));
      rr=sum[i]/((double)freq*Sqr(xs.grid));
      StaAdd(s,rr);
      if (0) prt_("%s=%g  %c",s,rr,"\0\n"[i==2]); }
  else {
    char *s="Xsection[cfg]";
    rr=sum[0]/(freq*Sqr(xs.grid));
    StaAdd(s,rr);
    if (0) prt("%s=%g",s,rr); }

  free(proj);

  return rr;
#endif
}
