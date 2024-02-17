static double Upot(vector *r,vector *f,species_t *spec,int toprint) /** Upot */
/***
    returns energy of molecule spec in configuration r
    forces are returned in f
    SHAKE-like constraint forces (not very good)
***/
{
  int i,j0,j1,j,k,l;
  double U,Ubond,Uangle,Udih,Uimp,Uaro;
#ifdef POLAR
  static double avrate=0;
  double Uself;
#endif /*# POLAR */
  int ns=spec->ns;
  site_t *site=spec->site;
  bond_t *b;
  angle_t *a;
  torsion_t *t;
  exception_t *exc;
  int cdyn=constraints;

  TRACE('.')

  nUpot++;
  /* test */
  if (testout) loop (i,0,ns) loop (j,0,3) if (fabs(r[i][j])>UNDEFATOM) {
    static int pass=0;

    if (pass++>20) ERROR(("too many atoms out of range"))
    WARNING(("%s: atom %d %c=%g out of range .. trying to fix", spec->fn,i,'x'+j,r[i][j]))
    /* trying to fix */
    VO(r[i],=71*rndgauss()); }

  makesym(r);

#ifdef CUTOFF
#  define SSPOT14 U+=sspot(r[j1], f[i],f[j1], &site[i],&site[j1],1)
#  define SSPOT   U+=sspot(r[j], f[i],f[j], &site[i],&site[j],0)
#  define PREPARE \
  loop (k,0,3) { \
    cutoff.r0[k]=r[i][k]; \
    cutoff.r0plusC2[k]=r[i][k]+cutoff.C2; \
    cutoff.r0minusC2[k]=r[i][k]-cutoff.C2; }

#  define IF_IN_CUTOFF_CUBE0(J) \
        if (cutoff.is_not || \
           (r[J][0]<cutoff.r0plusC2[0]  \
         && r[J][0]>cutoff.r0minusC2[0] \
         && r[J][1]<cutoff.r0plusC2[1]  \
         && r[J][1]>cutoff.r0minusC2[1] \
         && r[J][2]<cutoff.r0plusC2[2]  \
         && r[J][2]>cutoff.r0minusC2[2]) )

#  ifdef OMITKEPT
#    define IF_IN_CUTOFF_CUBE(J) \
        if ( !(keepmask & site[i].keep & site[J].keep) \
         && (cutoff.is_not || \
           (r[J][0]<cutoff.r0plusC2[0]  \
         && r[J][0]>cutoff.r0minusC2[0] \
         && r[J][1]<cutoff.r0plusC2[1]  \
         && r[J][1]>cutoff.r0minusC2[1] \
         && r[J][2]<cutoff.r0plusC2[2]  \
         && r[J][2]>cutoff.r0minusC2[2]) ) )
#  else /*# OMITKEPT */
#    define IF_IN_CUTOFF_CUBE IF_IN_CUTOFF_CUBE0
#  endif /*#!OMITKEPT */

  if (spec->C1>=spec->C2)
    ERROR(("C1=%g C2=%g bad cutoff range",spec->C1,spec->C2))
  if ( !(cutoff.is_not=spec->C1<2 || spec->C2<2) ) {
    cutoff.C1=spec->C1;
    cutoff.C1q=Sqr(cutoff.C1);
    cutoff.C2=spec->C2;
    cutoff.C2q=Sqr(cutoff.C2);
    cutoff.alpha=(cutoff.C2q+cutoff.C1q)/2;
    cutoff.beta=(cutoff.C2q-cutoff.C1q)/2;
    cutoff.f=-15./8/(Pow4(cutoff.beta)*cutoff.beta); /* =2A */
    cutoff.betaq=Sqr(cutoff.beta);
    cutoff.s0=4./15*cutoff.f*Sqr(cutoff.betaq);
    cutoff.s1=-2./15*cutoff.f*cutoff.betaq;
    cutoff.s2=cutoff.f/10; }
  else
    cutoff.C1q=cutoff.C2q=3e33;

#else /*# CUTOFF */

#  define PREPARE /* nothing */
#  ifndef POLAR /* POLAR version defined later */
#    define SSPOT14 U+=sspot(r[i],r[j1], f[i],f[j1], &site[i],&site[j1],1)
#    define SSPOT   U+=sspot(r[i],r[j], f[i],f[j], &site[i],&site[j],0)
#  endif /*? POLAR: POLAR version defined later */ /*# POLAR */
#  define IF_IN_CUTOFF_CUBE(J) /* nothing */
#  define IF_IN_CUTOFF_CUBE0(J) /* nothing */

#endif /*#!CUTOFF */

 AGAIN:
  
#ifdef POLAR
  Uself=
#endif /*# POLAR */
  U=Ubond=Uangle=Udih=Uimp=Uaro=ULJ=Uel=Uanglemax=0;
  verbose2=option('v')&8;

/* NEW: negative spec->opt_u = print of min energies */
  maxpairs=abs(spec->opt_u);
  signpairs=(double)(spec->opt_u/abs(spec->opt_u+!spec->opt_u));
  if ( (checkpairs = toprint*(maxpairs || verbose2)) )
    loop (i,0,maxpairs) spec->maxpair[i].U=-3e33;
  maxpair=spec->maxpair;

  measure=1;
  memset(f,0,ns*sizeof(vector));

  /* special probe patch */

  if (spec->probe.ns) {
#ifdef POLAR
    ERROR(("probe not implemented for POLAR"))
#else /*# POLAR */
    int source=ns-spec->probe.ns;

    if (source!=spec->probe.i)
      ERROR(("ns=%d probe.ns=%d probe.i=%d",ns,spec->probe.ns,spec->probe.i))

    loop (i,source,ns) {
      PREPARE
      loop (j,0,source)
        IF_IN_CUTOFF_CUBE0(j) SSPOT; }

    if (spec->probe.ns==3) {
      /* TIP3P water only -- note that there are two bonds and one angle
	 on the top of the the lists id this water has been added by a
	 `pw' statement (generally, the end of the probe should be marked) */
      i=b0->indx[0]; j=b0->indx[1];
      Ubond += bondpot(r[i],r[j], f[i],f[j], &b0->parm);
      b=b0->next;
      i=b->indx[0]; j=b->indx[1];
      Ubond += bondpot(r[i],r[j], f[i],f[j], &b->parm);
      i=a0->indx[0]; j=a0->indx[1]; k=a0->indx[2];
      Uangle += anglepot(r[i],r[j],r[k], f[i],f[j],f[k], &a0->parm); }
#endif /*#!POLAR */
    }
  else

    /* standard version (no patch probe):  */
    {
#ifdef POLAR
#  include "blendpol.c"
      /* NB: #includes also blendssi.c, blendbon.c */
      U+=Uself; /* changed in V2.3e */
#else /*# POLAR */
      /* site-site non-polar version -- see SSPOT and SSPOT14 above */
#  include "blendssi.c"
      /* bonded interactions */
#  include "blendbon.c"
#endif /*#!POLAR */
    }

  U+=Ubond+Uangle+Udih+Uimp+Uaro;

  if (toprint) {
    if (checkpairs) {
      if (maxpairs) prt("!=== %spair energies. * = 1-%d pairs, f = nbfixes ===",
	  "min \0\0\0\0max "+4+(int)signpairs*4,distance14);
      loop (i,0,maxpairs) {
	if (maxpair[i].U*signpairs<-2e33) break;
	prt("!%s  U=%g",maxpair[i].info,maxpair[i].U*signpairs); }
      _n }

#if defined(DEBUG) && defined(POLAR)
    /* ??? - cf. the same code in blendpol.c */
    prt_("! self-field: %d iterations, err=%g\n\
! rates=",
	 polnit,polerr);
    loop (i,0,NRATES) prt_("%8.5f", rate[i]);
    _n
#endif /*# defined(DEBUG) && defined(POLAR) */

    prts("! potential energy summary (in kcal/mol)");
#ifdef POLAR
    prt("! POLAR: Coulomb energy is without polarizable dipole self-energy");
    prt("!      SCF rate: %10.5f  self-energy:%11.5f   dihedrals: %11.5f",avrate,Uself,Udih);
#else /*# POLAR */
    prt("!                                                      dihedrals: %11.5f",Udih);
#endif /*#!POLAR */
    prt("! Lennard-Jones:%11.5f        bonds:%11.5f   impropers: %11.5f",ULJ,Ubond,Uimp);
    prt("!       Coulomb:%11.5f       angles:%11.5f   aromatics: %11.5f",Uel,Uangle,Uaro);
    prt("! sum nonbonded:%11.5f   sum bonded:%11.5f     total U: %11.5f",
                         ULJ+Uel,Ubond+Uangle+Udih+Uimp+Uaro,  U);
#ifdef CUTOFF
    if (!cutoff.is_not) prt("! cutoff range: [%g,%g]",cutoff.C1,cutoff.C2);
#endif /*# CUTOFF */
  }

  if (cdyn) {
    /* constraints and constraint forces - SHAKE-like */
    int MAXIT=spec->Xopt.S*spec->ns;
    fixdih_t *fd;

    bond_t *b;
    angle_t *a;
    double maxrr,xrr=1e-5;
    int icyc=0;

#ifdef OMITKEPT
    ERROR(("sorry, OMITKEPT version not compatible with constraints"))
#endif /*# OMITKEPT */

    do {
      maxrr=0;

/*****************************************************************************
   constraints - theory
   ^^^^^^^^^^^^^^^^^^^^

   Let omega({r_i}) be a function of positions {r_i}; e.g., omega=cos(phi)

   We want omega({r'_i})=omega_0. Linearization:
     r'=r+dr
   We will look for dr in the direction of constraint forces (gradients):
     dr_i = x grad_i omega
     omega({r'_i}) = omega({r_i}) + sum_i dr_i grad_i omega
                   = omega({r_i}) + x sum_i (grad_i omega)^2
   hence
     x = [omega_0-omega({r_i})]/[sum_i (grad_i omega)^2]
     r'_i = r_i + x grad_i omega

   Correcting forces - again in the direction of grad omega:
     f'_i = f_i + alpha grad_i omega
   omega should not change after a small diplacement in the direction of f'_i:
     omega({r_i + h f'_i}) = omega({r_i}) + h sum_i f'_i.grad_i omega
                           =! omega({r_i})
   So we need 
     0 = sum_i f'_i.grad_i omega 
       = sum_i [f_i + alpha grad_i omega].grad_i omega
   hence
     alpha = -[sum_i f_i.grad_i omega]/[sum_i (grad_i omega)^2]
     f'_i = f_i + alpha grad_i omega
*****************************************************************************/

      if (option('b')==0) for (b=b0; b; b=b->next) {
        /* constrained bonds - exact solution */
        int i=b->indx[0], j=b->indx[1];
        vector dr,df;
        double rr,x;

        VVV(dr,=r[i],-r[j])
        rr=SQR(dr);
        // x=(Sqr(b->parm.length)/SQR(dr)-1)/4; // approximate but faster
        x=(b->parm.length/sqrt(rr)-1)/2;
        VV(r[i],+=x*dr)
        VV(r[j],-=x*dr)
        Max(maxrr,fabs(x))
        VVV(df,=f[i],-f[j])

        x=SCAL(dr,df)/rr/2;
        VV(f[i],-=x*dr)
        VV(f[j],+=x*dr) }

      if (option('a')==0) for (a=a0; a; a=a->next) if (a->parm.angle>0) {
        /* constrained angles - linearized solution (see above) */
        /* BUG: angles 0 and 180 not supported */
        int i=a->indx[0], j=a->indx[1], k=a->indx[2];
        double ab,cf,norm,aa,bb,angle=a->parm.angle;
        double x,delta;
        vector a,b;
        vector grada,gradb,dra,drb;

        VVV(a,=r[i],-r[j])
        VVV(b,=r[k],-r[j])
        aa=SQR(a); bb=SQR(b);
        norm=sqrt(aa*bb);
        ab=SCAL(a,b);
        cf=ab/norm;

        if (fabs(cf)>=1)
          WARNING(("singular angle ignored: use randomization (e.g., hot key :)"))
        else {
          x=ab/aa; VVV(grada,=b,-x*a)
          x=ab/bb; VVV(gradb,=a,-x*b)

          x=cos(angle)-cf;
          while (fabs(x>0.1)) x*=0.6;
          delta=2*(SQR(grada)+SQR(gradb)+SCAL(grada,gradb));
          x/=norm*delta;
          Max(maxrr,fabs(x))
          VV(dra,=x*grada)
          VV(drb,=x*gradb)

          VV(r[i],+=dra)
          VVV(r[j],-=dra,+drb)
          VV(r[k],+=drb)

          x=(-SCAL(f[i],grada)+SCAL(f[j],grada)+SCAL(f[j],gradb)-SCAL(f[k],gradb))
            / delta;
          VV(dra,=x*grada)
          VV(drb,=x*gradb)
          VV(f[i],+=dra)
          VVV(f[j],-=dra,+drb)
          VV(f[k],+=drb) } }

      looplist (fd,fd0) {
        int i=fd->indx[0], j=fd->indx[1], k=fd->indx[2], l=fd->indx[3];

        vector a,b,c,ap,cp;
        vector grada,gradb,gradc;
        vector dra,drb,drc;
        /*        vector perpab;*/

        double ab,bb,bc,cf,norm,apq,cpq,x,delta;
        double abbb,bcbb;
#if 0
        double dphi, double sinphi;
#endif /*# 0 */
        VVV(a,=r[j],-r[i])
        VVV(b,=r[k],-r[j])
        VVV(c,=r[l],-r[k])
        ab=SCAL(a,b);
        bc=SCAL(b,c);
        bb=SQR(b);
        abbb=ab/bb; VVV(ap,= -a,+abbb*b)
        bcbb=bc/bb; VVV(cp,=  c,-bcbb*b)
        /* vectors ap and cp are perpendicular to b */

        apq=SQR(ap);
        cpq=SQR(cp);
        norm=1/sqrt(apq*cpq);
        /* cos phi */
        cf=SCAL(ap,cp)*norm;

        if (fabs(cf)>=1)
          WARNING(("singular dihedral ignored: use randomization (e.g., hot key :)"))
        else {

#if 0
          VECT(perpab,a,b) /* perpendicular to a,b */
          sinphi=SCAL(perpab,cp)/sqrt(SQR(perpab)*cpq);
          phi=atan2(sinphi,cf);
          dphi=phi-fd->phi;
          while (dphi>PI) dphi-=2*PI;
          while (dphi<-PI) dphi+=2*PI;
#endif /*# 0 */
          /* gradients: d cf/d r0 =      -grada
                        d cf/d r1 = grada-gradb
                        d cf/d r2 = gradb-gradc
                        d cf/d r3 = gradc       */
          x=cf/apq; VVV(grada,= -norm*cp, +x*ap)
          x=cf/cpq; VVV(gradc,=  norm*ap, -x*cp)
          VVV(gradb,= -abbb*grada, -bcbb*gradc)

          x=cos(fd->phi)-cf;
          while (fabs(x>0.1)) x*=0.6;
          delta=2*(SQR(grada)+SQR(gradb)+SQR(gradc)
                   -SCAL(grada,gradb)-SCAL(gradc,gradb));
          x/=-delta;
          Max(maxrr,fabs(x))
          VV(dra,=x*grada)
          VV(drb,=x*gradb)
          VV(drc,=x*gradc)

          VV(r[i],+=dra)
          VVV(r[j],-=dra,-drb)
          VVV(r[k],-=drb,-drc)
          VV(r[l],-=drc) 

#if 0
        VVV(a,=r[j],-r[i])
        VVV(b,=r[k],-r[j])
        VVV(c,=r[l],-r[k])
        ab=SCAL(a,b);
        bc=SCAL(b,c);
        bb=SQR(b);
        abbb=ab/bb; VVV(ap,= -a,+abbb*b)
        bcbb=bc/bb; VVV(cp,=  c,-bcbb*b)
        /* vectors ap and cp are perpendicular to b */

        apq=SQR(ap);
        cpq=SQR(cp);
        norm=1/sqrt(apq*cpq);
        /* cos phi */
        cf=SCAL(ap,cp)*norm;
        put2(acos(cf),acos(cf)-fd->phi)
            exit(0);
#endif /*# 0 */

        x=(-SCAL(f[i],grada)+SCAL(f[j],grada)-SCAL(f[j],gradb)+SCAL(f[k],gradb)-SCAL(f[k],gradc)+SCAL(f[l],gradc))
            / delta;
        VV(dra,=x*grada)
        VV(drb,=x*gradb)
        VV(drc,=x*gradc)
        VV(f[i],+=dra)
        VVV(f[j],-=dra,-drb)
        VVV(f[k],-=drb,-drc)
        VV(f[l],-=drc) 
      } }

      /* xrr=xrr*.8+maxrr*.2; */
      xrr=xrr*0.6+maxrr*0.4;
    } while (xrr>1e-14 && icyc++<MAXIT);

    if (xrr>1e-13) prt("WARNING: SHAKE not converged (%.2g) - results imprecise",xrr); }

  /* fixed sites */
  if (spec->opt_k)
    loop (i,0,ns) if (site[i].keep) memset(f[i],0,sizeof(vector));

  if (cdyn--) goto AGAIN;
  
  return U+spec->zero_energy;
}

static void prtU(int it,int by,double U,double h) /******************* prtU */
{
  const char fmt[]="%4i :  U=%-18.8g  h=%-g\n";
  if (it % by) return;
  if (out!=stdout) fprintf(stderr,fmt,it,U,h);
  TRACE('\n')
  prtc('!');
  prt_(fmt,it,U,h);
}
