/* DEBUG: table of eru, erd, and num.deriv of eru */
{
  double x,e,u1,u2,d=1e-4,dx=0.001,Eru,Erd,numder;
  double sumqu=0,sumqd=0,sumw=0,w=0;
  int ix;
  extern char *Fn(char*); /* bad habit */
  static int pass;

#ifdef QQTAB
  if (ga.qq==0)
    WARNING(("qq=0: omitted from dump"))
  else
    if (pass<100)
#else
  if (pass<5)
#endif
  {
    FILE *f=fopen(Fn(string("%d.ertest",pass)),"wt");

    pass++;
    fprintf(f,"#  grid=%d  alpha=%.15g  range=%.15g  cutoff=%.15g  sgrid=%g  shift=%d (not incl. in \"exact\")\n",
            Erfc.grid,Erfc.alpha,maxr,cutoff,Erfc.sgrid,shift);
#ifdef QQTAB
    fprintf(f,"#  QQTAB=%d qq=%g=%g e^2  qfactor_1=%g  sigma1=%g sigma2=%g  g12=%g\n",
            QQTAB,ga.qq,ga.qq/167100.75, ga.qfactor_1,
            sigma1,sigma2, ga.g12);
#endif
    fprintf(f,"#  deriv2 = 2nd-order numerical derivative, d=%g, r>0:symmetric formula; r=0:right\n",d);
    fprintf(f,"#  i=index of spline\n");
    fprintf(f,"#  r            eru                    erd             deriv2(eru)          exact_eru            exact_erd     i\n");

    loopto (ix,(int)(minr/dx+0.5),(int)(maxr/dx+0.5)) {
      x=ix*dx;
      if (x==0) {
        u1=eru(d*d)-eru(0);
        u2=eru(4*d*d)-eru(0);
        numder=(u1/7.5-u2/1.875)/(d*d); }                                                                           
      else {
        u1=eru(Sqr(x+d));
        u2=eru(Sqr(x-d));
        numder=(u2-u1)/(2*d*x); }
      e=erud(x*x);
#ifdef GAUSSIANCHARGES
      {
        double Erfc_h=Erfc.h; /* legacy workaround */
        Erfc.h=x*x;
        erfcs(1,&Eru,&Erd);
        Erfc.h=Erfc_h;
      }
#else
      Eru=erfc(x*alpha)/x;
      Erd=erfc(x*alpha)/Cub(x)+exp(-Sqr(x*alpha))*(2/SQRTPI)/Sqr(x)*alpha;
#endif
#if SUBGRID>0
      /* to test integration, put #if SUBGRID>0 above ... */
      /* if erfc() is not known, put #if SUBGRID==1 above and #define SUBGRID 4 */
      fprintf(f,"%.3f %20.15g %20.15g %20.15g %20.15g %20.15g %d\n",
              x,
              e,byerd,numder,
              Eru,Erd,
              (int)(x*x*Erfc.sgrid) );
      if (fabs(x-1)<0.5) {
        /* Gauus window, sig=0.1 */
        w=exp(-Sqr(x-1)*50);
        sumw+=w;
        sumqu+=w*Sqr(Eru-e);
        sumqd+=w*Sqr(Erd-byerd); }
#else
      fprintf(f,"%.3f %20.15g %20.15g %20.15g\n",x,e,byerd,(u2-u1)/(2*d*x) );
#endif
    }

    fprintf(f,"# Gauss-window-weighted error around (sig=0.1) distance r=1, incl. rshift\n\
# %.9g %.9g erru errd\n",sqrt(sumqu/w),sqrt(sumqd/w));
    fclose(f);
  }

#if 0
  {
    double s=0,x=1,xx;
    int i;
    extern double mytime(void);
    double t0=mytime();

    loop (i,0,100000000) {
      x+=0.3;
#  if 0
      xx=x;
      x=sqrt(x);
      s+=erfc(x*alpha)/(x*xx)+exp(-Sqr(x*alpha))*(2/SQRTPI)/xx*alpha;
#  else
      x=sqrt(x);
      s+=erud(x);
#  endif
    }
    prt("s=%g\ntime=%g s",s,mytime()-t0);
    exit(0);
  }
#endif

}

#if 0
    { /* point charges only */
      int i;
      double xx,y,e;

      loop (i,0,10000) {
        xx=i*0.001;
        e=exacterud_sqrt(xx,&y);

        prt("%g %.16g %.16g EEE1",xx,e,-erf(alpha*sqrt(xx))/sqrt(xx));
        e=exacterud_sqrt_1(xx,&y);
        prt("%g %.16g %.16g EEE2",xx,e,-erf(alpha*sqrt(xx))/sqrt(xx)+2/SQRTPI*alpha);
        prt("%g %.16g %.16g EEE3",xx,y,-erf(alpha*sqrt(xx))/sqrt(xx)/xx+alpha*2/SQRTPI/xx*exp(-alpha*alpha*xx));
      }
    }
#endif
