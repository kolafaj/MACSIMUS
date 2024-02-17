    if (method<-1 || method>4) return 0;

    if (err==-2 && !isstdev) err=0;

    dody(err,1);
    Minimize(method,N,A,ff,maxit,eps,D,par);

    sqsum=Val(ff(A));
    sigma=sqrt(sqsum/(n-N));

    prt("sigma=%g, adj. to deg.f.=%g",(double)(sqrt(sqsum/n)),(double)sigma);

    {
      FILE *f=fopen("fit.fit","wt");
      loop (i,0,N) fprintf(f,"# %.15g\n",Val(A[i]));
      fprintf(f,"# x   y   y_fit   dy  dy/sigma\n");
      loop (i,0,n)
	fprintf(f,"%.8g %.8g %.8g %.8g %.8g\n",
		x[i],y[i],Val(func(A,x[i])),dy[i],Val((func(A,x[i])-y[i])/dy[i]));
      fclose(f);
    }

    loop (i,0,N) {
      prt("%c=%.15g",i+'a',Val(A[i]));
      A0[i]=A[i];
      Asum[i]=Aqsum[i]=0; }
     
    loop (i,0,n) y_0[i]=y[i];

    if (nerr) {
      FILE *ferr=fopen("fit.err","wt");
    
      dody(err,sigma);

#if 0
      // debug
      Minimize(method,N,A,ff,maxit,eps,D,par);
      loop (i,0,N)
	prt("A[%d]=%.15g (debug)",i,A[i]);
#endif

      loop (ierr,0,nerr) {
	loop (i,0,N) A[i]=A0[i];
	loop (i,0,n) y[i]=y_0[i]+rndgauss()*dy[i];

	Minimize(method,N,A,ff,maxit,eps,D,par);
	loop (i,0,N) {
	  Asum[i]+=A[i];
	  Aqsum[i]+=Sqr(A[i]); } }
     
      loop (i,0,n) y[i]=y_0[i];

      loop (i,0,N) {
	A[i]=A0[i];
	if (nerr) {
          REAL Aerr=xsqrt((Aqsum[i]/nerr-Sqr(Asum[i]/nerr))/(1.-1./nerr));
 	  prt("%c=%14.11f %14.11f %14.11f",
	         i+'a',  Val(A[i]),   Val(Asum[i]/nerr), Val(Aerr));
          fprintf(ferr," %14.11f %14.11f%c",
                         Val(A[i]),   Val(Aerr),  i==N-1?'\n':' '); }
        Asum[i]=Aqsum[i]=0; }

    fclose(ferr); } }

  return 0;
}
